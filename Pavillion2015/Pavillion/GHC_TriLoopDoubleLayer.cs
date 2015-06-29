using System;
using System.Collections.Generic;
using System.Windows.Forms;

using Grasshopper.Kernel;
using Rhino.Geometry;

namespace Pavillion2015
{
    public class GHC_TriLoopDoubleLayer : GH_Component
    {
        SpringMesh iSpringMesh = null;
        double iMinThickness = double.NaN;
        double iMaxThickness = double.NaN;
        double iTangentScale = double.NaN;

        string oInfo = string.Empty;
        List<Point3d> oDebugPoints1 = null;
        List<Point3d> oDebugPoints2 = null;
        List<Curve> oDebugCurves1 = null;
        List<Curve> oDebugCurves2 = null;
        List<Vector3d> oDebugVectors1 = null;
        List<Vector3d> oDebugVectors2 = null;
        List<double> oDebugNumbers1 = null;
        List<double> oDebugNumbers2 = null;
        List<Brep> oDebugBreps1 = null;
        List<Brep> oDebugBreps2 = null;
        List<Brep> oStripes = null;
        List<Brep> oStructuralInteriorStripes = null;
        List<Brep> oStructuralExteriorStripes = null;

        List<Point3d> topCps= null;
        List<Point3d> bottomCps = null;

        ToolStripMenuItem menuItem_Standard = null;
        ToolStripMenuItem menuItem_StructuralAnalysis = null;
        ToolStripMenuItem menuItem_ExtendedSurfaces = null;

        bool menuItemValue_Standard = true;
        bool menuItemValue_StructuralAnalysis = false;
        bool menuItemValue_ExtendedSurfaces = false;

        double documentTolerance = DocumentTolerance();

        public GHC_TriLoopDoubleLayer()
            : base("Tri-Loop Double Layer", "Tri-Loop Double Layer",
                "Tri-Loop Double Layer",
                "Pavillion 2015", "Double Layer")
        {
        }


        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Spring Mesh", "Spring Mesh", "Spring Mesh", GH_ParamAccess.item);
            pManager.AddNumberParameter("Min. Thickness", "Min. Thickness", "Min. Thickness", GH_ParamAccess.item, 0.3);
            pManager.AddNumberParameter("Max. Thickness", "Max. Thickness", "Max. Thickness", GH_ParamAccess.item, 1.0);
            pManager.AddNumberParameter("Tangent Scale", "Tangent Scale", "Tangent Scale", GH_ParamAccess.item, 0.2);
        }


        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddTextParameter("Info", "00 - Info", "Information", GH_ParamAccess.item);
            pManager.AddGenericParameter("Debug 1", "01 - Debug 1", "This output is reserved for debugging", GH_ParamAccess.list);
            pManager.AddGenericParameter("Debug 2", "02 - Debug 2", "This output is reserved for debugging", GH_ParamAccess.list);
            pManager.AddGenericParameter("Debug 3", "03 - Debug 3", "This output is reserved for debugging", GH_ParamAccess.list);
            pManager.AddGenericParameter("Debug 4", "04 - Debug 4", "This output is reserved for debugging", GH_ParamAccess.list);
            pManager.AddGenericParameter("Debug 5", "05 - Debug 5", "This output is reserved for debugging", GH_ParamAccess.list);
            pManager.AddGenericParameter("Stripes", "06 - Stripes (Continous)", "Stripes (Continous)", GH_ParamAccess.list);
            pManager.AddGenericParameter("Stripes (Structural Interior)", "07 - Stripes (Structural Interior)", "Stripes (Structural Interior)", GH_ParamAccess.list);
            pManager.AddGenericParameter("Stripes (Structural Exterior)", "08 - Stripes (Structural Exterior)", "Stripes (Structural Exterior)", GH_ParamAccess.list);
        }


        protected override void BeforeSolveInstance()
        {
            oInfo = string.Empty;
            oDebugPoints1 = new List<Point3d>();
            oDebugPoints2 = new List<Point3d>();
            oDebugCurves1 = new List<Curve>();
            oDebugCurves2 = new List<Curve>();
            oDebugVectors1 = new List<Vector3d>();
            oDebugVectors2 = new List<Vector3d>();
            oDebugNumbers1 = new List<double>();
            oDebugNumbers2 = new List<double>();
            oDebugBreps1 = new List<Brep>();
            oDebugBreps2 = new List<Brep>();
        }


        protected override void SolveInstance(IGH_DataAccess DA)
        {
            if (!menuItemValue_Standard && !menuItemValue_StructuralAnalysis && !menuItemValue_ExtendedSurfaces)
                return;

            oStripes = new List<Brep>();
            oStructuralInteriorStripes = new List<Brep>();
            oStructuralExteriorStripes = new List<Brep>();

            DA.GetData<SpringMesh>("Spring Mesh", ref iSpringMesh);
            DA.GetData<double>("Min. Thickness", ref iMinThickness);
            DA.GetData<double>("Max. Thickness", ref iMaxThickness);
            DA.GetData<double>("Tangent Scale", ref iTangentScale);

            
            // =========================================================================================================

            List<Vector3d> vertexNormals = iSpringMesh.ComputeVertexNormals();

            foreach (Edge edge in iSpringMesh.Edges)
            {
                if (edge.SecondTriangleIndex >= 0) continue;

                Vector3d vertexNormal = vertexNormals[edge.FirstVertexIndex];
                vertexNormal.Z = 0.0;
                vertexNormal.Unitize();
                vertexNormals[edge.FirstVertexIndex] = vertexNormal;

                vertexNormal = vertexNormals[edge.SecondVertexIndex];
                vertexNormal.Z = 0.0;
                vertexNormal.Unitize();
                vertexNormals[edge.SecondVertexIndex] = vertexNormal;
            }


            // =========================================================================================================


            bottomCps = new List<Point3d>();
            topCps = new List<Point3d>();


            for (int i = 0; i < iSpringMesh.Vertices.Count; i++)
            {
                Point3d vertexPosition = iSpringMesh.Vertices[i].Position;
                double k = vertexPosition.Z * 0.125;
                double thickness = k * iMinThickness + (0.5 - k) * iMaxThickness;

                bottomCps.Add(vertexPosition - thickness * vertexNormals[i]);
                topCps.Add(vertexPosition + thickness * vertexNormals[i]);
            }


            // =========================================================================================================


            foreach (Triangle triangle in iSpringMesh.Triangles)
            {
                createStripe(triangle.FirstVertexIndex, triangle.SecondVertexIndex, triangle.ThirdVertexIndex);
                createStripe(triangle.SecondVertexIndex, triangle.ThirdVertexIndex, triangle.FirstVertexIndex);
                createStripe(triangle.ThirdVertexIndex, triangle.FirstVertexIndex, triangle.SecondVertexIndex);
                if (menuItemValue_StructuralAnalysis)
                    createSharedTriangles(triangle.FirstVertexIndex, triangle.SecondVertexIndex, triangle.ThirdVertexIndex);
            }


            // =========================================================================================================


            oInfo += "Triloop Segments: " + iSpringMesh.Triangles.Count;
            oInfo += "\nWood Stripes: " + iSpringMesh.Triangles.Count * 3;


            double totalArea = 0.0;
            if (menuItemValue_Standard)
                foreach (Brep brep in oStripes)
                    totalArea += brep.GetArea();
            else if (menuItemValue_StructuralAnalysis)
            {
                foreach (Brep brep in oStructuralInteriorStripes)
                    totalArea += brep.GetArea();
                foreach (Brep brep in oStructuralInteriorStripes)
                    totalArea += brep.GetArea();
            }
            oInfo += "\nTotal (Wood) Area: " + totalArea.ToString();


            // =========================================================================================================


            DA.SetData(0, oInfo);
            DA.SetDataList(1, oDebugCurves1);
            DA.SetDataList(2, oDebugCurves2);
            DA.SetDataList(6, oStripes);
            DA.SetDataList(7, oStructuralInteriorStripes);
            DA.SetDataList(8, oStructuralExteriorStripes);
        }


        protected override void AppendAdditionalComponentMenuItems(ToolStripDropDown menu)
        {
            base.AppendAdditionalComponentMenuItems(menu);

            Menu_AppendSeparator(menu);
            Menu_AppendItem(menu, "Recompute", MenuHandler_Recompute);

            Menu_AppendSeparator(menu);
            menuItem_Standard = Menu_AppendItem(menu, "Standard", MenuHandler_Standard, true, false);
            menuItem_StructuralAnalysis = Menu_AppendItem(menu, "Structural Analysis", MenuHandler_StructuralAnalysis, true, false);
            menuItem_ExtendedSurfaces = Menu_AppendItem(menu, "Extended Surfaces", MenuHandler_ExtendedSurfaces, true, false);

            //menuItem_Standard.CheckOnClick = true;
            //menuItem_StructuralAnalysis.CheckOnClick = true;
            //menuItem_ExtendedSurfaces.CheckOnClick = true;
        }


        private void MenuHandler_Recompute(Object sender, EventArgs e)
        {
            ExpireSolution(true);
        }


        private void MenuHandler_Standard(Object sender, EventArgs e)
        {
            menuItemValue_Standard = !menuItemValue_Standard;
            //menuItem_Standard.Checked = menuItemValue_Standard;          
        }


        private void MenuHandler_StructuralAnalysis(Object sender, EventArgs e)
        {
            menuItemValue_StructuralAnalysis = !menuItemValue_StructuralAnalysis;
            //menuItem_StructuralAnalysis.Checked = menuItemValue_StructuralAnalysis;
        }


        private void MenuHandler_ExtendedSurfaces(Object sender, EventArgs e)
        {
            menuItemValue_ExtendedSurfaces = !menuItemValue_ExtendedSurfaces;
            //menuItem_ExtendedSurfaces.Checked = menuItemValue_ExtendedSurfaces;
        }


        private void createSharedTriangles(int firstVertexIndex, int secondVertexIndex, int thirdVertexIndex)
        {
            Point3d a = bottomCps[firstVertexIndex];
            Point3d b = bottomCps[secondVertexIndex];
            Point3d c = bottomCps[thirdVertexIndex];

            Point3d ab = 0.5 * (a + b);
            Point3d bc = 0.5 * (b + c);
            Point3d ca = 0.5 * (c + a);         

            oStructuralInteriorStripes.Add(
                Brep.CreateFromLoft(
                       new List<Curve>() { new LineCurve(bc, ca), new LineCurve(bc, ab) },
                       Point3d.Unset, Point3d.Unset,
                       LoftType.Normal,
                       false
                       )[0]);
                

            Point3d A = topCps[firstVertexIndex];
            Point3d B = topCps[secondVertexIndex];
            Point3d C = topCps[thirdVertexIndex];

            Point3d AB = 0.5 * (A + B);
            Point3d BC = 0.5 * (B + C);
            Point3d CA = 0.5 * (C + A);

            oStructuralExteriorStripes.Add(
                Brep.CreateFromLoft(
                       new List<Curve>() { new LineCurve(BC, AB), new LineCurve(BC, CA) },
                       Point3d.Unset, Point3d.Unset,
                       LoftType.Normal,
                       false
                       )[0]);
        }

        private void createStripe(int firstVertexIndex,  int secondVertexIndex, int thirdVertexIndex)
        {           
            Point3d a = bottomCps[firstVertexIndex];
            Point3d A = topCps[firstVertexIndex];

            Point3d aA = 0.5 * (a + A);

            Point3d b = bottomCps[secondVertexIndex];
            Point3d B = topCps[secondVertexIndex];

            Point3d c = bottomCps[thirdVertexIndex];
            Point3d C = topCps[thirdVertexIndex];

            Point3d ab = 0.5 * (a + b);
            Point3d AB = 0.5 * (A + B);

            Point3d ac = 0.5 * (a + c);
            Point3d AC = 0.5 * (A + C);

            Point3d m = 0.33333 * (ab + ac + 0.5 * (b + c));
            Point3d M = 0.33333 * (AB + AC + 0.5 * (B + C));

            Point3d AAB = iTangentScale * A + (1.0 - iTangentScale) * AB;
            Point3d aab = iTangentScale * a + (1.0 - iTangentScale) * ab;

            Point3d AAC = iTangentScale * A + (1.0 - iTangentScale) * AC;
            Point3d aac = iTangentScale * a + (1.0 - iTangentScale) * ac;

            Curve profileCurve1 = Curve.CreateControlPointCurve(new List<Point3d>() { AB, AAB, aA, aab, ab });
            Curve profileCurve2 = Curve.CreateControlPointCurve(new List<Point3d>() { AC, AAC, aA, aac, ac });

            
         

            // =============================================
            // Stripe as one single (trimmed) surface
            // =============================================

            if (menuItemValue_Standard)
            {
                Curve profileCurveExtended1 = Curve.JoinCurves(
                    new List<Curve>() { new LineCurve(B, AB), profileCurve1, new LineCurve(ab, b) },
                    documentTolerance,
                    true)[0];

                Curve profileCurveExtended2 = Curve.JoinCurves(
                    new List<Curve>() { new LineCurve(C, AC), profileCurve2, new LineCurve(ac, c) },
                    documentTolerance,
                    true)[0];

                Brep brep = Brep.CreateFromLoft(
                       new List<Curve>() { profileCurveExtended1, profileCurveExtended2 },
                       Point3d.Unset, Point3d.Unset,
                       LoftType.Normal,
                       false
                       )[0];

                double offsetAmount = 0.2;
                Vector3d normal;
                Brep[] breps;

                normal = Vector3d.CrossProduct(B - A, C - A);

                Point3d A_ = A + offsetAmount * normal;

                breps = brep.Trim(new Plane(A_, AB, M), documentTolerance);
                if (breps.Length > 0) brep = breps[0];
                else { oDebugBreps2.Add(brep); return; }

                breps = brep.Trim(new Plane(A_, M, AC), documentTolerance);
                if (breps.Length > 0) brep = breps[0];
                else { oDebugBreps2.Add(brep); return; }

                normal = Vector3d.CrossProduct(b - a, c - a);

                Point3d a_ = a - offsetAmount * normal;

                breps = brep.Trim(new Plane(a_, m, ab), documentTolerance);
                if (breps.Length > 0) brep = breps[0];
                else { oDebugBreps2.Add(brep); return; }

                breps = brep.Trim(new Plane(a_, ac, m), documentTolerance);
                if (breps.Length > 0) brep = breps[0];
                else { oDebugBreps2.Add(brep); return; }

                oStripes.Add(brep);
            }

            // =============================================
            // Stripe as three separate surfaces
            // =============================================

            //profileCurve1.Domain = new Interval(0.0, 1.0);
            //profileCurve2.Domain = new Interval(0.0, 1.0);

            //Curve[] profileCurves1 = profileCurve1.Split(0.5);
            //Curve[] profileCurves2 = profileCurve2.Split(0.5);

            //oStructuralInteriorStripes.Add(
            //    Brep.CreateFromLoft(
            //       new List<Curve>() { new LineCurve(m, ac), new LineCurve(m, ab) },
            //       Point3d.Unset, Point3d.Unset,
            //       LoftType.Normal,
            //       false
            //       )[0]
            //   );

            //oStructuralInteriorStripes.Add(
            //    Brep.CreateFromLoft(
            //       new List<Curve>() { profileCurves1[1], profileCurves2[1] },
            //       Point3d.Unset, Point3d.Unset,
            //       LoftType.Normal,
            //       false
            //       )[0]
            //   );

            //oStructuralExteriorStripes.Add(
            //    Brep.CreateFromLoft(
            //       new List<Curve>() { profileCurves1[0], profileCurves2[0] },
            //       Point3d.Unset, Point3d.Unset,
            //       LoftType.Normal,
            //       false
            //       )[0]
            //   );

            //oStructuralExteriorStripes.Add(
            //    Brep.CreateFromLoft(
            //       new List<Curve>() { new LineCurve(M, AB), new LineCurve(M, AC) },
            //       Point3d.Unset, Point3d.Unset,
            //       LoftType.Normal,
            //       false
            //       )[0]
            //   );


            // =============================================
            // Stripe without triangles
            // =============================================

            if (menuItemValue_StructuralAnalysis)
            {
                profileCurve1.Domain = new Interval(0.0, 1.0);
                profileCurve2.Domain = new Interval(0.0, 1.0);

                Curve[] profileCurves1 = profileCurve1.Split(0.5);
                Curve[] profileCurves2 = profileCurve2.Split(0.5);

                oStructuralInteriorStripes.Add(
                    Brep.CreateFromLoft(
                       new List<Curve>() { profileCurves1[1], profileCurves2[1] },
                       Point3d.Unset, Point3d.Unset,
                       LoftType.Normal,
                       false
                       )[0]
                   );

                oStructuralExteriorStripes.Add(
                    Brep.CreateFromLoft(
                       new List<Curve>() { profileCurves1[0], profileCurves2[0] },
                       Point3d.Unset, Point3d.Unset,
                       LoftType.Normal,
                       false
                       )[0]
                   );
            }
        }


        protected override System.Drawing.Bitmap Icon {get{return null;}}
        public override Guid ComponentGuid {get { return new Guid("{cbca9974-644c-4fbd-a1ae-78d2576661ea}");}}
    }
}