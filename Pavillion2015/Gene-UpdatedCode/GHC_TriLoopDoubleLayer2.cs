using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;

namespace Pavillion2015
{
    public class GHC_TriLoopDoubleLayer2 : GH_Component
    {

        SpringMesh iSpringMesh = null;
        double iMinThickness = double.NaN;
        double iMaxThickness = double.NaN;
        double iTangentScale = double.NaN;
        bool iPolySrf = true;

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

        List<Point3d> topCps = null;
        List<Point3d> bottomCps = null;
        List<Point3d> topCentres = null;
        List<Point3d> bottomCentres = null;

        double documentTolerance = DocumentTolerance();

        public GHC_TriLoopDoubleLayer2()
            : base("Tri-Loop Double Layer2", "Tri-Loop Double Layer2",
                "Tri-Loop Double Layer2 with Switch",
                "Pavillion 2015", "Double Layer")
        {
        }


        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Spring Mesh", "Spring Mesh", "Spring Mesh", GH_ParamAccess.item);
            pManager.AddNumberParameter("Min. Thickness", "Min. Thickness", "Min. Thickness", GH_ParamAccess.item, 0.3);
            pManager.AddNumberParameter("Max. Thickness", "Max. Thickness", "Max. Thickness", GH_ParamAccess.item, 1.0);
            pManager.AddNumberParameter("Tangent Scale", "Tangent Scale", "Tangent Scale", GH_ParamAccess.item, 0.2);
            pManager.AddBooleanParameter("Allow PolySrf", "Allow PolySrf", "Allow PolySrf", GH_ParamAccess.item, true);
        }


        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddTextParameter("Info", "00 - Info", "Information", GH_ParamAccess.item);
            pManager.AddGenericParameter("Debug 1", "01 - Debug 1", "This output is reserved for debugging", GH_ParamAccess.list);
            pManager.AddGenericParameter("Debug 2", "02 - Debug 2", "This output is reserved for debugging", GH_ParamAccess.list);
            pManager.AddGenericParameter("Debug 3", "03 - Debug 3", "This output is reserved for debugging", GH_ParamAccess.list);
            pManager.AddGenericParameter("Debug 4", "04 - Debug 4", "This output is reserved for debugging", GH_ParamAccess.list);
            pManager.AddGenericParameter("Debug 5", "05 - Debug 5", "This output is reserved for debugging", GH_ParamAccess.list);
            pManager.AddGenericParameter("Components", "Components", "Components", GH_ParamAccess.list);
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
            DA.GetData<SpringMesh>("Spring Mesh", ref iSpringMesh);
            DA.GetData<double>("Min. Thickness", ref iMinThickness);
            DA.GetData<double>("Max. Thickness", ref iMaxThickness);
            DA.GetData<double>("Tangent Scale", ref iTangentScale);
            DA.GetData<bool>("Allow PolySrf", ref iPolySrf);


            // =========================================================================================================

            List<Vector3d> vertexNormals = iSpringMesh.ComputeVertexNormals();
            
            foreach (Edge edge in iSpringMesh.Edges)
            {
                if (edge.SecondTriangleIndex >= 0) continue;

                // Gene Added
                if ( iSpringMesh.Vertices[edge.FirstVertexIndex].Position.Z > 0.65 &&
                     iSpringMesh.Vertices[edge.SecondVertexIndex].Position.Z > 0.65) continue;

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
            }


            DA.SetDataList(1, oDebugBreps1);
            DA.SetDataList(2, oDebugBreps2);
            //DA.SetDataList(2, oDebugCurves1);
            //DA.SetDataList(3, oDebugPoints1);
        }


        private void createStripe(int firstVertexIndex, int secondVertexIndex, int thirdVertexIndex)
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

            Curve profileCurve1 = Curve.CreateControlPointCurve(new List<Point3d>() { ab, aab, aA, AAB, AB });
            Curve profileCurve2 = Curve.CreateControlPointCurve(new List<Point3d>() { ac, aac, aA, AAC, AC });

            if (iPolySrf)
            {

                PolyCurve polyCurve1 = new PolyCurve();
                polyCurve1.Append(new LineCurve(m, ab));
                polyCurve1.Append(profileCurve1);
                polyCurve1.Append(new LineCurve(AB, M));

                PolyCurve polyCurve2 = new PolyCurve();
                polyCurve2.Append(new LineCurve(m, ac));
                polyCurve2.Append(profileCurve2);
                polyCurve2.Append(new LineCurve(AC, M));

                oDebugBreps1.Add(
                    Brep.CreateFromLoft(
                       new List<Curve>() { polyCurve1, polyCurve2 },
                       Point3d.Unset, Point3d.Unset,
                       LoftType.Normal,
                       false
                       )[0]
                   );

                oDebugCurves1.Add(profileCurve1);
                oDebugCurves2.Add(profileCurve2);
            }
            else
            {
                profileCurve1 = Curve.JoinCurves(
                    new List<Curve>() { new LineCurve(b, ab), profileCurve1, new LineCurve(AB, B) },
                    documentTolerance,
                    true)[0];

                profileCurve2 = Curve.JoinCurves(
                    new List<Curve>() { new LineCurve(c, ac), profileCurve2, new LineCurve(AC, C) },
                    documentTolerance,
                    true)[0];

                Brep brep = Brep.CreateFromLoft(
                       new List<Curve>() { profileCurve1, profileCurve2 },
                       Point3d.Unset, Point3d.Unset,
                       LoftType.Normal,
                       false
                       )[0];

                // =============================================
                // Trim the brep using planes
                // =============================================

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

                //PolyCurve polyCurve1 = new PolyCurve();
                //polyCurve1.Append(new LineCurve(m, ab));
                //polyCurve1.Append(profileCurve1);
                //polyCurve1.Append(new LineCurve(AB, M));

                //PolyCurve polyCurve2 = new PolyCurve();
                //polyCurve2.Append(new LineCurve(m, ac));
                //polyCurve2.Append(profileCurve2);
                //polyCurve2.Append(new LineCurve(AC, M));

                //oDebugBreps1.Add(
                //    Brep.CreateFromLoft(
                //       new List<Curve>() { polyCurve1, polyCurve2 },
                //       Point3d.Unset, Point3d.Unset,
                //       LoftType.Normal,
                //       false
                //       )[0]
                //   );

                oDebugBreps1.Add(brep);
                oDebugCurves1.Add(profileCurve1);
                oDebugCurves2.Add(profileCurve2);
            }
        }

        
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                return null;
            }
        }

        
        public override Guid ComponentGuid
        {
            get { return new Guid("{8e903084-d394-4e52-b7d2-80c13f5b5429}"); }
        }
    }
}