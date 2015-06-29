using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;

namespace Pavillion2015
{
    public class GHC_PlanarPlateExperimental : GH_Component
    {

        SpringMesh iSpringMesh = null;
        double iMinThickness = double.NaN;
        double iMaxThickness = double.NaN;

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
        List<Line> oDebugLines1 = null;
        List<Line> oDebugLines2 = null;

        public GHC_PlanarPlateExperimental()
            : base("Plates Double Layer (Experimetial)", "Plates Double Layer (Experimetial)",
                "Planar Plates Double Layer",
                "Pavillion 2015", "Double Layer")
        {
        }

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Spring Mesh", "Spring Mesh", "Spring Mesh", GH_ParamAccess.item);
            pManager.AddNumberParameter("Min. Thickness", "Min. Thickness", "Min. Thickness", GH_ParamAccess.item, 0.3);
            pManager.AddNumberParameter("Max. Thickness", "Max. Thickness", "Max. Thickness", GH_ParamAccess.item, 1.0);
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
            oDebugLines1 = new List<Line>();
            oDebugLines2 = new List<Line>();
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            DA.GetData<SpringMesh>("Spring Mesh", ref iSpringMesh);
            DA.GetData<double>("Min. Thickness", ref iMinThickness);
            DA.GetData<double>("Max. Thickness", ref iMaxThickness);

            //---------------------------------------------------------------------------------------------------------------

            List<double> planarity = checkCoplanarity(iSpringMesh, iMinThickness, oDebugLines1, oDebugLines2);

            DA.SetDataList(1, oDebugLines1);
            DA.SetDataList(2, oDebugLines2);
            DA.SetDataList(3, planarity);

            plate();
            DA.SetDataList(4, oDebugCurves1);
        }

        private void plate()
        {

            //vertex thickness value
            List<double> thickness = new List<double>();
            foreach (Vertex vertex in iSpringMesh.Vertices)
                thickness.Add(0.5);

            List<List<LineCurve>> dualVertices = new List<List<LineCurve>>();

            for (int i = 0; i < iSpringMesh.Vertices.Count; i++)
                dualVertices.Add(new List<LineCurve>());

            foreach (Edge edge in iSpringMesh.Edges)
            {
                if (edge.SecondTriangleIndex >= 0)
                {
                    int firstVertexIndex = edge.FirstVertexIndex;
                    int secondVertexIndex = edge.SecondVertexIndex;
                    int firstTriangleIndex = edge.FirstTriangleIndex;
                    int secondTriangleIndex = edge.SecondTriangleIndex;

                    //Point3d firstTPI = iSpringMesh.ComputeDoubleLayerTPI(firstTriangleIndex, thickness);
                    //Point3d secondTPI = iSpringMesh.ComputeDoubleLayerTPI(secondTriangleIndex, thickness);

                    Point3d firstTPI = iSpringMesh.ComputeTPI(firstTriangleIndex);
                    Point3d secondTPI = iSpringMesh.ComputeTPI(secondTriangleIndex);

                    Point3d firstCenterPt = iSpringMesh.ComputeTriangleCentroid(firstTriangleIndex);
                    Point3d secondCenterPt = iSpringMesh.ComputeTriangleCentroid(secondTriangleIndex);

                    if (ICD.Utils.Distance(firstTPI, firstCenterPt) < iMaxThickness &&
                        ICD.Utils.Distance(secondTPI, secondCenterPt) < iMaxThickness)
                    {
                        dualVertices[firstVertexIndex].Add(new LineCurve(firstTPI, secondTPI));
                        dualVertices[secondVertexIndex].Add(new LineCurve(firstTPI, secondTPI));
                    }
                }
                else
                {

                }
            }


            for (int i = 0; i < iSpringMesh.Vertices.Count; i++)
            {
                if (dualVertices[i].Count == iSpringMesh.Vertices[i].NeighborVertexIndices.Count)
                {
                    Curve dualcurve = Curve.JoinCurves(dualVertices[i])[0];
                    oDebugCurves1.Add(dualcurve);
                }
            }


        }


        private List<double> checkCoplanarity(SpringMesh sMesh, double thickness, List<Line> linesUp, List<Line> linesDown)
        {

            List<double> coplanarity = new List<double>();

            foreach (Edge edge in sMesh.Edges)
            {
                if (edge.SecondTriangleIndex >= 0)
                {
                    int triLeft = edge.FirstTriangleIndex;
                    int triRight = edge.SecondTriangleIndex;

                    Point3d centreLeftO = sMesh.ComputeTriangleCentroid(triLeft);
                    Point3d centreRightO = sMesh.ComputeTriangleCentroid(triRight);

                    Point3d centreLeft = sMesh.ComputeTPI(triLeft);
                    Point3d centreRight = sMesh.ComputeTPI(triRight);

                    if (ICD.Utils.Distance(centreLeftO, centreLeft) < iMaxThickness && ICD.Utils.Distance(centreRightO, centreRight) < iMaxThickness)
                    {

                        Vector3d normalLeft = sMesh.ComputeTriangleNormal(triLeft);
                        Vector3d normalRight = sMesh.ComputeTriangleNormal(triRight);

                        Point3d centreLeftUp = centreLeft + normalLeft * thickness;
                        Point3d centreRightUp = centreRight + normalRight * thickness;

                        Point3d centreLeftDown = centreLeft - normalLeft * thickness;
                        Point3d centreRightDown = centreRight - normalRight * thickness;

                        linesUp.Add(new Line(centreLeftUp, centreRightUp));
                        linesDown.Add(new Line(centreLeftDown, centreRightDown));

                        Line planarCheckLine01 = new Line(centreLeftUp, centreRightDown);
                        Line planarCheckLine02 = new Line(centreRightUp, centreLeftDown);

                        double dist = planarCheckLine01.MinimumDistanceTo(planarCheckLine02);
                        coplanarity.Add(dist);
                    }
                }
            }
            return coplanarity;
        }


        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                //You can add image files to your project resources and access them like this:
                // return Resources.IconForThisComponent;
                return null;
            }
        }

        public override Guid ComponentGuid
        {
            get { return new Guid("{70d1a90f-1a4c-4680-9def-e1781912f0f1}"); }
        }
    }
}