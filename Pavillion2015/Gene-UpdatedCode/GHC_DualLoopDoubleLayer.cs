using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;

namespace Pavillion2015
{
    public class GHC_DoubleLayer : GH_Component
    {
        public SpringMesh sMesh;
        public double factor;
        public double thickness; 

        public List<double> vertexFactors;
        public double LayerThickness;

        /// <summary>
        /// Initializes a new instance of the GHC_DoubleLayer class.
        /// </summary>
        public GHC_DoubleLayer()
            : base("Dual-Loop Double Layer", "Dual-Loop Double Layer",
                "Dual-Loop Double Layer",
                "Pavillion 2015", "Double Layer")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Spring Mesh", "Spring Mesh", "Spring Mesh", GH_ParamAccess.item);
            pManager.AddNumberParameter("Factor", "Factor", "Factor", GH_ParamAccess.item, 0.9);
            pManager.AddNumberParameter("Thickness", "Thickness", "Thickness", GH_ParamAccess.item, 30.0);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddLineParameter("lines", "lines", "lines", GH_ParamAccess.item);
            //pManager.AddGenericParameter("surfaces", "surfaces", "surfaces", GH_ParamAccess.list);
            //pManager.AddBrepParameter("surfaces", "surfaces", "surfaces", GH_ParamAccess.list);
            pManager.AddCurveParameter("curveRight", "curveRight", "curveRight", GH_ParamAccess.list);
            pManager.AddCurveParameter("curveLeft", "curveLeft", "curveLeft", GH_ParamAccess.list);

            pManager.AddNumberParameter("Coplanarity", "Coplanarity", "Coplanarity", GH_ParamAccess.list);
            pManager.AddLineParameter("linesUp", "linesUp", "linesUp", GH_ParamAccess.list);
            pManager.AddLineParameter("linesDown", "linesDown", "linesDown", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            sMesh = null;
            DA.GetData<SpringMesh>(0, ref sMesh);

            factor = 0;
            DA.GetData<double>(1, ref factor);

            thickness = 0;
            DA.GetData<double>(2, ref thickness);

            List<Line> lines = SubdivideThreeLines(sMesh, factor);

            DA.SetDataList(0, lines);


            List<Curve> curveRight = new List<Curve>();
            List<Curve> curveLeft = new List<Curve>();

            vertexFactors = ComputeVertexFactor(sMesh);
            CalculateSurface(sMesh, thickness, curveRight, curveLeft);

            DA.SetDataList(1, curveRight);
            DA.SetDataList(2, curveLeft);

            List<Line> linesUp = new List<Line>();
            List<Line> linesDown = new List<Line>();

            List<double> planarity = checkCoplanarity(sMesh, thickness, linesUp, linesDown);

            DA.SetDataList(3, planarity);
            DA.SetDataList(4, linesUp);
            DA.SetDataList(5, linesDown);
        }

        /// <summary>
        /// Provides an Icon for the component.
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                //You can add image files to your project resources and access them like this:
                // return Resources.IconForThisComponent;
                return null;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("{7a0fed01-fc5d-4f6e-89b0-ec07471b1660}"); }
        }

        public List<Line> SubdivideThreeLines(SpringMesh sMesh, double factor)
        {
            List<Line> allLines = new List<Line>();

            foreach (Edge edge in sMesh.Edges)
            {
                if (edge.SecondTriangleIndex >= 0)
                {
                    int triangle1 = edge.FirstTriangleIndex;
                    int triangle2 = edge.SecondTriangleIndex;

                    //Point3d centreTri1 = sMesh.ComputeTriangleCenter(triangle1);
                    //Point3d centreTri2 = sMesh.ComputeTriangleCenter(triangle2);
                    Point3d centreTri1 = sMesh.ComputeCircumscribedCircleCenter(triangle1);
                    Point3d centreTri2 = sMesh.ComputeCircumscribedCircleCenter(triangle2);

                    List<Point3d> pointOnTri1 = sMesh.ComputeTriangleThreePts(triangle1, factor);
                    List<Point3d> pointOnTri2 = sMesh.ComputeTriangleThreePts(triangle2, factor);

                    foreach (Point3d IntPtr in pointOnTri1)
                        allLines.Add(new Line(centreTri1, IntPtr));
                    foreach (Point3d IntPtr in pointOnTri2)
                        allLines.Add(new Line(centreTri2, IntPtr));

                }
            }

            return allLines;
        }

        public List<double> checkCoplanarity(SpringMesh sMesh, double thickness, List<Line> linesUp, List<Line> linesDown)
        {

            List<double> coplanarity = new List<double>();

            foreach (Edge edge in sMesh.Edges)
            {
                if (edge.SecondTriangleIndex >= 0)
                {
                    int triLeft = edge.FirstTriangleIndex;
                    int triRight = edge.SecondTriangleIndex;

                    //Point3d centreLeft = sMesh.ComputeTriangleCenter(triLeft);
                    //Point3d centreRight = sMesh.ComputeTriangleCenter(triRight);

                    Point3d centreLeft = sMesh.ComputeCircumscribedCircleCenter(triLeft);
                    Point3d centreRight = sMesh.ComputeCircumscribedCircleCenter(triRight);

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
            return coplanarity;
        }

        public List<double> ComputeVertexFactor(SpringMesh sMesh)
        {
            List<double> vertexFactors = new List<double>();

            for (int i = 0; i < sMesh.Vertices.Count; i++)
            {
                vertexFactors.Add(0.5);
            }

            return vertexFactors;
        }

        public bool CalculateSurface(SpringMesh sMesh, double thickness, List<Curve> curvesRight, List<Curve> curvesLeft)
        {

            foreach (Edge edge in sMesh.Edges)
            {
                if (edge.SecondTriangleIndex >= 0)
                {

                    // first vertex
                    // first triangle 
                    #region curve01

                    int firstPtId = edge.FirstVertexIndex;
                    int firstTriangleIdRight = edge.FirstTriangleIndex;

                    Vertex firstPtVertex = sMesh.Vertices[firstPtId];
                    double firstPtVertexFactor = vertexFactors[firstPtId];
                    //Point3d firstTriCentreRight = sMesh.ComputeTriangleCenter(firstTriangleIdRight);
                    Point3d firstTriCentreRight = sMesh.ComputeCircumscribedCircleCenter(firstTriangleIdRight);

                    Vector3d firstNormalRight = thickness * sMesh.ComputeTriangleNormal(firstTriangleIdRight);

                    Point3d firstTriCentreUpRight = firstTriCentreRight + firstNormalRight;      // the first control point 
                    Point3d firstTriCentreDownRight = firstTriCentreRight - firstNormalRight;      // the seventh control point 

                    Point3d firstCPt02Right = (firstPtVertexFactor - 0.2) * firstPtVertex.Position + (1 - firstPtVertexFactor + 0.2) * firstTriCentreRight;

                    Point3d firstCPt02UpRight = firstCPt02Right + firstNormalRight;                // the second control point
                    Point3d firstCPt02DownRight = firstCPt02Right - firstNormalRight;              // the sixth control point

                    Point3d firstCPt03Right = firstPtVertexFactor * firstPtVertex.Position + (1 - firstPtVertexFactor) * firstTriCentreRight;

                    Point3d firstCPt03UpRight = firstCPt03Right + firstNormalRight;                  // the third control point
                    Point3d firstCPt03DownRight = firstCPt03Right - firstNormalRight;                // the fifth control point

                    Point3d firstUpPtRight = firstPtVertex.Position;                       // the fourth control point

                    List<Point3d> firstControlPointRight = new List<Point3d>() { 
                        firstTriCentreUpRight, firstCPt02UpRight, firstCPt03UpRight, firstUpPtRight, 
                        firstCPt03DownRight, firstCPt02DownRight, firstTriCentreDownRight
                    };
                    #endregion
                    NurbsCurve firstNurbCurveRight = NurbsCurve.Create(false, 5, firstControlPointRight);

                    // second triangle 
                    #region curve02

                    int firstTriangleIdLeft = edge.SecondTriangleIndex;

                    //Point3d firstTriCentreLeft = sMesh.ComputeTriangleCenter(firstTriangleIdLeft);
                    Point3d firstTriCentreLeft = sMesh.ComputeCircumscribedCircleCenter(firstTriangleIdLeft);

                    Vector3d firstNormalLeft = thickness * sMesh.ComputeTriangleNormal(firstTriangleIdLeft);

                    Point3d firstTriCentreUpLeft = firstTriCentreLeft + firstNormalLeft;      // the first control point 
                    Point3d firstTriCentreDownLeft = firstTriCentreLeft - firstNormalLeft;      // the seventh control point 

                    Point3d firstCPt02Left = (firstPtVertexFactor - 0.2) * firstPtVertex.Position + (1 - firstPtVertexFactor + 0.2) * firstTriCentreLeft;

                    Point3d firstCPt02UpLeft = firstCPt02Left + firstNormalLeft;                // the second control point
                    Point3d firstCPt02DownLeft = firstCPt02Left - firstNormalLeft;              // the sixth control point

                    Point3d firstCPt03Left = firstPtVertexFactor * firstPtVertex.Position + (1 - firstPtVertexFactor) * firstTriCentreLeft;

                    Point3d firstCPt03UpLeft = firstCPt03Left + firstNormalLeft;                  // the third control point
                    Point3d firstCPt03DownLeft = firstCPt03Left - firstNormalLeft;                // the fifth control point

                    Point3d firstUpPtLeft = firstPtVertex.Position;                       // the fourth control point

                    List<Point3d> firstControlPointLeft = new List<Point3d>() { 
                        firstTriCentreUpLeft, firstCPt02UpLeft, firstCPt03UpLeft, firstUpPtLeft, 
                        firstCPt03DownLeft, firstCPt02DownLeft, firstTriCentreDownLeft
                    };
                    #endregion
                    NurbsCurve firstNurbCurveLeft = NurbsCurve.Create(false, 5, firstControlPointLeft);

                    curvesRight.Add(firstNurbCurveRight);
                    curvesLeft.Add(firstNurbCurveLeft);

                    /*
                    Brep[] brep1 = Brep.CreateFromLoft(
                        new List<Curve>() { firstNurbCurveRight, firstNurbCurveLeft },
                        Point3d.Unset, Point3d.Unset, LoftType.Developable, false);
                    
                    if (brep1.Length > 0) breps.Add(brep1[0]); */

                    // second vertex
                    // first triangle 
                    #region curve03

                    int secondPtId = edge.SecondVertexIndex;
                    int secondTriangleIdRight = edge.FirstTriangleIndex;

                    Vertex secondPtVertex = sMesh.Vertices[secondPtId];
                    double secondPtVertexFactor = vertexFactors[secondPtId];
                    //Point3d secondTriCentreRight = sMesh.ComputeTriangleCenter(secondTriangleIdRight);
                    Point3d secondTriCentreRight = sMesh.ComputeCircumscribedCircleCenter(secondTriangleIdRight);

                    Vector3d secondNormalRight = thickness * sMesh.ComputeTriangleNormal(secondTriangleIdRight);

                    Point3d secondTriCentreUpRight = secondTriCentreRight + secondNormalRight;      // the second control point 
                    Point3d secondTriCentreDownRight = secondTriCentreRight - secondNormalRight;      // the seventh control point 

                    Point3d secondCPt02Right = (secondPtVertexFactor - 0.2) * secondPtVertex.Position + (1 - secondPtVertexFactor + 0.2) * secondTriCentreRight;

                    Point3d secondCPt02UpRight = secondCPt02Right + secondNormalRight;                // the second control point
                    Point3d secondCPt02DownRight = secondCPt02Right - secondNormalRight;              // the sixth control point

                    Point3d secondCPt03Right = secondPtVertexFactor * secondPtVertex.Position + (1 - secondPtVertexFactor) * secondTriCentreRight;

                    Point3d secondCPt03UpRight = secondCPt03Right + secondNormalRight;                  // the third control point
                    Point3d secondCPt03DownRight = secondCPt03Right - secondNormalRight;                // the fifth control point

                    Point3d secondUpPtRight = secondPtVertex.Position;                       // the fourth control point

                    List<Point3d> secondControlPointRight = new List<Point3d>() { 
                        secondTriCentreUpRight, secondCPt02UpRight, secondCPt03UpRight, secondUpPtRight, 
                        secondCPt03DownRight, secondCPt02DownRight, secondTriCentreDownRight
                    };

                    #endregion
                    NurbsCurve secondNurbCurveRight = NurbsCurve.Create(false, 5, secondControlPointRight);

                    // second triangle 
                    #region curve04

                    int secondTriangleIdLeft = edge.SecondTriangleIndex;

                    //Point3d secondTriCentreLeft = sMesh.ComputeTriangleCenter(secondTriangleIdLeft);
                    Point3d secondTriCentreLeft = sMesh.ComputeCircumscribedCircleCenter(secondTriangleIdLeft);

                    Vector3d secondNormalLeft = thickness * sMesh.ComputeTriangleNormal(secondTriangleIdLeft);

                    Point3d secondTriCentreUpLeft = secondTriCentreLeft + secondNormalLeft;      // the second control point 
                    Point3d secondTriCentreDownLeft = secondTriCentreLeft - secondNormalLeft;      // the seventh control point 

                    Point3d secondCPt02Left = (secondPtVertexFactor - 0.2) * secondPtVertex.Position + (1 - secondPtVertexFactor + 0.2) * secondTriCentreLeft;

                    Point3d secondCPt02UpLeft = secondCPt02Left + secondNormalLeft;                // the second control point
                    Point3d secondCPt02DownLeft = secondCPt02Left - secondNormalLeft;              // the sixth control point

                    Point3d secondCPt03Left = secondPtVertexFactor * secondPtVertex.Position + (1 - secondPtVertexFactor) * secondTriCentreLeft;

                    Point3d secondCPt03UpLeft = secondCPt03Left + secondNormalLeft;                  // the third control point
                    Point3d secondCPt03DownLeft = secondCPt03Left - secondNormalLeft;                // the fifth control point

                    Point3d secondUpPtLeft = secondPtVertex.Position;                       // the fourth control point

                    List<Point3d> secondControlPointLeft = new List<Point3d>() { 
                        secondTriCentreUpLeft, secondCPt02UpLeft, secondCPt03UpLeft, secondUpPtLeft, 
                        secondCPt03DownLeft, secondCPt02DownLeft, secondTriCentreDownLeft
                    };
                    #endregion
                    NurbsCurve secondNurbCurveLeft = NurbsCurve.Create(false, 5, secondControlPointLeft);

                    curvesRight.Add(secondNurbCurveRight);
                    curvesLeft.Add(secondNurbCurveLeft);

                    /*
                    Brep[] brep2 = Brep.CreateFromLoft(
                        new List<Curve>() { firstNurbCurveRight, firstNurbCurveLeft },
                        Point3d.Unset, Point3d.Unset, LoftType.Developable, false);

                    if (brep2.Length > 0) breps.Add(brep2[0]);  */

                }
            }
            return true;
        }


    }
}