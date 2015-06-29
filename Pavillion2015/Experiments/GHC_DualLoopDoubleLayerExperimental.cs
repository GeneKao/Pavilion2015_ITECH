using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;
using Rhino.Geometry.Intersect;

namespace Pavillion2015
{
    public class GHC_DualLoopDoubleLayerExperimental : GH_Component
    {

        SpringMesh iSpringMesh = null;
        bool iReset = false;
        bool iPlay = false;
        int iSubiterationCount = 0;
        int iPopulationSize = 0;
        double iMutationRate = double.NaN;
        double iMutationAmount = double.NaN;

        string oInfo = string.Empty;
        List<Point3d> oDebugPoints1 = null;
        List<Point3d> oDebugPoints2 = null;
        List<Curve> oDebugCurves1 = null;
        List<Curve> oDebugCurves2 = null;
        List<Vector3d> oDebugVectors1 = null;
        List<Vector3d> oDebugVectors2 = null;
        List<double> oDebugNumbers1 = null;
        List<double> oDebugNumbers2 = null;


        public GHC_DualLoopDoubleLayerExperimental()
            : base("Dual-Loop Double Layer (Experimetial)", "Dual-Loop Double Layer (Experimental)",
                "Double Layer (Experimetal)",
                "Pavillion 2015", "Double Layer")
        {
        }


        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Spring Mesh", "Spring Mesh", "Spring Mesh", GH_ParamAccess.item);
            pManager.AddBooleanParameter("Reset", "Reset", "Reset", GH_ParamAccess.item, false);
            pManager.AddBooleanParameter("Play", "Play", "Play", GH_ParamAccess.item, false);
            pManager.AddIntegerParameter("Subiteration Count", "Subiteration Count", "Subiteration Count", GH_ParamAccess.item, 1);
            pManager.AddIntegerParameter("Population Size", "Population Size", " Population Size", GH_ParamAccess.item, 50);
            pManager.AddNumberParameter("Mutation Rate", "Mutation Rate", "Mutation Rate", GH_ParamAccess.item, 0.1);
            pManager.AddNumberParameter("Mutation Amount", "Mutation Amount", "Mutation Amount", GH_ParamAccess.item, 90.0);
        }


        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddTextParameter("Info", "00 - Info", "Information", GH_ParamAccess.item);
            pManager.AddGenericParameter("Debug 1", "01 - Debug 1", "This output is reserved for debugging", GH_ParamAccess.list);
            pManager.AddGenericParameter("Debug 2", "02 - Debug 2", "This output is reserved for debugging", GH_ParamAccess.list);
            pManager.AddGenericParameter("Debug 3", "03 - Debug 3", "This output is reserved for debugging", GH_ParamAccess.list);
            pManager.AddGenericParameter("Debug 4", "04 - Debug 4", "This output is reserved for debugging", GH_ParamAccess.list);
            pManager.AddGenericParameter("Debug 5", "05 - Debug 5", "This output is reserved for debugging", GH_ParamAccess.list);
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
        }


        protected override void SolveInstance(IGH_DataAccess DA)
        {
            DA.GetData<SpringMesh>("Spring Mesh", ref iSpringMesh);
            DA.GetData<bool>("Reset", ref iReset);
            DA.GetData<bool>("Play", ref iPlay);
            DA.GetData<int>("Subiteration Count", ref iSubiterationCount);
            DA.GetData<int>("Population Size", ref iPopulationSize);
            DA.GetData<double>("Mutation Rate", ref iMutationRate);
            DA.GetData<double>("Mutation Amount", ref iMutationAmount);


            // -------------------------------------------------------------------------------------------------------

            List<Vector3d> triangleNormals = new List<Vector3d>();
            for (int i = 0; i < iSpringMesh.Triangles.Count; i++)
                triangleNormals.Add(iSpringMesh.ComputeTriangleNormal(i));

            List<List<Plane>> planeTriplets = new List<List<Plane>>();
            for (int i = 0; i < iSpringMesh.Triangles.Count; i++)
                planeTriplets.Add(new List<Plane>());

            List<Plane> edgePlanes = new List<Plane>();;

            foreach (Edge edge in iSpringMesh.Edges)
            {
                if (edge.SecondTriangleIndex >= 0)
                {
                    Vector3d planeNormal = triangleNormals[edge.FirstTriangleIndex] + triangleNormals[edge.SecondTriangleIndex];
                    Point3d planeOrigin = 0.5 * (iSpringMesh.Vertices[edge.FirstVertexIndex].Position + iSpringMesh.Vertices[edge.SecondVertexIndex].Position);
                    
                    
                    Vector3d localX = iSpringMesh.Vertices[edge.SecondVertexIndex].Position - iSpringMesh.Vertices[edge.FirstVertexIndex].Position;

                    localX.Unitize();
                    planeNormal.Unitize();

                    Vector3d localY = Vector3d.CrossProduct(localX, planeNormal);

                    Plane plane = new Plane(planeOrigin, localX, localY);

                    planeTriplets[edge.FirstTriangleIndex].Add(plane);
                    planeTriplets[edge.SecondTriangleIndex].Add(plane);
                    edgePlanes.Add(plane);
                }
                else
                {
                    /////////////////////
                }
            }

            List<Plane> nastyPlanes = new List<Plane>();
            List<Point3d> planeTripletIntersections = new List<Point3d>();

            foreach (List<Plane> planeTriplet in planeTriplets)
            {
                if (planeTriplet.Count != 3) continue;
                
                Point3d intersectionPoint;
                if (Intersection.PlanePlanePlane(planeTriplet[0], planeTriplet[1], planeTriplet[2], out intersectionPoint))
                    planeTripletIntersections.Add(intersectionPoint);
                else
                {
                    nastyPlanes.Add(planeTriplet[0]);
                    nastyPlanes.Add(planeTriplet[1]);
                    nastyPlanes.Add(planeTriplet[2]);
                }
            }

            

            DA.SetDataList(1, edgePlanes);
            DA.SetDataList(2, planeTripletIntersections);
            DA.SetDataList(3, nastyPlanes);
            

            //for (int i = 0; i < planeTriplets.Count; i++)
            //{
            //    if (planeTriplets[i].Count == 3) //////// remove later
            //    {
            //        Point3d intersectionPoint;

            //        if (Intersection.PlanePlanePlane(planeTriplets[i][0], planeTriplets[i][1], planeTriplets[i][2], out intersectionPoint))
            //            trianglePoints.Add(intersectionPoint);
            //        else
            //            trianglePoints.Add(iSpringMesh.ComputeTriangleCentroid(i));
            //    }
            //    else
            //    {
            //        Line intersectionLine;
            //        if (Intersection.PlanePlane(planeTriplets[i][0], planeTriplets[i][1], out intersectionLine))
            //        {
            //            trianglePoints.Add(
            //                Utils.ClosestPointOnLine(
            //                iSpringMesh.ComputeTriangleCentroid(i),
            //                intersectionLine.PointAt(0.0),
            //                intersectionLine.PointAt(1.0)));
            //        }

            //        else
            //            trianglePoints.Add(iSpringMesh.ComputeTriangleCentroid(i));
            //    }

            //}

            //Mesh oRhinoMesh = new Mesh();

            //List<PolylineCurve> quadPolylines = new List<PolylineCurve>();

            //foreach (Edge edge in iSpringMesh.Edges)
            //{
            //    if (edge.SecondTriangleIndex >= 0)
            //    {
            //        Point3d A = iSpringMesh.Vertices[edge.FirstVertexIndex].Position;
            //        Point3d C = iSpringMesh.Vertices[edge.SecondVertexIndex].Position;

            //        Point3d B = trianglePoints[edge.FirstTriangleIndex];
            //        Point3d D = trianglePoints[edge.SecondTriangleIndex];

            //        oRhinoMesh.Vertices.Add(A);
            //        oRhinoMesh.Vertices.Add(B);
            //        oRhinoMesh.Vertices.Add(C);
            //        oRhinoMesh.Vertices.Add(D);

            //        oRhinoMesh.Faces.AddFace(oRhinoMesh.Vertices.Count - 4, oRhinoMesh.Vertices.Count - 3, oRhinoMesh.Vertices.Count - 2, oRhinoMesh.Vertices.Count - 1);

            //        quadPolylines.Add(new PolylineCurve(new List<Point3d>() { A, B, C, D , A}));
            //    }
            //}
            
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
            get { return new Guid("{7a0fed01-fc5d-4f6e-89b0-ec07471b1666}"); }
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
                    //Point3d firstTriCentreRight = sMesh.ComputeTriangleCenter(firstTriangleIdRight);
                    Point3d firstTriCentreRight = sMesh.ComputeCircumscribedCircleCenter(firstTriangleIdRight);

                    Vector3d firstNormalRight = thickness * sMesh.ComputeTriangleNormal(firstTriangleIdRight);

                    Point3d firstTriCentreUpRight = firstTriCentreRight + firstNormalRight;      // the first control point 
                    Point3d firstTriCentreDownRight = firstTriCentreRight - firstNormalRight;      // the seventh control point 

                    Point3d firstCPt02Right = (firstPtVertex.Factor - 0.2) * firstPtVertex.Position + (1 - firstPtVertex.Factor + 0.2) * firstTriCentreRight;

                    Point3d firstCPt02UpRight = firstCPt02Right + firstNormalRight;                // the second control point
                    Point3d firstCPt02DownRight = firstCPt02Right - firstNormalRight;              // the sixth control point

                    Point3d firstCPt03Right = firstPtVertex.Factor * firstPtVertex.Position + (1 - firstPtVertex.Factor) * firstTriCentreRight;

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

                    Point3d firstCPt02Left = (firstPtVertex.Factor - 0.2) * firstPtVertex.Position + (1 - firstPtVertex.Factor + 0.2) * firstTriCentreLeft;

                    Point3d firstCPt02UpLeft = firstCPt02Left + firstNormalLeft;                // the second control point
                    Point3d firstCPt02DownLeft = firstCPt02Left - firstNormalLeft;              // the sixth control point

                    Point3d firstCPt03Left = firstPtVertex.Factor * firstPtVertex.Position + (1 - firstPtVertex.Factor) * firstTriCentreLeft;

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
                    //Point3d secondTriCentreRight = sMesh.ComputeTriangleCenter(secondTriangleIdRight);
                    Point3d secondTriCentreRight = sMesh.ComputeCircumscribedCircleCenter(secondTriangleIdRight);

                    Vector3d secondNormalRight = thickness * sMesh.ComputeTriangleNormal(secondTriangleIdRight);

                    Point3d secondTriCentreUpRight = secondTriCentreRight + secondNormalRight;      // the second control point 
                    Point3d secondTriCentreDownRight = secondTriCentreRight - secondNormalRight;      // the seventh control point 

                    Point3d secondCPt02Right = (secondPtVertex.Factor - 0.2) * secondPtVertex.Position + (1 - secondPtVertex.Factor + 0.2) * secondTriCentreRight;

                    Point3d secondCPt02UpRight = secondCPt02Right + secondNormalRight;                // the second control point
                    Point3d secondCPt02DownRight = secondCPt02Right - secondNormalRight;              // the sixth control point

                    Point3d secondCPt03Right = secondPtVertex.Factor * secondPtVertex.Position + (1 - secondPtVertex.Factor) * secondTriCentreRight;

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

                    Point3d secondCPt02Left = (secondPtVertex.Factor - 0.2) * secondPtVertex.Position + (1 - secondPtVertex.Factor + 0.2) * secondTriCentreLeft;

                    Point3d secondCPt02UpLeft = secondCPt02Left + secondNormalLeft;                // the second control point
                    Point3d secondCPt02DownLeft = secondCPt02Left - secondNormalLeft;              // the sixth control point

                    Point3d secondCPt03Left = secondPtVertex.Factor * secondPtVertex.Position + (1 - secondPtVertex.Factor) * secondTriCentreLeft;

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