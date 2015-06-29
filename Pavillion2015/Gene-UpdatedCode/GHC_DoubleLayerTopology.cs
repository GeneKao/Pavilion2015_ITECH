using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;

namespace Pavillion2015
{
    public class GHC_DoubleLayerTopology : GH_Component
    {
        SpringMesh iSpringMesh = null;
        double iMinThickness = double.NaN;
        double iMaxThickness = double.NaN;
        double iTangentScale = double.NaN;
        bool iPolySrf = true;
        double iPlatesOffset = double.NaN;
        double iPlatesThreads = double.NaN;
        double iPlatesThreadsP = double.NaN;

        List<int> triangleStatus = null;
        List<bool> vertexStatus = null;

        List<Point3d> topCps = null;
        List<Point3d> bottomCps = null;

        string oInfo = string.Empty;
        List<Point3d> oDebugPoints1 = null;
        List<Point3d> oDebugPoints2 = null;
        List<Curve> oDebugCurves1 = null;
        List<Curve> oDebugCurves2 = null;
        List<Curve> oDebugCurves3 = null;
        List<Curve> oDebugCurves4 = null;
        List<Vector3d> oDebugVectors1 = null;
        List<Vector3d> oDebugVectors2 = null;
        List<double> oDebugNumbers1 = null;
        List<double> oDebugNumbers2 = null;
        List<Brep> oDebugBreps1 = null;
        List<Brep> oDebugBreps2 = null;
        List<Brep> oDebugBreps3 = null;
        List<Brep> oDebugBreps4 = null;
        List<Line> oDebugLines1 = null;
        List<Line> oDebugLines2 = null;

        double documentTolerance = DocumentTolerance();


        public GHC_DoubleLayerTopology()
            : base("Double Layer Topology (Experimetial)", "Double Layer Topology (Experimetial)",
                "Double Layer Topology - combination plate and loops",
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
            pManager.AddNumberParameter("Plate Offset", "Plate Offset", "Plate Offset", GH_ParamAccess.item, 0.5);
            pManager.AddNumberParameter("Threads", "Threads", "Threads", GH_ParamAccess.item, 1.0);
            pManager.AddNumberParameter("ThreadsP", "ThreadsP", "ThreadsP", GH_ParamAccess.item, 1.0);
        }


        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddTextParameter("Info", "00 - Info", "Information", GH_ParamAccess.item);
            pManager.AddGenericParameter("Debug 1", "01 - Debug 1", "This output is reserved for debugging", GH_ParamAccess.list);
            pManager.AddGenericParameter("Debug 2", "02 - Debug 2", "This output is reserved for debugging", GH_ParamAccess.list);
            pManager.AddGenericParameter("Debug 3", "03 - Debug 3", "This output is reserved for debugging", GH_ParamAccess.list);
            pManager.AddGenericParameter("Debug 4", "04 - Debug 4", "This output is reserved for debugging", GH_ParamAccess.list);
            pManager.AddGenericParameter("Debug 5", "05 - Debug 5", "This output is reserved for debugging", GH_ParamAccess.list);
            pManager.AddGenericParameter("Debug 6", "06 - Debug 6", "This output is reserved for debugging", GH_ParamAccess.list);
            pManager.AddGenericParameter("Components", "Components", "Components", GH_ParamAccess.list);
        }


        protected override void BeforeSolveInstance()
        {
            oInfo = string.Empty;
            oDebugPoints1 = new List<Point3d>();
            oDebugPoints2 = new List<Point3d>();
            oDebugCurves1 = new List<Curve>();
            oDebugCurves2 = new List<Curve>();
            oDebugCurves3 = new List<Curve>();
            oDebugCurves4 = new List<Curve>();
            oDebugVectors1 = new List<Vector3d>();
            oDebugVectors2 = new List<Vector3d>();
            oDebugNumbers1 = new List<double>();
            oDebugNumbers2 = new List<double>();
            oDebugBreps1 = new List<Brep>();
            oDebugBreps2 = new List<Brep>();
            oDebugBreps3 = new List<Brep>();
            oDebugBreps4 = new List<Brep>();
            oDebugLines1 = new List<Line>();
            oDebugLines2 = new List<Line>();

            triangleStatus = new List<int>(); // value 3 is 3 plates, 2 is 2 plates, 1 is 1 plates, 0 is triLoop
            vertexStatus = new List<bool>(); // true is plate occupied, false should put triLoop

        }


        protected override void SolveInstance(IGH_DataAccess DA)
        {

            DA.GetData<SpringMesh>("Spring Mesh", ref iSpringMesh);
            DA.GetData<double>("Min. Thickness", ref iMinThickness);
            DA.GetData<double>("Max. Thickness", ref iMaxThickness);
            DA.GetData<double>("Tangent Scale", ref iTangentScale);
            DA.GetData<bool>("Allow PolySrf", ref iPolySrf);
            DA.GetData<double>("Plate Offset", ref iPlatesOffset);
            DA.GetData<double>("Threads", ref iPlatesThreads);
            DA.GetData<double>("ThreadsP", ref iPlatesThreadsP);

            //---------------------------------------------------------------------------------------------------------------
            // set each vertex not occupied by any plates
            foreach (Vertex vertex in iSpringMesh.Vertices)
                vertexStatus.Add(false);



            plate();

            //---------------------------------------------------------------------------------------------------------------

            DA.SetDataList(1, oDebugCurves1);
            DA.SetDataList(2, oDebugCurves2);

            //---------------------------------------------------------------------------------------------------------------

            triLoop();

            //---------------------------------------------------------------------------------------------------------------

            DA.SetDataList(3, oDebugBreps1);
            DA.SetDataList(4, oDebugBreps2);

            //DA.SetDataList(3, oDebugPoints1);
            //DA.SetDataList(4, oDebugPoints2);

            DA.SetDataList(5, oDebugBreps3);
            DA.SetDataList(6, oDebugBreps4);

        }

        private void plate()
        {
            //vertex thickness value
            List<double> thickness = new List<double>();
            foreach (Vertex vertex in iSpringMesh.Vertices)
                thickness.Add(iMinThickness);


            //List<List<LineCurve>> dualVerticesUp = new List<List<LineCurve>>();
            //List<List<LineCurve>> dualVerticesDown = new List<List<LineCurve>>();

            //List<Point3d> dualPtPlate = new List<Point3d>();
            /// this part half dual loop connect to plates
            List<List<Point3d>> curveControlPtsRight = new List<List<Point3d>>();
            List<List<Point3d>> curveControlPtsLeft = new List<List<Point3d>>();
            /// this part for transition half dual loop without connecting plates
            List<List<Point3d>> curveControlPtsPlate01 = new List<List<Point3d>>();
            List<List<Point3d>> curveControlPtsPlate02 = new List<List<Point3d>>();
            List<List<Point3d>> curveControlPtsTriLoop01 = new List<List<Point3d>>();
            List<List<Point3d>> curveControlPtsTriLoop02 = new List<List<Point3d>>();

            DualTPI(thickness, oDebugCurves1, oDebugPoints1, curveControlPtsRight, curveControlPtsLeft,
                    curveControlPtsPlate01, curveControlPtsPlate02, curveControlPtsTriLoop01, curveControlPtsTriLoop02, false);
            DualTPI(thickness, oDebugCurves2, oDebugPoints2, curveControlPtsRight, curveControlPtsLeft,
                    curveControlPtsPlate01, curveControlPtsPlate02, curveControlPtsTriLoop01, curveControlPtsTriLoop02, true);

            for (int i = 0; i < iSpringMesh.Edges.Count; i++)
            {
                if (iSpringMesh.Edges[i].SecondTriangleIndex >= 0)
                {
                    ///// dual Loops connected to plates
                    NurbsCurve curveRight = NurbsCurve.Create(false, 3, curveControlPtsRight[i]);
                    NurbsCurve curveLeft = NurbsCurve.Create(false, 3, curveControlPtsLeft[i]);

                    Brep[] loftHalfDualLoop = Brep.CreateFromLoft(new List<NurbsCurve>() { curveRight, curveLeft },
                           Point3d.Unset, Point3d.Unset, LoftType.Normal, false);
                    if (loftHalfDualLoop.Length == 1)
                        oDebugBreps2.Add(loftHalfDualLoop[0]);

                    //oDebugCurves3.Add(curveRight);
                    //oDebugCurves4.Add(curveLeft);

                    //// dual Loops without connecting to plates
                    NurbsCurve curvePlate01 = NurbsCurve.Create(false, 3, curveControlPtsPlate01[i]);
                    NurbsCurve curveTriLoop01 = NurbsCurve.Create(false, 3, curveControlPtsTriLoop01[i]);

                    oDebugCurves3.Add(curvePlate01);
                    oDebugCurves4.Add(curveTriLoop01);
                    
                    Brep[] loftHalfDualLoop01 = Brep.CreateFromLoft(new List<NurbsCurve>() { curvePlate01, curveTriLoop01 },
                           Point3d.Unset, Point3d.Unset, LoftType.Normal, false);
                    if (loftHalfDualLoop01.Length == 1)
                        oDebugBreps3.Add(loftHalfDualLoop01[0]);

                    NurbsCurve curvePlate02 = NurbsCurve.Create(false, 3, curveControlPtsPlate02[i]);
                    NurbsCurve curveTriLoop02 = NurbsCurve.Create(false, 3, curveControlPtsTriLoop02[i]);
                    
                    Brep[] loftHalfDualLoop02 = Brep.CreateFromLoft(new List<NurbsCurve>() { curvePlate02, curveTriLoop02 },
                           Point3d.Unset, Point3d.Unset, LoftType.Normal, false);
                    if (loftHalfDualLoop02.Length == 1)
                        oDebugBreps4.Add(loftHalfDualLoop02[0]);

                }
            }

        }

        private void DualTPI(List<double> thickness, List<Curve> oDebugCurves, List<Point3d> oDebugPoints,
            List<List<Point3d>> curveControlPtsRight, List<List<Point3d>> curveControlPtsLeft,
            List<List<Point3d>> curveControlPtsPlate01, List<List<Point3d>> curveControlPtsPlate02,
            List<List<Point3d>> curveControlPtsTriLoop01, List<List<Point3d>> curveControlPtsTriLoop02, 
            bool flip)
        {

            List<List<LineCurve>> dualVertices = new List<List<LineCurve>>();
            List<Point3d> dualPlatePoint = new List<Point3d>();


            for (int i = 0; i < iSpringMesh.Vertices.Count; i++)
                dualVertices.Add(new List<LineCurve>());

            for (int i = 0; i < iSpringMesh.Triangles.Count; i++)
                dualPlatePoint.Add(new Point3d());

            for (int i = 0; i < iSpringMesh.Edges.Count; i++)
            {
                curveControlPtsRight.Add(new List<Point3d>());
                curveControlPtsLeft.Add(new List<Point3d>());
                curveControlPtsPlate01.Add(new List<Point3d>());
                curveControlPtsPlate02.Add(new List<Point3d>());
                curveControlPtsTriLoop01.Add(new List<Point3d>());
                curveControlPtsTriLoop02.Add(new List<Point3d>());
            }

            ////// this is to create planar plate
            foreach (Edge edge in iSpringMesh.Edges)
            {
                if (edge.SecondTriangleIndex >= 0)
                {
                    int firstVertexIndex = edge.FirstVertexIndex;
                    int secondVertexIndex = edge.SecondVertexIndex;
                    int firstTriangleIndex = edge.FirstTriangleIndex;
                    int secondTriangleIndex = edge.SecondTriangleIndex;

                    Point3d firstTPI = iSpringMesh.ComputeDoubleLayerTPI(firstTriangleIndex, thickness, flip, iPlatesOffset);
                    Point3d secondTPI = iSpringMesh.ComputeDoubleLayerTPI(secondTriangleIndex, thickness, flip, iPlatesOffset);

                    Point3d firstTPIO = iSpringMesh.ComputeTPI(firstTriangleIndex);
                    Point3d secondTPIO = iSpringMesh.ComputeTPI(secondTriangleIndex);
                    Point3d firstCenterPtO = iSpringMesh.ComputeTriangleCentroid(firstTriangleIndex);
                    Point3d secondCenterPtO = iSpringMesh.ComputeTriangleCentroid(secondTriangleIndex);

                    Point3d firstCenterPt;
                    Point3d secondCenterPt;
                    if (!flip)
                    {
                        firstCenterPt = iSpringMesh.ComputeTriangleCentroid(firstTriangleIndex) +
                            iSpringMesh.ComputeTriangleNormal(firstTriangleIndex);
                        secondCenterPt = iSpringMesh.ComputeTriangleCentroid(secondTriangleIndex) +
                            iSpringMesh.ComputeTriangleNormal(secondTriangleIndex);
                    }
                    else
                    {
                        firstCenterPt = iSpringMesh.ComputeTriangleCentroid(firstTriangleIndex) -
                            iSpringMesh.ComputeTriangleNormal(firstTriangleIndex);
                        secondCenterPt = iSpringMesh.ComputeTriangleCentroid(secondTriangleIndex) -
                            iSpringMesh.ComputeTriangleNormal(secondTriangleIndex);
                    }

                    // ===========================================================================
                    // Added 14/06/2015 to make some tolerance for plates

                    Circle firstCircle = iSpringMesh.ComputeDoubleLayerIncircle(firstTriangleIndex, iPlatesThreadsP, thickness, flip, iPlatesOffset);
                    Circle secondCircle = iSpringMesh.ComputeDoubleLayerIncircle(secondTriangleIndex, iPlatesThreadsP, thickness, flip, iPlatesOffset);

                    Circle firstCircleO = iSpringMesh.ComputeIncircle(firstTriangleIndex, iPlatesThreadsP); // iPlatesThreads = 0.0-1.0
                    Circle secondCircleO = iSpringMesh.ComputeIncircle(secondTriangleIndex, iPlatesThreadsP);

                    //Plane firstP = iSpringMesh.ComputeDoubleLayerPlane(firstTriangleIndex, thickness, flip, iPlatesOffset);
                    //Plane secondP = iSpringMesh.ComputeDoubleLayerPlane(secondTriangleIndex, thickness, flip, iPlatesOffset);

                    Plane firstPO = iSpringMesh.ComputePlane(firstTriangleIndex);
                    Plane secondPO = iSpringMesh.ComputePlane(secondTriangleIndex);

                    Point3d projectFirstTPIO = firstPO.ClosestPoint(firstTPIO);
                    Point3d projectSecondTPIO = secondPO.ClosestPoint(secondTPIO);

                    if (ICD.Utils.Distance(projectFirstTPIO, firstCenterPtO) < iPlatesThreads &&
                        ICD.Utils.Distance(projectSecondTPIO, secondCenterPtO) < iPlatesThreads)
                    {
                        Point3d newFirstTPI = firstCircle.ClosestPoint(firstTPI);
                        Point3d newSecondTPI = secondCircle.ClosestPoint(secondTPI);

                        dualVertices[firstVertexIndex].Add(new LineCurve(newFirstTPI, newSecondTPI));
                        dualVertices[secondVertexIndex].Add(new LineCurve(newFirstTPI, newSecondTPI));

                        if (dualPlatePoint[firstTriangleIndex] == new Point3d())
                            dualPlatePoint[firstTriangleIndex] = newFirstTPI;
                        if (dualPlatePoint[secondTriangleIndex] == new Point3d())
                            dualPlatePoint[secondTriangleIndex] = newSecondTPI;
                    }

                    // ===========================================================================
                    /*
                    if (ICD.Utils.Distance(firstTPIO, firstCenterPtO) < iPlatesThreads &&
                        ICD.Utils.Distance(secondTPIO, secondCenterPtO) < iPlatesThreads)
                    {
                        dualVertices[firstVertexIndex].Add(new LineCurve(firstTPI, secondTPI));
                        dualVertices[secondVertexIndex].Add(new LineCurve(firstTPI, secondTPI));

                        if (dualPlatePoint[firstTriangleIndex] == new Point3d())
                            dualPlatePoint[firstTriangleIndex] = firstTPI;
                        if (dualPlatePoint[secondTriangleIndex] == new Point3d())
                            dualPlatePoint[secondTriangleIndex] = secondTPI;
                    }
                    */
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

                    //if (dualcurve.IsClosed) continue;

                    oDebugCurves.Add(dualcurve);

                    // change the vertex staus to be occupied 
                    vertexStatus[i] = true;

                    //----------------------------------------------------------------------------------
                    //dualLoop transition 
                    /*
                    List<double> factor = ComputeVertexFactor();

                    List<int> neighbourID = iSpringMesh.Vertices[i].NeighborVertexIndices;
                    
                    foreach (int id in neighbourID)
                    {
                        //Point3d ptOnPlate = 
                        Point3d pt01 = iSpringMesh.Vertices[id].Position;
                    }
                    */
                    //----------------------------------------------------------------------------------
                }
            }
            // show those TCPs which in the right condition
            for (int i = 0; i < iSpringMesh.Triangles.Count; i++)
            {
                int id01 = iSpringMesh.Triangles[i].FirstVertexIndex;
                int id02 = iSpringMesh.Triangles[i].SecondVertexIndex;
                int id03 = iSpringMesh.Triangles[i].ThirdVertexIndex;
                if (vertexStatus[id01] || vertexStatus[id02] || vertexStatus[id03] == true)
                    oDebugPoints.Add(dualPlatePoint[i]);
            }

            
            List<Vector3d> vertexNormals = iSpringMesh.ComputeVertexNormals();

            for (int i = 0; i < iSpringMesh.Edges.Count; i++)
            {
                if (iSpringMesh.Edges[i].SecondTriangleIndex >= 0)
                {
                    int firstVertexIndex = iSpringMesh.Edges[i].FirstVertexIndex;
                    int secondVertexIndex = iSpringMesh.Edges[i].SecondVertexIndex;
                    int firstTriangleIndex = iSpringMesh.Edges[i].FirstTriangleIndex;
                    int secondTriangleIndex = iSpringMesh.Edges[i].SecondTriangleIndex;
                    int firstAdjVIndex = iSpringMesh.Edges[i].FirstAdjacentVertexIndex;
                    int secondAdjVIndex = iSpringMesh.Edges[i].SecondAdjacentVertexIndex;

                    Vector3d vertexNormal = vertexNormals[iSpringMesh.Edges[i].FirstVertexIndex];
                    vertexNormal.Unitize();
                    vertexNormals[iSpringMesh.Edges[i].FirstVertexIndex] = vertexNormal;

                    vertexNormal = vertexNormals[iSpringMesh.Edges[i].SecondVertexIndex];
                    vertexNormal.Unitize();
                    vertexNormals[iSpringMesh.Edges[i].SecondVertexIndex] = vertexNormal;

                    
                    if (!flip)    // upper layer
                    {
                        ///// half dual loop which connect to the plates
                        Point3d platePtRightUp;
                        Point3d platePtLeftUp;
                        if (vertexStatus[firstVertexIndex] == true && vertexStatus[secondVertexIndex] == false)
                        {
                            platePtRightUp = dualPlatePoint[firstTriangleIndex];
                            platePtLeftUp = dualPlatePoint[secondTriangleIndex];

                            curveControlPtsRight[i].Add(platePtRightUp);     //plates points
                            curveControlPtsLeft[i].Add(platePtLeftUp);       
                            
                            //// here should be insert some control point so no kink in transition area
                            Point3d vertexPosition = iSpringMesh.Vertices[secondVertexIndex].Position;

                            Point3d ctPRight = (1 - iTangentScale) * platePtRightUp +
                                iTangentScale * (vertexPosition + iMinThickness * vertexNormals[secondVertexIndex]);
                            Point3d ctLeft = (1 - iTangentScale) * platePtLeftUp +
                                iTangentScale * (vertexPosition + iMinThickness * vertexNormals[secondVertexIndex]);

                            curveControlPtsRight[i].Add(ctPRight);
                            curveControlPtsLeft[i].Add(ctLeft);

                            curveControlPtsRight[i].Add(vertexPosition);
                            curveControlPtsLeft[i].Add(vertexPosition);
                        }
                        if (vertexStatus[firstVertexIndex] == false && vertexStatus[secondVertexIndex] == true)
                        {
                            platePtRightUp = dualPlatePoint[firstTriangleIndex];
                            platePtLeftUp = dualPlatePoint[secondTriangleIndex];

                            curveControlPtsRight[i].Add(platePtRightUp);
                            curveControlPtsLeft[i].Add(platePtLeftUp);

                            //// here should be insert some control point so no kink in transition area
                            Point3d vertexPosition = iSpringMesh.Vertices[firstVertexIndex].Position;

                            Point3d ctPRight = (1 - iTangentScale) * platePtRightUp +
                                iTangentScale * (vertexPosition + iMinThickness * vertexNormals[firstVertexIndex]);
                            Point3d ctLeft = (1 - iTangentScale) * platePtLeftUp +
                                iTangentScale * (vertexPosition + iMinThickness * vertexNormals[firstVertexIndex]);

                            curveControlPtsRight[i].Add(ctPRight);
                            curveControlPtsLeft[i].Add(ctLeft);

                            curveControlPtsRight[i].Add(vertexPosition);
                            curveControlPtsLeft[i].Add(vertexPosition);
                        }

                        ///------------------------------------------------------------------------------------------------
                        ///// half dual loop without connecting to the plates
                        if (vertexStatus[firstAdjVIndex] == true && vertexStatus[secondAdjVIndex] == false && 
                            vertexStatus[firstVertexIndex] == false && vertexStatus[secondVertexIndex] == false)
                        {
                            int firstTriV01 = iSpringMesh.Triangles[firstTriangleIndex].FirstVertexIndex;
                            int firstTriV02 = iSpringMesh.Triangles[firstTriangleIndex].SecondVertexIndex;
                            int firstTriV03 = iSpringMesh.Triangles[firstTriangleIndex].ThirdVertexIndex;
                            
                            // dualLoop connect to plate
                            Point3d platePtUp;

                            //// first adjacent vertex is belong to first triangle
                            if(firstTriV01 == firstAdjVIndex || firstTriV02 == firstAdjVIndex || firstTriV03 == firstAdjVIndex)
                                platePtUp = dualPlatePoint[firstTriangleIndex];
                            else //// first adjacent vertex is belong to second triangle
                                platePtUp = dualPlatePoint[secondTriangleIndex];

                            curveControlPtsPlate01[i].Add(platePtUp);
                            curveControlPtsPlate02[i].Add(platePtUp);

                            // half dual loop curve right close to plate
                            Point3d vertexPosition01 = iSpringMesh.Vertices[firstVertexIndex].Position;
                            
                            Point3d ctPPlate01 = (1 - iTangentScale) * platePtUp +
                                iTangentScale * (vertexPosition01 + iMinThickness * vertexNormals[firstVertexIndex]);

                            curveControlPtsPlate01[i].Add(ctPPlate01);
                            curveControlPtsPlate01[i].Add(vertexPosition01);

                            // half dual loop curve left close to plate
                            Point3d vertexPosition02 = iSpringMesh.Vertices[secondVertexIndex].Position;

                            Point3d ctPPlate02 = (1 - iTangentScale) * platePtUp +
                                iTangentScale * (vertexPosition02 + iMinThickness * vertexNormals[secondVertexIndex]);

                            curveControlPtsPlate02[i].Add(ctPPlate02);
                            curveControlPtsPlate02[i].Add(vertexPosition02);
                            
                            
                            // triLoop and dualLoop connection
                            Point3d midPt = 0.5 * (vertexPosition01 + iMinThickness * vertexNormals[firstVertexIndex] +
                                                   vertexPosition02 + iMinThickness * vertexNormals[secondVertexIndex]);

                            // half dual loop curve right close to triLoop
                            curveControlPtsTriLoop01[i].Add(midPt);
                            curveControlPtsTriLoop02[i].Add(midPt);

                            Point3d ctPTriLoop01 = (1 - iTangentScale) * midPt +
                                iTangentScale * (vertexPosition01 + iMinThickness * vertexNormals[firstVertexIndex]);
                            curveControlPtsTriLoop01[i].Add(ctPTriLoop01);
                            curveControlPtsTriLoop01[i].Add(vertexPosition01);

                            Point3d ctPTriLoop02 = (1 - iTangentScale) * midPt +
                                iTangentScale * (vertexPosition02 + iMinThickness * vertexNormals[secondVertexIndex]);
                            curveControlPtsTriLoop02[i].Add(ctPTriLoop02);
                            curveControlPtsTriLoop02[i].Add(vertexPosition02);

                        }

                        if (vertexStatus[firstAdjVIndex] == false && vertexStatus[secondAdjVIndex] == true &&
                            vertexStatus[firstVertexIndex] == false && vertexStatus[secondVertexIndex] == false)
                        {
                            int firstTriV01 = iSpringMesh.Triangles[firstTriangleIndex].FirstVertexIndex;
                            int firstTriV02 = iSpringMesh.Triangles[firstTriangleIndex].SecondVertexIndex;
                            int firstTriV03 = iSpringMesh.Triangles[firstTriangleIndex].ThirdVertexIndex;

                            // dualLoop connect to plate
                            Point3d platePtUp;

                            //// first adjacent vertex is belong to first triangle
                            if (firstTriV01 == firstAdjVIndex || firstTriV02 == firstAdjVIndex || firstTriV03 == firstAdjVIndex)
                                platePtUp = dualPlatePoint[secondTriangleIndex];
                            else //// first adjacent vertex is belong to second triangle
                                platePtUp = dualPlatePoint[firstTriangleIndex];

                            curveControlPtsPlate01[i].Add(platePtUp);
                            curveControlPtsPlate02[i].Add(platePtUp);

                            // half dual loop curve right close to plate
                            Point3d vertexPosition01 = iSpringMesh.Vertices[firstVertexIndex].Position;

                            Point3d ctPPlate01 = (1 - iTangentScale) * platePtUp +
                                iTangentScale * (vertexPosition01 + iMinThickness * vertexNormals[firstVertexIndex]);

                            curveControlPtsPlate01[i].Add(ctPPlate01);
                            curveControlPtsPlate01[i].Add(vertexPosition01);

                            // half dual loop curve left close to plate
                            Point3d vertexPosition02 = iSpringMesh.Vertices[secondVertexIndex].Position;

                            Point3d ctPPlate02 = (1 - iTangentScale) * platePtUp +
                                iTangentScale * (vertexPosition02 + iMinThickness * vertexNormals[secondVertexIndex]);

                            curveControlPtsPlate02[i].Add(ctPPlate02);
                            curveControlPtsPlate02[i].Add(vertexPosition02);


                            // triLoop and dualLoop connection
                            Point3d midPt = 0.5 * (vertexPosition01 + iMinThickness * vertexNormals[firstVertexIndex] +
                                                   vertexPosition02 + iMinThickness * vertexNormals[secondVertexIndex]);

                            // half dual loop curve right close to triLoop
                            curveControlPtsTriLoop01[i].Add(midPt);
                            curveControlPtsTriLoop02[i].Add(midPt);

                            Point3d ctPTriLoop01 = (1 - iTangentScale) * midPt +
                                iTangentScale * (vertexPosition01 + iMinThickness * vertexNormals[firstVertexIndex]);
                            curveControlPtsTriLoop01[i].Add(ctPTriLoop01);
                            curveControlPtsTriLoop01[i].Add(vertexPosition01);

                            Point3d ctPTriLoop02 = (1 - iTangentScale) * midPt +
                                iTangentScale * (vertexPosition02 + iMinThickness * vertexNormals[secondVertexIndex]);
                            curveControlPtsTriLoop02[i].Add(ctPTriLoop02);
                            curveControlPtsTriLoop02[i].Add(vertexPosition02);
                        }

                        // both plates
                        if (vertexStatus[firstAdjVIndex] == true && vertexStatus[secondAdjVIndex] == true &&
                            vertexStatus[firstVertexIndex] == false && vertexStatus[secondVertexIndex] == false)
                        {

                            Point3d platePtUp01 = dualPlatePoint[firstTriangleIndex];
                            Point3d platePtUp02 = dualPlatePoint[secondTriangleIndex];

                            curveControlPtsPlate01[i].Add(platePtUp01);
                            curveControlPtsTriLoop01[i].Add(platePtUp02);
                            //
                            curveControlPtsPlate02[i].Add(platePtUp01);
                            curveControlPtsTriLoop02[i].Add(platePtUp02);

                            Point3d vertexPosition01 = iSpringMesh.Vertices[firstVertexIndex].Position;
                            Point3d vertexPosition02 = iSpringMesh.Vertices[secondVertexIndex].Position;

                            Point3d ctPPlate01 = (1 - iTangentScale) * platePtUp01 +
                                iTangentScale * (vertexPosition01 + iMinThickness * vertexNormals[firstVertexIndex]);

                            curveControlPtsPlate01[i].Add(ctPPlate01);
                            curveControlPtsPlate01[i].Add(vertexPosition01);
                            //
                            Point3d ctPPlate11 = (1 - iTangentScale) * platePtUp01 +
                                iTangentScale * (vertexPosition02 + iMinThickness * vertexNormals[secondVertexIndex]);

                            curveControlPtsPlate02[i].Add(ctPPlate11);
                            curveControlPtsPlate02[i].Add(vertexPosition02);

                            Point3d ctPPlate02 = (1 - iTangentScale) * platePtUp02 +
                                iTangentScale * (vertexPosition01 + iMinThickness * vertexNormals[firstVertexIndex]);

                            curveControlPtsTriLoop01[i].Add(ctPPlate02);
                            curveControlPtsTriLoop01[i].Add(vertexPosition01);
                            //
                            Point3d ctPPlate12 = (1 - iTangentScale) * platePtUp02 +
                               iTangentScale * (vertexPosition02 + iMinThickness * vertexNormals[secondVertexIndex]);

                            curveControlPtsTriLoop02[i].Add(ctPPlate12);
                            curveControlPtsTriLoop02[i].Add(vertexPosition02);

                        }


                    }
                    else       // bottom layer
                    {
                        Point3d platePtRightDown;
                        Point3d platePtLeftDown;
                        if (vertexStatus[firstVertexIndex] == true && vertexStatus[secondVertexIndex] == false)
                        {
                            platePtRightDown = dualPlatePoint[firstTriangleIndex];
                            platePtLeftDown = dualPlatePoint[secondTriangleIndex];

                            //// here should be insert some control point so no kink in transition area
                            Point3d vertexPosition = iSpringMesh.Vertices[secondVertexIndex].Position;

                            Point3d ctPRight = (1 - iTangentScale) * platePtRightDown +
                                iTangentScale * (vertexPosition - iMinThickness * vertexNormals[secondVertexIndex]);
                            Point3d ctLeft = (1 - iTangentScale) * platePtLeftDown +
                                iTangentScale * (vertexPosition - iMinThickness * vertexNormals[secondVertexIndex]);

                            curveControlPtsRight[i].Add(ctPRight);
                            curveControlPtsLeft[i].Add(ctLeft);

                            curveControlPtsRight[i].Add(platePtRightDown);
                            curveControlPtsLeft[i].Add(platePtLeftDown);

                        }
                        if (vertexStatus[firstVertexIndex] == false && vertexStatus[secondVertexIndex] == true)
                        {
                            platePtRightDown = dualPlatePoint[firstTriangleIndex];
                            platePtLeftDown = dualPlatePoint[secondTriangleIndex];

                            //// here should be insert some control point so no kink in transition area
                            Point3d vertexPosition = iSpringMesh.Vertices[firstVertexIndex].Position;

                            Point3d ctPRight = (1 - iTangentScale) * platePtRightDown +
                                iTangentScale * (vertexPosition - iMinThickness * vertexNormals[firstVertexIndex]);
                            Point3d ctLeft = (1 - iTangentScale) * platePtLeftDown +
                                iTangentScale * (vertexPosition - iMinThickness * vertexNormals[firstVertexIndex]);

                            curveControlPtsRight[i].Add(ctPRight);
                            curveControlPtsLeft[i].Add(ctLeft);

                            curveControlPtsRight[i].Add(platePtRightDown);
                            curveControlPtsLeft[i].Add(platePtLeftDown);

                        }

                        ///------------------------------------------------------------------------------------------------
                        ///// half dual loop without connecting to the plates
                        if (vertexStatus[firstAdjVIndex] == true && vertexStatus[secondAdjVIndex] == false &&
                            vertexStatus[firstVertexIndex] == false && vertexStatus[secondVertexIndex] == false)
                        {
                            int firstTriV01 = iSpringMesh.Triangles[firstTriangleIndex].FirstVertexIndex;
                            int firstTriV02 = iSpringMesh.Triangles[firstTriangleIndex].SecondVertexIndex;
                            int firstTriV03 = iSpringMesh.Triangles[firstTriangleIndex].ThirdVertexIndex;

                            // dualLoop connect to plate
                            Point3d platePtDown;

                            //// first adjacent vertex is belong to first triangle
                            if (firstTriV01 == firstAdjVIndex || firstTriV02 == firstAdjVIndex || firstTriV03 == firstAdjVIndex)
                                platePtDown = dualPlatePoint[firstTriangleIndex];
                            else //// first adjacent vertex is belong to second triangle
                                platePtDown = dualPlatePoint[secondTriangleIndex];

                            // half dual loop curve right close to plate
                            Point3d vertexPosition01 = iSpringMesh.Vertices[firstVertexIndex].Position;

                            Point3d ctPPlate01 = (1 - iTangentScale) * platePtDown +
                                iTangentScale * (vertexPosition01 - iMinThickness * vertexNormals[firstVertexIndex]);

                            curveControlPtsPlate01[i].Add(ctPPlate01);

                            // half dual loop curve left close to plate
                            Point3d vertexPosition02 = iSpringMesh.Vertices[secondVertexIndex].Position;

                            Point3d ctPPlate02 = (1 - iTangentScale) * platePtDown +
                                iTangentScale * (vertexPosition02 - iMinThickness * vertexNormals[secondVertexIndex]);

                            curveControlPtsPlate02[i].Add(ctPPlate02);

                            curveControlPtsPlate01[i].Add(platePtDown);
                            curveControlPtsPlate02[i].Add(platePtDown);


                            // triLoop and dualLoop connection
                            Point3d midPt = 0.5 * (vertexPosition01 - iMinThickness * vertexNormals[firstVertexIndex] +
                                                   vertexPosition02 - iMinThickness * vertexNormals[secondVertexIndex]);

                            Point3d ctPTriLoop01 = (1 - iTangentScale) * midPt +
                                iTangentScale * (vertexPosition01 - iMinThickness * vertexNormals[firstVertexIndex]);
                            curveControlPtsTriLoop01[i].Add(ctPTriLoop01);

                            Point3d ctPTriLoop02 = (1 - iTangentScale) * midPt +
                                iTangentScale * (vertexPosition02 - iMinThickness * vertexNormals[secondVertexIndex]);
                            curveControlPtsTriLoop02[i].Add(ctPTriLoop02);

                            // half dual loop curve right close to triLoop
                            curveControlPtsTriLoop01[i].Add(midPt);
                            curveControlPtsTriLoop02[i].Add(midPt);
                        }

                        if (vertexStatus[firstAdjVIndex] == false && vertexStatus[secondAdjVIndex] == true &&
                            vertexStatus[firstVertexIndex] == false && vertexStatus[secondVertexIndex] == false)
                        {
                            int firstTriV01 = iSpringMesh.Triangles[firstTriangleIndex].FirstVertexIndex;
                            int firstTriV02 = iSpringMesh.Triangles[firstTriangleIndex].SecondVertexIndex;
                            int firstTriV03 = iSpringMesh.Triangles[firstTriangleIndex].ThirdVertexIndex;

                            // dualLoop connect to plate
                            Point3d platePtDown;

                            //// first adjacent vertex is belong to first triangle
                            if (firstTriV01 == firstAdjVIndex || firstTriV02 == firstAdjVIndex || firstTriV03 == firstAdjVIndex)
                                platePtDown = dualPlatePoint[secondTriangleIndex];
                            else //// first adjacent vertex is belong to second triangle
                                platePtDown = dualPlatePoint[firstTriangleIndex];

                            // half dual loop curve right close to plate
                            Point3d vertexPosition01 = iSpringMesh.Vertices[firstVertexIndex].Position;

                            Point3d ctPPlate01 = (1 - iTangentScale) * platePtDown +
                                iTangentScale * (vertexPosition01 - iMinThickness * vertexNormals[firstVertexIndex]);

                            curveControlPtsPlate01[i].Add(ctPPlate01);

                            // half dual loop curve left close to plate
                            Point3d vertexPosition02 = iSpringMesh.Vertices[secondVertexIndex].Position;

                            Point3d ctPPlate02 = (1 - iTangentScale) * platePtDown +
                                iTangentScale * (vertexPosition02 - iMinThickness * vertexNormals[secondVertexIndex]);

                            curveControlPtsPlate02[i].Add(ctPPlate02);

                            curveControlPtsPlate01[i].Add(platePtDown);
                            curveControlPtsPlate02[i].Add(platePtDown);


                            // triLoop and dualLoop connection
                            Point3d midPt = 0.5 * (vertexPosition01 - iMinThickness * vertexNormals[firstVertexIndex] +
                                                   vertexPosition02 - iMinThickness * vertexNormals[secondVertexIndex]);

                            Point3d ctPTriLoop01 = (1 - iTangentScale) * midPt +
                                iTangentScale * (vertexPosition01 - iMinThickness * vertexNormals[firstVertexIndex]);
                            curveControlPtsTriLoop01[i].Add(ctPTriLoop01);

                            Point3d ctPTriLoop02 = (1 - iTangentScale) * midPt +
                                iTangentScale * (vertexPosition02 - iMinThickness * vertexNormals[secondVertexIndex]);
                            curveControlPtsTriLoop02[i].Add(ctPTriLoop02);

                            // half dual loop curve right close to triLoop
                            curveControlPtsTriLoop01[i].Add(midPt);
                            curveControlPtsTriLoop02[i].Add(midPt);

                        }

                        // both plates
                        if (vertexStatus[firstAdjVIndex] == true && vertexStatus[secondAdjVIndex] == true &&
                            vertexStatus[firstVertexIndex] == false && vertexStatus[secondVertexIndex] == false)
                        {

                            Point3d platePtDown01 = dualPlatePoint[firstTriangleIndex];
                            Point3d platePtDown02 = dualPlatePoint[secondTriangleIndex];

                            Point3d vertexPosition01 = iSpringMesh.Vertices[firstVertexIndex].Position;
                            Point3d vertexPosition02 = iSpringMesh.Vertices[secondVertexIndex].Position;
                            

                            Point3d ctPPlate01 = (1 - iTangentScale) * platePtDown01 +
                                iTangentScale * (vertexPosition01 - iMinThickness * vertexNormals[firstVertexIndex]);

                            curveControlPtsPlate01[i].Add(ctPPlate01);
                            //
                            Point3d ctPPlate11 = (1 - iTangentScale) * platePtDown01 +
                                iTangentScale * (vertexPosition02 - iMinThickness * vertexNormals[secondVertexIndex]);

                            curveControlPtsPlate02[i].Add(ctPPlate11);

                            Point3d ctPPlate02 = (1 - iTangentScale) * platePtDown02 +
                                iTangentScale * (vertexPosition01 - iMinThickness * vertexNormals[firstVertexIndex]);

                            curveControlPtsTriLoop01[i].Add(ctPPlate02);
                            //
                            Point3d ctPPlate12 = (1 - iTangentScale) * platePtDown02 +
                               iTangentScale * (vertexPosition02 - iMinThickness * vertexNormals[secondVertexIndex]);

                            curveControlPtsTriLoop02[i].Add(ctPPlate12);

                            curveControlPtsPlate01[i].Add(platePtDown01);
                            //curveControlPtsPlate01[i].Add(platePtUp02);
                            //curveControlPtsPlate02[i].Add(platePtUp01);
                            curveControlPtsTriLoop01[i].Add(platePtDown02);
                            //
                            curveControlPtsPlate02[i].Add(platePtDown01);
                            curveControlPtsTriLoop02[i].Add(platePtDown02);

                        }
                        
                    }
                }
            }

        }


        private void dualLoop()
        {

        }


        public List<double> ComputeVertexFactor()
        {
            List<double> vertexFactors = new List<double>();

            for (int i = 0; i < iSpringMesh.Vertices.Count; i++)
            {
                vertexFactors.Add(0.5);
            }

            return vertexFactors;
        }

        private void triLoop()
        {
            List<Vector3d> vertexNormals = iSpringMesh.ComputeVertexNormals();

            foreach (Edge edge in iSpringMesh.Edges)
            {
                if (edge.SecondTriangleIndex >= 0) continue;

                // Gene Added
                if (iSpringMesh.Vertices[edge.FirstVertexIndex].Position.Z > 0.65 &&
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
            // here is to get all vertex normal and utilize them

            // =========================================================================================================


            bottomCps = new List<Point3d>();
            topCps = new List<Point3d>();


            for (int i = 0; i < iSpringMesh.Vertices.Count; i++)
            {
                Point3d vertexPosition = iSpringMesh.Vertices[i].Position;
                // different thickness
                double k = vertexPosition.Z * 0.125;
                double thickness = k * iMinThickness + (0.5 - k) * iMaxThickness;
                //bottomCps.Add(vertexPosition - thickness * vertexNormals[i]);
                //topCps.Add(vertexPosition + thickness * vertexNormals[i]);

                bottomCps.Add(vertexPosition - iMinThickness * vertexNormals[i]);
                topCps.Add(vertexPosition + iMinThickness * vertexNormals[i]);
            }


            // =========================================================================================================


            foreach (Triangle triangle in iSpringMesh.Triangles)
            {
                if (!vertexStatus[triangle.FirstVertexIndex] &&
                    !vertexStatus[triangle.SecondVertexIndex] &&
                    !vertexStatus[triangle.ThirdVertexIndex])
                {
                    createStripe(triangle.FirstVertexIndex, triangle.SecondVertexIndex, triangle.ThirdVertexIndex);
                    createStripe(triangle.SecondVertexIndex, triangle.ThirdVertexIndex, triangle.FirstVertexIndex);
                    createStripe(triangle.ThirdVertexIndex, triangle.FirstVertexIndex, triangle.SecondVertexIndex);
                }
            }

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

                //oDebugCurves1.Add(profileCurve1);
                //oDebugCurves2.Add(profileCurve2);
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
                //oDebugCurves1.Add(profileCurve1);
                //oDebugCurves2.Add(profileCurve2);
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
            get { return new Guid("{738b1c5c-8134-4658-a49b-f255bf4c4a2d}"); }
        }
    }
}