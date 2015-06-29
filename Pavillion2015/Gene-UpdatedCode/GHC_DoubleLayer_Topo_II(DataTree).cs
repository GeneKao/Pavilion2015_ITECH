using System;
using System.Collections.Generic;
using System.Linq;

using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;
using Rhino.Geometry;
using Rhino.Collections;

namespace Pavillion2015.Gene_UpdatedCode
{
    public class GHC_DoubleLayer_Topo_II_DataTree_ : GH_Component
    {

        // input
        SpringMesh iSpringMesh = null;
        List<Curve> iPolyLine = null;
        List<bool> iVertexStatus = null;
        List<int> iPolyLineID = null;
        List<double> iThickness = null;
        bool iPolySrf = true;
        double iTangentScaleMin = double.NaN;
        double iTangentScaleMax = double.NaN;
        List<Point3d> iAttractors = null;
        List<Point3d> iClosedPanelPts = null;
        double iClosePanelDist = 0.001;

        //=================================== EDITED BY JULIAN =========================================
        double iPlanarOffsetScaleMin = double.NaN;
        double iPlanarOffsetScaleMax = double.NaN;

        double CurvePointiness = double.NaN;
        //=================================== END EDITED BY JULIAN =====================================

        // output
        string oInfo = string.Empty;
        List<Point3d> oDebugList = null;
        DataTree<Brep> oTriLoop = null;
        DataTree<string> oTriLoopID = null;
        DataTree<Curve> oTriLoopCurves = null;
        DataTree<Brep> oDualLoop1 = null;
        DataTree<Brep> oDualLoop2 = null;
        DataTree<string> oDualLoop1ID = null;
        DataTree<string> oDualLoop2ID = null;
        DataTree<Curve> oDualLoop1Curves = null;
        DataTree<Curve> oDualLoop2Curves = null;
        DataTree<Brep> oClosedPanel = null;


        // internal usage
        List<Vector3d> vertexNormals;

        List<Point3d> topCps = null;
        List<Point3d> bottomCps = null;

        List<Point3d> topCenterPts = null;
        List<Point3d> bottomCenterPts = null;
        List<bool> topCenterPolygonOccupied = null;

        List<int[]> indexSortPolygon = null;

        //=================================== EDITED BY JULIAN =========================================
        // General Distance Relation to Attractor(s) ( before remaping to Min/Max - Range )
        List<double> verticesValues = null;

        // relative offset factors for individual opnening sizes (remaped verticesValues to iTangentScaleMin - iTangentScaleMax - range)
        List<double> curveVerticesValues = null;

        // offset distance in document unit for individual planar part sizes (remaped verticesValues to iPlanarOffsetScaleMin - iPlanarOffsetScaleMax - range)
        List<double> planarVerticesValues = null;

        // opens or closes openings (Plate "Type 2")
        List<bool> openClose = null;

        //=================================== END EDITED BY JULIAN =====================================

        double documentTolerance = DocumentTolerance();



        public GHC_DoubleLayer_Topo_II_DataTree_()
            : base("Double Layer Topo - II DT", "Double Layer Topo - II DT",
                "Double Layer Topology II DT - get plate and SpringMesh convert double layers - output data tree",
                "Pavillion 2015", "Double Layer")
        {
        }


        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            //=================================== EDITED BY JULIAN =========================================

            pManager.AddGenericParameter("Spring Mesh", "Spring Mesh", "Spring Mesh", GH_ParamAccess.item);
            pManager.AddCurveParameter("PolyLine", "PolyLine", "PolyLine from Plate", GH_ParamAccess.list);
            pManager.AddNumberParameter("Thickness", "Thickness", "Thickness of Component", GH_ParamAccess.list);
            pManager.AddBooleanParameter("Vertex status", "Vertex status", "Vertex status", GH_ParamAccess.list);
            pManager.AddIntegerParameter("PolyLine ID", "PolyLine ID", "PolyLine ID", GH_ParamAccess.list);
            pManager.AddBooleanParameter("Allow PolySrf", "Allow PolySrf", "Allow PolySrf", GH_ParamAccess.item, true);
            pManager.AddNumberParameter("TangentScale Min", "TangentScale Min", "TangentScale Min [relative]", GH_ParamAccess.item, 0.2);
            pManager.AddNumberParameter("TangentScale Max", "TangentScale Max", "TangentScale Max [relative]", GH_ParamAccess.item, 0.8);
            pManager.AddNumberParameter("Curve Pointiness", "Pointiness", "Pointiness of the bended Surfaces [relative]", GH_ParamAccess.item, 1.0);
            pManager.AddNumberParameter("PlanarOffset Min", "PlanarOffset Min", "Controlls minimal offset of planar parts [in doc. units]", GH_ParamAccess.item, 0.2);
            pManager.AddNumberParameter("PlanarOffset Max", "PlanarOffset Max", "Controlls maximal offset of planar parts [in doc. units]", GH_ParamAccess.item, 0.8);
            pManager.AddBooleanParameter("Open/Close", "O/C", "Open or close an opening", GH_ParamAccess.list);
            pManager.AddPointParameter("Attractors", "Attractors", "Attractors", GH_ParamAccess.list);

            //=================================== END EDITED BY JULIAN =====================================

            pManager.AddPointParameter("Closed Panel Area", "Closed Panel Area", "Closed Panel Area", GH_ParamAccess.list, new Point3d());
            pManager.AddNumberParameter("Panel Effect Area", "Panel Effect Area", "Panel Effect Area", GH_ParamAccess.item, 0.001);

        }


        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddTextParameter("Info", "Debug - Info", "Information", GH_ParamAccess.item);                  // 0
            pManager.AddGenericParameter("General List", "Degug - List ", "Debug", GH_ParamAccess.list);            // 1
            pManager.AddGenericParameter("TriLoop", "TriLoop", "TriLoop", GH_ParamAccess.list);                     // 2
            pManager.AddGenericParameter("TriLoopID", "TriLoopID", "TriLoopID", GH_ParamAccess.list);               // 3
            pManager.AddGenericParameter("TriLoopCurves", "TriLoopCurves", "TriLoopCurves", GH_ParamAccess.tree);   // 4
            pManager.AddGenericParameter("DualLoop1", "DualLoop1", "DualLoop1", GH_ParamAccess.tree);               // 5
            pManager.AddGenericParameter("DualLoop2", "DualLoop2", "DualLoop2", GH_ParamAccess.tree);               // 6
            pManager.AddGenericParameter("DualLoop1ID", "DualLoop1ID", "DualLoop1ID", GH_ParamAccess.tree);         // 7
            pManager.AddGenericParameter("DualLoop2ID", "DualLoop2ID", "DualLoop2ID", GH_ParamAccess.tree);         // 8
            pManager.AddGenericParameter("DualLoop1Curves", "DualLoop1Curves", "DualLoop1Curves", GH_ParamAccess.tree);         // 9
            pManager.AddGenericParameter("DualLoop2Curves", "DualLoop2Curves", "DualLoop2Curves", GH_ParamAccess.tree);         // 10
            pManager.AddGenericParameter("ClosedPanel", "ClosedPanel", "ClosedPanel", GH_ParamAccess.tree);         // 11
        }


        protected override void BeforeSolveInstance()
        {
            // input

            iSpringMesh = new SpringMesh();
            iPolyLine = new List<Curve>();
            iVertexStatus = new List<bool>(); // true is plate occupied, false should put triLoop
            iPolyLineID = new List<int>(); // for vertex status access the order of polyline
            iThickness = new List<double>();
            iAttractors = new List<Point3d>();
            iClosedPanelPts = new List<Point3d>();

            // output
            oInfo = string.Empty;
            oDebugList = new List<Point3d>();
            oTriLoop = new DataTree<Brep>();
            oTriLoopID = new DataTree<string>();
            oTriLoopCurves = new DataTree<Curve>();
            oDualLoop1 = new DataTree<Brep>();
            oDualLoop2 = new DataTree<Brep>();
            oDualLoop1ID = new DataTree<string>();
            oDualLoop2ID = new DataTree<string>();
            oDualLoop1Curves = new DataTree<Curve>();
            oDualLoop2Curves = new DataTree<Curve>();
            oClosedPanel = new DataTree<Brep>();

            // internal use data
            topCenterPts = new List<Point3d>();
            bottomCenterPts = new List<Point3d>();
            topCenterPolygonOccupied = new List<bool>(); // not be used at this moment 

            bottomCps = new List<Point3d>();
            topCps = new List<Point3d>();

            indexSortPolygon = new List<int[]>();   // array is reference type in csharp

            verticesValues = new List<double>();
        }


        protected override void SolveInstance(IGH_DataAccess DA)
        {
            // getting input
            DA.GetData<SpringMesh>("Spring Mesh", ref iSpringMesh);
            DA.GetDataList<Curve>("PolyLine", iPolyLine);
            DA.GetDataList<int>("PolyLine ID", iPolyLineID);
            DA.GetDataList<double>("Thickness", iThickness);
            DA.GetDataList<bool>("Vertex status", iVertexStatus);
            DA.GetData<bool>("Allow PolySrf", ref iPolySrf);
            DA.GetData<double>("TangentScale Min", ref iTangentScaleMin);
            DA.GetData<double>("TangentScale Max", ref iTangentScaleMax);
            DA.GetDataList<Point3d>("Attractors", iAttractors);
            DA.GetDataList<Point3d>("Closed Panel Area", iClosedPanelPts);
            DA.GetData<double>("Panel Effect Area", ref iClosePanelDist);
            //------------------------------------------------------------


            storePlatesTPI();

            calculateVertexNormals();
            calculateVertexCps();
            calculateVerticesValues();

            triLoop();

            half_dualLoop();

            //------------------------------------------------------------

            // setting output
            DA.SetData("Info", oInfo);
            DA.SetDataList("General List", oDebugList);
            DA.SetDataTree(2, oTriLoop);
            DA.SetDataTree(3, oTriLoopID);
            DA.SetDataTree(4, oTriLoopCurves);
            DA.SetDataTree(5, oDualLoop1);
            DA.SetDataTree(6, oDualLoop2);
            DA.SetDataTree(7, oDualLoop1ID);
            DA.SetDataTree(8, oDualLoop2ID);
            DA.SetDataTree(9, oDualLoop1Curves);
            DA.SetDataTree(10, oDualLoop2Curves);
            DA.SetDataTree(11, oClosedPanel);
            // -----------------------------------------------------------
        }

        private void storePlatesTPI()
        {
            // initial 
            // 3 for polygon id, 3 for 3 tri side id. use the number to compare if polygon is clockwise. 
            int[] triIndex = new int[6] { -1, -1, -1, -1, -1, -1 };

            // construct center pt
            for (int i = 0; i < iSpringMesh.Triangles.Count; i++)
            {
                topCenterPts.Add(iSpringMesh.ComputeTriangleCentroid(i));
                bottomCenterPts.Add(iSpringMesh.ComputeTriangleCentroid(i));
                topCenterPolygonOccupied.Add(false);    // to see if triangle be occupied by polygon

                int[] copy = (int[])triIndex.Clone();   // have to copy array otherwise it's reference type it will change all
                indexSortPolygon.Add(copy);
            }


            // compare center pt with TPI
            // get discontinuity points 
            for (int i = 0; i < iPolyLine.Count; i++)
            {
                Curve c = iPolyLine[i];
                double t;
                double start = c.Domain.Min;
                int counter = 0;
                while (c.GetNextDiscontinuity(Continuity.C1_locus_continuous, start, c.Domain.Max, out t))
                {

                    start = t;
                    Point3d pt = c.PointAt(t);
                    if (i < iPolyLine.Count / 2)
                    {
                        Point3d closestPt = Point3dList.ClosestPointInList(topCenterPts, pt);
                        int index = topCenterPts.IndexOf(closestPt);
                        topCenterPts[index] = pt;
                        topCenterPolygonOccupied[index] = true; // set occupied status

                        for (int j = 0; j < 3; j++)
                            if (indexSortPolygon[index][j] == -1)
                            {
                                indexSortPolygon[index][j] = i;
                                indexSortPolygon[index][j + 3] = counter;

                                break;
                            }
                        counter++;
                    }
                    else
                    {
                        Point3d closestPt = Point3dList.ClosestPointInList(bottomCenterPts, pt);
                        int index = bottomCenterPts.IndexOf(closestPt);
                        bottomCenterPts[index] = pt;
                    }
                }
            }
            foreach (Point3d p in topCenterPts)
                oDebugList.Add(p);
            foreach (Point3d p in bottomCenterPts)
                oDebugList.Add(p);

        }

        private void half_dualLoop()
        {
            for (int i = 0; i < iSpringMesh.Edges.Count; i++)
            {
                Edge edgeLocal = iSpringMesh.Edges[i];

                int tri01 = edgeLocal.FirstTriangleIndex;

                int vertex01 = edgeLocal.FirstVertexIndex;
                int vertex02 = edgeLocal.SecondVertexIndex;

                if (iSpringMesh.Edges[i].SecondTriangleIndex >= 0)
                {
                    int tri02 = edgeLocal.SecondTriangleIndex;

                    if ((iVertexStatus[vertex01] == true && iVertexStatus[vertex02] == false) ||
                        (iVertexStatus[vertex01] == false && iVertexStatus[vertex02] == true)) half_dualLoop_type_1(i);

                    if ((iVertexStatus[vertex01] == false && iVertexStatus[vertex02] == false))
                    {
                        half_dualLoop_type_2(i, tri01, tri02);
                        half_dualLoop_type_2(i, tri02, tri01);
                    }
                }
                else
                {
                    if ((iVertexStatus[vertex01] == false && iVertexStatus[vertex02] == false))
                    {
                        half_dualLoop_type_2(i, tri01, -1); // give it -1 if no neighbour
                    }
                }
            }
        }

        // dual loop directly connected to plates
        private void half_dualLoop_type_1(int edgeIndex)
        {
            Edge edge = iSpringMesh.Edges[edgeIndex];

            int rightTriIndex = edge.FirstTriangleIndex;
            int leftTriIndex = edge.SecondTriangleIndex;

            // vertex point
            int vertex1 = edge.FirstVertexIndex;
            int vertex2 = edge.SecondVertexIndex;

            int vertex = vertex1;
            int occupyiedVertex = vertex2;

            if (iVertexStatus[vertex1] == false && iVertexStatus[vertex2] == true) { vertex = vertex1; occupyiedVertex = vertex2; }
            else if (iVertexStatus[vertex2] == false && iVertexStatus[vertex1] == true) { vertex = vertex2; occupyiedVertex = vertex1; }
            else { oInfo += "got it..."; }

            Point3d vertexPt = iSpringMesh.Vertices[vertex].Position;

            // layer up
            Point3d rightUp01 = topCenterPts[rightTriIndex];
            Point3d leftUp01 = topCenterPts[leftTriIndex];

            Point3d up03 = topCps[vertex]; //reference point not for construct surfaces

            double scaler = verticesValues[vertex];

            Point3d rightUp02 = (1 - scaler) * rightUp01 + (scaler) * up03;
            Point3d leftUp02 = (1 - scaler) * leftUp01 + (scaler) * up03;

            // layer down
            Point3d rightDown01 = bottomCenterPts[rightTriIndex];
            Point3d leftDown01 = bottomCenterPts[leftTriIndex];

            Point3d down03 = bottomCps[vertex]; //reference point not for construct surfaces

            Point3d rightDown02 = (1 - scaler) * rightDown01 + (scaler) * down03;
            Point3d leftDown02 = (1 - scaler) * leftDown01 + (scaler) * down03;

            Curve right = Curve.CreateControlPointCurve(
                new List<Point3d>() { rightUp01, rightUp02, vertexPt, rightDown02, rightDown01 });

            Curve left = Curve.CreateControlPointCurve(
                new List<Point3d>() { leftUp01, leftUp02, vertexPt, leftDown02, leftDown01 });

            // dataTree 
            int plateID = iPolyLineID[occupyiedVertex];

            GH_Path path = new GH_Path(plateID);

            // compare their polygon index to sort the lofting sequence 
            int[] rightIndex = indexSortPolygon[rightTriIndex];
            int[] leftIndex = indexSortPolygon[leftTriIndex];

            for (int r = 0; r < 3; r++)
                for (int l = 0; l < 3; l++)
                    if (rightIndex[r] == leftIndex[l] && rightIndex[r] != -1 && leftIndex[l] != -1) // check to use right polygon
                        if ((rightIndex[r + 3] == 0 && leftIndex[l + 3] != 1) ||
                            (leftIndex[l + 3] == 0 && rightIndex[r + 3] != 1) && (rightIndex[r + 3] < leftIndex[l + 3]))
                        {
                            Brep[] brep = Brep.CreateFromLoft(
                                new List<Curve>() { right, left },
                                Point3d.Unset, Point3d.Unset,
                                LoftType.Normal,
                                false
                                );
                            if (brep.Length > 0)
                            {
                                oDualLoop1.Add(brep[0], path);
                                oDualLoop1ID.Add("H;" + rightTriIndex.ToString() + "-" + leftTriIndex.ToString() + ";" + plateID.ToString(), path);
                                
                                oDualLoop1Curves.Add(right, path);
                                oDualLoop1Curves.Add(left, path);
                            }
                        }
                        else if ((rightIndex[r + 3] == 0 && leftIndex[l + 3] != 1) ||
                                 (leftIndex[l + 3] == 0 && rightIndex[r + 3] != 1) && (rightIndex[r + 3] >= leftIndex[l + 3]))
                        {
                            Brep[] brep = Brep.CreateFromLoft(
                                new List<Curve>() { left, right },
                                Point3d.Unset, Point3d.Unset,
                                LoftType.Normal,
                                false
                                );
                            if (brep.Length > 0)
                            {
                                oDualLoop1.Add(brep[0], path);
                                oDualLoop1ID.Add("H;" + rightTriIndex.ToString() + "-" + leftTriIndex.ToString() + ";" + plateID.ToString(), path);

                                oDualLoop1Curves.Add(left, path);
                                oDualLoop1Curves.Add(right, path);
                            }
                        }
                        else if (rightIndex[r + 3] >= leftIndex[l + 3])
                        {
                            Brep[] brep = Brep.CreateFromLoft(
                                new List<Curve>() { right, left },
                                Point3d.Unset, Point3d.Unset,
                                LoftType.Normal,
                                false
                                );
                            if (brep.Length > 0)
                            {
                                oDualLoop1.Add(brep[0], path);
                                oDualLoop1ID.Add("H;" + rightTriIndex.ToString() + "-" + leftTriIndex.ToString() + ";" + plateID.ToString(), path);

                                oDualLoop1Curves.Add(right, path);
                                oDualLoop1Curves.Add(left, path);
                            }
                        }
                        else if (rightIndex[r + 3] < leftIndex[l + 3])
                        {
                            Brep[] brep = Brep.CreateFromLoft(
                                new List<Curve>() { left, right },
                                Point3d.Unset, Point3d.Unset,
                                LoftType.Normal,
                                false
                                );
                            if (brep.Length > 0)
                            {
                                oDualLoop1.Add(brep[0], path);
                                oDualLoop1ID.Add("H;" + rightTriIndex.ToString() + "-" + leftTriIndex.ToString() + ";" + plateID.ToString(), path);

                                oDualLoop1Curves.Add(left, path);
                                oDualLoop1Curves.Add(right, path);
                            }
                        }

        }

        // dual loop not directly connected to plates
        private void half_dualLoop_type_2(int edgeIndex, int triangleIndex, int neighbourTriIndex)
        {
            GH_Path path = new GH_Path(triangleIndex);

            Edge edge = iSpringMesh.Edges[edgeIndex];

            // vertex point
            int vertex1 = edge.FirstVertexIndex;
            int vertex2 = edge.SecondVertexIndex;

            Point3d vertexPt1 = iSpringMesh.Vertices[vertex1].Position;
            Point3d vertexPt2 = iSpringMesh.Vertices[vertex2].Position;

            double scaler1 = verticesValues[vertex1];
            double scaler2 = verticesValues[vertex2];

            Point3d vertexPtUp1 = topCps[vertex1];
            Point3d vertexPtUp2 = topCps[vertex2];
            Point3d vertexPtUpM = 0.5 * (vertexPtUp1 + vertexPtUp2);
            Point3d ctPtUp1 = (scaler1) * vertexPtUp1 + (1 - scaler1) * vertexPtUpM;
            Point3d ctPtUp2 = (scaler2) * vertexPtUp2 + (1 - scaler2) * vertexPtUpM;

            Point3d vertexPtDown1 = bottomCps[vertex1];
            Point3d vertexPtDown2 = bottomCps[vertex2];
            Point3d vertexPtDownM = 0.5 * (vertexPtDown1 + vertexPtDown2);
            Point3d ctPtDown1 = (scaler1) * vertexPtDown1 + (1 - scaler1) * vertexPtDownM;
            Point3d ctPtDown2 = (scaler2) * vertexPtDown2 + (1 - scaler2) * vertexPtDownM;

            Point3d upTPI = topCenterPts[triangleIndex];
            Point3d downTPI = bottomCenterPts[triangleIndex];
            Point3d ctUpTPI1 = (scaler1) * vertexPtUp1 + (1 - scaler1) * upTPI;
            Point3d ctDownTPI1 = (scaler1) * vertexPtDown1 + (1 - scaler1) * downTPI;
            Point3d ctUpTPI2 = (scaler2) * vertexPtUp2 + (1 - scaler2) * upTPI;
            Point3d ctDownTPI2 = (scaler2) * vertexPtDown2 + (1 - scaler2) * downTPI;

            Curve right1 = Curve.CreateControlPointCurve(
                new List<Point3d>() { upTPI, ctUpTPI1, vertexPt1, ctDownTPI1, downTPI });

            Curve left1 = Curve.CreateControlPointCurve(
                new List<Point3d>() { vertexPtUpM, ctPtUp1, vertexPt1, ctPtDown1, vertexPtDownM });

            Curve right2 = Curve.CreateControlPointCurve(
                new List<Point3d>() { upTPI, ctUpTPI2, vertexPt2, ctDownTPI2, downTPI });

            Curve left2 = Curve.CreateControlPointCurve(
                new List<Point3d>() { vertexPtUpM, ctPtUp2, vertexPt2, ctPtDown2, vertexPtDownM });

            // sorting sequence 

            Brep[] brep1 = Brep.CreateFromLoft(
                       new List<Curve>() { right1, left1 },
                       Point3d.Unset, Point3d.Unset,
                       LoftType.Normal,
                       false
                       );
            if (brep1.Length > 0)
            {
                //oDualLoop.Add(brep1[0]);
                Brep brepTestNormal = brep1[0];
                BrepFace brepF = brepTestNormal.Faces[0];

                // id 
                int neighbour2TriIndex = 0;
                if (vertex2 == iSpringMesh.Triangles[triangleIndex].FirstVertexIndex)
                    neighbour2TriIndex = iSpringMesh.Triangles[triangleIndex].FirstAdjTriIndex;
                if (vertex2 == iSpringMesh.Triangles[triangleIndex].SecondVertexIndex)
                    neighbour2TriIndex = iSpringMesh.Triangles[triangleIndex].SecondAdjTriIndex;
                if (vertex2 == iSpringMesh.Triangles[triangleIndex].ThirdVertexIndex)
                    neighbour2TriIndex = iSpringMesh.Triangles[triangleIndex].ThirdAdjTriIndex; 

                if (brepF.NormalAt(0, 0).Z > 0)
                {
                    Brep[] brepN1 = Brep.CreateFromLoft(
                       new List<Curve>() { left1, right1 },
                       Point3d.Unset, Point3d.Unset,
                       LoftType.Normal,
                       false
                       );
                    if (brep1.Length > 0)
                    {
                        oDualLoop2.Add(brepN1[0], path.AppendElement(0));
                        oDualLoop2ID.Add("J;" + triangleIndex.ToString() + ";" + neighbourTriIndex.ToString() + ";" + neighbour2TriIndex.ToString(), path.AppendElement(0));

                        oDualLoop2Curves.Add(left1, path.AppendElement(0));
                        oDualLoop2Curves.Add(right1, path.AppendElement(0));
                    }
                }
                else
                {
                    oDualLoop2.Add(brep1[0], path.AppendElement(0));
                    oDualLoop2ID.Add("J;" + triangleIndex.ToString() + ";" + neighbourTriIndex.ToString() + ";" + neighbour2TriIndex.ToString(), path.AppendElement(0));
                    
                    oDualLoop2Curves.Add(right1, path.AppendElement(0));
                    oDualLoop2Curves.Add(left1, path.AppendElement(0));
                }
            }


            Brep[] brep2 = Brep.CreateFromLoft(
                       new List<Curve>() { right2, left2 },
                       Point3d.Unset, Point3d.Unset,
                       LoftType.Normal,
                       false
                       );
            if (brep2.Length > 0)
            {
                //oDualLoop.Add(brep2[0]);
                Brep brepTestNormal = brep2[0];
                BrepFace brepF = brepTestNormal.Faces[0];

                // id 
                int neighbour2TriIndex = 0;
                if (vertex1 == iSpringMesh.Triangles[triangleIndex].FirstVertexIndex)
                    neighbour2TriIndex = iSpringMesh.Triangles[triangleIndex].FirstAdjTriIndex;
                if (vertex1 == iSpringMesh.Triangles[triangleIndex].SecondVertexIndex)
                    neighbour2TriIndex = iSpringMesh.Triangles[triangleIndex].SecondAdjTriIndex;
                if (vertex1 == iSpringMesh.Triangles[triangleIndex].ThirdVertexIndex)
                    neighbour2TriIndex = iSpringMesh.Triangles[triangleIndex].ThirdAdjTriIndex; 

                if (brepF.NormalAt(0, 0).Z > 0)
                {
                    Brep[] brepN2 = Brep.CreateFromLoft(
                       new List<Curve>() { left2, right2 },
                       Point3d.Unset, Point3d.Unset,
                       LoftType.Normal,
                       false
                       );
                    if (brep2.Length > 0)
                    {
                        oDualLoop2.Add(brepN2[0], path.AppendElement(1));
                        oDualLoop2ID.Add("J;" + triangleIndex.ToString() + ";" + neighbourTriIndex.ToString() + ";" + neighbour2TriIndex.ToString(), path.AppendElement(1));

                        oDualLoop2Curves.Add(left2, path.AppendElement(1));
                        oDualLoop2Curves.Add(right2, path.AppendElement(1));
                    }
                }
                else
                {
                    oDualLoop2.Add(brep2[0], path.AppendElement(1));
                    oDualLoop2ID.Add("J;" + triangleIndex.ToString() + ";" + neighbourTriIndex.ToString() + ";" + neighbour2TriIndex.ToString(), path.AppendElement(1));

                    oDualLoop2Curves.Add(right2, path.AppendElement(1));
                    oDualLoop2Curves.Add(left2, path.AppendElement(1));
                }
            }


        }

        // calculate normal on vertices and deal with boarder condition.
        private void calculateVertexNormals()
        {
            vertexNormals = iSpringMesh.ComputeVertexNormals();

            foreach (Edge edge in iSpringMesh.Edges)
            {
                if (edge.SecondTriangleIndex >= 0) continue;
                // Gene Added
                if (iSpringMesh.Vertices[edge.FirstVertexIndex].Position.Z > 0.05 &&
                     iSpringMesh.Vertices[edge.SecondVertexIndex].Position.Z > 0.05) continue;

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
        }

        private void calculateVertexCps()
        {
            for (int i = 0; i < iSpringMesh.Vertices.Count; i++)
            {
                Point3d vertexPosition = iSpringMesh.Vertices[i].Position;
                bottomCps.Add(vertexPosition - iThickness[i] * vertexNormals[i]);
                topCps.Add(vertexPosition + iThickness[i] * vertexNormals[i]);
            }
        }

        // caculate vertex value and use for opening
        private void calculateVerticesValues()
        {
            // set value
            foreach (Vertex v in iSpringMesh.Vertices)
            {
                // average method
                double distValue = 0;
                /*
                foreach (Point3d p in iAttractors)
                    distValue += ICD.Utils.Distance(v.Position, p);

                verticesValues.Add(distValue / iAttractors.Count);
                */
                // closest point method
                Point3d closestPt = Point3dList.ClosestPointInList(iAttractors, v.Position);
                distValue = ICD.Utils.Distance(v.Position, closestPt);
                verticesValues.Add(distValue);
            }
            // remap
            double max = verticesValues.Max();
            double min = verticesValues.Min();

            // map(value, low1, high1, low2, high2) = > low2 + (value - low1) * (high2 - low2) / (high1 - low1)
            for (int i = 0; i < verticesValues.Count; i++)
                verticesValues[i] = iTangentScaleMin + (verticesValues[i] - min) * (iTangentScaleMax - iTangentScaleMin) / (max - min);

        }

        private void triLoop()
        {
            for (int tri = 0; tri < iSpringMesh.Triangles.Count; tri++)
            {
                GH_Path path = new GH_Path(tri);
                //oInfo += "path = " + path.ToString() + "\n";

                Triangle triangle = iSpringMesh.Triangles[tri];

                if (!iVertexStatus[triangle.FirstVertexIndex] &&
                    !iVertexStatus[triangle.SecondVertexIndex] &&
                    !iVertexStatus[triangle.ThirdVertexIndex])
                {
                    // ids
                    oTriLoopID.Add("T;" + tri.ToString() + ";" + triangle.SecondAdjTriIndex.ToString() + ";" + triangle.ThirdAdjTriIndex.ToString(), path.AppendElement(0));
                    oTriLoopID.Add("T;" + tri.ToString() + ";" + triangle.ThirdAdjTriIndex.ToString() + ";" + triangle.FirstAdjTriIndex.ToString(), path.AppendElement(1));
                    oTriLoopID.Add("T;" + tri.ToString() + ";" + triangle.FirstAdjTriIndex.ToString() + ";" + triangle.SecondAdjTriIndex.ToString(), path.AppendElement(2));
                    
                    // stripes
                    createStripe(triangle.FirstVertexIndex, triangle.SecondVertexIndex, triangle.ThirdVertexIndex, path, 0);
                    createStripe(triangle.SecondVertexIndex, triangle.ThirdVertexIndex, triangle.FirstVertexIndex, path, 1);
                    createStripe(triangle.ThirdVertexIndex, triangle.FirstVertexIndex, triangle.SecondVertexIndex, path, 2);
                }
            }
        }

        private void createStripe(int firstVertexIndex, int secondVertexIndex, int thirdVertexIndex, GH_Path path, int item)
        {
            // first vertex in the stripe direction
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

            double scaler = verticesValues[firstVertexIndex];

            Point3d AAB = scaler * A + (1.0 - scaler) * AB;
            Point3d aab = scaler * a + (1.0 - scaler) * ab;

            Point3d AAC = scaler * A + (1.0 - scaler) * AC;
            Point3d aac = scaler * a + (1.0 - scaler) * ac;

            Curve profileCurve1 = Curve.CreateControlPointCurve(new List<Point3d>() { ab, aab, aA, AAB, AB });
            Curve profileCurve2 = Curve.CreateControlPointCurve(new List<Point3d>() { ac, aac, aA, AAC, AC });

            oTriLoopCurves.Add(profileCurve1, path.AppendElement(item));
            oTriLoopCurves.Add(profileCurve2, path.AppendElement(item));

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

                oTriLoop.Add(
                    Brep.CreateFromLoft(
                       new List<Curve>() { polyCurve1, polyCurve2 },
                       Point3d.Unset, Point3d.Unset,
                       LoftType.Normal,
                       false
                       )[0], path.AppendElement(item)
                   );
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
                else { oTriLoop.Add(brep, path.AppendElement(item)); return; }

                breps = brep.Trim(new Plane(A_, M, AC), documentTolerance);
                if (breps.Length > 0) brep = breps[0];
                else { oTriLoop.Add(brep, path.AppendElement(item)); return; }

                normal = Vector3d.CrossProduct(b - a, c - a);

                Point3d a_ = a - offsetAmount * normal;

                breps = brep.Trim(new Plane(a_, m, ab), documentTolerance);
                if (breps.Length > 0) brep = breps[0];
                else { oTriLoop.Add(brep, path.AppendElement(item)); return; }

                breps = brep.Trim(new Plane(a_, ac, m), documentTolerance);
                if (breps.Length > 0) brep = breps[0];
                else { oTriLoop.Add(brep, path.AppendElement(item)); return; }

                oTriLoop.Add(brep, path.AppendElement(item));
            }

            // close panel
            Point3d closeVertice = Point3dList.ClosestPointInList(iClosedPanelPts, iSpringMesh.Vertices[firstVertexIndex].Position);
            if (ICD.Utils.Distance(closeVertice, iSpringMesh.Vertices[firstVertexIndex].Position) < iClosePanelDist)
            {
                oClosedPanel.Add(
                    Brep.CreateFromLoft(
                       new List<Curve>() { new LineCurve(AB, A), new LineCurve(AC, A) },
                       Point3d.Unset, Point3d.Unset,
                       LoftType.Normal,
                       false
                       )[0], path.AppendElement(item)
                    );
                oClosedPanel.Add(
                    Brep.CreateFromLoft(
                       new List<Curve>() { new LineCurve(ab, a), new LineCurve(ac, a) },
                       Point3d.Unset, Point3d.Unset,
                       LoftType.Normal,
                       false
                       )[0], path.AppendElement(item)
                    );
            }

        }

        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                //You can add image files to your project resources and access them like this:
                // return Resources.IconForThisComponent;
                return Properties.Resources.Topo2;
            }
        }


        public override Guid ComponentGuid
        {
            get { return new Guid("{e52febc7-23b3-4589-ac59-210ec5c14f9a}"); }
        }
    }
}