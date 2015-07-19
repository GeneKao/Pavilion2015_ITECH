using System;
using System.Collections.Generic;
using System.Linq;

using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;
using Rhino.Geometry;
using Rhino.Collections;

using ICD;

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
        List<bool> iVertexPanel2 = null;
        double iPlanarOffsetScaleMin = double.NaN;
        double iPlanarOffsetScaleMax = double.NaN;
        double CurvePointinessMin = double.NaN;
        double CurvePointinessMax = double.NaN;

        double iOpeningWidthMin = double.NaN;
        double iOpeningWidthMax = double.NaN;

        GH_Structure<GH_Number> ManualValueTree = null;


        // output
        string oInfo = string.Empty;
        List<Vector3d> oDebugList = null;
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
        DataTree<Curve> oTriLoopPlanCrv = null;
        DataTree<Point3d> oTriLoopEffectorHoles = null;

        // internal usage
        List<Vector3d> vertexNormals;

        List<Point3d> topCps = null;
        List<Point3d> bottomCps = null;

        List<Point3d> topCenterPts = null;
        List<Point3d> bottomCenterPts = null;
        List<bool> topCenterPolygonOccupied = null;

        List<int[]> indexSortPolygon = null;

        List<double> verticesValues = null;  // General Distance Relation to Attractor(s) ( before remaping to Min/Max - Range )

        List<double> curveVerticesValues = null;  // relative offset factors for individual opnening sizes (remaped verticesValues to iTangentScaleMin - iTangentScaleMax - range)
        List<double> planarVerticesValues = null;  // offset distance in document unit for individual planar part sizes (remaped verticesValues to iPlanarOffsetScaleMin - iPlanarOffsetScaleMax - range)
        List<double> openingWidthVerticesValues = null;  // offset distance in document unit for individual opening widths
        List<double> pointinessValues = null;  // offset distance in document unit for individual opening widths
        List<int> ManualAdjustedVertexIndexies = null;

        double documentTolerance = DocumentTolerance();
        int curveDegree = 2;   // here change curve degree

        

        public GHC_DoubleLayer_Topo_II_DataTree_()
            : base("Double Layer Topo - II DT", "Double Layer Topo - II DT",
                "Double Layer Topology II DT - get plate and SpringMesh convert double layers - output data tree",
                "Pavillion 2015", "Double Layer")
        {
        }


        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Spring Mesh", "Spring Mesh", "Spring Mesh", GH_ParamAccess.item);
            pManager.AddCurveParameter("PolyLine", "PolyLine", "PolyLine from Plate", GH_ParamAccess.list);
            pManager.AddNumberParameter("Thickness", "Thickness", "Thickness of Component", GH_ParamAccess.list);
            pManager.AddBooleanParameter("Vertex status", "Vertex status", "Vertex status", GH_ParamAccess.list);
            pManager.AddIntegerParameter("PolyLine ID", "PolyLine ID", "PolyLine ID", GH_ParamAccess.list);
            pManager.AddBooleanParameter("Allow PolySrf", "Allow PolySrf", "Allow PolySrf", GH_ParamAccess.item, true);
            pManager.AddNumberParameter("TangentScale Min", "TangentScale Min", "TangentScale Min [relative]", GH_ParamAccess.item, 0.2);
            pManager.AddNumberParameter("TangentScale Max", "TangentScale Max", "TangentScale Max [relative]", GH_ParamAccess.item, 0.8);
            pManager.AddNumberParameter("Curve Pointiness Min", "Pointiness Min", "Pointiness of the bended Surfaces [relative]", GH_ParamAccess.item, 1.0);
            pManager.AddNumberParameter("Curve Pointiness Max", "Pointiness Max", "Pointiness of the bended Surfaces [relative]", GH_ParamAccess.item, 1.0);
            pManager.AddNumberParameter("PlanarOffset Min", "PlanarOffset Min", "Controlls minimal offset of planar parts [in doc. units]", GH_ParamAccess.item, 0.2);
            pManager.AddNumberParameter("PlanarOffset Max", "PlanarOffset Max", "Controlls maximal offset of planar parts [in doc. units]", GH_ParamAccess.item, 0.8);
            pManager.AddPointParameter("Attractors", "Attractors", "Attractors", GH_ParamAccess.list, new Point3d(10000, 10000, 10000));
            pManager.AddBooleanParameter("Vertex Panel2", "Vertex Panel2", "Vertex Panel2", GH_ParamAccess.list);
            pManager.AddNumberParameter("Opening Width Min", "OpeningWidth Min", "TangentScale Min [in doc. units]", GH_ParamAccess.item, 0.2);
            pManager.AddNumberParameter("Opening Width Max", "OpeningWidth Max", "TangentScale Max [in doc. units]", GH_ParamAccess.item, 0.8);

            pManager.AddNumberParameter("Manual Adjustments", "Manual Adjustments", "Tree of manual adjusted input data for selected verteices", GH_ParamAccess.tree);

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
            pManager.AddGenericParameter("Planar Parts of TriLoop", "TriLoopPlanarParts", "Planar Curves for Planar Parts of TriLoops", GH_ParamAccess.tree);       //12
            pManager.AddPointParameter("Holes for Effector", "EffectorHolePts", "Points for EffectorHoles", GH_ParamAccess.tree);       //13

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
            iVertexPanel2 = new List<bool>();
            //ManualValueTree = new GH_Structure<GH_Number>();


            // output
            oInfo = string.Empty;
            oDebugList = new List<Vector3d>();
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
            oTriLoopPlanCrv = new DataTree<Curve>();
            oTriLoopEffectorHoles = new DataTree<Point3d>();


            // internal use data
            topCenterPts = new List<Point3d>();
            bottomCenterPts = new List<Point3d>();
            topCenterPolygonOccupied = new List<bool>(); // not be used at this moment 

            bottomCps = new List<Point3d>();
            topCps = new List<Point3d>();

            indexSortPolygon = new List<int[]>();   // array is reference type in csharp

            verticesValues = new List<double>();
            curveVerticesValues = new List<double>();
            planarVerticesValues = new List<double>();
            openingWidthVerticesValues = new List<double>();
            pointinessValues = new List<double>();
            ManualAdjustedVertexIndexies = new List<int>();


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
            DA.GetDataList<bool>("Vertex Panel2", iVertexPanel2);
            DA.GetData<double>("PlanarOffset Min", ref iPlanarOffsetScaleMin);
            DA.GetData<double>("PlanarOffset Max", ref iPlanarOffsetScaleMax);
            DA.GetData<double>("Curve Pointiness Min", ref CurvePointinessMin);
            DA.GetData<double>("Curve Pointiness Max", ref CurvePointinessMax);

            DA.GetData<double>("Opening Width Min", ref iOpeningWidthMin);
            DA.GetData<double>("Opening Width Max", ref iOpeningWidthMax);
            DA.GetDataTree<GH_Number>("Manual Adjustments", out ManualValueTree);
            //------------------------------------------------------------

            // get all Vertex indexies that are manual adjusted
            for (int i = 0; i < ManualValueTree.Branches.Count; i++)
            { ManualAdjustedVertexIndexies.Add(ManualValueTree.Paths[i].Indices[0]); }


            storePlatesTPI();

            calculateVertexNormals();
            calculateVertexCps();
            calculateVerticesValues();

            triLoop();

            half_dualLoop();

            oDebugList = vertexNormals;

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
            DA.SetDataTree(12, oTriLoopPlanCrv);
            DA.SetDataTree(13, oTriLoopEffectorHoles);

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
            /*
            foreach (Point3d p in topCenterPts)
                oDebugList.Add(p);
            foreach (Point3d p in bottomCenterPts)
                oDebugList.Add(p);*/

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

                    if ((iVertexStatus[vertex01] == false && iVertexStatus[vertex02] == false) &&
                        (iVertexStatus[edgeLocal.FirstAdjacentVertexIndex] == true && iVertexStatus[edgeLocal.SecondAdjacentVertexIndex] == false))
                        half_dualLoop_type_2(i, tri01, tri02);

                    if ((iVertexStatus[vertex01] == false && iVertexStatus[vertex02] == false) &&
                        (iVertexStatus[edgeLocal.FirstAdjacentVertexIndex] == false && iVertexStatus[edgeLocal.SecondAdjacentVertexIndex] == true))
                        half_dualLoop_type_2(i, tri02, tri01);

                    if ((iVertexStatus[vertex01] == false && iVertexStatus[vertex02] == false) &&
                        (iVertexStatus[edgeLocal.FirstAdjacentVertexIndex] == true && iVertexStatus[edgeLocal.SecondAdjacentVertexIndex] == true))
                    {
                        half_dualLoop_type_2(i, tri01, tri02);
                        half_dualLoop_type_2(i, tri02, tri01);
                    }
                }
                else
                {
                    if ((iVertexStatus[vertex01] == false && iVertexStatus[vertex02] == false) &&
                        (iVertexStatus[edgeLocal.FirstAdjacentVertexIndex] == true))
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
            else { oInfo += "bug found..."; }


            // layer up
            Point3d rightUp01 = topCenterPts[rightTriIndex];
            Point3d leftUp01 = topCenterPts[leftTriIndex];

            Point3d up03 = topCps[vertex]; //reference point not for construct surfaces

            double planarOffset = planarVerticesValues[vertex];
            double scaler = curveVerticesValues[vertex];
            double openingScalar = openingWidthVerticesValues[vertex];
            double curvePointiness = pointinessValues[vertex];

            //Calculate offset planar points
            Vector3d v_oRightUp01 = up03 - rightUp01;
            Vector3d v_oLeftUp01 = up03 - leftUp01;
            v_oRightUp01.Unitize();
            v_oLeftUp01.Unitize();

            Point3d oRightUp01 = rightUp01 + (v_oRightUp01 * planarOffset);
            Point3d oLeftUp01 = leftUp01 + (v_oLeftUp01 * planarOffset);

            Point3d rightUp02 = (1 - scaler) * oRightUp01 + (scaler) * up03;
            Point3d leftUp02 = (1 - scaler) * oLeftUp01 + (scaler) * up03;

            // normalise hohles
            Point3d openingRightUp02 = up03 + (openingScalar * v_oRightUp01 * -1);
            Point3d openingLeftUp02 = up03 + (openingScalar * v_oLeftUp01 * -1);

            if (new Vector3d(rightUp02 - up03).Length > openingScalar)
                rightUp02 = openingRightUp02;

            if (new Vector3d(leftUp02 - up03).Length > openingScalar)
                leftUp02 = openingLeftUp02;
            // end normalise hohles

            // layer down
            Point3d rightDown01 = bottomCenterPts[rightTriIndex];
            Point3d leftDown01 = bottomCenterPts[leftTriIndex];

            Point3d down03 = bottomCps[vertex]; //reference point not for construct surfaces

            //Calculate offset planar points
            Vector3d v_oRightDownp01 = down03 - rightDown01;
            Vector3d v_oLeftDown01 = down03 - leftDown01;
            v_oRightDownp01.Unitize();
            v_oLeftDown01.Unitize();

            Point3d oRightDown01 = rightDown01 + (v_oRightDownp01 * planarOffset);
            Point3d oLeftDown01 = leftDown01 + (v_oLeftDown01 * planarOffset);

            Point3d rightDown02 = (1 - scaler) * oRightDown01 + (scaler) * down03;
            Point3d leftDown02 = (1 - scaler) * oLeftDown01 + (scaler) * down03;

            // normalise hohles
            Point3d openingRightDown02 = down03 + (openingScalar * v_oRightDownp01 * -1);
            Point3d openingLeftDown02 = down03 + (openingScalar * v_oLeftDown01 * -1);

            if (new Vector3d(rightDown02 - down03).Length > openingScalar)
                rightDown02 = openingRightDown02;

            if (new Vector3d(leftDown02 - down03).Length > openingScalar)
                leftDown02 = openingLeftDown02;
            // end normalise hohles

            //Pointiness
            Point3d maxVertexPt = iSpringMesh.Vertices[vertex].Position;
            Point3d min_rightVertexPt = 0.5 * (rightUp02 + rightDown02);
            Point3d min_leftVertexPt = 0.5 * (leftUp02 + leftDown02);
            Point3d rightVertexPt = (curvePointiness * maxVertexPt) + ((1 - curvePointiness) * min_rightVertexPt);
            Point3d leftVertexPt = (curvePointiness * maxVertexPt) + ((1 - curvePointiness) * min_leftVertexPt);


            // Right Curves (Curved Curve, Planar Part Curve at Top, Planar Part Curve at Bottom, Joined Curve)
            Curve rightLoop = Curve.CreateControlPointCurve(
                new List<Point3d>() { oRightDown01, rightDown02, rightVertexPt, rightUp02, oRightUp01 }, curveDegree);

            Curve rightPlanarUp = Curve.CreateControlPointCurve(
                new List<Point3d>() { oRightUp01, rightUp01 }, 1);

            Curve rightPlanarDown = Curve.CreateControlPointCurve(
                new List<Point3d>() { rightDown01, oRightDown01 }, 1);
            /*
            PolyCurve right = new PolyCurve();
            right.Append(rightPlanarUp);
            right.Append(rightLoop);
            right.Append(rightPlanarDown);*/

            Curve right = Curve.JoinCurves(new List<Curve>() { rightPlanarDown, rightLoop, rightPlanarUp }, documentTolerance, true)[0];

            if (iPolySrf)
            {
                right = rightLoop;
            }

            // Left Curves (Curved Curve, Planar Part Curve at Top, Planar Part Curve at Bottom, Joined Curve)
            Curve leftLoop = Curve.CreateControlPointCurve(
                new List<Point3d>() { oLeftDown01, leftDown02, leftVertexPt, leftUp02, oLeftUp01 }, curveDegree);

            Curve leftPlanarUp = Curve.CreateControlPointCurve(
                new List<Point3d>() { oLeftUp01, leftUp01 }, 1);

            Curve leftPlanarDown = Curve.CreateControlPointCurve(
                new List<Point3d>() { leftDown01, oLeftDown01 }, 1);

            /*
            PolyCurve left = new PolyCurve();
            left.Append(leftPlanarUp);
            left.Append(leftLoop);
            left.Append(leftPlanarDown);*/

            Curve left = Curve.JoinCurves(new List<Curve>() { leftPlanarDown, leftLoop, leftPlanarUp }, documentTolerance, true)[0];

            if (iPolySrf)
            {
                left = leftLoop;
            }

            // compare their polygon index to sort the lofting sequence 
            #region compair polygon index

            // dataTree 
            int plateID = iPolyLineID[occupyiedVertex];
            plateID = -10 - plateID;     // prevent -1 situation

            GH_Path path = new GH_Path(plateID);

            int[] rightIndex = indexSortPolygon[rightTriIndex];
            int[] leftIndex = indexSortPolygon[leftTriIndex];

            for (int r = 0; r < 3; r++)
                for (int l = 0; l < 3; l++)
                    if (rightIndex[r] == leftIndex[l] && rightIndex[r] != -1 && leftIndex[l] != -1)
                        if ((rightIndex[r + 3] == 0 && leftIndex[l + 3] != 1) ||
                            (leftIndex[l + 3] == 0 && rightIndex[r + 3] != 1) && (rightIndex[r + 3] < leftIndex[l + 3]))
                        {
                            Brep[] brep = Brep.CreateFromLoft(
                                new List<Curve>() { left, right },
                                Point3d.Unset, Point3d.Unset,
                                LoftType.Normal,
                                false
                                );

                            #region allow poly surface
                            if (iPolySrf)
                            {
                                Brep[] brepPlanarUp = Brep.CreateFromLoft(
                                    new List<Curve>() { leftPlanarUp, rightPlanarUp },
                                    Point3d.Unset, Point3d.Unset,
                                    LoftType.Normal, false
                                    );

                                Brep[] brepPlanarDown = Brep.CreateFromLoft(
                                    new List<Curve>() { leftPlanarDown, rightPlanarDown },
                                    Point3d.Unset, Point3d.Unset,
                                    LoftType.Normal, false
                                    );

                                brep = Brep.JoinBreps(new List<Brep>() { brepPlanarDown[0], brep[0], brepPlanarUp[0] }, documentTolerance);
                            }
                            #endregion allow poly surface

                            if (brep.Length > 0)
                            {
                                oDualLoop1.Add(brep[0], path.AppendElement(leftIndex[l + 3]));
                                oDualLoop1ID.Add("H;" + rightTriIndex.ToString() + "-" + leftTriIndex.ToString() + ";" + plateID.ToString(), path.AppendElement(leftIndex[l + 3]));

                                oDualLoop1Curves.Add(left, path.AppendElement(leftIndex[l + 3]));
                                oDualLoop1Curves.Add(right, path.AppendElement(leftIndex[l + 3]));
                            }
                        }
                        else if ((rightIndex[r + 3] == 0 && leftIndex[l + 3] != 1) ||
                                 (leftIndex[l + 3] == 0 && rightIndex[r + 3] != 1) && (rightIndex[r + 3] >= leftIndex[l + 3]))
                        {
                            Brep[] brep = Brep.CreateFromLoft(
                                new List<Curve>() { right, left },
                                Point3d.Unset, Point3d.Unset,
                                LoftType.Normal,
                                false
                                );

                            #region allow poly surface
                            if (iPolySrf)
                            {
                                Brep[] brepPlanarUp = Brep.CreateFromLoft(
                                    new List<Curve>() { rightPlanarUp, leftPlanarUp },
                                    Point3d.Unset, Point3d.Unset,
                                    LoftType.Normal, false
                                    );

                                Brep[] brepPlanarDown = Brep.CreateFromLoft(
                                    new List<Curve>() { rightPlanarDown, leftPlanarDown },
                                    Point3d.Unset, Point3d.Unset,
                                    LoftType.Normal, false
                                    );

                                brep = Brep.JoinBreps(new List<Brep>() { brepPlanarDown[0], brep[0], brepPlanarUp[0] }, documentTolerance);
                            }
                            #endregion allow poly surface

                            if (brep.Length > 0)
                            {
                                oDualLoop1.Add(brep[0], path.AppendElement(rightIndex[r + 3]));
                                oDualLoop1ID.Add("H;" + rightTriIndex.ToString() + "-" + leftTriIndex.ToString() + ";" + plateID.ToString(), path.AppendElement(rightIndex[r + 3]));

                                oDualLoop1Curves.Add(right, path.AppendElement(rightIndex[r + 3]));
                                oDualLoop1Curves.Add(left, path.AppendElement(rightIndex[r + 3]));
                            }
                        }
                        else if (rightIndex[r + 3] >= leftIndex[l + 3])
                        {
                            Brep[] brep = Brep.CreateFromLoft(
                                new List<Curve>() { left, right },
                                Point3d.Unset, Point3d.Unset,
                                LoftType.Normal,
                                false
                                );

                            #region allow poly surface
                            if (iPolySrf)
                            {
                                Brep[] brepPlanarUp = Brep.CreateFromLoft(
                                    new List<Curve>() { leftPlanarUp, rightPlanarUp },
                                    Point3d.Unset, Point3d.Unset,
                                    LoftType.Normal, false
                                    );

                                Brep[] brepPlanarDown = Brep.CreateFromLoft(
                                    new List<Curve>() { leftPlanarDown, rightPlanarDown },
                                    Point3d.Unset, Point3d.Unset,
                                    LoftType.Normal, false
                                    );

                                brep = Brep.JoinBreps(new List<Brep>() { brepPlanarDown[0], brep[0], brepPlanarUp[0] }, documentTolerance);
                            }
                            #endregion allow poly surface

                            if (brep.Length > 0)
                            {
                                oDualLoop1.Add(brep[0], path.AppendElement(leftIndex[l + 3]));
                                oDualLoop1ID.Add("H;" + rightTriIndex.ToString() + "-" + leftTriIndex.ToString() + ";" + plateID.ToString(), path.AppendElement(leftIndex[l + 3]));

                                oDualLoop1Curves.Add(left, path.AppendElement(leftIndex[l + 3]));
                                oDualLoop1Curves.Add(right, path.AppendElement(leftIndex[l + 3]));
                            }
                        }
                        else if (rightIndex[r + 3] < leftIndex[l + 3])
                        {
                            Brep[] brep = Brep.CreateFromLoft(
                                new List<Curve>() { right, left },
                                Point3d.Unset, Point3d.Unset,
                                LoftType.Normal,
                                false
                                );

                            #region allow poly surface
                            if (iPolySrf)
                            {
                                Brep[] brepPlanarUp = Brep.CreateFromLoft(
                                    new List<Curve>() { rightPlanarUp, leftPlanarUp },
                                    Point3d.Unset, Point3d.Unset,
                                    LoftType.Normal, false
                                    );

                                Brep[] brepPlanarDown = Brep.CreateFromLoft(
                                    new List<Curve>() { rightPlanarDown, leftPlanarDown },
                                    Point3d.Unset, Point3d.Unset,
                                    LoftType.Normal, false
                                    );

                                brep = Brep.JoinBreps(new List<Brep>() { brepPlanarDown[0], brep[0], brepPlanarUp[0] }, documentTolerance);
                            }
                            #endregion allow poly surface

                            if (brep.Length > 0)
                            {
                                oDualLoop1.Add(brep[0], path.AppendElement(rightIndex[r + 3]));
                                oDualLoop1ID.Add("H;" + rightTriIndex.ToString() + "-" + leftTriIndex.ToString() + ";" + plateID.ToString(), path.AppendElement(rightIndex[r + 3]));

                                oDualLoop1Curves.Add(right, path.AppendElement(rightIndex[r + 3]));
                                oDualLoop1Curves.Add(left, path.AppendElement(rightIndex[r + 3]));
                            }
                        }
            #endregion compair polygon index


        }

        // dual loop not directly connected to plates
        private void half_dualLoop_type_2(int edgeIndex, int triangleIndex, int neighbourTriIndex)
        {
            GH_Path path = new GH_Path(triangleIndex);

            Edge edge = iSpringMesh.Edges[edgeIndex];

            // added for recognizing where is the plate - plate id
            int occupyiedVertex = edge.FirstAdjacentVertexIndex;
            if (iVertexStatus[edge.FirstAdjacentVertexIndex]) { occupyiedVertex = edge.FirstAdjacentVertexIndex; }
            if (edge.SecondAdjacentVertexIndex > 0 && iVertexStatus[edge.SecondAdjacentVertexIndex]) { occupyiedVertex = edge.SecondAdjacentVertexIndex; }
            int plateID = iPolyLineID[occupyiedVertex];
            plateID = -10 - plateID;

            // vertex point
            int vertex1 = edge.FirstVertexIndex;
            int vertex2 = edge.SecondVertexIndex;

            Point3d vertexPt1 = iSpringMesh.Vertices[vertex1].Position;
            Point3d vertexPt2 = iSpringMesh.Vertices[vertex2].Position;

            double scaler1 = curveVerticesValues[vertex1];
            double scaler2 = curveVerticesValues[vertex2];

            double planarOffset1 = planarVerticesValues[vertex1];
            double planarOffset2 = planarVerticesValues[vertex2];

            double openingScalar1 = openingWidthVerticesValues[vertex1];
            double openingScalar2 = openingWidthVerticesValues[vertex2];

            double curvePointiness1 = pointinessValues[vertex1];
            double curvePointiness2 = pointinessValues[vertex2];


            //---- TriLoop Side UP -----------------------------------------------
            #region TriLoop Side UP

            Point3d vertexPtUp1 = topCps[vertex1];
            Point3d vertexPtUp2 = topCps[vertex2];
            Point3d vertexPtUpM = 0.5 * (vertexPtUp1 + vertexPtUp2);

            Vector3d v_vertexPtUp1 = vertexPtUp1 - vertexPtUpM;
            Vector3d v_vertexPtUp2 = vertexPtUp2 - vertexPtUpM;
            v_vertexPtUp1.Unitize();
            v_vertexPtUp2.Unitize();

            Point3d oVertexPtUpM1 = vertexPtUpM + (v_vertexPtUp1 * planarOffset1);
            Point3d oVertexPtUpM2 = vertexPtUpM + (v_vertexPtUp2 * planarOffset2);
            
            Point3d ctPtUp1 = (scaler1) * vertexPtUp1 + (1 - scaler1) * oVertexPtUpM1;
            Point3d ctPtUp2 = (scaler2) * vertexPtUp2 + (1 - scaler2) * oVertexPtUpM2;
            
            // normalise hohles
            Point3d openingCtPtUp1 = vertexPtUp1 + (openingScalar1 * v_vertexPtUp1 * -1);
            Point3d openingCtPtUp2 = vertexPtUp2 + (openingScalar2 * v_vertexPtUp2 * -1);

            if (new Vector3d(ctPtUp1 - vertexPtUp1).Length > openingScalar1)
                ctPtUp1 = openingCtPtUp1;

            if (new Vector3d(ctPtUp2 - vertexPtUp2).Length > openingScalar2)
                ctPtUp2 = openingCtPtUp2;
            // end normalise hohles

            #endregion TriLoop Side UP

            //---- TriLoop Side Down -----------------------------------------------
            #region TriLoop Side Down
            Point3d vertexPtDown1 = bottomCps[vertex1];
            Point3d vertexPtDown2 = bottomCps[vertex2];
            Point3d vertexPtDownM = 0.5 * (vertexPtDown1 + vertexPtDown2);

            Vector3d v_vertexPtDown1 = vertexPtDown1 - vertexPtDownM;
            Vector3d v_vertexPtDown2 = vertexPtDown2 - vertexPtDownM;
            v_vertexPtDown1.Unitize();
            v_vertexPtDown2.Unitize();

            Point3d oVertexPtDownM1 = vertexPtDownM + (v_vertexPtDown1 * planarOffset1);
            Point3d oVertexPtDownM2 = vertexPtDownM + (v_vertexPtDown2 * planarOffset2);

            Point3d ctPtDown1 = (scaler1) * vertexPtDown1 + (1 - scaler1) * oVertexPtDownM1;
            Point3d ctPtDown2 = (scaler2) * vertexPtDown2 + (1 - scaler2) * oVertexPtDownM2;

            // normalise hohles
            Point3d openingCtPtDown1 = vertexPtDown1 + (openingScalar1 * v_vertexPtDown1 * -1);
            Point3d openingCtPtDown2 = vertexPtDown2 + (openingScalar2 * v_vertexPtDown2 * -1);

            if (new Vector3d(ctPtDown1 - vertexPtDown1).Length > openingScalar1)
                ctPtDown1 = openingCtPtDown1;

            if (new Vector3d(ctPtDown2 - vertexPtDown2).Length > openingScalar2)
                ctPtDown2 = openingCtPtDown2;
            // end normalise hohles

            #endregion TriLoop Side Down

            //---- Plate Side UP -------------------------------------------------
            #region Plate Side UP
            Point3d upTPI = topCenterPts[triangleIndex];

            Vector3d v_vertexPtUp1_upTPI = vertexPtUp1 - upTPI;
            Vector3d v_vertexPtUp2_upTPI = vertexPtUp2 - upTPI;
            v_vertexPtUp1_upTPI.Unitize();
            v_vertexPtUp2_upTPI.Unitize();

            Point3d oUpTPI01 = upTPI + (v_vertexPtUp1_upTPI * planarOffset1);
            Point3d oUpTPI02 = upTPI + (v_vertexPtUp2_upTPI * planarOffset2);

            Point3d ctUpTPI1 = (scaler1) * vertexPtUp1 + (1 - scaler1) * oUpTPI01;
            Point3d ctUpTPI2 = (scaler2) * vertexPtUp2 + (1 - scaler2) * oUpTPI02;

            // normalise hohles
            Point3d openingCtUpTPI1 = vertexPtUp1 + (openingScalar1 * v_vertexPtUp1_upTPI * -1);
            Point3d openingCtUpTPI2 = vertexPtUp2 + (openingScalar2 * v_vertexPtUp2_upTPI * -1);

            if (new Vector3d(ctUpTPI1 - vertexPtUp1).Length > openingScalar1)
                ctUpTPI1 = openingCtUpTPI1;

            if (new Vector3d(ctUpTPI2 - vertexPtUp2).Length > openingScalar2)
                ctUpTPI2 = openingCtUpTPI2;
            // end normalise hohles

            #endregion Plate Side UP

            //---- Plate Side DOWN -----------------------------------------------
            #region Plate Side DOWN
            Point3d downTPI = bottomCenterPts[triangleIndex];

            Vector3d v_vertexPtUp1_downTPI = vertexPtDown1 - downTPI;
            Vector3d v_vertexPtUp2_downTPI = vertexPtDown2 - downTPI;
            v_vertexPtUp1_downTPI.Unitize();
            v_vertexPtUp2_downTPI.Unitize();

            Point3d oDownTPI01 = downTPI + (v_vertexPtUp1_downTPI * planarOffset1);
            Point3d oDownTPI02 = downTPI + (v_vertexPtUp2_downTPI * planarOffset2);

            Point3d ctDownTPI1 = (scaler1) * vertexPtDown1 + (1 - scaler1) * oDownTPI01;
            Point3d ctDownTPI2 = (scaler2) * vertexPtDown2 + (1 - scaler2) * oDownTPI02;

            // normalise hohles
            Point3d openingCtDownTPI1 = vertexPtDown1 + (openingScalar1 * v_vertexPtUp1_downTPI * -1);
            Point3d openingCtDownTPI2 = vertexPtDown2 + (openingScalar2 * v_vertexPtUp2_downTPI * -1);

            if (new Vector3d(ctDownTPI1 - vertexPtDown1).Length > openingScalar1)
                ctDownTPI1 = openingCtDownTPI1;

            if (new Vector3d(ctDownTPI2 - vertexPtDown2).Length > openingScalar2)
                ctDownTPI2 = openingCtDownTPI2;
            // end normalise hohles

            #endregion Plate Side DOWN

            //---- Curves Right 1 -----------------------------------------------------------
            #region Curves Right 1

            // Pointiness of the Surface
            Point3d min_vertexPtRight1 = 0.5 * (ctUpTPI1 + ctDownTPI1);
            Point3d vertexPtRight1 = (curvePointiness1 * vertexPt1) + ((1 - curvePointiness1) * min_vertexPtRight1);


            Curve right1Loop = Curve.CreateControlPointCurve(
                new List<Point3d>() { oDownTPI01, ctDownTPI1, vertexPtRight1, ctUpTPI1, oUpTPI01 }, curveDegree);

            Curve right1PlanarUp = Curve.CreateControlPointCurve(
                new List<Point3d>() { oUpTPI01, upTPI }, 1);

            Curve right1PlanarDown = Curve.CreateControlPointCurve(
                new List<Point3d>() { downTPI, oDownTPI01 }, 1);
            /*
            PolyCurve right1 = new PolyCurve();
            right1.Append(right1PlanarUp);
            right1.Append(right1Loop);
            right1.Append(right1PlanarDown);*/
            Curve right1 = Curve.JoinCurves(new List<Curve>() { right1PlanarDown, right1Loop, right1PlanarUp }, documentTolerance, true)[0];

            if (iPolySrf)
            {
                right1 = right1Loop;
            }

            #endregion Curves Right 1

            //---- Curves Left 1 -----------------------------------------------------------
            #region Curves Left 1

            // Pointiness of the Surface
            Point3d min_vertexPtLeft1 = 0.5 * (ctPtUp1 + ctPtDown1);
            Point3d vertexPtLeft1 = (curvePointiness1 * vertexPt1) + ((1 - curvePointiness1) * min_vertexPtLeft1);

            Curve left1Loop = Curve.CreateControlPointCurve(
                new List<Point3d>() { oVertexPtDownM1, ctPtDown1, vertexPtLeft1, ctPtUp1, oVertexPtUpM1 }, curveDegree);

            Curve left1PlanarUp = Curve.CreateControlPointCurve(
                new List<Point3d>() { oVertexPtUpM1, vertexPtUpM }, 1);

            Curve left1PlanarDown = Curve.CreateControlPointCurve(
                new List<Point3d>() { vertexPtDownM, oVertexPtDownM1 }, 1);

            /*
            PolyCurve left1 = new PolyCurve();
            left1.Append(left1PlanarUp);
            left1.Append(left1Loop);
            left1.Append(left1PlanarDown);*/
            Curve left1 = Curve.JoinCurves(new List<Curve>() { left1PlanarDown, left1Loop, left1PlanarUp }, documentTolerance, true)[0];

            if (iPolySrf)
            {
                left1 = left1Loop;
            }

            #endregion Curves Left 1


            //---- Curves Right 2 -----------------------------------------------------------
            #region Curves Right 2

            // Pointiness of the Surface
            Point3d min_vertexPtRight2 = 0.5 * (ctUpTPI2 + ctDownTPI2);
            Point3d vertexPtRight2 = (curvePointiness2 * vertexPt2) + ((1 - curvePointiness2) * min_vertexPtRight2);

            Curve right2Loop = Curve.CreateControlPointCurve(
                new List<Point3d>() { oDownTPI02, ctDownTPI2, vertexPtRight2, ctUpTPI2, oUpTPI02 }, curveDegree);

            Curve right2PlanarUp = Curve.CreateControlPointCurve(
                new List<Point3d>() { oUpTPI02, upTPI }, 1);

            Curve right2PlanarDown = Curve.CreateControlPointCurve(
                new List<Point3d>() { downTPI, oDownTPI02 }, 1);

            /*
            PolyCurve right2 = new PolyCurve();
            right2.Append(right2PlanarUp);
            right2.Append(right2Loop);
            right2.Append(right2PlanarDown);*/
            Curve right2 = Curve.JoinCurves(new List<Curve>() { right2PlanarDown, right2Loop, right2PlanarUp }, documentTolerance, true)[0];

            if (iPolySrf)
            {
                right2 = right2Loop;
            }

            #endregion Curves Right 2

            //---- Curves Left 2 -----------------------------------------------------------
            #region Curves Left 2

            // Pointiness of the Surface
            Point3d min_vertexPtLeft2 = 0.5 * (ctPtUp2 + ctPtDown2);
            Point3d vertexPtLeft2 = (curvePointiness2 * vertexPt2) + ((1 - curvePointiness2) * min_vertexPtLeft2);

            Curve left2Loop = Curve.CreateControlPointCurve(
                new List<Point3d>() { oVertexPtDownM2, ctPtDown2, vertexPtLeft2, ctPtUp2, oVertexPtUpM2 }, curveDegree);

            Curve left2PlanarUp = Curve.CreateControlPointCurve(
                new List<Point3d>() { oVertexPtUpM2, vertexPtUpM }, 1);

            Curve left2PlanarDown = Curve.CreateControlPointCurve(
                new List<Point3d>() { vertexPtDownM, oVertexPtDownM2 }, 1);

            /*
            PolyCurve left2 = new PolyCurve();
            left2.Append(left2PlanarUp);
            left2.Append(left2Loop);
            left2.Append(left2PlanarDown);*/
            Curve left2 = Curve.JoinCurves(new List<Curve>() { left2PlanarDown, left2Loop, left2PlanarUp }, documentTolerance, true)[0];

            if (iPolySrf)
            {
                left2 = left2Loop;
            }

            #endregion Curves Left 2


            // sorting sequence 

            Brep[] brep1 = Brep.CreateFromLoft(
                       new List<Curve>() { left1, right1 },
                       Point3d.Unset, Point3d.Unset,
                       LoftType.Normal,
                       false
                       );

            #region allow poly surface
            if (iPolySrf)
            {
                Brep[] brep1PlanarUp = Brep.CreateFromLoft(
                    new List<Curve>() { left1PlanarUp, right1PlanarUp },
                    Point3d.Unset, Point3d.Unset,
                    LoftType.Normal, false
                    );

                Brep[] brep1PlanarDown = Brep.CreateFromLoft(
                    new List<Curve>() { left1PlanarDown, right1PlanarDown },
                    Point3d.Unset, Point3d.Unset,
                    LoftType.Normal, false
                    );

                brep1 = Brep.JoinBreps(new List<Brep>(){brep1PlanarDown[0], brep1[0], brep1PlanarUp[0]}, documentTolerance);
            }
            #endregion allow poly surface


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

                if (brepF.NormalAt(0, 0).Z < 0)  // conditional flipping surface, according face normal. special case for our pavilion mesh
                {
                    Brep[] brepN1 = Brep.CreateFromLoft(
                       new List<Curve>() { right1, left1 },
                       Point3d.Unset, Point3d.Unset,
                       LoftType.Normal,
                       false
                       );

                    #region allow poly surface
                    if (iPolySrf)
                    {
                        Brep[] brepN1PlanarUp = Brep.CreateFromLoft(
                            new List<Curve>() { right1PlanarUp, left1PlanarUp },
                            Point3d.Unset, Point3d.Unset,
                            LoftType.Normal, false
                            );

                        Brep[] brepN1PlanarDown = Brep.CreateFromLoft(
                            new List<Curve>() { right1PlanarDown, left1PlanarDown },
                            Point3d.Unset, Point3d.Unset,
                            LoftType.Normal, false
                            );

                        brepN1 = Brep.JoinBreps(new List<Brep>() { brepN1PlanarDown[0], brepN1[0], brepN1PlanarUp[0] }, documentTolerance);
                    }
                    #endregion allow poly surface

                    if (brep1.Length > 0)
                    {
                        oDualLoop2.Add(brepN1[0], path.AppendElement(0));
                        oDualLoop2ID.Add("J;" + triangleIndex.ToString() + ";" + neighbourTriIndex.ToString() + ";" + neighbour2TriIndex.ToString() + ";" + plateID.ToString(), path.AppendElement(0));

                        oDualLoop2Curves.Add(right1, path.AppendElement(0));
                        oDualLoop2Curves.Add(left1, path.AppendElement(0));
                    }
                }
                else
                {
                    oDualLoop2.Add(brep1[0], path.AppendElement(0));
                    oDualLoop2ID.Add("J;" + triangleIndex.ToString() + ";" + neighbourTriIndex.ToString() + ";" + neighbour2TriIndex.ToString() + ";" + plateID.ToString(), path.AppendElement(0));

                    oDualLoop2Curves.Add(left1, path.AppendElement(0));
                    oDualLoop2Curves.Add(right1, path.AppendElement(0));
                }
            }


            Brep[] brep2 = Brep.CreateFromLoft(
                       new List<Curve>() { left2, right2 },
                       Point3d.Unset, Point3d.Unset,
                       LoftType.Normal,
                       false
                       );

            #region allow poly surface
            if (iPolySrf)
            {
                Brep[] brep2PlanarUp = Brep.CreateFromLoft(
                    new List<Curve>() { left2PlanarUp, right2PlanarUp },
                    Point3d.Unset, Point3d.Unset,
                    LoftType.Normal, false
                    );

                Brep[] brep2PlanarDown = Brep.CreateFromLoft(
                    new List<Curve>() { left2PlanarDown, right2PlanarDown },
                    Point3d.Unset, Point3d.Unset,
                    LoftType.Normal, false
                    );

                brep2 = Brep.JoinBreps(new List<Brep>() { brep2PlanarDown[0], brep2[0], brep2PlanarUp[0] }, documentTolerance);
            }
            #endregion allow poly surface

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

                if (brepF.NormalAt(0, 0).Z < 0)
                {
                    Brep[] brepN2 = Brep.CreateFromLoft(
                       new List<Curve>() { right2, left2 },
                       Point3d.Unset, Point3d.Unset,
                       LoftType.Normal,
                       false
                       );

                    #region allow poly surface
                    if (iPolySrf)
                    {
                        Brep[] brepN2PlanarUp = Brep.CreateFromLoft(
                            new List<Curve>() { right2PlanarUp, left2PlanarUp },
                            Point3d.Unset, Point3d.Unset,
                            LoftType.Normal, false
                            );

                        Brep[] brepN2PlanarDown = Brep.CreateFromLoft(
                            new List<Curve>() { right2PlanarDown, left2PlanarDown },
                            Point3d.Unset, Point3d.Unset,
                            LoftType.Normal, false
                            );

                        brepN2 = Brep.JoinBreps(new List<Brep>() { brepN2PlanarDown[0], brepN2[0], brepN2PlanarUp[0] }, documentTolerance);
                    }
                    #endregion allow poly surface

                    if (brep2.Length > 0)
                    {
                        oDualLoop2.Add(brepN2[0], path.AppendElement(1));
                        oDualLoop2ID.Add("J;" + triangleIndex.ToString() + ";" + neighbourTriIndex.ToString() + ";" + neighbour2TriIndex.ToString() + ";" + plateID.ToString(), path.AppendElement(1));

                        oDualLoop2Curves.Add(right2, path.AppendElement(1));
                        oDualLoop2Curves.Add(left2, path.AppendElement(1));
                    }
                }
                else
                {
                    oDualLoop2.Add(brep2[0], path.AppendElement(1));
                    oDualLoop2ID.Add("J;" + triangleIndex.ToString() + ";" + neighbourTriIndex.ToString() + ";" + neighbour2TriIndex.ToString() + ";" + plateID.ToString(), path.AppendElement(1));

                    oDualLoop2Curves.Add(left2, path.AppendElement(1));
                    oDualLoop2Curves.Add(right2, path.AppendElement(1));
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

            // get the panel 2 condition and equalize the vertex normal around it. 
            for (int i = 0; i < iSpringMesh.Vertices.Count; i++)
            {
                Vertex vt = iSpringMesh.Vertices[i];
                foreach (int n in vt.NeighborVertexIndices)
                {
                    if (iVertexPanel2[n] == true)
                    {
                        Vector3d temp = vertexNormals[n];
                        vertexNormals[i] = temp;
                    }
                }
            }
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
            // ---- remap curveVerticesValues ------------------------------------------------------------------
            for (int i = 0; i < verticesValues.Count; i++)
            {
                // ---- remap curveVerticesValues ------------------------------------------------------------------
                double remapedValue0 = iTangentScaleMin + (verticesValues[i] - min) * (iTangentScaleMax - iTangentScaleMin) / (max - min);
                curveVerticesValues.Add(remapedValue0);

                // ---- remap planarVerticesValues ------------------------------------------------------------------
                double remapedValue1 = iPlanarOffsetScaleMin + (verticesValues[i] - min) * (iPlanarOffsetScaleMax - iPlanarOffsetScaleMin) / (max - min);
                planarVerticesValues.Add(remapedValue1);

                // ---- remap openingWidthVerticesValues ------------------------------------------------------------------
                double remapedValue2 = iOpeningWidthMin + (verticesValues[i] - min) * (iOpeningWidthMax - iOpeningWidthMin) / (max - min);
                openingWidthVerticesValues.Add(remapedValue2);

                // ---- remap openingWidthVerticesValues ------------------------------------------------------------------
                double remapedValue3 = CurvePointinessMin + (verticesValues[i] - min) * (CurvePointinessMax - CurvePointinessMin) / (max - min);
                pointinessValues.Add(remapedValue3);
            }

            /*
            // ---- remap planarVerticesValues ------------------------------------------------------------------
            for (int i = 0; i < verticesValues.Count; i++)
            {
                double remapedValue = iPlanarOffsetScaleMin + (verticesValues[i] - min) * (iPlanarOffsetScaleMax - iPlanarOffsetScaleMin) / (max - min);
                planarVerticesValues.Add(remapedValue);
            }

            // ---- remap openingWidthVerticesValues ------------------------------------------------------------------
            for (int i = 0; i < verticesValues.Count; i++)
            {
                double remapedValue = iOpeningWidthMin + (verticesValues[i] - min) * (iOpeningWidthMax - iOpeningWidthMin) / (max - min);
                openingWidthVerticesValues.Add(remapedValue);
            }

            // ---- remap openingWidthVerticesValues ------------------------------------------------------------------
            for (int i = 0; i < verticesValues.Count; i++)
            {
                double remapedValue = CurvePointinessMin + (verticesValues[i] - min) * (CurvePointinessMax - CurvePointinessMin) / (max - min);
                pointinessValues.Add(remapedValue);
            }
            */
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


            Point3d b = bottomCps[secondVertexIndex];
            Point3d B = topCps[secondVertexIndex];

            Point3d c = bottomCps[thirdVertexIndex];
            Point3d C = topCps[thirdVertexIndex];



            // Centre of Edge AB
            Point3d ab = 0.5 * (a + b);
            Point3d AB = 0.5 * (A + B);

            // Centre of Edge AC
            Point3d ac = 0.5 * (a + c);
            Point3d AC = 0.5 * (A + C);


            // Centre of Mesh Face
            Point3d m = (ab + ac + 0.5 * (b + c)) / 3;
            Point3d M = (AB + AC + 0.5 * (B + C)) / 3;


            // perform planar Offset
            double planarOffset = planarVerticesValues[firstVertexIndex];
            double curveScaler = curveVerticesValues[firstVertexIndex];
            double openingScaler = openingWidthVerticesValues[firstVertexIndex];
            double curvePointiness = pointinessValues[firstVertexIndex];

            if (ManualAdjustedVertexIndexies.Contains(firstVertexIndex))
            {
                GH_Path pth = new GH_Path(firstVertexIndex);

                planarOffset = ManualValueTree[new GH_Path(firstVertexIndex)][2].Value;
                curveScaler = ManualValueTree[new GH_Path(firstVertexIndex)][0].Value;
                openingScaler = ManualValueTree[new GH_Path(firstVertexIndex)][3].Value;
                curvePointiness = ManualValueTree[new GH_Path(firstVertexIndex)][1].Value;

            }

            // Offset : Centre of Edge AB
            Vector3d v_ab = a - b;
            Vector3d v_AB = A - B;
            v_ab.Unitize();
            v_AB.Unitize();

            // Offset : Centre of Edge AC
            Vector3d v_ac = a - c;
            Vector3d v_AC = A - C;
            v_ac.Unitize();
            v_AC.Unitize();

            // not used
            #region Correct Offset
            /*
           
            // Angles
            Vector3d ac_ab = ac - ab; ac_ab.Unitize();
            Vector3d AC_AB = AC - AB; AC_AB.Unitize();

            double dotProduct_b = ac_ab * v_ab;
            double dotProduct_B = AC_AB * v_AB;
            double dotProduct_c = -ac_ab * v_ac;
            double dotProduct_C = -AC_AB * v_AC;
            
            
            Utils.ToDegree(dotProduct_b);
            Utils.ToDegree(dotProduct_B);
            Utils.ToDegree(dotProduct_c);
            Utils.ToDegree(dotProduct_C);
            

            v_ab *= planarOffset / Math.Tan(Math.Acos(dotProduct_b));
            v_AB *= planarOffset / Math.Tan(Math.Acos(dotProduct_B));
            v_ac *= planarOffset / Math.Tan(Math.Acos(dotProduct_c));
            v_AC *= planarOffset / Math.Tan(Math.Acos(dotProduct_C));
                          
            */
            #endregion correct offset



            // create offset point 
            Point3d oab = ab + (v_ab * planarOffset);
            Point3d oAB = AB + (v_AB * planarOffset);

            // create offset point 
            Point3d oac = ac + (v_ac * planarOffset);
            Point3d oAC = AC + (v_AC * planarOffset);



            Point3d AAB = curveScaler * A + (1.0 - curveScaler) * oAB;
            Point3d aab = curveScaler * a + (1.0 - curveScaler) * oab;

            Point3d AAC = curveScaler * A + (1.0 - curveScaler) * oAC;
            Point3d aac = curveScaler * a + (1.0 - curveScaler) * oac;



            Point3d openingAAB = (-1 * v_AB * openingScaler) + A;
            Point3d openingaab = (-1 * v_ab * openingScaler) + a;

            Point3d openingAAC = (-1 * v_AC * openingScaler) + A;
            Point3d openingaac = (-1 * v_ac * openingScaler) + a;

            // Condition, choose between Opening Width and tangentscale

            if (new Vector3d(AAB - A).Length > openingScaler)
                AAB = openingAAB;

            if (new Vector3d(aab - a).Length > openingScaler)
                aab = openingaab;

            if (new Vector3d(AAC - A).Length > openingScaler)
                AAC = openingAAC;

            if (new Vector3d(aac - a).Length > openingScaler)
                aac = openingaac;


            // Pointiness of the Surface


            Point3d max_aA = 0.5 * (a + A);
            Point3d min_aA1 = 0.5 * (aab + AAB);
            Point3d min_aA2 = 0.5 * (aac + AAC);
            Point3d aA1 = (curvePointiness * max_aA) + ((1 - curvePointiness) * min_aA1);
            Point3d aA2 = (curvePointiness * max_aA) + ((1 - curvePointiness) * min_aA2);

            // Create Curves
            Curve profileCurve1 = Curve.CreateControlPointCurve(new List<Point3d>() { oab, aab, aA1, AAB, oAB }, curveDegree);
            Curve profileCurve2 = Curve.CreateControlPointCurve(new List<Point3d>() { oac, aac, aA2, AAC, oAC }, curveDegree);
            Curve profileCrv11 = Curve.CreateControlPointCurve(new List<Point3d>() { ab, oab }, 1);
            Curve profileCrv12 = Curve.CreateControlPointCurve(new List<Point3d>() { oAB, AB }, 1);
            Curve profileCrv21 = Curve.CreateControlPointCurve(new List<Point3d>() { ac, oac }, 1);
            Curve profileCrv22 = Curve.CreateControlPointCurve(new List<Point3d>() { oAC, AC }, 1);

            Curve profileCurveNew1 = Curve.JoinCurves(new List<Curve>() { profileCrv11, profileCurve1, profileCrv12 }, documentTolerance, true)[0];
            Curve profileCurveNew2 = Curve.JoinCurves(new List<Curve>() { profileCrv21, profileCurve2, profileCrv22 }, documentTolerance, true)[0];

            oTriLoopCurves.Add(profileCurveNew1, path.AppendElement(item));
            oTriLoopCurves.Add(profileCurveNew2, path.AppendElement(item));

            if (iPolySrf)
            {
                //not used method : if triloop need to be polySurface. 
                #region Genes Method
                /*
                PolyCurve polyCurve1 = new PolyCurve();
                polyCurve1.Append(new LineCurve(m, ab));
                polyCurve1.Append(new LineCurve(ab, oab));
                polyCurve1.Append(profileCurve1);
                polyCurve1.Append(new LineCurve(oAB, AB));
                polyCurve1.Append(new LineCurve(AB, M));

                PolyCurve polyCurve2 = new PolyCurve();
                polyCurve2.Append(new LineCurve(m, ac));
                polyCurve2.Append(new LineCurve(ac, oac));
                polyCurve2.Append(profileCurve2);
                polyCurve2.Append(new LineCurve(oAC, AC));
                polyCurve2.Append(new LineCurve(AC, M));

                oTriLoop.Add(
                    Brep.CreateFromLoft(
                       new List<Curve>() { polyCurve1, polyCurve2 },
                       Point3d.Unset, Point3d.Unset,
                       LoftType.Normal,
                       false
                       )[0], path.AppendElement(item)
                   );

                */
                #endregion Genes Method

                //used method
                #region Seperate Planar from Curved Method
                PolyCurve polyCurve1 = new PolyCurve();
                polyCurve1.Append(profileCurve1);

                PolyCurve polyCurve2 = new PolyCurve();
                polyCurve2.Append(profileCurve2);

                // Planar Part Curves
                Curve planarCurveBottom = Curve.CreateControlPointCurve(new List<Point3d>() { m, ac, oac, oab, ab, m }, 1);
                Curve planarCurveTop = Curve.CreateControlPointCurve(new List<Point3d>() { M, AB, oAB, oAC, AC, M }, 1);

                // Add Planar Crv
                oTriLoopPlanCrv.Add(planarCurveBottom, path.AppendElement(item));
                oTriLoopPlanCrv.Add(planarCurveTop, path.AppendElement(item));

                // Add Bended Brep
                oTriLoop.Add(
                    Brep.CreateFromLoft(
                       new List<Curve>() { polyCurve1, polyCurve2 },
                       Point3d.Unset, Point3d.Unset,
                       LoftType.Normal,
                       false
                       )[0], path.AppendElement(item)
                   );
                #endregion Seperate Planar from Curved Method

            }
            else
            #region single surface
            {
                profileCurve1 = Curve.JoinCurves(
                    new List<Curve>() { new LineCurve(b, oab), profileCurve1, new LineCurve(oAB, B) },
                    documentTolerance,
                    true)[0];

                profileCurve2 = Curve.JoinCurves(
                    new List<Curve>() { new LineCurve(c, oac), profileCurve2, new LineCurve(oAC, C) },
                    documentTolerance,
                    true)[0];

                // Planar Part Curves
                Curve planarCurveBottom = Curve.CreateControlPointCurve(new List<Point3d>() { m, ac, oac, oab, ab, m }, 1);
                Curve planarCurveTop = Curve.CreateControlPointCurve(new List<Point3d>() { M, AB, oAB, oAC, AC, M }, 1);
                // Add Planar Crvs
                oTriLoopPlanCrv.Add(planarCurveBottom, path.AppendElement(item));
                oTriLoopPlanCrv.Add(planarCurveTop, path.AppendElement(item));

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
            #endregion single surface

            #region close panel

            // close panel
            //Point3d closeVertice = Point3dList.ClosestPointInList(iClosedPanelPts, iSpringMesh.Vertices[firstVertexIndex].Position);
            //if (ICD.Utils.Distance(closeVertice, iSpringMesh.Vertices[firstVertexIndex].Position) < iClosePanelDist)
            if (iVertexPanel2[firstVertexIndex])
            {
                oClosedPanel.Add(
                    Brep.CreateFromLoft(
                       new List<Curve>() { new LineCurve(A, oAB), new LineCurve(A, oAC) },
                       Point3d.Unset, Point3d.Unset,
                       LoftType.Normal,
                       false
                       )[0], path.AppendElement(item)
                    );
                oClosedPanel.Add(
                    Brep.CreateFromLoft(
                       new List<Curve>() { new LineCurve(oab, a), new LineCurve(oac, a) },
                       Point3d.Unset, Point3d.Unset,
                       LoftType.Normal,
                       false
                       )[0], path.AppendElement(item)
                    );
            }

            #endregion close panel


            // do EffectorHoles
            List<Point3d> EffectorHoleTop = EffectorHoles(A, B, C, M, oAB, oAC, item);
            List<Point3d> EffectorHoleBottom = EffectorHoles(a, b, c, m, oab, oac, item);

            foreach (Point3d pt in EffectorHoleTop)
                oTriLoopEffectorHoles.Add(pt, path.AppendElement(item).AppendElement(1));

            foreach (Point3d pt in EffectorHoleBottom)
                oTriLoopEffectorHoles.Add(pt, path.AppendElement(item).AppendElement(0));
        }

        private List<Point3d> EffectorHoles(Point3d a, Point3d b, Point3d c, Point3d m, Point3d oab, Point3d oac, int idx)
        {
            double distance1 = 0.075;
            double distance2 = 0.145;
            Vector3d normal = Vector3d.CrossProduct(new Vector3d(a - b), new Vector3d(a - c));
            double angleRadians = Utils.ToRadian(120);

            List<Point3d> effectorHoles = new List<Point3d>();

            // check distance to end-of-planarity-line
            Vector3d v_toPlanarLine = (0.5 * (oac + oab)) - m;

            if (v_toPlanarLine.Length < 0.145)
            {
                distance1 = 0.035;
                distance2 = 0.075;
            }

            // for 1st stipe (main stripe)
            Vector3d v_am = a - m;
            v_am.Unitize();

            // for 2nd Stripe
            if (idx == 1)
            {
                v_am = c - m;
                v_am.Unitize();
                v_am.Rotate(angleRadians, normal);
            }

            // for 3rd Stripe
            else if (idx == 2)
            {
                v_am = b - m;
                v_am.Unitize();
                v_am.Rotate((2 * angleRadians), normal);
            }

            Point3d pt1 = m + (v_am * distance1);
            Point3d pt2 = m + (v_am * distance2);

            effectorHoles.Add(pt1);
            effectorHoles.Add(pt2);

            return effectorHoles;

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