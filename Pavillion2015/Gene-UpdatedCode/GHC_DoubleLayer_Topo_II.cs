using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;
using Rhino.Collections;

namespace Pavillion2015
{
    public class GHC_DoubleLayer_Topo_II : GH_Component
    {

        // input
        SpringMesh iSpringMesh = null;
        List<Curve> iPolyLine = null;
        List<bool> iVertexStatus = null;
        List<int> iPolyLineID = null;
        List<double> iThickness = null;
        double iTangentScale = double.NaN;
        bool iPolySrf = true;

        // output
        string oInfo = string.Empty;
        List<Point3d> oDebugList = null;
        List<Brep> oTriLoop = null;
        List<Brep> oDualLoop = null;

        // internal usage
        List<Vector3d> vertexNormals;

        List<Point3d> topCps = null;
        List<Point3d> bottomCps = null;

        List<Point3d> topCenterPts = null;
        List<Point3d> bottomCenterPts = null;
        List<bool> topCenterPolygonOccupied = null;

        List<int[]> indexSortPolygon = null;

        double documentTolerance = DocumentTolerance();


        public GHC_DoubleLayer_Topo_II()
            : base("Double Layer Topo - II", "Double Layer Topo - II",
                "Double Layer Topology II - get plate and SpringMesh convert double layers",
                "Pavillion 2015", "Double Layer")
        {
        }


        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Spring Mesh", "Spring Mesh", "Spring Mesh", GH_ParamAccess.item);
            pManager.AddCurveParameter("PolyLine", "PolyLine", "PolyLine", GH_ParamAccess.list);
            pManager.AddNumberParameter("Thickness", "Thickness", "Thickness", GH_ParamAccess.list);
            pManager.AddBooleanParameter("Vertex status", "Vertex status", "Vertex status", GH_ParamAccess.list);
            pManager.AddIntegerParameter("PolyLine ID", "PolyLine ID", "PolyLine ID", GH_ParamAccess.list);
            pManager.AddNumberParameter("Tangent Scale", "Tangent Scale", "Tangent Scale", GH_ParamAccess.item, 0.2);
            pManager.AddBooleanParameter("Allow PolySrf", "Allow PolySrf", "Allow PolySrf", GH_ParamAccess.item, true);
        }


        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddTextParameter("Info", "Debug - Info", "Information", GH_ParamAccess.item);
            pManager.AddGenericParameter("General List", "Degug - List ", "Debug", GH_ParamAccess.list);
            pManager.AddGenericParameter("TriLoop", "TriLoop", "TriLoop", GH_ParamAccess.list);
            pManager.AddGenericParameter("DualLoop", "DualLoop", "DualLoop", GH_ParamAccess.list);
        }


        protected override void BeforeSolveInstance()
        {
            // input

            iSpringMesh = new SpringMesh();
            iPolyLine = new List<Curve>();
            iVertexStatus = new List<bool>(); // true is plate occupied, false should put triLoop
            iPolyLineID = new List<int>(); // for vertex status access the order of polyline
            iThickness = new List<double>();

            // output
            oInfo = string.Empty;
            oDebugList = new List<Point3d>();
            oTriLoop = new List<Brep>();
            oDualLoop = new List<Brep>();

            // internal use data
            topCenterPts = new List<Point3d>();
            bottomCenterPts = new List<Point3d>();
            topCenterPolygonOccupied = new List<bool>(); // not use at this moment 

            bottomCps = new List<Point3d>();
            topCps = new List<Point3d>();

            indexSortPolygon = new List<int[]>();


        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            // getting input
            DA.GetData<SpringMesh>("Spring Mesh", ref iSpringMesh);
            DA.GetDataList<Curve>("PolyLine", iPolyLine);
            DA.GetDataList<int>("PolyLine ID", iPolyLineID);
            DA.GetDataList<double>("Thickness", iThickness);
            DA.GetDataList<bool>("Vertex status", iVertexStatus);
            DA.GetData<double>("Tangent Scale", ref iTangentScale);
            DA.GetData<bool>("Allow PolySrf", ref iPolySrf);
            //------------------------------------------------------------


            storePlatesTPI();

            calculateVertexNormals();
            calculateVertexCps();

            triLoop();

            half_dualLoop();

            //------------------------------------------------------------

            // setting output
            DA.SetData("Info", oInfo);
            DA.SetDataList("General List", oDebugList);
            DA.SetDataList("TriLoop", oTriLoop);
            DA.SetDataList("DualLoop", oDualLoop);
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
                    //oInfo += "outsider i = " + i.ToString() + "\n";
                    start = t;
                    Point3d pt = c.PointAt(t);
                    if (i < iPolyLine.Count / 2)
                    {
                        Point3d closestPt = Point3dList.ClosestPointInList(topCenterPts, pt);
                        int index = topCenterPts.IndexOf(closestPt);
                        topCenterPts[index] = pt;
                        topCenterPolygonOccupied[index] = true; // set occupied status

                        // sort triangle id ///////////////////////////////////////
                        //oInfo += "index = " + index.ToString() + "\n";
                        for (int j = 0; j < 3; j++)
                        {
                            //oInfo += "in for loop \n";
                            if (indexSortPolygon[index][j] == -1)
                            {
                                //oInfo += "indexP " + indexSortPolygon[index][j] + "\n";
                                //oInfo += "indexP " + indexSortPolygon[index][j + 3] + "\n";

                                //oInfo += "j = " + j.ToString() + "\n";
                                //oInfo += "i = " + i.ToString() + "\n";

                                indexSortPolygon[index][j] = i;
                                indexSortPolygon[index][j + 3] = counter;

                                //oInfo += "indexP " + indexSortPolygon[index][j] + "\n";
                                //oInfo += "indexP " + indexSortPolygon[index][j + 3] + "\n";

                                break;
                            }
                        }
                        counter++;
                        //oInfo += "counter = " + counter.ToString() + "\n";
                        //oInfo += "================= \n";
                        /////////////////////////////////////////////////////////////
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

            // debug
            /*
            foreach (int[] ip in indexSortPolygon)
                oInfo += ip[0].ToString() + ", " +
                         ip[1].ToString() + ", " +
                         ip[2].ToString() + ", " +
                         ip[3].ToString() + ", " +
                         ip[4].ToString() + ", " +
                         ip[5].ToString() + "\n";*/
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

                    // the one directly connected to plates
                    if ((iVertexStatus[vertex01] == true && iVertexStatus[vertex02] == false) ||
                        (iVertexStatus[vertex01] == false && iVertexStatus[vertex02] == true)) 
                        half_dualLoop_type_1(i);

                    // point connected dual loop
                    if ((iVertexStatus[vertex01] == false && iVertexStatus[vertex02] == false))
                    {
                        half_dualLoop_type_2(i, tri01);
                        half_dualLoop_type_2(i, tri02);
                    }
                }
                else // boarder condition
                {
                    if ((iVertexStatus[vertex01] == false && iVertexStatus[vertex02] == false))
                    {
                        half_dualLoop_type_2(i, tri01);
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

            if (iVertexStatus[vertex1] == false && iVertexStatus[vertex2] == true) { vertex = vertex1; }
            else if (iVertexStatus[vertex2] == false && iVertexStatus[vertex1] == true) { vertex = vertex2; }
            else { oInfo += "got it..."; }

            Point3d vertexPt = iSpringMesh.Vertices[vertex].Position;

            // layer up
            Point3d rightUp01 = topCenterPts[rightTriIndex];
            Point3d leftUp01 = topCenterPts[leftTriIndex];

            Point3d up03 = topCps[vertex]; //reference point not for construct surfaces

            Point3d rightUp02 = (1 - iTangentScale) * rightUp01 + (iTangentScale) * up03;
            Point3d leftUp02 = (1 - iTangentScale) * leftUp01 + (iTangentScale) * up03;

            // layer down
            Point3d rightDown01 = bottomCenterPts[rightTriIndex];
            Point3d leftDown01 = bottomCenterPts[leftTriIndex];

            Point3d down03 = bottomCps[vertex]; //reference point not for construct surfaces

            Point3d rightDown02 = (1 - iTangentScale) * rightDown01 + (iTangentScale) * down03;
            Point3d leftDown02 = (1 - iTangentScale) * leftDown01 + (iTangentScale) * down03;

            Curve right = Curve.CreateControlPointCurve(
                new List<Point3d>() { rightUp01, rightUp02, vertexPt, rightDown02, rightDown01 });

            Curve left = Curve.CreateControlPointCurve(
                new List<Point3d>() { leftUp01, leftUp02, vertexPt, leftDown02, leftDown01 });


            // compare their polygon index to sort the lofting sequence 
            int[] rightIndex = indexSortPolygon[rightTriIndex];
            int[] leftIndex = indexSortPolygon[leftTriIndex];

            for (int r = 0; r < 3; r++)
                for (int l = 0; l < 3; l++)
                    if (rightIndex[r] == leftIndex[l] && rightIndex[r] != -1 && leftIndex[l] != -1)
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
                                oDualLoop.Add(brep[0]);
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
                                oDualLoop.Add(brep[0]);
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
                                oDualLoop.Add(brep[0]);
                            //break;
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
                                oDualLoop.Add(brep[0]);
                            //break;
                        }

            /*
            Brep[] brep = Brep.CreateFromLoft(
                      new List<Curve>() { right, left },
                      Point3d.Unset, Point3d.Unset,
                      LoftType.Normal,
                      false
                      );
           
            if (brep.Length > 0)
            {
                //oInfo += "got it...";
                oDualLoop.Add(brep[0]);
            }*/
        }

        // dual loop not directly connected to plates
        private void half_dualLoop_type_2(int edgeIndex, int triangleIndex)
        {
            Edge edge = iSpringMesh.Edges[edgeIndex];

            // vertex point
            int vertex1 = edge.FirstVertexIndex;
            int vertex2 = edge.SecondVertexIndex;

            Point3d vertexPt1 = iSpringMesh.Vertices[vertex1].Position;
            Point3d vertexPt2 = iSpringMesh.Vertices[vertex2].Position;

            Point3d vertexPtUp1 = topCps[vertex1];
            Point3d vertexPtUp2 = topCps[vertex2];
            Point3d vertexPtUpM = 0.5 * (vertexPtUp1 + vertexPtUp2);
            Point3d ctPtUp1 = (iTangentScale) * vertexPtUp1 + (1 - iTangentScale) * vertexPtUpM;
            Point3d ctPtUp2 = (iTangentScale) * vertexPtUp2 + (1 - iTangentScale) * vertexPtUpM;

            Point3d vertexPtDown1 = bottomCps[vertex1];
            Point3d vertexPtDown2 = bottomCps[vertex2];
            Point3d vertexPtDownM = 0.5 * (vertexPtDown1 + vertexPtDown2);
            Point3d ctPtDown1 = (iTangentScale) * vertexPtDown1 + (1 - iTangentScale) * vertexPtDownM;
            Point3d ctPtDown2 = (iTangentScale) * vertexPtDown2 + (1 - iTangentScale) * vertexPtDownM;

            Point3d upTPI = topCenterPts[triangleIndex];
            Point3d downTPI = bottomCenterPts[triangleIndex];
            Point3d ctUpTPI1 = (iTangentScale) * vertexPtUp1 + (1 - iTangentScale) * upTPI;
            Point3d ctDownTPI1 = (iTangentScale) * vertexPtDown1 + (1 - iTangentScale) * downTPI;
            Point3d ctUpTPI2 = (iTangentScale) * vertexPtUp2 + (1 - iTangentScale) * upTPI;
            Point3d ctDownTPI2 = (iTangentScale) * vertexPtDown2 + (1 - iTangentScale) * downTPI;

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
                if (brepF.NormalAt(0, 0).Z > 0)
                {
                    Brep[] brepN1 = Brep.CreateFromLoft(
                       new List<Curve>() { left1, right1 },
                       Point3d.Unset, Point3d.Unset,
                       LoftType.Normal,
                       false
                       );
                    if (brep1.Length > 0)
                        oDualLoop.Add(brepN1[0]);
                }
                else
                    oDualLoop.Add(brep1[0]);
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
                if (brepF.NormalAt(0, 0).Z > 0)
                {
                    Brep[] brepN2 = Brep.CreateFromLoft(
                       new List<Curve>() { left2, right2 },
                       Point3d.Unset, Point3d.Unset,
                       LoftType.Normal,
                       false
                       );
                    if (brep2.Length > 0)
                        oDualLoop.Add(brepN2[0]);
                }
                else
                    oDualLoop.Add(brep2[0]);
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

        private void triLoop()
        {
            foreach (Triangle triangle in iSpringMesh.Triangles)
            {
                if (!iVertexStatus[triangle.FirstVertexIndex] &&
                    !iVertexStatus[triangle.SecondVertexIndex] &&
                    !iVertexStatus[triangle.ThirdVertexIndex])
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

                oTriLoop.Add(
                    Brep.CreateFromLoft(
                       new List<Curve>() { polyCurve1, polyCurve2 },
                       Point3d.Unset, Point3d.Unset,
                       LoftType.Normal,
                       false
                       )[0]
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
                else { oTriLoop.Add(brep); return; }

                breps = brep.Trim(new Plane(A_, M, AC), documentTolerance);
                if (breps.Length > 0) brep = breps[0];
                else { oTriLoop.Add(brep); return; }

                normal = Vector3d.CrossProduct(b - a, c - a);

                Point3d a_ = a - offsetAmount * normal;

                breps = brep.Trim(new Plane(a_, m, ab), documentTolerance);
                if (breps.Length > 0) brep = breps[0];
                else { oTriLoop.Add(brep); return; }

                breps = brep.Trim(new Plane(a_, ac, m), documentTolerance);
                if (breps.Length > 0) brep = breps[0];
                else { oTriLoop.Add(brep); return; }

                oTriLoop.Add(brep);
            }

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
            get { return new Guid("{f3b6e70b-11ab-40c3-9429-7976bb521289}"); }
        }
    }
}