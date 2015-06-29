using System;
using System.Collections.Generic;
using System.Linq;

using Grasshopper.Kernel;
using Rhino.Geometry;

using Pavillion2015.DelaunayTriangulation;

using ICD;

namespace Pavillion2015
{
    public class GHC_PointToMesh : GH_Component
    {
        List<Point3d> iPoints = null;
        List<Circle> iApicals = null;
        List<double> iMasses = null;
        Curve iBoundaryCurve = null;
        double iRestLengthScale = double.NaN;
        double iRestLengthOffset = double.NaN;
        double iStiffness = double.NaN;
        double iBendingStiffness = double.NaN;
        bool iInOutSwitch = false;
        double iInOutThreshold = double.NaN;
        double iFixedPointThreshold = double.NaN;
        List<Curve> iFixedPointExclusion = null;
        List<Curve> iFixedPointInclusion = null;
        double iInitialCurvingBias = double.NaN;
        bool iBiasApicalRegion = false;

        String oInfo = string.Empty;
        List<Point3d> oDebugPoints1 = null;
        List<Point3d> oDebugPoints2 = null;
        List<Curve> oDebugCurves1 = null;
        List<Curve> oDebugCurves2 = null;
        List<Vector3d> oDebugVectors1 = null;
        List<Vector3d> oDebugVectors2 = null;
        List<double> oDebugNumbers1 = null;
        List<double> oDebugNumbers2 = null;

        SpringMesh oSpringMesh = null;


        //===============================================================================================


        public GHC_PointToMesh()
            : base("Points to Mesh", "Points to Mesh",
                "Points to Mesh",
                "Pavillion 2015", "Pavillion 2015")
        {
        }


        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddPointParameter("Points", "Points - 00", "Points", GH_ParamAccess.list);
            pManager.AddCircleParameter("Apicals", "Apicals - 01", "Apicals", GH_ParamAccess.list);
            pManager.AddNumberParameter("Masses", "Masses - 02", "Masses", GH_ParamAccess.list);
            pManager.AddCurveParameter("Boundary Curve", "Boundary Curve - 03", "Boundary Curve", GH_ParamAccess.item);
            pManager.AddNumberParameter("Rest Length Scale", "Rest Length Scale - 04", "Rest Length Scale", GH_ParamAccess.item, 1.0);
            pManager.AddNumberParameter("Rest Length Offset", "Rest Length Offset - 05", "Rest Length Offset", GH_ParamAccess.item, 0.0);
            pManager.AddNumberParameter("Stiffness", "Stiffness- 06", "Stiffness Scale", GH_ParamAccess.item, 1.0);
            pManager.AddNumberParameter("Bending Stiffness", "Bending Stiffness - 07", "Stiffness Offset", GH_ParamAccess.item, 0.0);
            pManager.AddBooleanParameter("In/Out Switch", "In/Out Switch - 08", "In/Out Switch", GH_ParamAccess.item, false);
            pManager.AddNumberParameter("In/Out Threshold", "In/Out Threshold - 09", "In/Out Threshold", GH_ParamAccess.item, 0.1);
            pManager.AddNumberParameter("Fixed Point Threshold", "Fixed Point Threshold - 10", "Fixed Point Threshold", GH_ParamAccess.item, 0.05);
            pManager.AddCurveParameter("Fixed Point Exclusion", "Fixed Point Exclusion - 11", "Fixed Point Exclusion", GH_ParamAccess.list);
            pManager.AddCurveParameter("Fixed Point Inclusion", "Fixed Point Inclusion - 12", "Fixed Point Inclusion", GH_ParamAccess.list);
            pManager.AddNumberParameter("Initial Curving Bias", "Initial Curving Bias - 13", "Initial Curving Bias", GH_ParamAccess.item, 0.0);
            pManager.AddBooleanParameter("Bias Apical Regions", "Bias Apical Regions - 14", "Bias Apical Regions", GH_ParamAccess.item, false);
        }


        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddTextParameter("Info", "00 - Info", "Information", GH_ParamAccess.item);
            pManager.AddGenericParameter("Debug 1", "01 - Debug 1", "This output is reserved for debugging", GH_ParamAccess.list);
            pManager.AddGenericParameter("Debug 2", "02 - Debug 2", "This output is reserved for debugging", GH_ParamAccess.list);
            pManager.AddGenericParameter("Debug 3", "03 - Debug 3", "This output is reserved for debugging", GH_ParamAccess.list);
            pManager.AddGenericParameter("Debug 4", "04 - Debug 4", "This output is reserved for debugging", GH_ParamAccess.list);
            pManager.AddGenericParameter("Debug 5", "05 - Debug 5", "This output is reserved for debugging", GH_ParamAccess.list);
            pManager.AddGenericParameter("Mesh", "06 - Mesh", "Mesh", GH_ParamAccess.item);
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
            iPoints = new List<Point3d>();
            DA.GetDataList<Point3d>("Points", iPoints);
            iApicals = new List<Circle>();
            DA.GetDataList<Circle>("Apicals", iApicals);
            iMasses = new List<double>();
            DA.GetDataList<double>("Masses", iMasses);
            DA.GetData<Curve>("Boundary Curve", ref iBoundaryCurve);
            DA.GetData<double>("Rest Length Scale", ref iRestLengthScale);
            DA.GetData<double>("Rest Length Offset", ref iRestLengthOffset);
            DA.GetData<double>("Stiffness", ref iStiffness);
            DA.GetData<double>("Bending Stiffness", ref iBendingStiffness);
            DA.GetData<bool>("In/Out Switch", ref iInOutSwitch);
            DA.GetData<double>("In/Out Threshold", ref iInOutThreshold);
            DA.GetData<double>("Fixed Point Threshold", ref iFixedPointThreshold);
            iFixedPointExclusion = new List<Curve>();
            DA.GetDataList<Curve>("Fixed Point Exclusion", iFixedPointExclusion);
            iFixedPointInclusion = new List<Curve>();
            DA.GetDataList<Curve>("Fixed Point Inclusion", iFixedPointInclusion);


            // ===========================================================================================
            // Compute Delaunay Triangulation
            // ===========================================================================================

            List<DelaunayVertex> delaunayVertices = new List<DelaunayVertex>();
            foreach (Point3d point in iPoints)
                delaunayVertices.Add(new DelaunayVertex(point.X, point.Y));

            List<Triad> triads = new DelaunayTriangulator().Triangulation(delaunayVertices);

            HashSet<Tuple<int, int>> delaunayEdgeTuples = new HashSet<Tuple<int, int>>();


            for (int i = 0; i < triads.Count; i++)
            {
                Triad triad = triads[i];
                delaunayEdgeTuples.Add(triad.a < triad.b ? new Tuple<int, int>(triad.a, triad.b) : new Tuple<int, int>(triad.b, triad.a));
                delaunayEdgeTuples.Add(triad.b < triad.c ? new Tuple<int, int>(triad.b, triad.c) : new Tuple<int, int>(triad.c, triad.b));
                delaunayEdgeTuples.Add(triad.c < triad.a ? new Tuple<int, int>(triad.c, triad.a) : new Tuple<int, int>(triad.a, triad.c));
            }


            // ===========================================================================================
            // Convert Delaunay mesh to particle-spring mesh
            // ===========================================================================================

            oSpringMesh = new SpringMesh();

            // Create edge list -----------------------------------------------------------------------------------------------

            foreach (Tuple<int, int> delaunayEdgeTuple in delaunayEdgeTuples)
            {
                Point3d A = iPoints[delaunayEdgeTuple.Item1];
                Point3d B = iPoints[delaunayEdgeTuple.Item2];
                Point3d M = 0.5 * (A + B);

                // Skip if the edge lies outside of the boundary
                double t;
                iBoundaryCurve.ClosestPoint(M, out t);
                Point3d N = iBoundaryCurve.PointAt(t);
                if (Vector3d.CrossProduct(iBoundaryCurve.TangentAt(t), M - N).Z * (iInOutSwitch ? -1.0 : 1.0) < 0.0
                    && Utils.DistanceSquared(M, N) > iInOutThreshold * iInOutThreshold)
                    continue;

                double edgeLength = Utils.Distance(A, B);
                double restLength = iRestLengthScale * edgeLength + iRestLengthOffset;
                oSpringMesh.Edges.Add(new Edge(delaunayEdgeTuple.Item1, delaunayEdgeTuple.Item2, restLength, iStiffness, Math.PI, iBendingStiffness));
            }

            // Create vertex list -----------------------------------------------------------------------------------------------

            List<HashSet<int>> neighborVerticesSets = new List<HashSet<int>>();
            for (int i = 0; i < iPoints.Count; i++)
                neighborVerticesSets.Add(new HashSet<int>());

            foreach (Edge edge in oSpringMesh.Edges)
            {
                neighborVerticesSets[edge.FirstVertexIndex].Add(edge.SecondVertexIndex);
                neighborVerticesSets[edge.SecondVertexIndex].Add(edge.FirstVertexIndex);
            }

            for (int i = 0; i < iPoints.Count; i++)
            {
                Point3d p = iPoints[i];
                double t;
                iBoundaryCurve.ClosestPoint(p, out t);

                bool vertexFixedness = true;

                if (Utils.Distance(p, iBoundaryCurve.PointAt(t)) > iFixedPointThreshold)
                    vertexFixedness = false;
                else
                    foreach (Curve curve in iFixedPointExclusion)
                        if (curve.Contains(p) == PointContainment.Inside)
                        {
                            vertexFixedness = false;
                            break;
                        }

                foreach (Curve curve in iFixedPointInclusion)
                    if (curve.Contains(p) == PointContainment.Inside)
                    {
                        vertexFixedness = true;
                        break;
                    }

                Vertex vertex = new Vertex(p, Vector3d.Zero, neighborVerticesSets[i].ToList<int>(), iMasses.Count == 1 ? iMasses[0] : iMasses[i], vertexFixedness);
                oSpringMesh.Vertices.Add(vertex);
            }


            // Set boundary edge -----------------------------------------------------------------------------------------------

            foreach (Edge edge in oSpringMesh.Edges)
            {
                if (oSpringMesh.Vertices[edge.FirstVertexIndex].IsFixed && oSpringMesh.Vertices[edge.SecondVertexIndex].IsFixed)
                {
                    edge.IsBoundaryEdge = true;                
                }
            }


            // Create triangle list ------------------------------------------------------------------------------------------------

            Dictionary<Tuple<int, int, int>, int> tripletDict = new Dictionary<Tuple<int, int, int>, int>();

            for (int k = 0; k < oSpringMesh.Edges.Count; k++)
            {
                Edge edge = oSpringMesh.Edges[k];
                Vertex A = oSpringMesh.Vertices[edge.FirstVertexIndex];
                Vertex B = oSpringMesh.Vertices[edge.SecondVertexIndex];

                for (int i = 0; i < A.NeighborVertexIndices.Count; i++)
                {
                    for (int j = 0; j < B.NeighborVertexIndices.Count; j++)
                    {
                        if (A.NeighborVertexIndices[i] == B.NeighborVertexIndices[j])
                        {
                            Tuple<int, int, int> triplet = sortTriplet(edge.FirstVertexIndex, edge.SecondVertexIndex, A.NeighborVertexIndices[i]);

                            if (tripletDict.ContainsKey(triplet))
                            {
                                if (edge.FirstTriangleIndex < 0)
                                {
                                    edge.FirstTriangleIndex = tripletDict[triplet];
                                    edge.FirstAdjacentVertexIndex = A.NeighborVertexIndices[i];
                                }
                                else
                                {
                                    edge.SecondTriangleIndex = tripletDict[triplet];
                                    edge.SecondAdjacentVertexIndex = A.NeighborVertexIndices[i];
                                }
                            }
                            else
                            {
                                oSpringMesh.Triangles.Add(new Triangle(triplet.Item1, triplet.Item2, triplet.Item3));

                                int triangleIndex = oSpringMesh.Triangles.Count - 1;

                                if (edge.FirstTriangleIndex < 0)
                                {
                                    edge.FirstTriangleIndex = triangleIndex;
                                    edge.FirstAdjacentVertexIndex = A.NeighborVertexIndices[i];
                                }
                                else
                                {
                                    edge.SecondTriangleIndex = triangleIndex;
                                    edge.SecondAdjacentVertexIndex = A.NeighborVertexIndices[i];
                                }

                                tripletDict.Add(triplet, triangleIndex);
                            }
                        }
                    }
                }
            }   

            // ===========================================================================================
            // Compute edge indices for each triangle
            // ===========================================================================================

            for (int i = 0; i < oSpringMesh.Edges.Count; i++)
            {
                Edge edge = oSpringMesh.Edges[i];

                Triangle triangle = oSpringMesh.Triangles[edge.FirstTriangleIndex];
                if (triangle.FirstEdgeIndex == -1) triangle.FirstEdgeIndex = i;
                else if (triangle.SecondEdgeIndex == -1) triangle.SecondEdgeIndex = i;
                else triangle.ThirdEdgeIndex = i;

                if (edge.SecondTriangleIndex == -1) continue;

                triangle = oSpringMesh.Triangles[edge.SecondTriangleIndex];
                if (triangle.FirstEdgeIndex == -1) triangle.FirstEdgeIndex = i;
                else if (triangle.SecondEdgeIndex == -1) triangle.SecondEdgeIndex = i;
                else triangle.ThirdEdgeIndex = i;
            }


            // ===========================================================================================
            // Rearange the vertex order in each triangle so the normal calculation is consistent
            // ===========================================================================================

            for (int i = 0; i < oSpringMesh.Triangles.Count; i++)
            {
                if (oSpringMesh.ComputeTriangleNormal(i).Z < 0.0)
                {
                    int temp = oSpringMesh.Triangles[i].SecondVertexIndex;
                    oSpringMesh.Triangles[i].SecondVertexIndex = oSpringMesh.Triangles[i].ThirdVertexIndex;
                    oSpringMesh.Triangles[i].ThirdVertexIndex = temp;
                }
            }


            // ===========================================================================================
            // Rearranging edge adjacent vertex indices for consitency
            // ===========================================================================================

            foreach (Edge edge in oSpringMesh.Edges)
            {
                if (edge.SecondAdjacentVertexIndex == -1)
                {
                    Point3d A = oSpringMesh.Vertices[edge.FirstVertexIndex].Position;
                    Point3d B = oSpringMesh.Vertices[edge.SecondVertexIndex].Position;
                    Point3d M = oSpringMesh.Vertices[edge.FirstAdjacentVertexIndex].Position;

                    if (Vector3d.CrossProduct(B - A, M - A) * oSpringMesh.ComputeTriangleNormal(edge.FirstTriangleIndex) < 0.0)
                    {
                        Point3d temp = A;
                        A = B;
                        B = temp;
                    }
                }
                else
                {
                    Point3d A = oSpringMesh.Vertices[edge.FirstVertexIndex].Position;
                    Point3d B = oSpringMesh.Vertices[edge.SecondVertexIndex].Position;
                    Point3d M = oSpringMesh.Vertices[edge.FirstAdjacentVertexIndex].Position;
                    Point3d N = oSpringMesh.Vertices[edge.SecondAdjacentVertexIndex].Position;

                    if (Vector3d.CrossProduct(B - A, M - A) * oSpringMesh.ComputeTriangleNormal(edge.FirstAdjacentVertexIndex) < 0.0)
                    {
                        int temp = edge.FirstAdjacentVertexIndex;
                        edge.FirstAdjacentVertexIndex = edge.SecondAdjacentVertexIndex;
                        edge.SecondAdjacentVertexIndex = temp;

                        temp = edge.FirstTriangleIndex;
                        edge.FirstTriangleIndex = edge.SecondTriangleIndex;
                        edge.SecondTriangleIndex = temp;
                    }
                }
            }

            // ===========================================================================================
            // Compute adjacent vertex index for each triangle
            // ===========================================================================================

            //foreach (Triangle triangle in oSpringMesh.Triangles)
            //{
            //    Vertex firstVertex = oSpringMesh.Vertices[triangle.FirstVertexIndex];
            //    Vertex secondVertex = oSpringMesh.Vertices[triangle.SecondVertexIndex];
            //    Vertex thirdVertex = oSpringMesh.Vertices[triangle.ThirdVertexIndex];

            //    foreach (int firstNeighbourIndex in firstVertex.NeighborVertexIndices)
            //        foreach (int secondNeighbourIndex in secondVertex.NeighborVertexIndices)
            //            if (firstNeighbourIndex == secondNeighbourIndex && firstNeighbourIndex != triangle.ThirdVertexIndex)
            //                triangle.FirstSecondAdjacentVertexIndex = firstNeighbourIndex;

            //    foreach (int secondNeighbourIndex in secondVertex.NeighborVertexIndices)
            //        foreach (int thirdNeighbourIndex in thirdVertex.NeighborVertexIndices)
            //            if (secondNeighbourIndex == thirdNeighbourIndex && secondNeighbourIndex != triangle.FirstVertexIndex)
            //                triangle.SecondThirdAdjacentVertexIndex = secondNeighbourIndex;

            //    foreach (int thirdNeighbourIndex in thirdVertex.NeighborVertexIndices)
            //        foreach (int firstNeighbourIndex in firstVertex.NeighborVertexIndices)
            //            if (thirdNeighbourIndex == firstNeighbourIndex && thirdNeighbourIndex != triangle.SecondVertexIndex)
            //                triangle.ThirdFirstAdjacentVertexIndex = thirdNeighbourIndex;
            //}


            // ===========================================================================================
            // Initial curving bias
            // ===========================================================================================

            DA.GetData<double>("Initial Curving Bias", ref iInitialCurvingBias);
            DA.GetData<bool>("Bias Apical Regions", ref iBiasApicalRegion);

            foreach (Vertex vertex in oSpringMesh.Vertices)
            {
                vertex.Position.Z = computeBiasHeight(vertex.Position, iBiasApicalRegion);
            }


            // ===========================================================================================
            // Conclusion
            // ===========================================================================================

            foreach (Edge edge in oSpringMesh.Edges)
                oDebugCurves1.Add(new LineCurve(oSpringMesh.Vertices[edge.FirstVertexIndex].Position, oSpringMesh.Vertices[edge.SecondVertexIndex].Position));

            DA.SetData(0, oInfo);
            DA.SetDataList(1, oDebugCurves1);
            DA.SetData(6, oSpringMesh);
        }


        private double computeBiasHeight(Point3d P, bool biasApicalRegions)
        {
            double t;
            iBoundaryCurve.ClosestPoint(P, out t);
            double minDist = Utils.Distance(P, iBoundaryCurve.PointAt(t));

            if (biasApicalRegions)
                foreach (Circle apical in iApicals)
                {
                    double dist = Utils.Distance(P, apical.Center) - (apical.Radius + 0.1);
                    if (dist < 0.0) dist = 0.0;

                    if (dist < minDist)
                        minDist = dist;
                }

            return iInitialCurvingBias * (1.0 - Math.Pow(0.9, minDist * 15.0));
        }


        private Tuple<int, int, int> sortTriplet(int a, int b, int c)
        {
            if (a < b)
            {
                if (b < c) return new Tuple<int, int, int>(a, b, c);
                else
                {
                    if (a < c) return new Tuple<int, int, int>(a, c, b);
                    else return new Tuple<int, int, int>(c, a, b);
                }
            }
            else // b < a
            {
                if (a < c) return new Tuple<int, int, int>(b, a, c);
                else // c < a
                {
                    if (b < c) return new Tuple<int, int, int>(b, c, a);
                    else return new Tuple<int, int, int>(c, b, a);
                }
            }
        }


        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                // return Resources.IconForThisComponent;
                return null;
            }
        }


        public override Guid ComponentGuid
        {
            get { return new Guid("{85855199-ac9b-4a6c-af34-cb72a518bfac}"); }
        }
    }
}