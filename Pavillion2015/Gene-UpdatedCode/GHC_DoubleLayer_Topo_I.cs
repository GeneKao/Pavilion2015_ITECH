using System;
using System.Collections.Generic;
using System.Linq;

using Grasshopper.Kernel;
using Rhino.Geometry;
using Rhino;
using Rhino.Collections;

namespace Pavillion2015
{
    public class GHC_DoubleLayer_Topo_I : GH_Component
    {
        // input
        Mesh iRhinoMesh = null;
        List<int> iID = null;
        double iMinThickness = double.NaN;
        double iMaxThickness = double.NaN;
        List<Point3d> iAttrThicknessPts = null;
        bool iHeightThickness = false;
        // ------------------------------------------
        List<Point3d> iClosedPanelPts = null;  // close panel second type
        double iClosePanelDist = 0.001;
        // ---------------------------------------------
        bool iAutoGenPlates = true;
        double iPlatesOffset = double.NaN;
        double iPlatesThreads = double.NaN;
        double iPlatesThreadsP = double.NaN;

        // output
        string oInfo = string.Empty;
        List<Curve> oDebugList = null;
        List<Curve> oSMeshLine = null;
        SpringMesh oSpringMesh = null;
        List<Curve> oPolyLine = null;
        List<bool> oVertexStatus = null;
        List<int> oPolyLineID = null;
        List<double> othickness = null;
        List<bool> oVertexPanel2 = null;

        // internal use data




        public GHC_DoubleLayer_Topo_I()
            : base("Double Layer Topo - I", "Double Layer Topo - I",
                "Double Layer Topology I - get plate ID and convert to SpringMesh",
                "Pavillion 2015", "Double Layer")
        {
        }


        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddMeshParameter("RhinoMesh", "RhinoMesh", "RhinoMesh", GH_ParamAccess.item);
            pManager.AddIntegerParameter("Plates ID List", "Plates ID List", "Plates ID List", GH_ParamAccess.list);
            pManager.AddNumberParameter("Min. Thickness", "Min. Thickness", "Min. Thickness", GH_ParamAccess.item, 0.3);
            pManager.AddNumberParameter("Max. Thickness", "Max. Thickness", "Max. Thickness", GH_ParamAccess.item, 1.0);
            pManager.AddPointParameter("Attr. ThicknessPts", "Attr. ThicknessPts", "Attr. ThicknessPts", GH_ParamAccess.list, new Point3d(10000, 10000, 10000));
            pManager.AddBooleanParameter("H. EffectThickness", "H. EffectThickness", "H. EffectThickness", GH_ParamAccess.item, false);
            // ------------------------------------------------
            pManager.AddPointParameter("Closed Panel Area", "Closed Panel Area", "Closed Panel Area", GH_ParamAccess.list, new Point3d());
            pManager.AddNumberParameter("Panel Effect Area", "Panel Effect Area", "Panel Effect Area", GH_ParamAccess.item, 0.001);
            // ---------------------------------------------------
            pManager.AddBooleanParameter("Auto Gen Plates", "Auto Gen Plates", "Auto Gen Plates", GH_ParamAccess.item, true);
            pManager.AddNumberParameter("Plate Offset", "Plate Offset", "Plate Offset", GH_ParamAccess.item, 0.5);
            pManager.AddNumberParameter("Threads", "Threads", "Threads", GH_ParamAccess.item, 1.0);
            pManager.AddNumberParameter("ThreadsP", "ThreadsP", "ThreadsP", GH_ParamAccess.item, 1.0);
        }


        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddTextParameter("Info", "Debug - Info", "Information", GH_ParamAccess.item);
            pManager.AddGenericParameter("General List", "Degug - List ", "Debug", GH_ParamAccess.list);
            pManager.AddGenericParameter("SMesh Line", "SMesh Line", "SMesh Line", GH_ParamAccess.list);
            pManager.AddGenericParameter("Spring Mesh", "Spring Mesh", "Spring Mesh", GH_ParamAccess.item);
            pManager.AddCurveParameter("PolyLines", "PolyLines", "PolyLines", GH_ParamAccess.list);
            pManager.AddNumberParameter("Thickness", "Thickness", "Thickness", GH_ParamAccess.list);
            pManager.AddBooleanParameter("Vertex Status", "Vertex Status", "Vertex Status", GH_ParamAccess.list);
            pManager.AddIntegerParameter("PolyLineID", "PolyLineID", "PolyLineID", GH_ParamAccess.list);
            pManager.AddBooleanParameter("Vertex Panel2", "Vertex Panel2", "Vertex Panel2", GH_ParamAccess.list);
        }


        protected override void BeforeSolveInstance()
        {
            // input
            iID = new List<int>();  
            iAttrThicknessPts = new List<Point3d>();
            //--------------------------------------------------
            iClosedPanelPts = new List<Point3d>();

            // output
            oInfo = string.Empty;
            oDebugList = new List<Curve>();
            oSMeshLine = new List<Curve>();
            oSpringMesh = new SpringMesh();
            oPolyLine = new List<Curve>();
            oVertexStatus = new List<bool>(); // true is plate occupied, false should put triLoop
            oPolyLineID = new List<int>(); // for vertex status access the order of polyline
            othickness = new List<double>();
            oVertexPanel2 = new List<bool>();

            // internal use data

        }


        protected override void SolveInstance(IGH_DataAccess DA)
        {
            // getting input
            DA.GetData<Mesh>("RhinoMesh", ref iRhinoMesh);
            DA.GetDataList<int>("Plates ID List", iID);
            DA.GetData<double>("Min. Thickness", ref iMinThickness);
            DA.GetData<double>("Max. Thickness", ref iMaxThickness);
            DA.GetDataList<Point3d>("Attr. ThicknessPts", iAttrThicknessPts);
            DA.GetData<bool>("H. EffectThickness", ref iHeightThickness);
            // -----------------------------------------------------------------
            DA.GetDataList<Point3d>("Closed Panel Area", iClosedPanelPts);
            DA.GetData<double>("Panel Effect Area", ref iClosePanelDist);
            // -------------------------------------------------------------------
            DA.GetData<bool>("Auto Gen Plates", ref iAutoGenPlates);
            DA.GetData<double>("Plate Offset", ref iPlatesOffset);
            DA.GetData<double>("Threads", ref iPlatesThreads);
            DA.GetData<double>("ThreadsP", ref iPlatesThreadsP);

            //------------------------------------------------------------
            // function start 

            // convert rhino mesh to spring mesh
            ConvertToSpringMesh();

            // initial some spring mesh properties
            // set each vertex not occupied by any plates
            foreach (Vertex vertex in oSpringMesh.Vertices)
            {
                oVertexStatus.Add(false);
                oPolyLineID.Add(-1);
                oVertexPanel2.Add(false);
            }

            ComputePlates();

            oInfo += "Debuging...";
            //------------------------------------------------------------

            // setting output
            DA.SetData("Info", oInfo);
            DA.SetDataList("General List", oDebugList);
            DA.SetData("Spring Mesh", oSpringMesh);
            DA.SetDataList("SMesh Line", oSMeshLine);
            DA.SetDataList("PolyLines", oPolyLine);
            DA.SetDataList("Thickness", othickness);
            DA.SetDataList("Vertex Status", oVertexStatus);
            DA.SetDataList("PolyLineID", oPolyLineID);
            DA.SetDataList("Vertex Panel2", oVertexPanel2);
            // -----------------------------------------------------------
        }


        // compute plate and output polyline
        public void ComputePlates()
        {
            // set up second type of plate 
            for (int i = 0; i < oSpringMesh.Vertices.Count; i++)
            {
                Point3d vertexPosition = oSpringMesh.Vertices[i].Position;
                Point3d closeVertice = Point3dList.ClosestPointInList(iClosedPanelPts, vertexPosition);
                if (ICD.Utils.Distance(closeVertice, vertexPosition) < iClosePanelDist)
                    oVertexPanel2[i] = true;
            }

            // attractor thickness
            for (int i = 0; i < oSpringMesh.Vertices.Count; i++)
            {
                Point3d vertexPosition = oSpringMesh.Vertices[i].Position;

                // different thickness
                Point3d closestPt = Point3dList.ClosestPointInList(iAttrThicknessPts, vertexPosition);
                double distValue = ICD.Utils.Distance(vertexPosition, closestPt);

                if (iHeightThickness)  // adding z coordinate as account
                {
                    double k = vertexPosition.Z * 0.125; 
                    distValue *= k * iMinThickness + (0.5 - k) * iMaxThickness;
                }
                othickness.Add(distValue);
            }

            // fix the second type plate area thickness
            for (int i = 0; i < oSpringMesh.Vertices.Count; i++)
            {
                Vertex v = oSpringMesh.Vertices[i];
                List<int> Neighbours = v.NeighborVertexIndices.ToList();
                oInfo += "i = " + i.ToString() + ", value = " + oVertexPanel2[i] + "\n";
                foreach ( int n in Neighbours)
                {
                    //if (oVertexPanel2[n] == true)
                    //{
                    oInfo += "n = " + n.ToString() + ", value = " + oVertexPanel2[n] + "\n";
                        //oInfo += n.ToString() + " get it \n";
                        //double temp = othickness[n];
                        //othickness[i] = temp;
                    //}
                }
                oInfo += "==============================\n";
            }

            // remap
            double max = othickness.Max();
            double min = othickness.Min();

            // map(value, low1, high1, low2, high2) = > low2 + (value - low1) * (high2 - low2) / (high1 - low1)
            for (int i = 0; i < othickness.Count; i++)
                othickness[i] = iMinThickness + (othickness[i] - min) * (iMaxThickness - iMinThickness) / (max - min);
            // ------------------------- end attractors --------------------------------------------------------------------------

            // test
            foreach (int id in iID)
            {
                Circle c = new Circle(oSpringMesh.Vertices[id].Position, 0.2);
                oDebugList.Add(NurbsCurve.CreateFromCircle(c));
            }

            // create plate using TPI
            DualTPI(false);
            DualTPI(true);

        }

        public void DualTPI(bool flip)
        {

            List<List<LineCurve>> dualVertices = new List<List<LineCurve>>();

            for (int i = 0; i < oSpringMesh.Vertices.Count; i++)
                dualVertices.Add(new List<LineCurve>());

            ////// this is to create planar plate
            foreach (Edge edge in oSpringMesh.Edges)
            {
                if (edge.SecondTriangleIndex >= 0)
                {
                    int firstVertexIndex = edge.FirstVertexIndex;
                    int secondVertexIndex = edge.SecondVertexIndex;
                    int firstTriangleIndex = edge.FirstTriangleIndex;
                    int secondTriangleIndex = edge.SecondTriangleIndex;

                    Point3d firstTPI = oSpringMesh.ComputeDoubleLayerTPI(firstTriangleIndex, othickness, flip, iPlatesOffset);
                    Point3d secondTPI = oSpringMesh.ComputeDoubleLayerTPI(secondTriangleIndex, othickness, flip, iPlatesOffset);

                    Point3d firstTPIO = oSpringMesh.ComputeTPI(firstTriangleIndex);
                    Point3d secondTPIO = oSpringMesh.ComputeTPI(secondTriangleIndex);
                    Point3d firstCenterPtO = oSpringMesh.ComputeTriangleCentroid(firstTriangleIndex);
                    Point3d secondCenterPtO = oSpringMesh.ComputeTriangleCentroid(secondTriangleIndex);

                    Point3d firstCenterPt;
                    Point3d secondCenterPt;

                    if (!flip)
                    {
                        firstCenterPt = oSpringMesh.ComputeTriangleCentroid(firstTriangleIndex) +
                            oSpringMesh.ComputeTriangleNormal(firstTriangleIndex);
                        secondCenterPt = oSpringMesh.ComputeTriangleCentroid(secondTriangleIndex) +
                            oSpringMesh.ComputeTriangleNormal(secondTriangleIndex);
                    }
                    else
                    {
                        firstCenterPt = oSpringMesh.ComputeTriangleCentroid(firstTriangleIndex) -
                            oSpringMesh.ComputeTriangleNormal(firstTriangleIndex);
                        secondCenterPt = oSpringMesh.ComputeTriangleCentroid(secondTriangleIndex) -
                            oSpringMesh.ComputeTriangleNormal(secondTriangleIndex);
                    }

                    // ===========================================================================
                    // Added 14/06/2015 to make some tolerance for plates

                    Circle firstCircle = oSpringMesh.ComputeDoubleLayerIncircle(firstTriangleIndex, iPlatesThreadsP, othickness, flip, iPlatesOffset);
                    Circle secondCircle = oSpringMesh.ComputeDoubleLayerIncircle(secondTriangleIndex, iPlatesThreadsP, othickness, flip, iPlatesOffset);

                    //Circle firstCircle = oSpringMesh.ComputeDoubleLayerCircumscribedCircle(firstTriangleIndex, iPlatesThreadsP, othickness, flip, iPlatesOffset);
                    //Circle secondCircle = oSpringMesh.ComputeDoubleLayerCircumscribedCircle(secondTriangleIndex, iPlatesThreadsP, othickness, flip, iPlatesOffset);

                    Circle firstCircleO = oSpringMesh.ComputeIncircle(firstTriangleIndex, iPlatesThreadsP); // iPlatesThreads = 0.0-1.0
                    Circle secondCircleO = oSpringMesh.ComputeIncircle(secondTriangleIndex, iPlatesThreadsP);

                    Plane firstPO = oSpringMesh.ComputePlane(firstTriangleIndex);
                    Plane secondPO = oSpringMesh.ComputePlane(secondTriangleIndex);

                    Point3d projectFirstTPIO = firstPO.ClosestPoint(firstTPIO);
                    Point3d projectSecondTPIO = secondPO.ClosestPoint(secondTPIO);


                    Point3d newFirstTPI = firstCircle.ClosestPoint(firstTPI);
                    Point3d newSecondTPI = secondCircle.ClosestPoint(secondTPI);

                    // manuel selected plates
                    foreach (int id in iID)
                    {
                        if (id == firstVertexIndex)
                            dualVertices[firstVertexIndex].Add(new LineCurve(newFirstTPI, newSecondTPI));
                        if (id == secondVertexIndex)
                            dualVertices[secondVertexIndex].Add(new LineCurve(newFirstTPI, newSecondTPI));
                    }
                    // auto generated plates
                    if (ICD.Utils.Distance(projectFirstTPIO, firstCenterPtO) < iPlatesThreads &&
                        ICD.Utils.Distance(projectSecondTPIO, secondCenterPtO) < iPlatesThreads &&
                        iAutoGenPlates && !iID.Contains(firstVertexIndex) && !iID.Contains(secondVertexIndex))
                    {
                        dualVertices[firstVertexIndex].Add(new LineCurve(newFirstTPI, newSecondTPI));
                        dualVertices[secondVertexIndex].Add(new LineCurve(newFirstTPI, newSecondTPI));
                    }
                }
                else
                {
                }
            }

            int counter = -1;

            for (int i = 0; i < oSpringMesh.Vertices.Count; i++)
            {
                // if is close and all right
                if (dualVertices[i].Count == oSpringMesh.Vertices[i].NeighborVertexIndices.Count
                    && !iID.Contains(i))
                {
                    Curve dualcurve = Curve.JoinCurves(dualVertices[i])[0];

                    oInfo += dualcurve.IsClosed.ToString() + "\n";

                    if (dualcurve.IsClosed)
                    {
                        oPolyLine.Add(dualcurve);
                        counter++;
                        // change the vertex staus to be occupied 
                        oVertexStatus[i] = true;
                        oPolyLineID[i] = counter;
                    }
                }
                // manuel selected plate 
                foreach (int id in iID)
                {
                    if (id == i)
                    {
                        Curve dualcurve = Curve.JoinCurves(dualVertices[i])[0];
                        if (dualcurve.IsClosed)
                        {
                            oPolyLine.Add(dualcurve);
                            counter++;
                            // change the vertex staus to be occupied 
                            oVertexStatus[i] = true;
                            oPolyLineID[i] = counter;
                        }
                    }
                }
            }
        }

        // function to convert rhino mesh to spring mesh
        public void ConvertToSpringMesh()
        {

            bool[] naked = iRhinoMesh.GetNakedEdgePointStatus();

            for (int i = 0; i < iRhinoMesh.Vertices.Count; i++)
            {
                Point3d vertex = iRhinoMesh.Vertices[i];
                oSpringMesh.Vertices.Add(new Vertex(vertex, Vector3d.Zero, new List<int>()));

                // boundary and fixed condition
                if (naked[i] == true)
                {
                    oSpringMesh.Vertices[i].IsBoundaryVertex = true;
                    oSpringMesh.Vertices[i].IsFixed = true;
                }

                // vertex neighbours
                int topoIndex = iRhinoMesh.TopologyVertices.TopologyVertexIndex(i);
                int[] connectIndex = iRhinoMesh.TopologyVertices.ConnectedTopologyVertices(topoIndex);
                for (int j = 0; j < connectIndex.Length; j++)
                    oSpringMesh.Vertices[i].NeighborVertexIndices.Add(iRhinoMesh.TopologyVertices.MeshVertexIndices(j)[0]);
            }

            foreach (MeshFace face in iRhinoMesh.Faces)
                oSpringMesh.Triangles.Add(new Triangle(face.A, face.B, face.C));

            for (int i = 0; i < iRhinoMesh.TopologyEdges.Count; i++)
            {
                IndexPair indexPair = iRhinoMesh.TopologyEdges.GetTopologyVertices(i);
                // should convert TopologyVertices to mesh vertices first
                int firstIndex = iRhinoMesh.TopologyVertices.MeshVertexIndices(indexPair.I)[0];
                int secondIndex = iRhinoMesh.TopologyVertices.MeshVertexIndices(indexPair.J)[0];

                double len = iRhinoMesh.TopologyEdges.EdgeLine(i).Length;

                oSpringMesh.Edges.Add(new Edge(firstIndex, secondIndex, len, 1.0, Math.PI * 0.8, 0.0));

                // edge vertex
                oSpringMesh.Edges[i].FirstVertexIndex = firstIndex;
                oSpringMesh.Edges[i].SecondVertexIndex = secondIndex;

            }

            for (int i = 0; i < iRhinoMesh.TopologyEdges.Count; i++)
            {
                // connected faces
                int[] connectedFacesIndex = iRhinoMesh.TopologyEdges.GetConnectedFaces(i);
                int firstTriIndex = -1;
                int secondTriIndex = -1;

                IndexPair indexPair = iRhinoMesh.TopologyEdges.GetTopologyVertices(i);
                int firstIndex = iRhinoMesh.TopologyVertices.MeshVertexIndices(indexPair.I)[0];
                int secondIndex = iRhinoMesh.TopologyVertices.MeshVertexIndices(indexPair.J)[0];

                
                if (connectedFacesIndex.Length == 2)
                {
                    firstTriIndex = connectedFacesIndex[0];
                    oSpringMesh.Edges[i].FirstTriangleIndex = firstTriIndex;
                    secondTriIndex = connectedFacesIndex[1];
                    oSpringMesh.Edges[i].SecondTriangleIndex = secondTriIndex;

                    int triangleAdj01 = 0;
                    int triangleAdj02 = 0;
                    int vertex01 = oSpringMesh.Triangles[firstTriIndex].FirstVertexIndex;
                    int vertex02 = oSpringMesh.Triangles[firstTriIndex].SecondVertexIndex;
                    int vertex03 = oSpringMesh.Triangles[firstTriIndex].ThirdVertexIndex;
                    int vertex11 = oSpringMesh.Triangles[secondTriIndex].FirstVertexIndex;
                    int vertex12 = oSpringMesh.Triangles[secondTriIndex].SecondVertexIndex;
                    int vertex13 = oSpringMesh.Triangles[secondTriIndex].ThirdVertexIndex;

                    int[] triangleVertexList01 = { vertex01, vertex02, vertex03 };
                    int[] triangleVertexList02 = { vertex11, vertex12, vertex13 };

                    for (int j = 0; j < triangleVertexList01.Length; j++)
                        if (triangleVertexList01[j] != firstIndex && triangleVertexList01[j] != secondIndex)
                            triangleAdj01 = triangleVertexList01[j];
                    for (int j = 0; j < triangleVertexList02.Length; j++)
                        if (triangleVertexList02[j] != firstIndex && triangleVertexList02[j] != secondIndex)
                            triangleAdj02 = triangleVertexList02[j];

                    oSpringMesh.Edges[i].FirstAdjacentVertexIndex = triangleAdj01;
                    oSpringMesh.Edges[i].SecondAdjacentVertexIndex = triangleAdj02;

                    oSpringMesh.Edges[i].IsBoundaryEdge = false;

                    // triangle edge and triangle
                    if (triangleAdj01 == vertex01)
                    {
                        oSpringMesh.Triangles[firstTriIndex].FirstEdgeIndex = i;
                        oSpringMesh.Triangles[firstTriIndex].FirstAdjTriIndex = secondTriIndex;
                    }
                    if (triangleAdj01 == vertex02)
                    {
                        oSpringMesh.Triangles[firstTriIndex].SecondEdgeIndex = i;
                        oSpringMesh.Triangles[firstTriIndex].SecondAdjTriIndex = secondTriIndex;
                    }
                    if (triangleAdj01 == vertex03)
                    {
                        oSpringMesh.Triangles[firstTriIndex].ThirdEdgeIndex = i;
                        oSpringMesh.Triangles[firstTriIndex].ThirdAdjTriIndex = secondTriIndex;
                    }
                    if (triangleAdj02 == vertex11)
                    {
                        oSpringMesh.Triangles[secondTriIndex].FirstEdgeIndex = i;
                        oSpringMesh.Triangles[secondTriIndex].FirstAdjTriIndex = firstTriIndex;
                    }
                    if (triangleAdj02 == vertex12)
                    { 
                        oSpringMesh.Triangles[secondTriIndex].SecondEdgeIndex = i;
                        oSpringMesh.Triangles[secondTriIndex].SecondAdjTriIndex = firstTriIndex;
                    }
                    if (triangleAdj02 == vertex13)
                    {
                        oSpringMesh.Triangles[secondTriIndex].ThirdEdgeIndex = i;
                        oSpringMesh.Triangles[secondTriIndex].ThirdAdjTriIndex = firstTriIndex;
                    }

                }
                if (connectedFacesIndex.Length == 1)
                {
                    firstTriIndex = connectedFacesIndex[0];
                    oSpringMesh.Edges[i].FirstTriangleIndex = firstTriIndex;

                    int triangleAdj01 = 0;
                    int vertex01 = oSpringMesh.Triangles[firstTriIndex].FirstVertexIndex;
                    int vertex02 = oSpringMesh.Triangles[firstTriIndex].SecondVertexIndex;
                    int vertex03 = oSpringMesh.Triangles[firstTriIndex].ThirdVertexIndex;

                    int[] triangleVertexList01 = { vertex01, vertex02, vertex03 };

                    for (int j = 0; j < triangleVertexList01.Length; j++)
                        if (triangleVertexList01[j] != firstIndex && triangleVertexList01[j] != secondTriIndex)
                            triangleAdj01 = triangleVertexList01[j];

                    oSpringMesh.Edges[i].FirstAdjacentVertexIndex = triangleAdj01;

                    oSpringMesh.Edges[i].IsBoundaryEdge = true;

                    // triangle edge
                    if (triangleAdj01 == vertex01)
                        oSpringMesh.Triangles[firstTriIndex].FirstEdgeIndex = i;
                    if (triangleAdj01 == vertex02)
                        oSpringMesh.Triangles[firstTriIndex].SecondEdgeIndex = i;
                    if (triangleAdj01 == vertex03)
                        oSpringMesh.Triangles[firstTriIndex].ThirdEdgeIndex = i;
                }
                
            }

            // visualize springMesh
            foreach (Edge edge in oSpringMesh.Edges)
                oSMeshLine.Add(new LineCurve(oSpringMesh.Vertices[edge.FirstVertexIndex].Position, oSpringMesh.Vertices[edge.SecondVertexIndex].Position));
        }

        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                //You can add image files to your project resources and access them like this:
                // return Resources.IconForThisComponent;
                return Properties.Resources.Topo1;
            }
        }


        public override Guid ComponentGuid
        {
            get { return new Guid("{a001484b-717e-412a-a815-d9a72c3e86b6}"); }
        }
    }

}