using System;
using System.Collections.Generic;

using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Rhino.Geometry;

namespace Pavillion2015
{
    public class GHC_01_ConvertMesh : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the GHC_01_ConvertMesh class.
        /// </summary>
        public GHC_01_ConvertMesh()
            : base("GHC_01_ConvertMesh", "01-ConvertMesh",
                "Setup a mesh for the Triloop",
                "Pavillion 2015", "DoubleLayer")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddMeshParameter("Mesh", "M", "Input the base mesh here", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddPointParameter("Face Points", "FacePt", "Points of the face", GH_ParamAccess.tree); //0
            pManager.AddIntegerParameter("Topology Points", "TopoPtIdx", "Topology points index of the face", GH_ParamAccess.tree);  //1
            pManager.AddVectorParameter("Face Normals", "FaceN", "Normal of the face", GH_ParamAccess.list);    //2
            pManager.AddIntegerParameter("Neighbour Face Index", "NFIdx", "Neighbouring faces", GH_ParamAccess.tree); //3

            pManager.AddPointParameter("All Topology Points", "ATPC", "all topology points coordinates", GH_ParamAccess.list);  //4
            pManager.AddVectorParameter("Topology vertex normals", "TopNormals", "Normals of each topology vertex", GH_ParamAccess.list);   //5

            pManager.AddIntegerParameter("Face Topology Edge Indexies", "FaceTopEdgeIdx", "Indexies of topology edges foreach face", GH_ParamAccess.tree);  //6
            pManager.AddIntegerParameter("Topology Edge Topology Point Index", "TopEdgeTopPtsIdx", "Indexies of Topology Points foreach Topology Edge", GH_ParamAccess.tree);   //7
            pManager.AddIntegerParameter("Topology Point Connected Topology Edges Index", "TopPtsConnectedTopEdgeIdx", "Indexies of connected Topology Edges foreach Topology Point", GH_ParamAccess.tree);   //8

        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {

            //----Declareing--------------------------------------------------------------------------

            // contains the 3 points of each face
            DataTree<Point3f> facePoints = new DataTree<Point3f>();
            // contains the coresponding topology points of each face
            DataTree<int> faceTopologyPoints = new DataTree<int>();
            // contains the face normals of each face
            List<Vector3f> faceNormals = new List<Vector3f>();

            // contains the 3 topology edges of each face
            DataTree<int> faceTopoEdgesIdx = new DataTree<int>();
            // contains the points of each topology edge
            DataTree<int> topologyEdgesTopPtsIdx = new DataTree<int>();

            // Contains the coordinates of each topology point
            List<Point3d> topologyPoints = new List<Point3d>();
            // Contains the index of neighbouring faces for each face
            DataTree<int> faceNeighbours = new DataTree<int>();
            // Contains Normals of topology vertices
            List<Vector3d> topologyNormals = new List<Vector3d>();

            // Contains the index of topology Edges for each Topology Point
            DataTree<int> TopPt_Connected_TopEdges = new DataTree<int>();

            //get Mesh from input
            Mesh M = new Mesh();
            DA.GetData<Mesh>("Mesh", ref M);

            //----End Declareing-----------------------------------------------------------------------


            //----Functions------------------------------------------------------------------------------

            // get List with sublist of 3 points per face
            for (int face_id = 0; face_id < M.Faces.Count; face_id++)
            {
                // set up the branch index
                GH_Path pth = new GH_Path(face_id);

                # region FacePoints

                //---- Face Points (Point3f)
                Point3f A, B, C, D;
                M.Faces.GetFaceVertices(face_id, out A, out B, out C, out D);

                facePoints.Add(A, pth);
                facePoints.Add(B, pth);
                facePoints.Add(C, pth);

                #endregion FacePoints

                #region FaceNormals
                //---- Face Normals (Vector3f)
                M.FaceNormals.ComputeFaceNormals();
                faceNormals.Add(M.FaceNormals[face_id]);
                #endregion FaceNormals

                #region faceTopologyPoints

                //---- Topology Points of the face (int)
                int TA = M.Faces.GetTopologicalVertices(face_id)[0];
                int TB = M.Faces.GetTopologicalVertices(face_id)[1];
                int TC = M.Faces.GetTopologicalVertices(face_id)[2];

                faceTopologyPoints.Add(TA, pth);
                faceTopologyPoints.Add(TB, pth);
                faceTopologyPoints.Add(TC, pth);

                #endregion faceTopologyPoints

                #region faceNeighbours

                //---- Neighbours of face (int)  

                foreach (int i in M.TopologyEdges.GetEdgesForFace(face_id))
                {
                    if (M.TopologyEdges.GetConnectedFaces(i).Length > 1)
                    {
                        foreach (int j in M.TopologyEdges.GetConnectedFaces(i))
                        {
                            if (j != face_id)
                            { faceNeighbours.Add(j, pth); }
                        }
                    }
                    else
                    { faceNeighbours.Add(-1, pth); }
                }


                #endregion faceNeighbours

                #region Face Topology Edges

                //---- Topology Edges (int)  

                foreach (int i in M.TopologyEdges.GetEdgesForFace(face_id))
                { faceTopoEdgesIdx.Add(i, pth); }

                #endregion Face Topology Edges

            }



            for (int i = 0; i < M.TopologyVertices.Count; i++)
            {
                #region topologyPoints
                //---- Topology Points (point3f)
                int[] vertIdx = M.TopologyVertices.MeshVertexIndices(i);
                topologyPoints.Add(M.Vertices[vertIdx[0]]);
                #endregion topologyPoints

                #region topologyNormals
                //---- Topology Normals
                M.FaceNormals.ComputeFaceNormals();
                Vector3d normal = new Vector3d(0, 0, 0);
                int count = 0;

                foreach (int face in M.TopologyVertices.ConnectedFaces(i))
                {
                    Vector3f temp = new Vector3f();
                    temp = M.FaceNormals[face];
                    Vector3d temp2 = new Vector3d(temp.X, temp.Y, temp.Z);
                    normal += temp2;
                    count++;
                }

                normal /= count;
                topologyNormals.Add(normal);
                #endregion topologyNormals

            }

            #region Topology Edges

            for (int i = 0; i < M.TopologyEdges.Count; i++)
            {
                topologyEdgesTopPtsIdx.Add(M.TopologyEdges.GetTopologyVertices(i).I, new GH_Path(i));
                topologyEdgesTopPtsIdx.Add(M.TopologyEdges.GetTopologyVertices(i).J, new GH_Path(i));
            }

            #endregion Topology Edges

            #region Topology Vertex connected Topology Edge

            for (int i = 0; i < topologyPoints.Count; i++)
            {
                // i = index of Topology Point
                GH_Path pth = new GH_Path(i);

                for (int j = 0; j < topologyEdgesTopPtsIdx.BranchCount; j++)
                {
                    // j = index of Topology Edge
                    foreach (int k in topologyEdgesTopPtsIdx.Branch(j))
                    {
                        if (k == i)
                        // add multiple Topology Edges to the branch index, which is representing the topology point index
                        { TopPt_Connected_TopEdges.Add(j, pth); }
                    }
                }
            }

            #endregion Topology Vertex connected Topology Edge

            //----End Functions--------------------------------------------------------------------------

            //----Set Output-----------------------------------------------------------------------------

            DA.SetDataTree(0, facePoints);
            DA.SetDataTree(1, faceTopologyPoints);
            DA.SetDataList(2, faceNormals);
            DA.SetDataTree(3, faceNeighbours);
            DA.SetDataList(4, topologyPoints);
            DA.SetDataList(5, topologyNormals);
            DA.SetDataTree(6, faceTopoEdgesIdx);
            DA.SetDataTree(7, topologyEdgesTopPtsIdx);
            DA.SetDataTree(8, TopPt_Connected_TopEdges);

            //----End Set Output-------------------------------------------------------------------------



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
                return Properties.Resources.ConvertMesh_Icon1;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("{143574c7-bca8-4af7-9812-b1b8faacc6a2}"); }
        }
    }
}