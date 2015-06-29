using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino;
using Rhino.Geometry;

namespace Pavillion2015
{
    public class GHC_RhinoMeshToSpringMesh : GH_Component
    {

        public GHC_RhinoMeshToSpringMesh()
            : base("Rhino Mesh To Spring Mesh", "Spring Mesh",
                "Description",
                "Pavillion 2015", "Utilities")
        {
        }


        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddMeshParameter("Rhino Mesh", "Rhino Mesh", "Rhino Mesh", GH_ParamAccess.item);
        }


        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Debug 1", "Debug 1", "Debug 1", GH_ParamAccess.item);
            pManager.AddGenericParameter("Spring Mesh", "Spring Mesh", "Spring Mesh", GH_ParamAccess.item);
        }


        protected override void SolveInstance(IGH_DataAccess DA)
        {
            Mesh iRhinoMesh = null;
            DA.GetData<Mesh>(0, ref iRhinoMesh);

            SpringMesh oSpringMesh = new SpringMesh();

            bool [] naked = iRhinoMesh.GetNakedEdgePointStatus();

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

                    int[] triangleVertexList01 = {vertex01, vertex02, vertex03}; 
                    int[] triangleVertexList02 = {vertex11, vertex12, vertex13};

                    for (int j = 0; j < triangleVertexList01.Length; j++)
                        if (triangleVertexList01[j] != firstIndex && triangleVertexList01[j] != secondIndex)
                            triangleAdj01 = triangleVertexList01[j];
                    for (int j = 0; j < triangleVertexList02.Length; j++)
                        if (triangleVertexList02[j] != firstIndex && triangleVertexList02[j] != secondIndex)
                            triangleAdj02 = triangleVertexList02[j];

                    oSpringMesh.Edges[i].FirstAdjacentVertexIndex = triangleAdj01;
                    oSpringMesh.Edges[i].SecondAdjacentVertexIndex = triangleAdj02;

                    oSpringMesh.Edges[i].IsBoundaryEdge = false;
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
                }
            }

            DA.SetData(1, oSpringMesh);


            List<LineCurve> oDebugCurves = new List<LineCurve>();

            foreach (Edge edge in oSpringMesh.Edges)
                oDebugCurves.Add(new LineCurve(oSpringMesh.Vertices[edge.FirstVertexIndex].Position, oSpringMesh.Vertices[edge.SecondVertexIndex].Position));

            DA.SetDataList(0, oDebugCurves);
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
            get { return new Guid("{b372f8f2-7007-4b5d-99fc-0db34b7675b9}"); }
        }
    }
}