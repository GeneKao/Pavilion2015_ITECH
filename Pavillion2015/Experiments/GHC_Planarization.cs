using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;

using ICD;

namespace Pavillion2015.Experiments
{
    public class GHC_Planarization : GH_Component
    {
        List<Sphere> iSphere = new List<Sphere>();
        List<Vertex> flatPlateCenter = new List<Vertex>();
        SpringMesh iSpringMesh = new SpringMesh();

        public GHC_Planarization()
            : base("Planarization", "Planarization",
                "Planarization groups of triangles surrounding chosen vertices",
                "Pavillion 2015", "Experiments")
        {
        }


        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Spring Mesh", "Spring Mesh", "Spring Mesh", GH_ParamAccess.item);
            pManager.AddCircleParameter("Planarization Regions", "Planarization Regions", "Planarization Regions", GH_ParamAccess.list);
            pManager.AddNumberParameter("Bending Stiffness", "Bending Stiffness", "Bending Stiffness", GH_ParamAccess.item, 0.0);
            pManager.AddNumberParameter("Stiffness", "Stiffness", "Stiffness", GH_ParamAccess.item, 0.0);
            pManager.AddNumberParameter("Bending Stiffness (Planar Regions)", "Bending Stiffness (Planar Regions)", "Bending Stiffness (Planar Regions)", GH_ParamAccess.list);
            pManager.AddNumberParameter("Stiffness (Planar Regions)", "Stiffness (Planar Regions)", "Stiffness (Planar Regions)", GH_ParamAccess.list);
        }


        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Spring Mesh", "Spring Mesh", "Spring Mesh", GH_ParamAccess.item);
        }


        protected override void SolveInstance(IGH_DataAccess DA)
        {
            SpringMesh iSpringMesh = new SpringMesh();
            List<Circle> iCircles = new List<Circle>();
            double iBendingStiffness = double.NaN;
            double iStiffness = double.NaN;
            List<double> iPlanarRegionBendingStiffness = new List<double>();
            List<double> iPlanarRegionStiffness = new List<double>();

            DA.GetData<SpringMesh>(0, ref iSpringMesh);
            DA.GetDataList<Circle>(1, iCircles);
            DA.GetData<double>(2, ref iBendingStiffness);
            DA.GetData<double>(3, ref iStiffness);
            DA.GetDataList<double>(4, iPlanarRegionBendingStiffness);
            DA.GetDataList<double>(5, iPlanarRegionStiffness);

            SpringMesh oSpringMesh = new SpringMesh(iSpringMesh);

            foreach (Edge edge in iSpringMesh.Edges)
            {
                edge.RestLength = Utils.Distance(oSpringMesh.Vertices[edge.FirstVertexIndex].Position, oSpringMesh.Vertices[edge.SecondVertexIndex].Position);
                edge.Stiffness = iStiffness;
                if (edge.SecondTriangleIndex >= 0)
                {
                    edge.RestAngle = Utils.AngleBetweenTwoTriangles(
                        iSpringMesh.Vertices[edge.FirstVertexIndex].Position,
                        iSpringMesh.Vertices[edge.SecondVertexIndex].Position,
                        iSpringMesh.Vertices[edge.FirstAdjacentVertexIndex].Position,
                        iSpringMesh.Vertices[edge.SecondAdjacentVertexIndex].Position
                        );
                }
                edge.BendingStiffness = iBendingStiffness;
            }

            for (int i = 0; i < oSpringMesh.Vertices.Count; i++)
                for (int j = 0; j < iCircles.Count; j++ )
                    if (Utils.Distance(iCircles[j].Center, oSpringMesh.Vertices[i].Position) <= iCircles[j].Radius)
                    {
                        foreach (Edge edge in oSpringMesh.Edges)
                            if (edge.FirstVertexIndex == i || edge.SecondVertexIndex == i)
                            {
                                edge.Stiffness = iPlanarRegionStiffness[j];
                                edge.RestAngle = Math.PI;
                                edge.BendingStiffness = iPlanarRegionBendingStiffness[j];
                            }
                        break;
                    }
            DA.SetData(0, oSpringMesh);
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
            get { return new Guid("{c7a47126-50a6-430f-8733-916e6f88ce41}"); }
        }
    }
}