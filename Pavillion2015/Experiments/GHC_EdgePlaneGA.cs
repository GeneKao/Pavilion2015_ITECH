using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;
using Rhino.Geometry.Intersect;

using ICD;
using ICD.Beagle;


namespace Pavillion2015
{
    public class GHC_EdgePlaneGA : GH_Component
    {

        SpringMesh iSpringMesh = null;
        bool iReset = false;
        bool iPlay = false;
        int iSubiterationCount = 0;
        int iPopulationSize = 0;
        double iMaxCouplingDistance = double.NaN;
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

        List<Point3d> triangleCentres = new List<Point3d>();
        List<Plane> trianglePlanes = new List<Plane>();
        List<Vector3d> triangleNormals = new List<Vector3d>();
        List<Vector3d> edgeXAxes = new List<Vector3d>();
        List<Vector3d> edgeYAxes = new List<Vector3d>();
        List<Vector3d> edgeZAxes = new List<Vector3d>();
        List<Plane> edgePlanes = new List<Plane>();
        List<Point3d> edgeOrigins = new List<Point3d>();

        Beagle beagle = null;

        List<double> temps = new List<double>();
        
        public GHC_EdgePlaneGA()
            : base("Edge Plane GA", "Edge Plane GA",
                "Edge Plane GA",
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
            pManager.AddNumberParameter("Max. Coupling Distance", "Max. Coupling Distance", "Max. Coupling Distance", GH_ParamAccess.item, 1.0);
            pManager.AddNumberParameter("Mutation Rate", "Mutation Rate", "Mutation Rate", GH_ParamAccess.item, 0.1);
            pManager.AddNumberParameter("Mutation Amount", "Mutation Amount", "Mutation Amount", GH_ParamAccess.item, 0.5);
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
            // =============================================================================================================
            // Reset
            // =============================================================================================================

            DA.GetData<bool>("Reset", ref iReset);

            if (iReset)
            {
                DA.GetData<SpringMesh>("Spring Mesh", ref iSpringMesh);
                DA.GetData<int>("Population Size", ref iPopulationSize);
              
                triangleCentres = new List<Point3d>();
                trianglePlanes = new List<Plane>();
                triangleNormals = new List<Vector3d>();
                for (int i = 0; i < iSpringMesh.Triangles.Count; i++)
                {                    
                    triangleCentres.Add(iSpringMesh.ComputeCircumscribedCircleCenter(i));
                    trianglePlanes.Add(new Plane(
                        iSpringMesh.Vertices[iSpringMesh.Triangles[i].FirstVertexIndex].Position,
                        iSpringMesh.Vertices[iSpringMesh.Triangles[i].SecondVertexIndex].Position,
                        iSpringMesh.Vertices[iSpringMesh.Triangles[i].ThirdVertexIndex].Position
                        ));
                    triangleNormals.Add(iSpringMesh.ComputeTriangleNormal(i));
                }

                edgeXAxes = new List<Vector3d>();
                edgeYAxes = new List<Vector3d>();
                edgeZAxes = new List<Vector3d>();
                edgePlanes = new List<Plane>();
                edgeOrigins = new List<Point3d>();

                List<double> geneMinValues = new List<double>();
                List<double> geneMaxValues = new List<double>();

                foreach (Edge edge in iSpringMesh.Edges)
                {
                    if (edge.SecondTriangleIndex >= 0)
                    {
                        Vector3d edgeXAxis = iSpringMesh.Vertices[edge.SecondVertexIndex].Position - iSpringMesh.Vertices[edge.FirstVertexIndex].Position;
                        edgeXAxes.Add(edgeXAxis);

                        Vector3d edgeZAxis = triangleNormals[edge.FirstTriangleIndex] + triangleNormals[edge.SecondTriangleIndex];
                        edgeZAxis.Unitize();
                        edgeZAxes.Add(edgeZAxis);

                        Vector3d edgeYAxis = Vector3d.CrossProduct(edgeZAxis, edgeXAxis);
                        edgeYAxes.Add(edgeYAxis);

                        Point3d edgePlaneOrigin = 0.5 * (iSpringMesh.Vertices[edge.FirstVertexIndex].Position + iSpringMesh.Vertices[edge.SecondVertexIndex].Position);

                        edgeOrigins.Add(edgePlaneOrigin);


                        edgePlanes.Add(new Plane(edgePlaneOrigin, edgeXAxis, edgeYAxis));

                        double angle = Utils.AngleBetweenTwoUnitVectors(edgeZAxis, triangleNormals[edge.FirstTriangleIndex]);

                        geneMinValues.Add(-angle);
                        geneMaxValues.Add(angle);
                    }
                    else
                    {
                        Vector3d edgeXAxis = iSpringMesh.Vertices[edge.SecondVertexIndex].Position - iSpringMesh.Vertices[edge.FirstVertexIndex].Position;
                        edgeXAxes.Add(edgeXAxis);

                        Vector3d edgeZAxis = triangleNormals[edge.FirstTriangleIndex];
                        if (edgeZAxis.Z < 0.0)
                        {
                            edgeZAxis += -Vector3d.ZAxis;
                        }
                        else
                        {
                            Vector3d temp = new Vector3d(-edgeZAxis.X, -edgeZAxis.Y, 0.0);
                            temp.Unitize();
                            edgeZAxis += temp;
                        }

                        edgeZAxis.Unitize();
                        edgeZAxes.Add(edgeZAxis);

                        Vector3d edgeYAxis = Vector3d.CrossProduct(triangleNormals[edge.FirstTriangleIndex], edgeXAxis);                    
                        edgeYAxes.Add(edgeYAxis);

                        Point3d edgePlaneOrigin = 0.5 * (iSpringMesh.Vertices[edge.FirstVertexIndex].Position + iSpringMesh.Vertices[edge.SecondVertexIndex].Position);

                        edgeOrigins.Add(edgePlaneOrigin);

                        edgePlanes.Add(new Plane(edgePlaneOrigin, edgeXAxis, edgeYAxis));

                        double angle = Utils.AngleBetweenTwoUnitVectors(edgeZAxis, triangleNormals[edge.FirstTriangleIndex]);

                        geneMinValues.Add(-angle);
                        geneMaxValues.Add(angle);
                    }
                }

                beagle = new Beagle(iPopulationSize, geneMinValues, geneMaxValues);

                foreach (Genome genome in beagle.Genomes)
                    computeFitness(genome);
              
                goto Conclusion;
            }


            DA.GetData<bool>("Play", ref iPlay);

            if (iPlay) ExpireSolution(true);
            else
            {
                List<double> values = new List<double>();
                foreach (Gene gene in beagle.Genomes[0].Genes)
                    values.Add(gene.Value);
                DA.SetDataList(4, values);
                goto Conclusion;
            }
            //if (!iPlay) goto Conclusion;


            // ====================================================================================================
            // Evolution
            // ====================================================================================================

            DA.GetData<int>("Subiteration Count", ref iSubiterationCount);
            DA.GetData<double>("Max. Coupling Distance", ref iMaxCouplingDistance);
            DA.GetData<double>("Mutation Rate", ref iMutationRate);
            DA.GetData<double>("Mutation Amount", ref iMutationAmount);

            beagle.MaxCouplingDistance = iMaxCouplingDistance;
            beagle.MutationRate = iMutationRate;
            beagle.MutationAmount = iMutationAmount;

            for (int i = 0; i < iSubiterationCount; i++)
            {
                beagle.SelectParents();
                beagle.Reproduce();

                for (int j = beagle.PopulationSize; j < beagle.Genomes.Count; j++)
                    computeFitness(beagle.Genomes[j]);

                beagle.SelectFittests();
            }

            DA.SetDataList(3, new List<double> {beagle.Genomes[0].Fitness});

            // ====================================================================================================
            // Conclusion
            // ====================================================================================================

            Conclusion:

            Genome bestGenome = beagle.Genomes[0];

            for (int i = 0; i < bestGenome.Genes.Count; i++)
            {
                double angle = bestGenome.Genes[i].Value;

                Vector3d planeNormal = Math.Sin(angle) * edgeYAxes[i] + Math.Cos(angle) * edgeZAxes[i];

                edgePlanes[i] = new Plane(
                    edgePlanes[i].Origin,
                    edgePlanes[i].XAxis,
                    Vector3d.CrossProduct(planeNormal, edgePlanes[i].XAxis)
                    );
            }

            //List<Plane> nastyPlanes = new List<Plane>();
            List<Point3d> planeTripletIntersections = new List<Point3d>();

            foreach (Triangle triangle in iSpringMesh.Triangles)
            {
                Point3d intersectionPoint;
                if (Intersection.PlanePlanePlane(edgePlanes[triangle.FirstEdgeIndex], edgePlanes[triangle.SecondEdgeIndex], edgePlanes[triangle.ThirdEdgeIndex], out intersectionPoint))
                    planeTripletIntersections.Add(intersectionPoint);
                //else
                //{
                //    nastyPlanes.Add(edgePlanes[triangle.FirstEdgeIndex]);
                //    nastyPlanes.Add(edgePlanes[triangle.FirstEdgeIndex]);
                //    nastyPlanes.Add(edgePlanes[triangle.FirstEdgeIndex]);
                //}
            }

            DA.SetDataList(1, edgePlanes);
            DA.SetDataList(2, planeTripletIntersections);



            

            //DA.SetDataList(1, edgePlanes);
            //DA.SetDataList(2, edgeOrigins);
            //DA.SetDataList(3, edgeXAxes);
            //DA.SetDataList(4, edgeYAxes);
            //DA.SetDataList(5, edgeZAxes);
        }


        private void computeFitness(Genome genome)
        {
            for (int i = 0; i < genome.Genes.Count; i++)
            {
                double angle = genome.Genes[i].Value;

                edgePlanes[i] = new Plane(
                    edgePlanes[i].Origin,
                    edgePlanes[i].XAxis,
                    Vector3d.CrossProduct(
                        edgeYAxes[i] * Math.Sin(angle) + edgeZAxes[i] * Math.Cos(angle),
                        edgePlanes[i].XAxis)
                    );
            }

            genome.Fitness = 0.0;

            for (int i = 0; i < iSpringMesh.Triangles.Count; i++)
            {
                Triangle triangle = iSpringMesh.Triangles[i];
                Point3d intersectionPoint;
                if (Intersection.PlanePlanePlane(
                    edgePlanes[triangle.FirstEdgeIndex],
                    edgePlanes[triangle.SecondEdgeIndex],
                    edgePlanes[triangle.ThirdEdgeIndex],
                    out intersectionPoint))
                {
                    genome.Fitness += Utils.Distance(
                        trianglePlanes[i].ClosestPoint(intersectionPoint),
                        triangleCentres[i]);
                }
            }

            if (double.IsNaN(genome.Fitness))
            {
                int a = 2;
                int b = a + 1;
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
            get { return new Guid("{7a0fed01-fc5d-4f6e-89b0-ec07471b1669}"); }
        }
    }
}