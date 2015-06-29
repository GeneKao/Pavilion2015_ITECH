using System;
using System.Collections.Generic;
using System.Diagnostics;

using Grasshopper.Kernel;
using Rhino.Geometry;

namespace Pavillion2015
{
    public class GHC_MeshRelaxer : GH_Component
    {
        bool iReset = true;
        int iIntegrator = 0;
        bool iPlay = true;
        int iIterationCount = 0;
        SpringMesh iSpringMesh = null;
        double iDamping = double.NaN;
        double iTimeStep = double.NaN;
        public List<Circle> iForceGenerators = null;
        double iGravityScale = double.NaN;
        double iRestLengthScale = double.NaN;
        double iRestLengthOffset = double.NaN;
        double iStiffnessOffset = double.NaN;
        double iBendingStiffness = double.NaN;
        bool iEnableBendingStiffnessOffset = false;
        double iBoundaryBendingStiffness = double.NaN;
        bool iEnableColumnStiffnessOffset = false;
        double iColumnBendingStiffness = double.NaN;

        //added by Gene
        int iEnableEquilateralTriangle = 0;
        double iEquilateralStrength = double.NaN;

        String oInfo = string.Empty;
        List<Point3d> oDebugPoints1 = null;
        List<Point3d> oDebugPoints2 = null;
        List<Curve> oDebugCurves1 = null;
        List<Curve> oDebugCurves2 = null;
        List<Vector3d> oDebugVectors1 = null;
        List<Vector3d> oDebugVectors2 = null;
        List<double> oDebugNumbers1 = null;
        List<double> oDebugNumbers2 = null;

        SpringMesh oSpringMesh;

        int iteration = 0;
        Stopwatch stopwatch = null;

        public GHC_MeshRelaxer()
            : base("Mesh Relaxer", "Mesh Relaxer",
                "Mesh Relaxer",
                "Pavillion 2015", "Pavillion 2015")
        {
        }


        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddBooleanParameter("Reset", "Reset - 00", "Reset", GH_ParamAccess.item);
            pManager.AddIntegerParameter("Integrator", "Integrator - 01", "Select numerical integration method. 1 = Euler, 2 = Velocity Verlet", GH_ParamAccess.item);
            pManager.AddBooleanParameter("Play", "Play - 02", "Play", GH_ParamAccess.item);
            pManager.AddIntegerParameter("Iteration Count", "Iteration Count - 03", "Iteration Count", GH_ParamAccess.item, 800);
            pManager.AddGenericParameter("Spring Mesh", "Spring Mesh - 04", "Spring Mesh", GH_ParamAccess.item);
            pManager.AddNumberParameter("Damping", "Damping - 05", "Damping", GH_ParamAccess.item, 0.2);
            pManager.AddNumberParameter("Time Step", "Time Step - 06", "Time Step", GH_ParamAccess.item, 0.005);
            pManager.AddCircleParameter("Force Generators", "Force Generators - 07", "Force Generators", GH_ParamAccess.list);
            pManager.AddNumberParameter("Gravity Scale", "Gravity Scale - 08", "Gravity Scale", GH_ParamAccess.item, 0.0);
            pManager.AddNumberParameter("Rest Length Scale", "Rest Length Scale - 09", "Rest Length Scale", GH_ParamAccess.item, 1.0);
            pManager.AddNumberParameter("Rest Length Offset", "Rest Length Offset - 10", "Rest Length Offset", GH_ParamAccess.item, 0.0);
            pManager.AddNumberParameter("Stiffness Offset", "Stiffness Offset - 11", "Stiffness Offset", GH_ParamAccess.item, 0.0);
            pManager.AddNumberParameter("Bending Stiffness Offset", "Bending Stiffness Offset - 12", "Bending Stiffness Offset", GH_ParamAccess.item, 0.0);
            pManager.AddBooleanParameter("Enable Boundary Bending Stiffness", "Enable Boundary Bending Stiffness - 13", "Enable Boundary Bending Stiffness", GH_ParamAccess.item, false);
            pManager.AddNumberParameter("Boundary Bending Stiffness", "Boundary Bending Stiffness - 14", "Boundary Bending Stiffness", GH_ParamAccess.item, 0.0);
            pManager.AddBooleanParameter("Enable Column Bending Stiffness", "Enable Column Bending Stiffness - 15", "Enable Column Bending Stiffness", GH_ParamAccess.item, false);
            pManager.AddNumberParameter("Column Bending Stiffness", "Column Bending Stiffness - 16", "Column Bending Stiffness", GH_ParamAccess.item, 0.0);

            //added by Gene
            pManager.AddIntegerParameter("Enable Equilateral Triangle", "Enable Equilateral Triangle - 15", "Enable Equilateral Triangle", GH_ParamAccess.item, 0);
            pManager.AddNumberParameter("Equilateral Strength", "Equilateral Strength - 16", "Equilateral Strength", GH_ParamAccess.item, 1.0);        
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
            // ==================================================================================
            // Reset
            // ==================================================================================

            DA.GetData<bool>("Reset", ref iReset);

            if (iReset)
            {
                iteration = 0;

                DA.GetData<SpringMesh>("Spring Mesh", ref iSpringMesh);
                oSpringMesh = new SpringMesh(iSpringMesh);

                DA.GetData<int>(1, ref iIntegrator);
                oSpringMesh.Integrator = iIntegrator;

                goto Conclusion;
            }

            DA.GetData<bool>("Play", ref iPlay);
            if (!iPlay) goto Conclusion;
            if ((iIterationCount != 0) && (iteration >= iIterationCount)) goto Conclusion;
            iteration++;


            // ==================================================================================
            // Update the mesh
            // ==================================================================================          

            DA.GetData<double>("Damping", ref iDamping);
            oSpringMesh.Damping = iDamping;

            DA.GetData<double>("Time Step", ref iTimeStep);
            oSpringMesh.DeltaT = iTimeStep;

            DA.GetData<double>("Gravity Scale", ref iGravityScale);
            oSpringMesh.GravityScale = iGravityScale;

            DA.GetData<double>("Rest Length Scale", ref iRestLengthScale);
            oSpringMesh.RestLengthScale = iRestLengthScale;

            DA.GetData<double>("Rest Length Offset", ref iRestLengthOffset);
            oSpringMesh.RestLengthOffset = iRestLengthOffset;

            DA.GetData<double>("Stiffness Offset", ref iStiffnessOffset);
            oSpringMesh.Stiffness = iStiffnessOffset;

            DA.GetData<double>("Bending Stiffness Offset", ref iBendingStiffness);
            oSpringMesh.BendingStiffness = iBendingStiffness;

            DA.GetData<bool>("Enable Boundary Bending Stiffness", ref iEnableBendingStiffnessOffset);
            oSpringMesh.EnableBoundaryBendingStiffness = iEnableBendingStiffnessOffset;

            DA.GetData<double>("Boundary Bending Stiffness", ref iBoundaryBendingStiffness);
            oSpringMesh.BoundaryBendingStiffness = iBoundaryBendingStiffness;

            DA.GetData<bool>("Enable Column Bending Stiffness", ref iEnableColumnStiffnessOffset);
            oSpringMesh.EnableColumnBendingStiffness = iEnableColumnStiffnessOffset;

            DA.GetData<double>("Column Bending Stiffness", ref iColumnBendingStiffness);
            oSpringMesh.ColumnnBendingStiffness = iColumnBendingStiffness;

            // Added by Gene   ================================================================
            DA.GetData<int>("Enable Equilateral Triangle", ref iEnableEquilateralTriangle);
            oSpringMesh.EnableEquilateralTriangle = iEnableEquilateralTriangle;

            DA.GetData<double>("Equilateral Strength", ref iEquilateralStrength);
            oSpringMesh.EquilateralStrength = iEquilateralStrength;
            // ================================================================================


            oSpringMesh.Update();


            // ==================================================================================
        // Conclusion
        // ==================================================================================

            Conclusion:

            DA.GetData<int>("Iteration Count", ref iIterationCount);
            if (iPlay && (iteration < iIterationCount || iIterationCount == 0)) ExpireSolution(true);

            oInfo += "Iteration: " + iteration.ToString();

            foreach (Edge edge in oSpringMesh.Edges)
                oDebugCurves1.Add(new LineCurve(oSpringMesh.Vertices[edge.FirstVertexIndex].Position, oSpringMesh.Vertices[edge.SecondVertexIndex].Position));

            DA.SetData(0, oInfo);
            DA.SetDataList(1, oDebugCurves1);
            //DA.SetDataList(2, new List<Mesh> { oSpringMesh.ConvertToRhinoMesh() });
            DA.SetData(6, oSpringMesh);
        }


        protected override System.Drawing.Bitmap Icon { get { return null; } }

        public override Guid ComponentGuid { get { return new Guid("{9e9ffa9d-6f43-4d50-89c2-987a61a25298}"); } }
    }
}