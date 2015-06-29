//using System;
//using System.Collections.Generic;

//using Grasshopper.Kernel;
//using Rhino.Geometry;

//using EvoluteRhinoCommon;

//using DevelopableLofterDevelopabilityLevel = EvoluteRhinoCommon.EvoluteDevelopableLofting.DevelopableLofterDevelopabilityLevel;
//using DevelopableLofterSideExtend = EvoluteRhinoCommon.EvoluteDevelopableLofting.DevelopableLofterSideExtend;
//using DevelopableLofterError = EvoluteRhinoCommon.EvoluteDevelopableLofting.DevelopableLofterError;

//namespace Pavillion2015.Experiments
//{
//    public class GHC_EvoluteLoft : GH_Component
//    {
//        public GHC_EvoluteLoft()
//            : base("Evolute Loft", "Evolute Loft",
//                "Create developable loft surface using Evolute D.LOFT tool",
//                "Pavillion 2015", "Experiments")
//        {
//        }


//        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
//        {
//            pManager.AddCurveParameter("First Curve", "First Curve", "First Curve", GH_ParamAccess.item);
//            pManager.AddCurveParameter("Second Curve", "Second Curve", "Second Curve", GH_ParamAccess.item);
//            pManager.AddNumberParameter("Closeness to Tolerance", "Closeness to Tolerance", "Closeness to Tolerance", GH_ParamAccess.item, 0.001);
//            pManager.AddIntegerParameter("Number of Control Points", "Number of Control Points", "Number of Control Points", GH_ParamAccess.item, 0);
//        }


//        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
//        {
//            pManager.AddTextParameter("Info", "Info", "Info", GH_ParamAccess.item);
//            pManager.AddBrepParameter("Loft", "Loft", "Loft", GH_ParamAccess.item);
//        }


//        protected override void SolveInstance(IGH_DataAccess DA)
//        {
//            Curve iFirstCurve = null;
//            Curve iSecondCurve = null;
//            double iClosenessToTolerance = double.NaN;
//            int iNumberOfControlPoints = 0;
//            DA.GetData<Curve>("First Curve", ref iFirstCurve);
//            DA.GetData<Curve>("Second Curve", ref iSecondCurve);
//            DA.GetData<double>("Closeness to Tolerance", ref iClosenessToTolerance);
//            DA.GetData<int>("Number of Control Points", ref iNumberOfControlPoints);

//            //DevelopableLofterDevelopabilityLevel developabilityLevel, DevelopableLofterSideExtend sideExtend, int nbControlPts, ref Brep brep, ref DevelopableLofterError lastError)

//            Brep oBrep = new Brep();
//            DevelopableLofterError developableLofterError = DevelopableLofterError.CURVE_EXTENSION_FAILED;

//            EvoluteDevelopableLofting.DevelopableLofter(
//                ref iFirstCurve, 
//                ref iSecondCurve, 
//                iClosenessToTolerance, 
//                DevelopableLofterDevelopabilityLevel.VERYTIGHT, 
//                DevelopableLofterSideExtend.EXTEND_BOTH, 
//                iNumberOfControlPoints, 
//                ref oBrep, 
//                ref developableLofterError);

//            DA.SetData("Info", developableLofterError.ToString());
//            DA.SetData("Loft", oBrep);
//        }


//        protected override System.Drawing.Bitmap Icon
//        {
//            get
//            {
//                return null;
//            }
//        }


//        public override Guid ComponentGuid
//        {
//            get { return new Guid("{b12acdc9-c5fd-41dd-8696-f93499d6f37a}"); }
//        }
//    }
//}