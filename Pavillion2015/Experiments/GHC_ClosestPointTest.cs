using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;

using ICD;

namespace Pavillion2015.Experiments
{
    public class GHC_ClosestPointTest : GH_Component
    {        
        public GHC_ClosestPointTest()
            : base("GHC_ClosestPointTest", "Nickname",
                "Description",
                "Pavillion 2015", "Debug")
        {
        }

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddPointParameter("M", "M", "M", GH_ParamAccess.item);
            pManager.AddPointParameter("A", "A", "A", GH_ParamAccess.item);
            pManager.AddPointParameter("B", "B", "B", GH_ParamAccess.item);
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddPointParameter("N1", "N1", "N1", GH_ParamAccess.item);
            pManager.AddPointParameter("N2", "N2", "N2", GH_ParamAccess.item);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            Point3d M = Point3d.Origin;
            Point3d A = Point3d.Origin;
            Point3d B = Point3d.Origin;

            DA.GetData<Point3d>(0, ref M);
            DA.GetData<Point3d>(1, ref A);
            DA.GetData<Point3d>(2, ref B);

            LineCurve AB = new LineCurve(A, B);
            double t;
            AB.ClosestPoint(M, out t);

            DA.SetData(0, Utils.ClosestPointOnLine(M, A, B));
            DA.SetData(1, AB.PointAt(t));
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
            get { return new Guid("{73117719-ffa4-4ae1-bce1-d7a910022186}"); }
        }
    }
}