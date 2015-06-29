using System;
using System.Collections.Generic;

using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;
using Rhino.Geometry;

namespace Pavillion2015
{
    public class GHC_04_Explode : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the GHC_04_Explode class.
        /// </summary>
        public GHC_04_Explode()
            : base("GHC_04_Explode", "04_Explode",
                "Explodes a Component into its elements",
                "Pavillion 2015", "Extract Data")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Components", "Components", "All instances of Tri-Loop-Components-Class", GH_ParamAccess.list);

        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddPointParameter("Centre Points L01", "Centers01", "Centre Points L01", GH_ParamAccess.tree); //0
            pManager.AddPointParameter("Centre Points L02", "Centers02", "Centre Points L02", GH_ParamAccess.tree); //1
            pManager.AddCurveParameter("Curves", "Curves", "Curves", GH_ParamAccess.tree); //2
            pManager.AddCurveParameter("PlateCurves", "PlateCrvs", "Curves", GH_ParamAccess.tree); //3
            pManager.AddPointParameter("Pentagon Points L01", "PentPtsL01", "planar pentagon like shape at L01", GH_ParamAccess.tree);//4
            pManager.AddPointParameter("Pentagon Points L02", "PentPtsL02", "planar pentagon like shape at L02", GH_ParamAccess.tree);//5
            pManager.AddIntegerParameter("Stripe Indentification", "StripeIdx", "Identification Numbers for each Stripe of a Component", GH_ParamAccess.tree);//6

            pManager.AddBrepParameter("Stipe Geometry", "StripeGeo", "Geometry for each Stripe of a Component", GH_ParamAccess.tree);//7
            pManager.AddPlaneParameter("Stipe Intersection Planes", "StripXPlanes", "Stipe planes for each Stripe of a Component", GH_ParamAccess.tree);//8


            //// =======================================================================================
            // Added by Gene 
            pManager.AddCurveParameter("ExtendedCurves", "ExtendedCurves", "ExtendedCurves", GH_ParamAccess.tree); //9
            pManager.AddBrepParameter("SingleStripeBreps", "SingleStripeBreps", "SingleStripeBreps", GH_ParamAccess.tree); //10
            //// =======================================================================================

        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            //---- Declareing ---------------------------------------------------------------------------

            List<Tri_Loop_Component> All_Components = new List<Tri_Loop_Component>();
            DA.GetDataList<Tri_Loop_Component>("Components", All_Components);


            DataTree<Curve> All_Curves = new DataTree<Curve>();

            //// =======================================================================================
            // Added by Gene 
            DataTree<Curve> All_ExtendedCurves = new DataTree<Curve>();
            DataTree<Brep> All_SingleStripeBreps = new DataTree<Brep>();

            //// =======================================================================================           

            DataTree<Point3d> All_CentrePointsL01 = new DataTree<Point3d>();
            DataTree<Point3d> All_CentrePointsL02 = new DataTree<Point3d>();
            DataTree<Curve> All_PlateCrv = new DataTree<Curve>();

            DataTree<Point3d> All_PentagonPts_L01 = new DataTree<Point3d>();
            DataTree<Point3d> All_PentagonPts_L02 = new DataTree<Point3d>();

            DataTree<int> All_StripeIds = new DataTree<int>();
            DataTree<Brep> All_StripeBreps = new DataTree<Brep>();
            DataTree<Curve> All_StripeIntersectCrvs = new DataTree<Curve>();
            DataTree<Plane> All_StripeIntersectPlane = new DataTree<Plane>();

            //---- End Declareing -----------------------------------------------------------------------
            //---- Functions ----------------------------------------------------------------------------

            int componentIdx = 0;

            foreach (Tri_Loop_Component component in All_Components)
            {
                GH_Path pth1 = new GH_Path(componentIdx);

                // run / generate Informations for Sofistik
                // create output information
                for (int p = 0; p < component.CentersL01.Count; p++)
                {
                    All_CentrePointsL01.Add(component.CentersL01[p], pth1);
                    All_CentrePointsL02.Add(component.CentersL02[p], pth1);
                }

                for (int c = 0; c < component.Curves.BranchCount; c++)
                {
                    All_Curves.Add(component.Curves.Branch(c)[0], pth1.AppendElement(c));
                    All_Curves.Add(component.Curves.Branch(c)[1], pth1.AppendElement(c));

                    //// =======================================================================================
                    // Added by Gene 
                    All_ExtendedCurves.Add(component.ExtendedCurves.Branch(c)[0], pth1.AppendElement(c));
                    All_ExtendedCurves.Add(component.ExtendedCurves.Branch(c)[1], pth1.AppendElement(c));
                    //// =======================================================================================
                }

                for (int sid = 0; sid < component.StripeID.BranchCount; sid++)
                {
                    All_StripeIds.Add(component.StripeID.Branch(sid)[0], pth1.AppendElement(sid));
                    All_StripeIds.Add(component.StripeID.Branch(sid)[1], pth1.AppendElement(sid));
                    All_StripeIds.Add(component.StripeID.Branch(sid)[2], pth1.AppendElement(sid));
                }

                for (int p = 0; p < component.PlateCrv.BranchCount; p++)
                {
                    All_PlateCrv.Add(component.PlateCrv.Branch(p)[0], pth1.AppendElement(p));
                    All_PlateCrv.Add(component.PlateCrv.Branch(p)[1], pth1.AppendElement(p));
                }

                for (int pent = 0; pent < component.PlanarPentagonL01.BranchCount; pent++)
                {
                    All_PentagonPts_L01.Add(component.PlanarPentagonL01.Branch(pent)[0], pth1.AppendElement(pent));
                    All_PentagonPts_L01.Add(component.PlanarPentagonL01.Branch(pent)[1], pth1.AppendElement(pent));
                    All_PentagonPts_L01.Add(component.PlanarPentagonL01.Branch(pent)[2], pth1.AppendElement(pent));
                    All_PentagonPts_L01.Add(component.PlanarPentagonL01.Branch(pent)[3], pth1.AppendElement(pent));
                    All_PentagonPts_L01.Add(component.PlanarPentagonL01.Branch(pent)[4], pth1.AppendElement(pent));
                }

                for (int pent = 0; pent < component.PlanarPentagonL02.BranchCount; pent++)
                {
                    All_PentagonPts_L02.Add(component.PlanarPentagonL02.Branch(pent)[0], pth1.AppendElement(pent));
                    All_PentagonPts_L02.Add(component.PlanarPentagonL02.Branch(pent)[1], pth1.AppendElement(pent));
                    All_PentagonPts_L02.Add(component.PlanarPentagonL02.Branch(pent)[2], pth1.AppendElement(pent));
                    All_PentagonPts_L02.Add(component.PlanarPentagonL02.Branch(pent)[3], pth1.AppendElement(pent));
                    All_PentagonPts_L02.Add(component.PlanarPentagonL02.Branch(pent)[4], pth1.AppendElement(pent));
                }

                // Information from Stripe Class
                for (int s = 0; s < component.Stripes.Count; s++)
                {
                    for (int t = 0; t < component.Stripes[s].SurfaceAtLoop.Length; t++)
                    { All_StripeBreps.Add(component.Stripes[s].SurfaceAtLoop[t], pth1.AppendElement(s)); }

                    //for (int u = 0; u < component.Stripes[s].IntersectionCurve.Length; u++)
                    //{ All_StripeIntersectCrvs.Add(component.Stripes[s].IntersectionCurve[u], pth1.AppendElement(s));}

                    All_StripeIntersectPlane.Add(component.Stripes[s].IntersectionPlane, pth1.AppendElement(s));
                }

                // ======================================================================================================
                // Gene Added

                for (int s = 0; s < component.StripesSingle.Count; s++)
                {
                    All_SingleStripeBreps.Add(component.StripesSingle[s].loop, pth1.AppendElement(s));

                    //for (int u = 0; u < component.Stripes[s].IntersectionCurve.Length; u++)
                    //{ All_StripeIntersectCrvs.Add(component.Stripes[s].IntersectionCurve[u], pth1.AppendElement(s));}

                    //All_StripeIntersectPlane.Add(component.StripesSingle[s].IntersectionPlane, pth1.AppendElement(s));
                }

                // ======================================================================================================

                componentIdx++;
            }

            //---- End Functions ------------------------------------------------------------------------
            //---- Set Output ---------------------------------------------------------------------------

            DA.SetDataTree(0, All_CentrePointsL01);
            DA.SetDataTree(1, All_CentrePointsL02);
            DA.SetDataTree(2, All_Curves);
            DA.SetDataTree(3, All_PlateCrv);
            DA.SetDataTree(4, All_PentagonPts_L01);
            DA.SetDataTree(5, All_PentagonPts_L02);
            DA.SetDataTree(6, All_StripeIds);
            DA.SetDataTree(7, All_StripeBreps);
            DA.SetDataTree(8, All_StripeIntersectPlane);

            //// =======================================================================================
            // Added by Gene 

            DA.SetDataTree(9, All_ExtendedCurves);
            DA.SetDataTree(10, All_SingleStripeBreps);
            //// =======================================================================================

            //---- End Set Output -----------------------------------------------------------------------

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
                return Properties.Resources.Explode;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("{3ffaa054-978d-49bf-a3ed-e948843e1c53}"); }
        }
    }
}