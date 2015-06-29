using System;
using System.Collections.Generic;

using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;
using Rhino.Geometry;


namespace Pavillion2015
{
    public class GHC_Unroller : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the GHC_Unroller class.
        /// </summary>
        public GHC_Unroller()
            : base("GHC_Unroller", "05_Unroller",
                "Unrolls Something",
                "Pavillion 2015", "Extract Data")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddBrepParameter("Triloop Stipes", "Stripes", "Input the 3 Stripes of each Triloop", GH_ParamAccess.tree);
            pManager.AddBooleanParameter("Merge Stripes", "Y", "Reorients unrolled surfaces to one connected surface", GH_ParamAccess.item);
            pManager.AddBooleanParameter("Switch", "Switch", "Switch", GH_ParamAccess.item);
            pManager.AddIntegerParameter("Seam", "Seam", "0: Seam on Inside, 1: Seam on Outside", GH_ParamAccess.item);
            pManager.AddPointParameter("Points", "P", "Points to unroll with Brep", GH_ParamAccess.tree);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddBrepParameter("Unrolled Breps", "UnrolledBreps", "Unrolled Breps", GH_ParamAccess.tree); //0
            pManager.AddPlaneParameter("Planes for Orientation", "OrientPlanes", "1st item origin plane, 2nd item target plane", GH_ParamAccess.tree); //1
            pManager.AddCurveParameter("Shared Curves", "SharedCrv", "Shared Curves between Stripes", GH_ParamAccess.tree); //2
            pManager.AddPointParameter("Points", "P", "Points unrolled with Brep", GH_ParamAccess.tree);//3
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            //---- Declareing ---------------------------------------------------------------------------
            GH_Structure<GH_Brep> AllStripes;
            DA.GetDataTree("Triloop Stipes", out AllStripes);

            GH_Structure<GH_Point> AllPoints;
            DA.GetDataTree("Points", out AllPoints);

            bool Reorient = false;
            DA.GetData<bool>("Merge Stripes", ref Reorient);

            bool Switch = false;
            DA.GetData<bool>("Switch", ref Switch);

            int Seam = 0;
            DA.GetData<int>("Seam", ref Seam);



            DataTree<Brep> AllUnrolledBreps = new DataTree<Brep>();
            DataTree<Brep> ForReorientBreps = new DataTree<Brep>();
            DataTree<Plane> AllOrientPlanes = new DataTree<Plane>();
            DataTree<Curve> AllSharedCurves = new DataTree<Curve>();
            DataTree<Point3d> AllUnrolledPoints = new DataTree<Point3d>();

            //---- End Declareing -----------------------------------------------------------------------
            //---- Functions ----------------------------------------------------------------------------

            #region Unroll

            for (int i = 0; i < AllStripes.Branches.Count; i++)
            {
                GH_Path pth = new GH_Path(i);
                GH_Path originalPath = AllStripes.Paths[i];
                int stripecounter = 0;

                foreach (GH_Brep gbrep in AllStripes[i])
                {
                    Unroller unroll = new Unroller(gbrep.Value);
                    // Add points to unroll with
                    if (AllPoints.Branches.Count != 0)
                    {
                        foreach (GH_Point pt in AllPoints[i])
                        { unroll.AddFollowingGeometry(pt.Value); }
                    }


                    unroll.ExplodeOutput = false;

                    Curve[] curves;
                    Point3d[] unrolledPoints;
                    TextDot[] dots;
                    Brep[] unrolledBreps = unroll.PerformUnroll(out curves, out unrolledPoints, out dots);

                    if (Reorient == false)
                    {
                        foreach (Brep b in unrolledBreps)
                        { AllUnrolledBreps.Add(b, originalPath); }

                        foreach (Point3d p in unrolledPoints)
                        { AllUnrolledPoints.Add(p, originalPath); }
                    }

                    else
                    {
                        foreach (Brep b in unrolledBreps)
                        { AllUnrolledBreps.Add(b, pth.AppendElement(stripecounter)); }
                    }

                    // For reorientation
                    if (Reorient == true)
                    { ForReorientBreps.Add(unrolledBreps[Seam], pth); }

                    stripecounter++;
                }
            }
            #endregion unroll

            if (Reorient == true)
            {
                //ForReorientBreps.Branch(0)[0].Curves3D[0].PointAtEnd;

                for (int i = 0; i < ForReorientBreps.BranchCount; i++)
                {
                    GH_Path pth = new GH_Path(i);

                    foreach (Curve crv0 in ForReorientBreps.Branch(i)[0].Curves3D)
                    {
                        foreach (Curve crv1 in ForReorientBreps.Branch(i)[1].Curves3D)
                        {
                            double l0 = crv0.GetLength();
                            double l1 = crv1.GetLength();

                            if (Math.Abs(l0 - l1) < 0.00001)
                            {

                                // orient crv0
                                Plane origin0 = new Plane(new Point3d(0, 0, 0), new Vector3d(1, 0, 0), new Vector3d(0, 1, 0));
                                Plane target0 = new Plane(origin0);
                                AllOrientPlanes.Add(origin0, pth.AppendElement(0));
                                AllOrientPlanes.Add(target0, pth.AppendElement(0));


                                // orient crv1
                                Vector3d vect0 = crv1.TangentAtStart;
                                Vector3d vect1 = Vector3d.CrossProduct(vect0, new Vector3d(0.0, 0.0, 1.0));
                                Plane origin1 = new Plane(crv1.PointAtStart, vect0, vect1);

                                Vector3d vect2 = new Vector3d();
                                Vector3d vect3 = new Vector3d();
                                Plane target1 = new Plane();

                                if (Switch == true)
                                {
                                    vect2 = crv0.TangentAtStart;
                                    vect3 = Vector3d.CrossProduct(vect2, new Vector3d(0.0, 0.0, 1.0));
                                    target1 = new Plane(crv0.PointAtStart, vect2, vect3);
                                }

                                else
                                {
                                    vect2 = crv0.TangentAtStart;
                                    vect3 = Vector3d.CrossProduct(-vect2, new Vector3d(0.0, 0.0, 1.0));
                                    target1 = new Plane(crv0.PointAtEnd, -vect2, vect3);
                                }


                                AllOrientPlanes.Add(origin1, pth.AppendElement(1));
                                AllOrientPlanes.Add(target1, pth.AppendElement(1));
                                // shared curve of stripe0 and stripe 1
                                AllSharedCurves.Add(crv0, pth.AppendElement(0));
                                AllSharedCurves.Add(crv1, pth.AppendElement(0));

                            }
                        }

                        // orient crv2
                        foreach (Curve crv2 in ForReorientBreps.Branch(i)[2].Curves3D)
                        {
                            double l0 = crv0.GetLength();
                            double l1 = crv2.GetLength();

                            if (Math.Abs(l0 - l1) < 0.00001)
                            {
                                Vector3d vect0 = crv2.TangentAtStart;
                                Vector3d vect1 = Vector3d.CrossProduct(vect0, new Vector3d(0.0, 0.0, 1.0));
                                Plane origin2 = new Plane(crv2.PointAtStart, vect0, vect1);

                                Vector3d vect2 = new Vector3d();
                                Vector3d vect3 = new Vector3d();
                                Plane target2 = new Plane();

                                if (Switch == true)
                                {
                                    vect2 = crv0.TangentAtStart;
                                    vect3 = Vector3d.CrossProduct(vect2, new Vector3d(0.0, 0.0, 1.0));
                                    target2 = new Plane(crv0.PointAtStart, vect2, vect3);
                                }

                                else
                                {
                                    vect2 = crv0.TangentAtStart;
                                    vect3 = Vector3d.CrossProduct(-vect2, new Vector3d(0.0, 0.0, 1.0));
                                    target2 = new Plane(crv0.PointAtEnd, -vect2, vect3);
                                }

                                AllOrientPlanes.Add(origin2, pth.AppendElement(2));
                                AllOrientPlanes.Add(target2, pth.AppendElement(2));
                                // shared curve of stripe0 and stripe 2
                                AllSharedCurves.Add(crv2, pth.AppendElement(2));
                                AllSharedCurves.Add(crv0, pth.AppendElement(2));
                            }
                        }

                    }
                    // find the shared curve oft stripe 1 and stripe 2
                    foreach (Curve crv1 in ForReorientBreps.Branch(i)[1].Curves3D)
                    {
                        foreach (Curve crv2 in ForReorientBreps.Branch(i)[2].Curves3D)
                        {
                            double l1 = crv1.GetLength();
                            double l2 = crv2.GetLength();

                            // shared curve of stripe1 and stripe 2
                            if (Math.Abs(l1 - l2) < 0.00001)
                            {
                                AllSharedCurves.Add(crv1, pth.AppendElement(1));
                                AllSharedCurves.Add(crv2, pth.AppendElement(1));
                            }

                        }

                    }
                }


            }



            //---- End Functions --------------------------------------------------------------------------
            //----Set Output-------------------------------------------------------------------------------

            DA.SetDataTree(0, AllUnrolledBreps);
            DA.SetDataTree(1, AllOrientPlanes);
            DA.SetDataTree(2, AllSharedCurves);
            DA.SetDataTree(3, AllUnrolledPoints);

            //----End Set Output---------------------------------------------------------------------------



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
                return Properties.Resources.Unroller;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("{89277a4a-4bac-4b81-b684-9c24a56a3f7c}"); }
        }
    }
}