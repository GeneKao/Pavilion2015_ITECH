using System;
using System.Collections.Generic;

using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;
using Rhino.Geometry;

namespace Pavillion2015
{
    public class GHC_03_Sofistik : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the GHC_03_Sofistik class.
        /// </summary>
        public GHC_03_Sofistik()
            : base("GHC_03_Sofistik", "03_Sofistik",
                "Extract Data for Sofistik",
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
            pManager.AddPointParameter("Sofistik Points", "Points", "Points for Sofistik", GH_ParamAccess.tree);//0
            pManager.AddIntegerParameter("Sofistik Indexies", "PointIdxs", "Indexies for Sofistik", GH_ParamAccess.tree);//1

            pManager.AddPointParameter("Sofistik Curve Points", "CrvPoints", "Curve Points for Sofistik", GH_ParamAccess.tree);//2
            pManager.AddIntegerParameter("Sofistik Curve Indexies", "CrvPointIdxs", "Curve Indexies for Sofistik", GH_ParamAccess.tree);//3

            pManager.AddBrepParameter("Sofistik Breps", "Breps", "Surfaces from each Component", GH_ParamAccess.tree); //4
            pManager.AddIntegerParameter("Sofistik Breps Indexies", "BrepsIdx", "Breps Indexies for Sofistik", GH_ParamAccess.tree);//5
            pManager.AddCurveParameter("Sofistik Curves", "Curves", "All Curves for Sofistik", GH_ParamAccess.tree);//6

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

            DataTree<Point3d> SofistikPoints = new DataTree<Point3d>();
            DataTree<int> SofistikIndexies = new DataTree<int>();


            DataTree<int> SofistikCrvPtIdx = new DataTree<int>();
            DataTree<Point3d> SofistikCrvPtCoo = new DataTree<Point3d>();

            DataTree<Brep> SofistikBreps = new DataTree<Brep>();
            DataTree<int> SofistikBrepsIdx = new DataTree<int>();

            DataTree<Curve> SofistikCurves = new DataTree<Curve>();

            //---- End Declareing -----------------------------------------------------------------------

            //---- Functions ----------------------------------------------------------------------------

            int componentIdx = 0;

            foreach (Tri_Loop_Component component in All_Components)
            {
                GH_Path pth1 = new GH_Path(componentIdx);

                // run / generate Informations for Sofistik
                component.SofistikInformation();
                component.SofistikCreateSurfaces();

                // Points
                for (int j = 0; j < component.SofistikPlatePoints.BranchCount; j++)
                {

                    foreach (Point3d p in component.SofistikPlatePoints.Branch(j))
                    { SofistikPoints.Add(p, pth1.AppendElement(j)); }

                    foreach (int idx in component.SofistikPlateIndexies.Branch(j))
                    { SofistikIndexies.Add(idx, pth1.AppendElement(j)); }
                }

                // CurvePoints
                for (int j = 0; j < component.SofistikCrvPtCoo.BranchCount; j++)
                {
                    GH_Path pth0 = component.SofistikCrvPtCoo.Path(j);

                    foreach (Point3d pt in component.SofistikCrvPtCoo.Branch(j))
                    {
                        //pth1.AppendElement(pth0[0]);
                        //pth1.AppendElement(pth0[1]);
                        //SofistikCrvPtCoo.Add(pt, pth1.AppendElement(pth0[0]));

                        // 2 level Tree
                        SofistikCrvPtCoo.Add(pt, pth1.AppendElement(pth0[0]));

                    }

                    foreach (int idx in component.SofistikCrvPtIdx.Branch(j))
                    {
                        //pth1.AppendElement(pth0[0]);
                        //pth1.AppendElement(pth0[1]);
                        //SofistikCrvPtIdx.Add(idx, pth1.AppendElement(pth0[0]));

                        // 2 level Tree
                        SofistikCrvPtIdx.Add(idx, pth1.AppendElement(pth0[0]));

                    }
                }

                // Curves
                for (int i = 0; i < component.SofistikCurves.BranchCount; i++)
                {
                    GH_Path pth0 = component.SofistikCurves.Path(i);

                    foreach (Curve c in component.SofistikCurves.Branch(i))
                    {
                        SofistikCurves.Add(c, pth1.AppendElement(pth0[0]));
                    }

                }

                // Surfaces
                for (int j = 0; j < component.SofistikSurfaces.BranchCount; j++)
                {

                    foreach (Brep[] p in component.SofistikSurfaces.Branch(j))
                    { SofistikBreps.Add(p[0], pth1.AppendElement(j)); }

                    foreach (int idx in component.SofistikSurfacesIdx.Branch(j))
                    { SofistikBrepsIdx.Add(idx, pth1.AppendElement(j)); }
                }

                componentIdx++;
            }

            //---- End Functions ------------------------------------------------------------------------
            //---- Set Output ---------------------------------------------------------------------------

            DA.SetDataTree(0, SofistikPoints);
            DA.SetDataTree(1, SofistikIndexies);

            DA.SetDataTree(2, SofistikCrvPtCoo);
            DA.SetDataTree(3, SofistikCrvPtIdx);

            DA.SetDataTree(4, SofistikBreps);
            DA.SetDataTree(5, SofistikBrepsIdx);

            DA.SetDataTree(6, SofistikCurves);



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
                return Properties.Resources.Sofistik;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("{fd65da52-4037-4c55-9e88-79389efef229}"); }
        }
    }
}