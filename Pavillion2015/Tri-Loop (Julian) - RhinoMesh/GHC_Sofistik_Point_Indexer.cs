using System;
using System.Collections.Generic;

using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;
using Rhino.Geometry;

namespace Pavillion2015
{
    public class GHC_Sofistik_Point_Indexer : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the GHC_Sofistik_Point_Indexer class.
        /// </summary>
        public GHC_Sofistik_Point_Indexer()
            : base("GHC_Sofistik_Point_Indexer", "Do not Use",
                "Creates a tree of indexies",
                "Pavillion 2015", "DoubleLayer")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddPointParameter("SofistikPoints", "SofistikPts", "Tree of Points for Sofistik", GH_ParamAccess.tree);
            pManager.AddPointParameter("Flatten Points", "flatclean", "List of Points (doubles reduced)", GH_ParamAccess.list);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddIntegerParameter("Tree of indexies", "IndexTree", "Tree of indexies", GH_ParamAccess.tree);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            //---- Declareing -------------------------------------------------------------------------------
            
            GH_Structure<GH_Point> SofistikPts;
            DA.GetDataTree("SofistikPoints", out SofistikPts);

            List<Point3d> flatclean = new List<Point3d>();
            DA.GetDataList<Point3d>("Flatten Points", flatclean);

            DataTree<int> SofistikInt = new DataTree<int>();

            //---- End Declareing ---------------------------------------------------------------------------

            //---- Functions --------------------------------------------------------------------------------

            for (int i = 0; i < flatclean.Count; i++)
            {
                Point3d fpoint = flatclean[i];

                for (int j = 0; j < SofistikPts.Paths.Count; j++)
                {
                    GH_Path pth = SofistikPts.Paths[j];
                    for (int k = 0; k < SofistikPts[pth].Count; k++)
                    {
                        GH_Point ghp = SofistikPts[pth][k];

                        if ((Math.Abs(fpoint.X - ghp.Value.X) <= 0.01) &&
                                (Math.Abs(fpoint.Y - ghp.Value.Y) <= 0.01) &&
                                (Math.Abs(fpoint.Z - ghp.Value.Z) <= 0.01))
                        {
                            SofistikInt.Add(i, pth.AppendElement(k));
                        }
                    }
                       
                }
            }

            
            //---- End Functions ----------------------------------------------------------------------------
            
            //---- Set Output -----------------------------------------------------------------------------

            DA.SetDataTree(0, SofistikInt);


            //---- End Set Output -----------------------------------------------------------------------------

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
            get { return new Guid("{2dff1906-a5f8-46ed-adc8-3eb6f2d33ca2}"); }
        }
    }
}