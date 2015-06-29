using System;
using System.Collections.Generic;

using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;
using Rhino.Geometry;

namespace Pavillion2015
{
    public class GHC_02_DoubleLayer : GH_Component
    {

        // ==================================================================================================
        // Gene Added 

        double documentTolerance = DocumentTolerance();

        bool iAllowPolySrf = true;

        // ==================================================================================================

        public GHC_02_DoubleLayer()
            : base("GHC_02_DoubleLayer", "02_DoubleLayer",
                "Creates the double layer",
                "Pavillion 2015", "DoubleLayer")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddPointParameter("Layer 01 Points", "PtL01", "[Tree<Point3d>input tree of points of layer 01", GH_ParamAccess.tree);
            pManager.AddPointParameter("Layer 02 Points", "PtL02", "[Tree<Point3d>input tree of points of layer 02", GH_ParamAccess.tree);
            pManager.AddIntegerParameter("Topology Points", "TopoPtIdx", "[Tree<int>]Topology points index of the face", GH_ParamAccess.tree);

            pManager.AddIntegerParameter("Neighbour Face Index", "NFIdx", "Neighbouring faces", GH_ParamAccess.tree);

            pManager.AddIntegerParameter("Face Topology Edge Indexies", "FaceTopEdgeIdx", "Indexies of topology edges foreach face", GH_ParamAccess.tree);
            pManager.AddPointParameter("Edge Centre L01", "EdgeCentreL01", "[Tree<Point3d>] Edge Centre Points foreach Face, L01", GH_ParamAccess.tree);
            pManager.AddPointParameter("Edge Centre L02", "EdgeCentreL02", "[Tree<Point3d>] Edge Centre Points foreach Face, L02", GH_ParamAccess.tree);

            pManager.AddNumberParameter("Curve Control Point Tangent Scalar", "CtrlPScalar", "[List<double>]Controls the size of the wholes per topological point", GH_ParamAccess.list);
            pManager.AddNumberParameter("Planar Pentagon Point Tangent Scalar", "PentPScalar", "[List<double>]Controls the size of the wholes per topological point", GH_ParamAccess.list);
            pManager.AddBooleanParameter("Plate?", "O/C", "[List<boolean>]Is the module closed or open", GH_ParamAccess.list);

            // ==================================================================================================
            // Gene Added 
            pManager.AddBooleanParameter("Allow PolySrf", "Allow PolySrf", "Allow PolySrf", GH_ParamAccess.item, true);

            // ==================================================================================================


        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Components", "Components", "All instances of Tri-Loop-Components-Class", GH_ParamAccess.list);//0
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            //---- Declareing ---------------------------------------------------------------------------
            GH_Structure<GH_Point> PtL01;
            DA.GetDataTree("Layer 01 Points", out PtL01);

            GH_Structure<GH_Point> PtL02;
            DA.GetDataTree("Layer 02 Points", out PtL02);

            GH_Structure<GH_Integer> TopoPtIdx;
            DA.GetDataTree("Topology Points", out TopoPtIdx);

            GH_Structure<GH_Integer> NFIdx;
            DA.GetDataTree("Neighbour Face Index", out NFIdx);

            GH_Structure<GH_Integer> All_TopoEdges;
            DA.GetDataTree("Face Topology Edge Indexies", out All_TopoEdges);

            GH_Structure<GH_Point> EdgeCentreL01;
            DA.GetDataTree("Edge Centre L01", out EdgeCentreL01);

            GH_Structure<GH_Point> EdgeCentreL02;
            DA.GetDataTree("Edge Centre L02", out EdgeCentreL02);


            List<double> CtrlPSc = new List<double>();
            DA.GetDataList<double>("Curve Control Point Tangent Scalar", CtrlPSc);

            List<double> PentPSc = new List<double>();
            DA.GetDataList<double>("Planar Pentagon Point Tangent Scalar", PentPSc);

            List<bool> openclosed = new List<bool>();
            DA.GetDataList<bool>("Plate?", openclosed);

            // ==================================================================================================
            // Gene Added 

            DA.GetData<bool>("Allow PolySrf", ref iAllowPolySrf);
            // ==================================================================================================

            List<Tri_Loop_Component> All_Components = new List<Tri_Loop_Component>();

            //---- End Declareing -----------------------------------------------------------------------

            //---- Functions ----------------------------------------------------------------------------

            for (int i = 0; i < PtL01.Branches.Count; i++)
            {
                // List 01
                List<Point3d> _ptL01 = new List<Point3d>();
                foreach (GH_Point pt in PtL01[i])
                    _ptL01.Add(pt.Value);

                // List 02
                List<Point3d> _ptL02 = new List<Point3d>();
                foreach (GH_Point pt in PtL02[i])
                    _ptL02.Add(pt.Value);

                // List TopoPtIdx
                List<int> _topoVidx = new List<int>();
                foreach (GH_Integer integer in TopoPtIdx[i])
                    _topoVidx.Add(integer.Value);

                // List TopoPtIdx
                List<int> _nFIdx = new List<int>();
                foreach (GH_Integer integer in NFIdx[i])
                    _nFIdx.Add(integer.Value);

                // List TopoEdges
                List<int> _topoEdges = new List<int>();
                foreach (GH_Integer integer in All_TopoEdges[i])
                    _topoEdges.Add(integer.Value);

                // EdgeCentreL01
                List<Point3d> _edgeCentreL01 = new List<Point3d>();
                foreach (GH_Point cntr in EdgeCentreL01[i])
                    _edgeCentreL01.Add(cntr.Value);

                // EdgeCentreL02
                List<Point3d> _edgeCentreL02 = new List<Point3d>();
                foreach (GH_Point cntr in EdgeCentreL02[i])
                    _edgeCentreL02.Add(cntr.Value);

                // Create Components
                //Tri_Loop_Component component = new Tri_Loop_Component(_ptL01, _ptL02, _topoVidx, _nFIdx, _topoEdges, _edgeCentreL01, _edgeCentreL02, i);

                // ======================================================================================================================
                // Gene Added with documentTolerance
                Tri_Loop_Component component = new Tri_Loop_Component(_ptL01, _ptL02, _topoVidx, _nFIdx, _topoEdges, _edgeCentreL01, _edgeCentreL02, i, documentTolerance);
                // ======================================================================================================================



                component.Run(PentPSc, CtrlPSc, openclosed);




                // ======================================================================================================================
                // Gene Added with documentTolerance
                if (!iAllowPolySrf)
                    component.RunSingle();
                // ======================================================================================================================

                // add to component list
                All_Components.Add(component);
            }
            //---- End Functions --------------------------------------------------------------------------

            //---- Set Output -----------------------------------------------------------------------------

            DA.SetDataList(0, All_Components);

            //---- End Set Output -------------------------------------------------------------------------

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
                return Properties.Resources.DoubleLayer_Icon;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("{b672aa5c-2bf6-4bc7-9cff-2c69e041a81c}"); }
        }
    }
}