using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;

namespace Pavillion2015
{
    public class GHC_SpringMeshToRhinoMesh : GH_Component
    {
        public GHC_SpringMeshToRhinoMesh()
            : base("Spring Mesh To Rhino Mesh", "Rhino Mesh",
                "Convert Spring Mesh to Rhino Mesh",
                "Pavillion 2015", "Utilities")
        {
        }


        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Spring Mesh", "Spring Mesh", "Spring Mesh", GH_ParamAccess.item);
        }


        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Rhino Mesh", "Rhino Mesh", "Rhino Mesh", GH_ParamAccess.item);
        }


        protected override void SolveInstance(IGH_DataAccess DA)
        {
            SpringMesh iSpringMesh = null;
            DA.GetData<SpringMesh>("Spring Mesh", ref iSpringMesh);

            DA.SetData("Rhino Mesh", iSpringMesh.ConvertToRhinoMesh());
        }


        protected override System.Drawing.Bitmap Icon{get{return null;}}


        public override Guid ComponentGuid{ get { return new Guid("{1cc33dce-6ce8-46d0-bbde-f0ea7995f03b}"); }}
    }
}