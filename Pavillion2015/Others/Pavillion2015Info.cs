using System;
using System.Drawing;
using Grasshopper.Kernel;

namespace Pavillion2015
{
    public class Pavillion2015Info : GH_AssemblyInfo
    {
        public override string Name
        {
            get
            {
                return "Pavillion 2015";
            }
        }
        public override Bitmap Icon
        {
            get { return null; }
        }
        public override string Description
        {
            get
            {
                //Return a short string describing the purpose of this GHA library.
                return @"ICD/ITKE Research Pavillion 2015";
            }
        }
        public override Guid Id
        {
            get
            {
                return new Guid("77def88f-829a-4e46-a365-25439455f0ac");
            }
        }

        public override string AuthorName
        {
            get
            {
                //Return a string identifying you or your company.
                return "Institute for Computational Design (ICD), University of Stuttgart";
            }
        }
        public override string AuthorContact
        {
            get
            {
                //Return a string representing your preferred contact details.
                return "long.nguyen@icd.uni-stuttgart.de";
            }
        }
    }
}
