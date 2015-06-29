using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

using ICD;

namespace ICD.Beagle
{
    class Gene
    {
        public double Min { get; set; }
        public double Max { get; set; }

        public double Value{ get; set; }

        public Gene()
        {
            this.Min = double.NaN;
            this.Max = double.NaN;
            this.Value = double.NaN;
        }

        public Gene(double min, double max)
        {
            this.Min = min;
            this.Max = max;
            this.Value = Utils.GetRandomDouble(min, max);
        }

        public Gene(double min, double max, double value)
        {
            this.Min = min;
            this.Max = max;
            this.Value = value;
        }

        public void MakeRandom()
        {
            Value = Utils.GetRandomDouble(Min, Max);
        }

        public void MakeValueValid()
        {
            Value = (Value < Min ? Min : (Value > Max ? Max : Value));
        }
    }
}
