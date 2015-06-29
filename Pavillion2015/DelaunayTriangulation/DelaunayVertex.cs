using System;
using System.Collections.Generic;
using System.Text;

/*
  copyright s-hull.org 2011
  released under the contributors beerware license

  contributors: Phil Atkin, Dr Sinclair.
*/
namespace Pavillion2015.DelaunayTriangulation
{
    public class DelaunayVertex
    {
        public double x, y;

        protected DelaunayVertex() { }

        public DelaunayVertex(double x, double y) 
        {
            this.x = x; this.y = y;
        }

        public double distance2To(DelaunayVertex other)
        {
            double dx = x - other.x;
            double dy = y - other.y;
            return dx * dx + dy * dy;
        }

        public double distanceTo(DelaunayVertex other)
        {
            return (double)Math.Sqrt(distance2To(other));
        }

        public override string ToString()
        {
            return string.Format("({0},{1})", x, y);
        }
    }

}
