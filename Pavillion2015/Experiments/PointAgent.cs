using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

using Rhino.Geometry;

namespace Pavillion2015
{
    class PointAgent
    {
        public Point3d Position = Point3d.Origin;
        public double Separation = 0.0;
        public int age = 0;

        public PointAgent(Point3d position, double separation)
        {
            Position = position;
            Separation = separation;
        }

        public PointAgent(Circle circle)
        {
            Position = circle.Center;
            Separation = circle.Radius;
        }
    }
}
