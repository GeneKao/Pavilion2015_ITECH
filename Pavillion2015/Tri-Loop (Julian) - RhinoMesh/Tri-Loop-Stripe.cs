using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

using Grasshopper;
using Grasshopper.Kernel.Data;
using Rhino.Geometry;
using Rhino.Geometry.Intersect;

namespace Pavillion2015
{
    class Tri_Loop_Stripe
    {

        // top curves
        public List<Point3d> PointsAtTop = new List<Point3d>();
        public List<Curve> CurvesAtTop = new List<Curve>();
        // bottom curves
        public List<Point3d> PointsAtBottom = new List<Point3d>();
        public List<Curve> CurvesAtBottom = new List<Curve>();
        // loop curves
        public DataTree<Point3d> PointsAtLoop = new DataTree<Point3d>();
        public List<Curve> CurvesAtLoop = new List<Curve>();
        // top plate
        public List<Point3d> PointsAtTopPlate = new List<Point3d>();
        public List<Curve> CurvesAtTopPlate = new List<Curve>();
        // bottom plate
        public List<Point3d> PointsAtBottomPlate = new List<Point3d>();
        public List<Curve> CurvesAtBottomPlate = new List<Curve>();

        // ====================================================================
        // Added by Gene

        public Curve curvesLeft;
        public Curve curvesRight;

        List<Point3d> cuttingPts;

        public double documentTolerance;

        public Brep loop;
        public Brep finalLoop;

        // ====================================================================

        // planar top Surfaces
        public Brep[] SurfaceAtTop = null;
        // planar bottom Surfaces
        public Brep[] SurfaceAtBottom = null;
        // curved Surfaces
        public Brep[] SurfaceAtLoop = null;
        // top plate surfaces
        public Brep[] SurfaceAtTopPlate = null;
        // bottom plate surfaces
        public Brep[] SurfaceAtBottomPlate = null;

        //StripeID
        public List<int> StripeID = new List<int>();

        //Intersection Plane
        public Plane IntersectionPlane;
        public Curve[] IntersectionCurve;
        public Point3d[] IntersectionPoint;

        //Material BuildUp
        public List<string> BuildUp;
        public List<Point3d> UVPoints;


        //---- CONSTRUCTOR-----------------------------------------------------
        public Tri_Loop_Stripe(
            List<Point3d> _pointsAtTop, List<Point3d> _pointsAtBottom,
            List<Point3d> _pointsAtTopPlate, List<Point3d> _pointsAtBottomPlate,
            List<Point3d> _pointsAtLoopLeft, List<Point3d> _pointsAtLoopRight,
            List<int> _stripeID
            )
        {
            PointsAtTop = _pointsAtTop;
            PointsAtBottom = _pointsAtBottom;
            PointsAtTopPlate = _pointsAtTopPlate;
            PointsAtBottomPlate = _pointsAtBottomPlate;
            StripeID = _stripeID;

            foreach (Point3d p in _pointsAtLoopLeft)
            { PointsAtLoop.Add(p, new GH_Path(0)); }

            foreach (Point3d p in _pointsAtLoopRight)
            { PointsAtLoop.Add(p, new GH_Path(1)); }



            // IntersectionPlane
            #region CreateIntersectionPlane

            Point3d c0 = (PointsAtTop[0] + PointsAtTop[PointsAtTop.Count - 1]) * 0.5;
            Vector3d v0 = (PointsAtTop[PointsAtTop.Count - 1] - PointsAtTop[0]);
            Vector3d v1 = (PointsAtBottom[PointsAtBottom.Count - 1] - PointsAtBottom[0]);
            Vector3d normal = (v0 + v1) * 0.5;

            IntersectionPlane = new Plane(c0, normal);
            #endregion CreateIntersectionPlane


            // Run Functions
            Run();
        }

        // ====================================================================
        // Added by Gene

        public Tri_Loop_Stripe(Curve curvesLeft, Curve curvesRight, List<Point3d> cuttingPts, double documentTolerance)
        {
            this.curvesLeft = curvesLeft;
            this.curvesRight = curvesRight;
            this.cuttingPts = cuttingPts;
            this.documentTolerance = documentTolerance;

            RunSingleSrf();
        }

        // ====================================================================

        //---- END CONSTRUCTOR-------------------------------------------------

        // ====================================================================
        // Added by Gene

        public void RunSingleSrf()
        {
            Point3d pLeft = curvesLeft.PointAtStart;
            Point3d pRight = curvesRight.PointAtStart;
            LineCurve rail = new LineCurve(pLeft, pRight);

            Point3d pLeft2 = curvesLeft.PointAtEnd;
            Point3d pRight2 = curvesRight.PointAtEnd;
            LineCurve rail2 = new LineCurve(pLeft2, pRight2);

            loop = Brep.CreateFromSweep(curvesLeft, curvesRight, new List<Curve>() { rail, rail2 }, false, documentTolerance)[0];
            //loop = Brep.CreateFromLoftRebuild(new List<Curve>() { curvesLeft, curvesRight }, Point3d.Unset, Point3d.Unset, LoftType.Normal, false, 50)[0];
            //loop = Brep.CreateFromLoft(new List<Curve>() { curvesLeft, curvesRight }, Point3d.Unset, Point3d.Unset, LoftType.Normal, false )[0];

            Brep[] loops;
            double offsetAmount = 0.002;
            //double offsetAmount2 = 0.02;
            Vector3d normal1, normal2;

            normal1 = Vector3d.CrossProduct(cuttingPts[0] - cuttingPts[4], cuttingPts[0] - cuttingPts[2]);
            normal2 = Vector3d.CrossProduct(cuttingPts[1] - cuttingPts[3], cuttingPts[1] - cuttingPts[5]);

            Point3d A_ = cuttingPts[7] + offsetAmount * normal2;
            Point3d B_ = cuttingPts[6] + offsetAmount * normal1;

            //Brep cutter1 = Brep.CreateFromCornerPoints(cuttingPts[0], cuttingPts[1], cuttingPts[2], cuttingPts[3], documentTolerance);
            //Brep cutter2 = Brep.CreateFromCornerPoints(cuttingPts[0], cuttingPts[1], cuttingPts[4], cuttingPts[5], documentTolerance);

            loops = loop.Trim(new Plane(A_, cuttingPts[1], cuttingPts[3]), documentTolerance);
            //loops = loop.Trim(cutter1, documentTolerance);

            if (loops.Length > 0) loop = loops[0];
            //else { return; }

            loops = loop.Trim(new Plane(A_, cuttingPts[5], cuttingPts[1]), documentTolerance);

            if (loops.Length > 0) loop = loops[0];
            //else { return; }

            loops = loop.Trim(new Plane(B_, cuttingPts[2], cuttingPts[0]), documentTolerance);
            if (loops.Length > 0) loop = loops[0];

            loops = loop.Trim(new Plane(B_, cuttingPts[0], cuttingPts[4]), documentTolerance);
            if (loops.Length > 0) loop = loops[0];

        }

        // ====================================================================

        public void Run()
        {
            CreateCurves();
            CreateBreps();
            AnalyseCurvature();
        }

        public void CreateCurves()
        {
            //Curve At Top
            #region PointsAtTop

            List<Point3d> _PointsAtTop = new List<Point3d>(PointsAtTop);
            _PointsAtTop.Add(PointsAtTop[0]);
            CurvesAtTop.Add(Curve.CreateControlPointCurve(_PointsAtTop));

            #endregion PointsAtTop

            //Curve AT Bottom
            #region PointsAtBottom

            List<Point3d> _PointsAtBottom = new List<Point3d>(PointsAtBottom);
            _PointsAtBottom.Add(PointsAtBottom[0]);
            CurvesAtBottom.Add(Curve.CreateControlPointCurve(_PointsAtBottom));

            #endregion PointsAtBottom

            //Curve AT Top Plate
            #region PointsAtTopPlate
            if (PointsAtTopPlate != null)
            {
                List<Point3d> _PointsAtTopPlate = new List<Point3d>(PointsAtTopPlate);
                _PointsAtTopPlate.Add(PointsAtTopPlate[0]);
                CurvesAtTopPlate.Add(Curve.CreateControlPointCurve(_PointsAtTopPlate));
            }
            #endregion PointsAtTopPlate

            //Curves At Bottom Plate
            #region PointsAtBottomPlate
            if (PointsAtBottomPlate != null)
            {
                List<Point3d> _PointsAtBottomPlate = new List<Point3d>(PointsAtBottomPlate);
                _PointsAtBottomPlate.Add(PointsAtBottomPlate[0]);
                CurvesAtBottomPlate.Add(Curve.CreateControlPointCurve(_PointsAtBottomPlate));
            }
            #endregion PointsAtBottomPlate

            //Curves At Loop
            #region PointsAtLoop

            for (int i = 0; i < PointsAtLoop.BranchCount; i++)
            {
                CurvesAtLoop.Add(Curve.CreateControlPointCurve(PointsAtLoop.Branch(i), 3));
            }

            #endregion PointsAtLoop

        }

        public void CreateBreps()
        {
            // Surface At Top
            SurfaceAtTop = Brep.CreatePlanarBreps(CurvesAtTop);

            // Surface At Bottom
            SurfaceAtBottom = Brep.CreatePlanarBreps(CurvesAtBottom);

            // Surface At Top Plate
            if (PointsAtTopPlate != null)
            { SurfaceAtTopPlate = Brep.CreatePlanarBreps(CurvesAtTopPlate); }

            // Surface At Bottom Plate
            if (PointsAtBottomPlate != null)
            { SurfaceAtBottomPlate = Brep.CreatePlanarBreps(CurvesAtBottomPlate); }

            // Surface At Top Plate
            SurfaceAtLoop = Brep.CreateFromLoft(CurvesAtLoop, Point3d.Unset, Point3d.Unset, LoftType.Normal, false);

        }


        public void AnalyseCurvature()
        {
            // create IntersectionCurve[] and IntersectionPoint[]
            Intersection.BrepPlane(SurfaceAtLoop[0], IntersectionPlane, 0.001, out IntersectionCurve, out IntersectionPoint);


        }

        public void SetBuildUp(List<string> buildUp, List<Point3d> uvPoints)
        {
            BuildUp = buildUp;
            UVPoints = uvPoints;
        }


    }
}
