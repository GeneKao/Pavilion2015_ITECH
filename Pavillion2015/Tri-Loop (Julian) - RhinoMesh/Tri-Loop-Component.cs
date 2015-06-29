using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

using Grasshopper;
using Grasshopper.Kernel.Data;
using Rhino.Geometry;

namespace Pavillion2015
{
    class Tri_Loop_Component
    {
        // Input parameters
        public List<Point3d> PointsL01 = new List<Point3d>();
        public List<Point3d> PointsL02 = new List<Point3d>();
        public List<int> TopoPts = new List<int>();
        public List<int> TopoEdges = new List<int>();
        public List<Point3d> EdgeCentreL01 = new List<Point3d>();
        public List<Point3d> EdgeCentreL02 = new List<Point3d>();

        public List<int> NFIdx = new List<int>();
        public int ID = -1;

        // calculated parameters
        public List<Point3d> CentersL01 = new List<Point3d>();
        public List<Point3d> CentersL02 = new List<Point3d>();

        public Point3d AreaCentreL01 = new Point3d();
        public Point3d AreaCentreL02 = new Point3d();

        public DataTree<Point3d> PlanarPentagonL01 = new DataTree<Point3d>();       // contains all points for planar part L01
        public DataTree<Point3d> PlanarPentagonL02 = new DataTree<Point3d>();       // contains all points for planar part L00
        public DataTree<Curve> Curves = new DataTree<Curve>();                      // contains all the curves at loop
        public DataTree<List<Point3d>> CurvesPts = new DataTree<List<Point3d>>();   // contains all control points of curves



        ////=======================================================
        // Gene Added
        public DataTree<Curve> ExtendedCurves = new DataTree<Curve>();              // contains all the extended curves at loop to trim
        public DataTree<List<Point3d>> cuttingPoints = new DataTree<List<Point3d>>();
        // contains all cutting points
        public double documentTolerance;

        public List<Tri_Loop_Stripe> StripesSingle = new List<Tri_Loop_Stripe>();

        ////=======================================================

        public DataTree<Curve> PlateCrv = new DataTree<Curve>();                    // contains all curves for plates 
        public DataTree<List<Point3d>> PlateCrvPts = new DataTree<List<Point3d>>(); // contains all control points for ülates
        public DataTree<int> StripeID = new DataTree<int>();                        // contains all the StripeIds (Stripe A, Stripe B, Stripe C)

        public List<Tri_Loop_Stripe> Stripes = new List<Tri_Loop_Stripe>();         // [0] Stripe A; [1] Stipe B; [2] Stripe C


        //Output for Sofistik
        #region Output for Sofistik
        //---------------------------------------------------------------------
        public DataTree<Point3d> AllSofistikPoints = new DataTree<Point3d>();
        public DataTree<int> SofistikIndexies = new DataTree<int>();
        public List<Point3d> LeftfromCenterL01 = new List<Point3d>();
        public List<Point3d> RightfromCenterL01 = new List<Point3d>();
        public List<Point3d> LeftfromCenterL02 = new List<Point3d>();
        public List<Point3d> RightfromCenterL02 = new List<Point3d>();

        public DataTree<int> SofistikCrvPtIdx = new DataTree<int>();
        public DataTree<Point3d> SofistikCrvPtCoo = new DataTree<Point3d>();

        public DataTree<Curve> SofistikCurves = new DataTree<Curve>();

        public DataTree<Brep[]> SofistikSurfaces = new DataTree<Brep[]>();
        public DataTree<int> SofistikSurfacesIdx = new DataTree<int>();

        public DataTree<Point3d> SofistikPlatePoints = new DataTree<Point3d>();
        public DataTree<int> SofistikPlateIndexies = new DataTree<int>();
        //---------------------------------------------------------------------
        #endregion Output for Sofistik

        // Constructor
        public Tri_Loop_Component(
            List<Point3d> pointsL01, List<Point3d> pointsL02,
            List<int> topoPts, List<int> nFIdx, List<int> topoEdges,
            List<Point3d> edgeCentreL01, List<Point3d> edgeCentreL02,
            int id)
        {
            PointsL01 = pointsL01;
            PointsL02 = pointsL02;
            TopoPts = topoPts;
            TopoEdges = topoEdges;
            EdgeCentreL01 = edgeCentreL01;
            EdgeCentreL02 = edgeCentreL02;
            NFIdx = nFIdx;
            ID = id;
        }

        // ================================================================================================
        // Gene Added Constructor with documentTolerance
        public Tri_Loop_Component(
            List<Point3d> pointsL01, List<Point3d> pointsL02,
            List<int> topoPts, List<int> nFIdx, List<int> topoEdges,
            List<Point3d> edgeCentreL01, List<Point3d> edgeCentreL02,
            int id, double documentTolerance)
            : this(pointsL01, pointsL02, topoPts, nFIdx, topoEdges, edgeCentreL01, edgeCentreL02, id)
        {
            this.documentTolerance = documentTolerance;
        }
        // ================================================================================================

        public void Run(List<double> _offsetvalue, List<double> _controlpointvalue, List<bool> _openclosed)
        {
            Find_Midpoints();
            Create_Curves(_offsetvalue, _controlpointvalue);
            Create_Plate(_openclosed);
            Create_Stripe();
        }

        // =====================================
        // Added by Gene
        public void RunSingle()
        {
            Create_Stipe_SinSrf();
        }
        // =====================================

        // Geometry
        public void Find_Midpoints()
        {
            //---- calculate Midpoints
            // Layer 01
            CentersL01.Add(EdgeCentreL01[1]);
            CentersL01.Add(EdgeCentreL01[2]);
            CentersL01.Add(EdgeCentreL01[0]);


            // Layer 02
            CentersL02.Add(EdgeCentreL02[1]);
            CentersL02.Add(EdgeCentreL02[2]);
            CentersL02.Add(EdgeCentreL02[0]);


            //---- calculate AreaCentres
            // Layer 01
            AreaCentreL01 = (CentersL01[0] + CentersL01[1] + CentersL01[2]) / 3;
            // Layer 02
            AreaCentreL02 = (CentersL02[0] + CentersL02[1] + CentersL02[2]) / 3;

        }

        public void Create_Curves(List<double> offsetvalue, List<double> controlpointvalue)
        {
            // Data tree path
            GH_Path pthA = new GH_Path(0);
            GH_Path pthB = new GH_Path(1);
            GH_Path pthC = new GH_Path(2);


            #region planar offset points

            // calculate the planar offset points (points offsetet from the centre points)

            // Side A left
            Point3d pOAL0 = offsetvalue[TopoPts[0]] * PointsL01[0] + CentersL01[0] * (1.0 - offsetvalue[TopoPts[0]]);
            Point3d pOAL1 = offsetvalue[TopoPts[0]] * PointsL02[0] + CentersL02[0] * (1.0 - offsetvalue[TopoPts[0]]);
            LeftfromCenterL01.Add(pOAL0);
            LeftfromCenterL02.Add(pOAL1);

            // Side A right
            Point3d pOAR0 = offsetvalue[TopoPts[1]] * PointsL01[1] + CentersL01[0] * (1.0 - offsetvalue[TopoPts[1]]);
            Point3d pOAR1 = offsetvalue[TopoPts[1]] * PointsL02[1] + CentersL02[0] * (1.0 - offsetvalue[TopoPts[1]]);
            RightfromCenterL01.Add(pOAR0);
            RightfromCenterL02.Add(pOAR1);

            // Side B left
            Point3d pOBL0 = offsetvalue[TopoPts[1]] * PointsL01[1] + CentersL01[1] * (1.0 - offsetvalue[TopoPts[1]]);
            Point3d pOBL1 = offsetvalue[TopoPts[1]] * PointsL02[1] + CentersL02[1] * (1.0 - offsetvalue[TopoPts[1]]);
            LeftfromCenterL01.Add(pOBL0);
            LeftfromCenterL02.Add(pOBL1);

            // Side B right
            Point3d pOBR0 = offsetvalue[TopoPts[2]] * PointsL01[2] + CentersL01[1] * (1.0 - offsetvalue[TopoPts[2]]);
            Point3d pOBR1 = offsetvalue[TopoPts[2]] * PointsL02[2] + CentersL02[1] * (1.0 - offsetvalue[TopoPts[2]]);
            RightfromCenterL01.Add(pOBR0);
            RightfromCenterL02.Add(pOBR1);

            // Side C left
            Point3d pOCL0 = offsetvalue[TopoPts[2]] * PointsL01[2] + CentersL01[2] * (1.0 - offsetvalue[TopoPts[2]]);
            Point3d pOCL1 = offsetvalue[TopoPts[2]] * PointsL02[2] + CentersL02[2] * (1.0 - offsetvalue[TopoPts[2]]);
            LeftfromCenterL01.Add(pOCL0);
            LeftfromCenterL02.Add(pOCL1);

            // Side C right
            Point3d pOCR0 = offsetvalue[TopoPts[0]] * PointsL01[0] + CentersL01[2] * (1.0 - offsetvalue[TopoPts[0]]);
            Point3d pOCR1 = offsetvalue[TopoPts[0]] * PointsL02[0] + CentersL02[2] * (1.0 - offsetvalue[TopoPts[0]]);
            RightfromCenterL01.Add(pOCR0);
            RightfromCenterL02.Add(pOCR1);

            #endregion planar offset points

            #region curve control points

            // calculate control points for curves

            // Side A left
            Point3d ctrlAL0 = controlpointvalue[TopoPts[0]] * PointsL01[0] + pOAL0 * (1.0 - controlpointvalue[TopoPts[0]]);
            Point3d ctrlAL1 = controlpointvalue[TopoPts[0]] * PointsL02[0] + pOAL1 * (1.0 - controlpointvalue[TopoPts[0]]);
            // Side A right
            Point3d ctrlAR0 = controlpointvalue[TopoPts[1]] * PointsL01[1] + pOAR0 * (1.0 - controlpointvalue[TopoPts[1]]);
            Point3d ctrlAR1 = controlpointvalue[TopoPts[1]] * PointsL02[1] + pOAR1 * (1.0 - controlpointvalue[TopoPts[1]]);

            // Side B left
            Point3d ctrlBL0 = controlpointvalue[TopoPts[1]] * PointsL01[1] + pOBL0 * (1.0 - controlpointvalue[TopoPts[1]]);
            Point3d ctrlBL1 = controlpointvalue[TopoPts[1]] * PointsL02[1] + pOBL1 * (1.0 - controlpointvalue[TopoPts[1]]);
            // Side B right
            Point3d ctrlBR0 = controlpointvalue[TopoPts[2]] * PointsL01[2] + pOBR0 * (1.0 - controlpointvalue[TopoPts[2]]);
            Point3d ctrlBR1 = controlpointvalue[TopoPts[2]] * PointsL02[2] + pOBR1 * (1.0 - controlpointvalue[TopoPts[2]]);

            // Side C left
            Point3d ctrlCL0 = controlpointvalue[TopoPts[2]] * PointsL01[2] + pOCL0 * (1.0 - controlpointvalue[TopoPts[2]]);
            Point3d ctrlCL1 = controlpointvalue[TopoPts[2]] * PointsL02[2] + pOCL1 * (1.0 - controlpointvalue[TopoPts[2]]);
            // Side C right
            Point3d ctrlCR0 = controlpointvalue[TopoPts[0]] * PointsL01[0] + pOCR0 * (1.0 - controlpointvalue[TopoPts[0]]);
            Point3d ctrlCR1 = controlpointvalue[TopoPts[0]] * PointsL02[0] + pOCR1 * (1.0 - controlpointvalue[TopoPts[0]]);

            #endregion curve control points

            #region add points to data tree

            // A Level 01
            foreach (Point3d p in new List<Point3d>() { pOAL0, pOBR0, CentersL01[1], AreaCentreL01, CentersL01[0] })
            { PlanarPentagonL01.Add(p, pthA); }
            // A Level 02
            foreach (Point3d p in new List<Point3d>() { pOAL1, pOBR1, CentersL02[1], AreaCentreL02, CentersL02[0] })
            { PlanarPentagonL02.Add(p, pthA); }

            // B Level 01
            foreach (Point3d p in new List<Point3d>() { pOBL0, pOCR0, CentersL01[2], AreaCentreL01, CentersL01[1] })
            { PlanarPentagonL01.Add(p, pthB); }
            // B Level 02
            foreach (Point3d p in new List<Point3d>() { pOBL1, pOCR1, CentersL02[2], AreaCentreL02, CentersL02[1] })
            { PlanarPentagonL02.Add(p, pthB); }

            // C Level 01
            foreach (Point3d p in new List<Point3d>() { pOCL0, pOAR0, CentersL01[0], AreaCentreL01, CentersL01[2] })
            { PlanarPentagonL01.Add(p, pthC); }
            // C Level 02
            foreach (Point3d p in new List<Point3d>() { pOCL1, pOAR1, CentersL02[0], AreaCentreL02, CentersL02[2] })
            { PlanarPentagonL02.Add(p, pthC); }

            #endregion add points to data tree

            #region curves
            //---- Create Curves ---------------------------------------------------------

            // Side A left
            List<Point3d> cALPts = new List<Point3d>() { pOAL0, ctrlAL0, ctrlAL1, pOAL1 };
            Curve cAL = Curve.CreateControlPointCurve(cALPts);
            // Side A right
            List<Point3d> cARPts = new List<Point3d>() { pOAR0, ctrlAR0, ctrlAR1, pOAR1 };
            Curve cAR = Curve.CreateControlPointCurve(cARPts);

            // Side B left
            List<Point3d> cBLPts = new List<Point3d>() { pOBL0, ctrlBL0, ctrlBL1, pOBL1 };
            Curve cBL = Curve.CreateControlPointCurve(cBLPts);
            // Side B right
            List<Point3d> cBRPts = new List<Point3d>() { pOBR0, ctrlBR0, ctrlBR1, pOBR1 };
            Curve cBR = Curve.CreateControlPointCurve(cBRPts);

            // Side C left
            List<Point3d> cCLPts = new List<Point3d>() { pOCL0, ctrlCL0, ctrlCL1, pOCL1 };
            Curve cCL = Curve.CreateControlPointCurve(cCLPts);
            // Side C right
            List<Point3d> cCRPts = new List<Point3d>() { pOCR0, ctrlCR0, ctrlCR1, pOCR1 };
            Curve cCR = Curve.CreateControlPointCurve(cCRPts);


            // Add to data trees
            #region Add to Trees

            // A
            Curves.Add(cAL, pthA);
            Curves.Add(cBR, pthA);

            CurvesPts.Add(cALPts, pthA);
            CurvesPts.Add(cBRPts, pthA);

            // B
            Curves.Add(cBL, pthB);
            Curves.Add(cCR, pthB);

            CurvesPts.Add(cBLPts, pthB);
            CurvesPts.Add(cCRPts, pthB);

            // C
            Curves.Add(cCL, pthC);
            Curves.Add(cAR, pthC);

            CurvesPts.Add(cCLPts, pthC);
            CurvesPts.Add(cARPts, pthC);

            // Create and Add StripeID
            StripeID.Add(ID, pthA); StripeID.Add(NFIdx[1], pthA); StripeID.Add(NFIdx[2], pthA);

            StripeID.Add(ID, pthB); StripeID.Add(NFIdx[2], pthB); StripeID.Add(NFIdx[0], pthB);

            StripeID.Add(ID, pthC); StripeID.Add(NFIdx[0], pthC); StripeID.Add(NFIdx[1], pthC);

            #endregion Add to Trees


            #endregion curves


            //// =======================================================================================
            // Added by Gene 
            #region extended Curves
            //---- Create Extended Curves --- Join Curve with Line -------------------------------------

            // Side A left
            Curve ecAL = Curve.JoinCurves(new List<Curve>() { 
                new LineCurve(PointsL01[0], pOAL0), cAL, 
                new LineCurve(pOAL1, PointsL02[0]) }, documentTolerance, true)[0];

            // Side A right
            Curve ecAR = Curve.JoinCurves(new List<Curve>() { 
                new LineCurve(PointsL01[1], pOAR0), cAR, 
                new LineCurve(pOAR1, PointsL02[1]) }, documentTolerance, true)[0];

            // Side B left
            Curve ecBL = Curve.JoinCurves(new List<Curve>() { 
                new LineCurve(PointsL01[1], pOBL0), cBL, 
                new LineCurve(pOBL1, PointsL02[1]) }, documentTolerance, true)[0];

            // Side B right
            Curve ecBR = Curve.JoinCurves(new List<Curve>() { 
                new LineCurve(PointsL01[2], pOBR0), cBR, 
                new LineCurve(pOBR1, PointsL02[2]) }, documentTolerance, true)[0];

            // Side C left
            Curve ecCL = Curve.JoinCurves(new List<Curve>() { 
                new LineCurve(PointsL01[2], pOCL0), cCL, 
                new LineCurve(pOCL1, PointsL02[2]) }, documentTolerance, true)[0];

            // Side C right
            Curve ecCR = Curve.JoinCurves(new List<Curve>() { 
                new LineCurve(PointsL01[0], pOCR0), cCR, 
                new LineCurve(pOCR1, PointsL02[0]) }, documentTolerance, true)[0];


            // A' cutting points
            List<Point3d> cuttingA = new List<Point3d>() { 
                AreaCentreL01, AreaCentreL02, 
                (pOAR0+pOAL0)*0.5, (pOAR1+pOAL1)*0.5,
                (pOBR0+pOBL0)*0.5, (pOBR1+pOBL1)*0.5,
                PointsL01[1], PointsL02[1]
            };

            // B' cutting points
            List<Point3d> cuttingB = new List<Point3d>() { 
                AreaCentreL01, AreaCentreL02, 
                (pOBR0+pOBL0)*0.5, (pOBR1+pOBL1)*0.5,
                (pOCR0+pOCL0)*0.5, (pOCR1+pOCL1)*0.5,
                PointsL01[2], PointsL02[2]
            };

            // C' cutting points
            List<Point3d> cuttingC = new List<Point3d>() { 
                AreaCentreL01, AreaCentreL02, 
                (pOCR0+pOCL0)*0.5, (pOCR1+pOCL1)*0.5,
                (pOAR0+pOAL0)*0.5, (pOAR1+pOAL1)*0.5,
                PointsL01[0], PointsL02[0]
            };


            // Add to data trees
            #region Add to Trees

            // A'
            ExtendedCurves.Add(ecAL, pthA);
            ExtendedCurves.Add(ecBR, pthA);

            cuttingPoints.Add(cuttingA, pthA);

            // B'
            ExtendedCurves.Add(ecBL, pthB);
            ExtendedCurves.Add(ecCR, pthB);

            cuttingPoints.Add(cuttingB, pthB);

            // C'
            ExtendedCurves.Add(ecCL, pthC);
            ExtendedCurves.Add(ecAR, pthC);

            cuttingPoints.Add(cuttingC, pthC);


            #endregion Add to Trees


            #endregion  extended Curves

            //// =======================================================================================

        }

        public void Create_Plate(List<bool> openclosed)
        {
            GH_Path pthA = new GH_Path(0);
            GH_Path pthB = new GH_Path(1);
            GH_Path pthC = new GH_Path(2);

            //---- Stripe A -----------------------------------------------------------------------------------------------------
            if (openclosed[TopoPts[0]])
            {
                List<Point3d> PlateL01Pts = new List<Point3d>() { CentersL01[0], PointsL01[0], CentersL01[2], CentersL01[0] };
                Curve PlateL01 = Curve.CreateControlPointCurve(PlateL01Pts, 1);

                List<Point3d> PlateL02Pts = new List<Point3d>() { CentersL02[0], PointsL02[0], CentersL02[2], CentersL02[0] };
                Curve PlateL02 = Curve.CreateControlPointCurve(PlateL02Pts, 1);

                // Curves
                PlateCrv.Add(PlateL01, pthA);
                PlateCrv.Add(PlateL02, pthA);
                // Control Points
                PlateCrvPts.Add(PlateL01Pts, pthA);
                PlateCrvPts.Add(PlateL02Pts, pthA);
            }

            else
            {
                // Curves
                PlateCrv.Add(null, pthA);
                PlateCrv.Add(null, pthA);
                // Control Points
                PlateCrvPts.Add(null, pthA);
                PlateCrvPts.Add(null, pthA);
            }

            //---- Stripe B -----------------------------------------------------------------------------------------------------            
            if (openclosed[TopoPts[1]])
            {

                List<Point3d> PlateL01Pts = new List<Point3d>() { CentersL01[1], PointsL01[1], CentersL01[0], CentersL01[1] };
                Curve PlateL01 = Curve.CreateControlPointCurve(PlateL01Pts, 1);

                List<Point3d> PlateL02Pts = new List<Point3d>() { CentersL02[1], PointsL02[1], CentersL02[0], CentersL02[1] };
                Curve PlateL02 = Curve.CreateControlPointCurve(PlateL02Pts, 1);

                // Curves
                PlateCrv.Add(PlateL01, pthB);
                PlateCrv.Add(PlateL02, pthB);
                // Control Points
                PlateCrvPts.Add(PlateL01Pts, pthB);
                PlateCrvPts.Add(PlateL02Pts, pthB);
            }

            else
            {
                // Curves
                PlateCrv.Add(null, pthB);
                PlateCrv.Add(null, pthB);
                // Control Points
                PlateCrvPts.Add(null, pthB);
                PlateCrvPts.Add(null, pthB);
            }

            //---- Stripe C -----------------------------------------------------------------------------------------------------            
            if (openclosed[TopoPts[2]])
            {
                List<Point3d> PlateL01Pts = new List<Point3d>() { CentersL01[2], PointsL01[2], CentersL01[1], CentersL01[2] };
                Curve PlateL01 = Curve.CreateControlPointCurve(PlateL01Pts, 1);

                List<Point3d> PlateL02Pts = new List<Point3d>() { CentersL02[2], PointsL02[2], CentersL02[1], CentersL02[2] };
                Curve PlateL02 = Curve.CreateControlPointCurve(PlateL02Pts, 1);

                // Curves
                PlateCrv.Add(PlateL01, pthC);
                PlateCrv.Add(PlateL02, pthC);
                // Control Points
                PlateCrvPts.Add(PlateL01Pts, pthC);
                PlateCrvPts.Add(PlateL02Pts, pthC);
            }

            else
            {
                // Curves
                PlateCrv.Add(null, pthC);
                PlateCrv.Add(null, pthC);
                // Control Points
                PlateCrvPts.Add(null, pthC);
                PlateCrvPts.Add(null, pthC);
            }

        }

        //// =======================================================================================
        // Added by Gene 
        public void Create_Stipe_SinSrf()
        {
            GH_Path pthA = new GH_Path(0);
            GH_Path pthB = new GH_Path(1);
            GH_Path pthC = new GH_Path(2);

            Tri_Loop_Stripe StripeAS = new Tri_Loop_Stripe(
                ExtendedCurves.Branch(pthA)[0], ExtendedCurves.Branch(pthA)[1], cuttingPoints.Branch(pthA)[0], documentTolerance);
            Tri_Loop_Stripe StripeBS = new Tri_Loop_Stripe(
                ExtendedCurves.Branch(pthB)[0], ExtendedCurves.Branch(pthB)[1], cuttingPoints.Branch(pthB)[0], documentTolerance);
            Tri_Loop_Stripe StripeCS = new Tri_Loop_Stripe(
                ExtendedCurves.Branch(pthC)[0], ExtendedCurves.Branch(pthC)[1], cuttingPoints.Branch(pthC)[0], documentTolerance);

            StripesSingle.Add(StripeAS);
            StripesSingle.Add(StripeBS);
            StripesSingle.Add(StripeCS);

        }

        //// =======================================================================================


        public void Create_Stripe()
        {
            GH_Path pthA = new GH_Path(0);
            GH_Path pthB = new GH_Path(1);
            GH_Path pthC = new GH_Path(2);

            Tri_Loop_Stripe StripeA = new Tri_Loop_Stripe
                (
                PlanarPentagonL01.Branch(pthA), PlanarPentagonL02.Branch(pthA),
                PlateCrvPts.Branch(pthA)[0], PlateCrvPts.Branch(pthA)[1],
                CurvesPts.Branch(pthA)[0], CurvesPts.Branch(pthA)[1],
                StripeID.Branch(pthA)
                );

            Tri_Loop_Stripe StripeB = new Tri_Loop_Stripe
              (
              PlanarPentagonL01.Branch(pthB), PlanarPentagonL02.Branch(pthB),
              PlateCrvPts.Branch(pthB)[0], PlateCrvPts.Branch(pthB)[1],
              CurvesPts.Branch(pthB)[0], CurvesPts.Branch(pthB)[1],
              StripeID.Branch(pthB)
              );

            Tri_Loop_Stripe StripeC = new Tri_Loop_Stripe
              (
              PlanarPentagonL01.Branch(pthC), PlanarPentagonL02.Branch(pthC),
              PlateCrvPts.Branch(pthC)[0], PlateCrvPts.Branch(pthC)[1],
              CurvesPts.Branch(pthC)[0], CurvesPts.Branch(pthC)[1],
              StripeID.Branch(pthC)
              );

            Stripes.Add(StripeA); // [0]
            Stripes.Add(StripeB); // [1]
            Stripes.Add(StripeC); // [2]
        }

        // Sofistik Information
        public void SofistikInformation()
        {
            GH_Path pthPlanarL01 = new GH_Path(0);
            GH_Path pthPlanarL02 = new GH_Path(1);
            GH_Path pthLoop00 = new GH_Path(2);
            GH_Path pthLoop01 = new GH_Path(3);
            GH_Path pthLoop02 = new GH_Path(4);

            #region Sofistik Points

            for (int i = 0; i < 3; i++)
            {
                //planar Part Points L01
                AllSofistikPoints.Add(RightfromCenterL01[i], pthPlanarL01);
                SofistikIndexies.Add(i + 0 + (i * 2), pthPlanarL01);

                AllSofistikPoints.Add(CentersL01[i], pthPlanarL01);
                SofistikIndexies.Add(i + 1 + (i * 2), pthPlanarL01);

                AllSofistikPoints.Add(LeftfromCenterL01[i], pthPlanarL01);
                SofistikIndexies.Add(i + 2 + (i * 2), pthPlanarL01);


                //planar Part Points L02
                AllSofistikPoints.Add(RightfromCenterL02[i], pthPlanarL02);
                SofistikIndexies.Add(i + 10 + (i * 2), pthPlanarL02);

                AllSofistikPoints.Add(CentersL02[i], pthPlanarL02);
                SofistikIndexies.Add(i + 11 + (i * 2), pthPlanarL02);

                AllSofistikPoints.Add(LeftfromCenterL02[i], pthPlanarL02);
                SofistikIndexies.Add(i + 12 + (i * 2), pthPlanarL02);

                #region Sofistik Plate

                //---- Sofistik Plate Curves ---------------------------------------------------------------------------
                //planar Part Points L01
                SofistikPlatePoints.Add(RightfromCenterL01[i], pthPlanarL01);
                SofistikPlateIndexies.Add(i + 0 + (i * 2), pthPlanarL01);

                SofistikPlatePoints.Add(CentersL01[i], pthPlanarL01);
                SofistikPlateIndexies.Add(i + 1 + (i * 2), pthPlanarL01);

                SofistikPlatePoints.Add(LeftfromCenterL01[i], pthPlanarL01);
                SofistikPlateIndexies.Add(i + 2 + (i * 2), pthPlanarL01);


                //planar Part Points L02
                SofistikPlatePoints.Add(RightfromCenterL02[i], pthPlanarL02);
                SofistikPlateIndexies.Add(i + 10 + (i * 2), pthPlanarL02);

                SofistikPlatePoints.Add(CentersL02[i], pthPlanarL02);
                SofistikPlateIndexies.Add(i + 11 + (i * 2), pthPlanarL02);

                SofistikPlatePoints.Add(LeftfromCenterL02[i], pthPlanarL02);
                SofistikPlateIndexies.Add(i + 12 + (i * 2), pthPlanarL02);
                #endregion Sofistik Plate

            }


            for (int i = 0; i < 3; i++)
            {
                int counter = 0;
                foreach (Curve c in Curves.Branch(i))
                {
                    AllSofistikPoints.Add(c.PointAtNormalizedLength(0.0), new GH_Path(i + 2));

                    // check for corresponding point index
                    for (int k = 0; k < AllSofistikPoints.Branch(0).Count; k++)
                    {
                        if ((Math.Abs(c.PointAtNormalizedLength(0.0).X - AllSofistikPoints.Branch(0)[k].X) < 0.001) &&
                           (Math.Abs(c.PointAtNormalizedLength(0.0).Y - AllSofistikPoints.Branch(0)[k].Y) < 0.001))
                        { SofistikIndexies.Add(SofistikIndexies.Branch(0)[k], new GH_Path(i + 2)); }
                    }

                    //AllSofistikPoints.Add(c.PointAtNormalizedLength(0.5), new GH_Path(i + 2));

                    // check for corresponding point index
                    //SofistikIndexies.Add(i * 2 + counter + 20, new GH_Path(i + 2));

                    AllSofistikPoints.Add(c.PointAtNormalizedLength(1.0), new GH_Path(i + 2));

                    // check for corresponding point index
                    for (int k = 0; k < AllSofistikPoints.Branch(1).Count; k++)
                    {
                        if ((Math.Abs(c.PointAtNormalizedLength(1.0).X - AllSofistikPoints.Branch(1)[k].X) < 0.001) &&
                            (Math.Abs(c.PointAtNormalizedLength(1.0).Y - AllSofistikPoints.Branch(1)[k].Y) < 0.001))
                        { SofistikIndexies.Add(SofistikIndexies.Branch(1)[k], new GH_Path(i + 2)); }
                    }

                    counter++;
                }
            }
            #endregion Sofistik Points

            #region Sofistik Curves

            // planar Curves

            //L01
            #region L01
            for (int i = 0; i < AllSofistikPoints.Branch(0).Count; i++)
            {
                if (i < AllSofistikPoints.Branch(0).Count - 1)
                {
                    int p = i + 0;



                    SofistikCrvPtIdx.Add(SofistikIndexies.Branch(0)[i], new GH_Path(p));
                    SofistikCrvPtIdx.Add(SofistikIndexies.Branch(0)[i + 1], new GH_Path(p));

                    SofistikCrvPtCoo.Add(AllSofistikPoints.Branch(0)[i], new GH_Path(p));
                    SofistikCrvPtCoo.Add(AllSofistikPoints.Branch(0)[i + 1], new GH_Path(p));

                    SofistikCurves.Add(
                        Curve.CreateControlPointCurve(new List<Point3d>() { AllSofistikPoints.Branch(0)[i], AllSofistikPoints.Branch(0)[i + 1] }),
                        new GH_Path(p));

                }
                else
                {
                    int p = i + 0;

                    SofistikCrvPtIdx.Add(SofistikIndexies.Branch(0)[i], new GH_Path(p));
                    SofistikCrvPtIdx.Add(SofistikIndexies.Branch(0)[0], new GH_Path(p));

                    SofistikCrvPtCoo.Add(AllSofistikPoints.Branch(0)[i], new GH_Path(p));
                    SofistikCrvPtCoo.Add(AllSofistikPoints.Branch(0)[0], new GH_Path(p));

                    SofistikCurves.Add(
                        Curve.CreateControlPointCurve(new List<Point3d>() { AllSofistikPoints.Branch(0)[i], AllSofistikPoints.Branch(0)[0] }),
                        new GH_Path(p));

                }
            }
            #endregion L01

            //L02
            #region L02
            for (int i = 0; i < AllSofistikPoints.Branch(1).Count; i++)
            {
                if (i < AllSofistikPoints.Branch(1).Count - 1)
                {
                    int p = i + 10;

                    SofistikCrvPtIdx.Add(SofistikIndexies.Branch(1)[i], new GH_Path(p));
                    SofistikCrvPtIdx.Add(SofistikIndexies.Branch(1)[i + 1], new GH_Path(p));

                    SofistikCrvPtCoo.Add(AllSofistikPoints.Branch(1)[i], new GH_Path(p));
                    SofistikCrvPtCoo.Add(AllSofistikPoints.Branch(1)[i + 1], new GH_Path(p));

                    SofistikCurves.Add(
                        Curve.CreateControlPointCurve(new List<Point3d>() { AllSofistikPoints.Branch(1)[i], AllSofistikPoints.Branch(1)[i + 1] }),
                        new GH_Path(p));
                }
                else
                {
                    int p = i + 10;

                    SofistikCrvPtIdx.Add(SofistikIndexies.Branch(1)[i], new GH_Path(p));
                    SofistikCrvPtIdx.Add(SofistikIndexies.Branch(1)[0], new GH_Path(p));

                    SofistikCrvPtCoo.Add(AllSofistikPoints.Branch(1)[i], new GH_Path(p));
                    SofistikCrvPtCoo.Add(AllSofistikPoints.Branch(1)[0], new GH_Path(p));

                    SofistikCurves.Add(
                        Curve.CreateControlPointCurve(new List<Point3d>() { AllSofistikPoints.Branch(1)[i], AllSofistikPoints.Branch(1)[0] }),
                        new GH_Path(p));
                }
            }
            #endregion L02

            //Curves A
            #region Crv A

            SofistikCurves.Add(Curves.Branch(0)[0], new GH_Path(20));
            SofistikCurves.Add(Curves.Branch(0)[1], new GH_Path(21));

            for (int i = 0; i < AllSofistikPoints.Branch(2).Count; i++)
            {
                if (i < AllSofistikPoints.Branch(2).Count * 0.5)
                {
                    // Left Curve
                    int p = 20;

                    SofistikCrvPtIdx.Add(SofistikIndexies.Branch(2)[i], new GH_Path(p));
                    SofistikCrvPtCoo.Add(AllSofistikPoints.Branch(2)[i], new GH_Path(p));

                }
                else
                {
                    // Right Curve
                    int p = 21;

                    SofistikCrvPtIdx.Add(SofistikIndexies.Branch(2)[i], new GH_Path(p));
                    SofistikCrvPtCoo.Add(AllSofistikPoints.Branch(2)[i], new GH_Path(p));

                }
            }
            #endregion Crv A

            //Curves B
            #region Crv B

            SofistikCurves.Add(Curves.Branch(1)[0], new GH_Path(30));
            SofistikCurves.Add(Curves.Branch(1)[1], new GH_Path(31));

            for (int i = 0; i < AllSofistikPoints.Branch(3).Count; i++)
            {
                if (i < AllSofistikPoints.Branch(3).Count * 0.5)
                {
                    // Left Curve
                    int p = 30;

                    SofistikCrvPtIdx.Add(SofistikIndexies.Branch(3)[i], new GH_Path(p));
                    SofistikCrvPtCoo.Add(AllSofistikPoints.Branch(3)[i], new GH_Path(p));

                }
                else
                {
                    // Right Curve
                    int p = 31;

                    SofistikCrvPtIdx.Add(SofistikIndexies.Branch(3)[i], new GH_Path(p));
                    SofistikCrvPtCoo.Add(AllSofistikPoints.Branch(3)[i], new GH_Path(p));

                }
            }
            #endregion Crv B


            //Curves C
            #region Crv C

            SofistikCurves.Add(Curves.Branch(2)[0], new GH_Path(40));
            SofistikCurves.Add(Curves.Branch(2)[1], new GH_Path(41));

            for (int i = 0; i < AllSofistikPoints.Branch(4).Count; i++)
            {
                if (i < AllSofistikPoints.Branch(4).Count * 0.5)
                {
                    // Left Curve
                    int p = 40;

                    SofistikCrvPtIdx.Add(SofistikIndexies.Branch(4)[i], new GH_Path(p));
                    SofistikCrvPtCoo.Add(AllSofistikPoints.Branch(4)[i], new GH_Path(p));

                }
                else
                {
                    // Right Curve
                    int p = 41;

                    SofistikCrvPtIdx.Add(SofistikIndexies.Branch(4)[i], new GH_Path(p));
                    SofistikCrvPtCoo.Add(AllSofistikPoints.Branch(4)[i], new GH_Path(p));

                }
            }
            #endregion Crv C

            //Curves AL-BR
            #region  AL-BR
            //L01
            SofistikCrvPtCoo.Add(SofistikCrvPtCoo.Branch(new GH_Path(20))[0], new GH_Path(22));
            SofistikCrvPtCoo.Add(SofistikCrvPtCoo.Branch(new GH_Path(21))[0], new GH_Path(22));

            SofistikCrvPtIdx.Add(SofistikCrvPtIdx.Branch(new GH_Path(20))[0], new GH_Path(22));
            SofistikCrvPtIdx.Add(SofistikCrvPtIdx.Branch(new GH_Path(21))[0], new GH_Path(22));

            SofistikCurves.Add(
                    Curve.CreateControlPointCurve(new List<Point3d>() { SofistikCrvPtCoo.Branch(new GH_Path(20))[0], SofistikCrvPtCoo.Branch(new GH_Path(21))[0] }),
                    new GH_Path(22));
            //L02
            SofistikCrvPtCoo.Add(SofistikCrvPtCoo.Branch(new GH_Path(20))[1], new GH_Path(23));
            SofistikCrvPtCoo.Add(SofistikCrvPtCoo.Branch(new GH_Path(21))[1], new GH_Path(23));

            SofistikCrvPtIdx.Add(SofistikCrvPtIdx.Branch(new GH_Path(20))[1], new GH_Path(23));
            SofistikCrvPtIdx.Add(SofistikCrvPtIdx.Branch(new GH_Path(21))[1], new GH_Path(23));

            SofistikCurves.Add(
                Curve.CreateControlPointCurve(new List<Point3d>() { SofistikCrvPtCoo.Branch(new GH_Path(20))[1], SofistikCrvPtCoo.Branch(new GH_Path(21))[1] }),
                new GH_Path(23));
            #endregion  AL-BR

            //Curves BL-CR
            #region  BL-CR
            //L01
            SofistikCrvPtCoo.Add(SofistikCrvPtCoo.Branch(new GH_Path(30))[0], new GH_Path(32));
            SofistikCrvPtCoo.Add(SofistikCrvPtCoo.Branch(new GH_Path(31))[0], new GH_Path(32));

            SofistikCrvPtIdx.Add(SofistikCrvPtIdx.Branch(new GH_Path(30))[0], new GH_Path(32));
            SofistikCrvPtIdx.Add(SofistikCrvPtIdx.Branch(new GH_Path(31))[0], new GH_Path(32));

            SofistikCurves.Add(
                Curve.CreateControlPointCurve(new List<Point3d>() { SofistikCrvPtCoo.Branch(new GH_Path(30))[0], SofistikCrvPtCoo.Branch(new GH_Path(31))[0] }),
                new GH_Path(32));
            //L02
            SofistikCrvPtCoo.Add(SofistikCrvPtCoo.Branch(new GH_Path(30))[1], new GH_Path(33));
            SofistikCrvPtCoo.Add(SofistikCrvPtCoo.Branch(new GH_Path(31))[1], new GH_Path(33));

            SofistikCrvPtIdx.Add(SofistikCrvPtIdx.Branch(new GH_Path(30))[1], new GH_Path(33));
            SofistikCrvPtIdx.Add(SofistikCrvPtIdx.Branch(new GH_Path(31))[1], new GH_Path(33));

            SofistikCurves.Add(
                Curve.CreateControlPointCurve(new List<Point3d>() { SofistikCrvPtCoo.Branch(new GH_Path(30))[1], SofistikCrvPtCoo.Branch(new GH_Path(31))[1] }),
                new GH_Path(33));
            #endregion  BL-CR

            //Curves CL-AR
            #region  CL-AR
            //L01
            SofistikCrvPtCoo.Add(SofistikCrvPtCoo.Branch(new GH_Path(40))[0], new GH_Path(42));
            SofistikCrvPtCoo.Add(SofistikCrvPtCoo.Branch(new GH_Path(41))[0], new GH_Path(42));

            SofistikCrvPtIdx.Add(SofistikCrvPtIdx.Branch(new GH_Path(40))[0], new GH_Path(42));
            SofistikCrvPtIdx.Add(SofistikCrvPtIdx.Branch(new GH_Path(41))[0], new GH_Path(42));

            SofistikCurves.Add(
                Curve.CreateControlPointCurve(new List<Point3d>() { SofistikCrvPtCoo.Branch(new GH_Path(40))[0], SofistikCrvPtCoo.Branch(new GH_Path(41))[0] }),
                new GH_Path(42));

            //L02
            SofistikCrvPtCoo.Add(SofistikCrvPtCoo.Branch(new GH_Path(40))[1], new GH_Path(43));
            SofistikCrvPtCoo.Add(SofistikCrvPtCoo.Branch(new GH_Path(41))[1], new GH_Path(43));

            SofistikCrvPtIdx.Add(SofistikCrvPtIdx.Branch(new GH_Path(40))[1], new GH_Path(43));
            SofistikCrvPtIdx.Add(SofistikCrvPtIdx.Branch(new GH_Path(41))[1], new GH_Path(43));

            SofistikCurves.Add(
                Curve.CreateControlPointCurve(new List<Point3d>() { SofistikCrvPtCoo.Branch(new GH_Path(40))[1], SofistikCrvPtCoo.Branch(new GH_Path(41))[1] }),
                new GH_Path(43));
            #endregion  CL-AR

            #endregion Sofistik Curves

        }

        public void SofistikInformation_Version2()
        {
            GH_Path pthPlanarL01 = new GH_Path(0);
            GH_Path pthPlanarL02 = new GH_Path(1);
            GH_Path pthLoop00 = new GH_Path(2);
            GH_Path pthLoop01 = new GH_Path(3);
            GH_Path pthLoop02 = new GH_Path(4);

            #region Sofistik Points

            for (int i = 0; i < 3; i++)
            {
                //planar Part Points L01
                AllSofistikPoints.Add(RightfromCenterL01[i], pthPlanarL01);
                SofistikIndexies.Add(i + 0 + (i * 2), pthPlanarL01);

                AllSofistikPoints.Add(CentersL01[i], pthPlanarL01);
                SofistikIndexies.Add(i + 1 + (i * 2), pthPlanarL01);

                AllSofistikPoints.Add(LeftfromCenterL01[i], pthPlanarL01);
                SofistikIndexies.Add(i + 2 + (i * 2), pthPlanarL01);


                //planar Part Points L02
                AllSofistikPoints.Add(RightfromCenterL02[i], pthPlanarL02);
                SofistikIndexies.Add(i + 10 + (i * 2), pthPlanarL02);

                AllSofistikPoints.Add(CentersL02[i], pthPlanarL02);
                SofistikIndexies.Add(i + 11 + (i * 2), pthPlanarL02);

                AllSofistikPoints.Add(LeftfromCenterL02[i], pthPlanarL02);
                SofistikIndexies.Add(i + 12 + (i * 2), pthPlanarL02);

                #region Sofistik Plate

                //---- Sofistik Plate Curves ---------------------------------------------------------------------------
                //planar Part Points L01
                SofistikPlatePoints.Add(RightfromCenterL01[i], pthPlanarL01);
                SofistikPlateIndexies.Add(i + 0 + (i * 2), pthPlanarL01);

                SofistikPlatePoints.Add(CentersL01[i], pthPlanarL01);
                SofistikPlateIndexies.Add(i + 1 + (i * 2), pthPlanarL01);

                SofistikPlatePoints.Add(LeftfromCenterL01[i], pthPlanarL01);
                SofistikPlateIndexies.Add(i + 2 + (i * 2), pthPlanarL01);


                //planar Part Points L02
                SofistikPlatePoints.Add(RightfromCenterL02[i], pthPlanarL02);
                SofistikPlateIndexies.Add(i + 10 + (i * 2), pthPlanarL02);

                SofistikPlatePoints.Add(CentersL02[i], pthPlanarL02);
                SofistikPlateIndexies.Add(i + 11 + (i * 2), pthPlanarL02);

                SofistikPlatePoints.Add(LeftfromCenterL02[i], pthPlanarL02);
                SofistikPlateIndexies.Add(i + 12 + (i * 2), pthPlanarL02);
                #endregion Sofistik Plate

            }


            for (int i = 0; i < 3; i++)
            {
                int counter = 0;
                foreach (Curve c in Curves.Branch(i))
                {
                    AllSofistikPoints.Add(c.PointAtNormalizedLength(0.0), new GH_Path(i + 2));

                    // check for corresponding point index
                    for (int k = 0; k < AllSofistikPoints.Branch(0).Count; k++)
                    {
                        if ((Math.Abs(c.PointAtNormalizedLength(0.0).X - AllSofistikPoints.Branch(0)[k].X) < 0.001) &&
                           (Math.Abs(c.PointAtNormalizedLength(0.0).Y - AllSofistikPoints.Branch(0)[k].Y) < 0.001))
                        { SofistikIndexies.Add(SofistikIndexies.Branch(0)[k], new GH_Path(i + 2)); }
                    }

                    //AllSofistikPoints.Add(c.PointAtNormalizedLength(0.5), new GH_Path(i + 2));

                    // check for corresponding point index
                    //SofistikIndexies.Add(i * 2 + counter + 20, new GH_Path(i + 2));

                    AllSofistikPoints.Add(c.PointAtNormalizedLength(1.0), new GH_Path(i + 2));

                    // check for corresponding point index
                    for (int k = 0; k < AllSofistikPoints.Branch(1).Count; k++)
                    {
                        if ((Math.Abs(c.PointAtNormalizedLength(1.0).X - AllSofistikPoints.Branch(1)[k].X) < 0.001) &&
                            (Math.Abs(c.PointAtNormalizedLength(1.0).Y - AllSofistikPoints.Branch(1)[k].Y) < 0.001))
                        { SofistikIndexies.Add(SofistikIndexies.Branch(1)[k], new GH_Path(i + 2)); }
                    }

                    counter++;
                }
            }
            #endregion Sofistik Points

            #region Sofistik Curves

            // planar Curves

            //L01
            #region L01
            for (int i = 0; i < AllSofistikPoints.Branch(0).Count; i++)
            {
                if (i < AllSofistikPoints.Branch(0).Count - 1)
                {
                    GH_Path p = new GH_Path(0);
                    p.AppendElement(i);

                    SofistikCrvPtIdx.Add(SofistikIndexies.Branch(0)[i], p);
                    SofistikCrvPtIdx.Add(SofistikIndexies.Branch(0)[i + 1], p);

                    SofistikCrvPtCoo.Add(AllSofistikPoints.Branch(0)[i], p);
                    SofistikCrvPtCoo.Add(AllSofistikPoints.Branch(0)[i + 1], p);
                }
                else
                {
                    GH_Path p = new GH_Path(0);
                    p.AppendElement(i);

                    SofistikCrvPtIdx.Add(SofistikIndexies.Branch(0)[i], p);
                    SofistikCrvPtIdx.Add(SofistikIndexies.Branch(0)[0], p);

                    SofistikCrvPtCoo.Add(AllSofistikPoints.Branch(0)[i], p);
                    SofistikCrvPtCoo.Add(AllSofistikPoints.Branch(0)[0], p);


                }
            }
            #endregion L01

            //L02
            #region L02
            for (int i = 0; i < AllSofistikPoints.Branch(1).Count; i++)
            {
                if (i < AllSofistikPoints.Branch(1).Count - 1)
                {
                    GH_Path p = new GH_Path(1);
                    p.AppendElement(i);

                    SofistikCrvPtIdx.Add(SofistikIndexies.Branch(1)[i], p);
                    SofistikCrvPtIdx.Add(SofistikIndexies.Branch(1)[i + 1], p);

                    SofistikCrvPtCoo.Add(AllSofistikPoints.Branch(1)[i], p);
                    SofistikCrvPtCoo.Add(AllSofistikPoints.Branch(1)[i + 1], p);

                }
                else
                {
                    GH_Path p = new GH_Path(1);
                    p.AppendElement(i);

                    SofistikCrvPtIdx.Add(SofistikIndexies.Branch(1)[i], p);
                    SofistikCrvPtIdx.Add(SofistikIndexies.Branch(1)[0], p);

                    SofistikCrvPtCoo.Add(AllSofistikPoints.Branch(1)[i], p);
                    SofistikCrvPtCoo.Add(AllSofistikPoints.Branch(1)[0], p);

                }
            }
            #endregion L02

            //Curves A
            #region Crv A
            for (int i = 0; i < AllSofistikPoints.Branch(2).Count; i++)
            {
                if (i < AllSofistikPoints.Branch(2).Count * 0.5)
                {
                    // Left Curve
                    GH_Path p = new GH_Path(2);
                    p.AppendElement(i);

                    SofistikCrvPtIdx.Add(SofistikIndexies.Branch(2)[i], p);
                    SofistikCrvPtCoo.Add(AllSofistikPoints.Branch(2)[i], p);

                }
                else
                {
                    // Right Curve
                    GH_Path p = new GH_Path(2);
                    p.AppendElement(i);

                    SofistikCrvPtIdx.Add(SofistikIndexies.Branch(2)[i], p);
                    SofistikCrvPtCoo.Add(AllSofistikPoints.Branch(2)[i], p);
                }
            }
            #endregion Crv A

            //Curves B
            #region Crv B
            for (int i = 0; i < AllSofistikPoints.Branch(3).Count; i++)
            {
                if (i < AllSofistikPoints.Branch(3).Count * 0.5)
                {
                    // Left Curve
                    GH_Path p = new GH_Path(3);
                    p.AppendElement(i);

                    SofistikCrvPtIdx.Add(SofistikIndexies.Branch(3)[i], p);
                    SofistikCrvPtCoo.Add(AllSofistikPoints.Branch(3)[i], p);
                }
                else
                {
                    // Right Curve
                    GH_Path p = new GH_Path(3);
                    p.AppendElement(i);

                    SofistikCrvPtIdx.Add(SofistikIndexies.Branch(3)[i], p);
                    SofistikCrvPtCoo.Add(AllSofistikPoints.Branch(3)[i], p);
                }
            }
            #endregion Crv B


            //Curves C
            #region Crv C
            for (int i = 0; i < AllSofistikPoints.Branch(4).Count; i++)
            {
                if (i < AllSofistikPoints.Branch(4).Count * 0.5)
                {
                    // Left Curve
                    GH_Path p = new GH_Path(4);
                    p.AppendElement(i);

                    SofistikCrvPtIdx.Add(SofistikIndexies.Branch(4)[i], new GH_Path(p));
                    SofistikCrvPtCoo.Add(AllSofistikPoints.Branch(4)[i], new GH_Path(p));
                }
                else
                {
                    // Right Curve
                    GH_Path p = new GH_Path(4);
                    p.AppendElement(i);

                    SofistikCrvPtIdx.Add(SofistikIndexies.Branch(4)[i], new GH_Path(p));
                    SofistikCrvPtCoo.Add(AllSofistikPoints.Branch(4)[i], new GH_Path(p));
                }
            }
            #endregion Crv C

            //Curves AL-BR
            #region  AL-BR
            // adds te planar curves at the top and bottum of the bended curves in order to close the boundary curves of th surface
            GH_Path albr = new GH_Path(2);
            // i takes the start points and the endpoints of the bended curves and puts them in a branch 
            //L01            
            SofistikCrvPtCoo.Add(SofistikCrvPtCoo.Branch(albr.AppendElement(0))[0], albr.AppendElement(2)); // start point of first bended curve in path {2;0} element idx [0]
            SofistikCrvPtCoo.Add(SofistikCrvPtCoo.Branch(albr.AppendElement(1))[0], albr.AppendElement(2)); // start point of second bended curve in path {2;1} element idx [0]
            // stores it in path {2;2}
            SofistikCrvPtIdx.Add(SofistikCrvPtIdx.Branch(albr.AppendElement(0))[0], albr.AppendElement(2));
            SofistikCrvPtIdx.Add(SofistikCrvPtIdx.Branch(albr.AppendElement(1))[0], albr.AppendElement(2));
            //L02
            SofistikCrvPtCoo.Add(SofistikCrvPtCoo.Branch(albr.AppendElement(0))[1], albr.AppendElement(3));
            SofistikCrvPtCoo.Add(SofistikCrvPtCoo.Branch(albr.AppendElement(1))[1], albr.AppendElement(3));

            SofistikCrvPtIdx.Add(SofistikCrvPtIdx.Branch(albr.AppendElement(0))[1], albr.AppendElement(3));
            SofistikCrvPtIdx.Add(SofistikCrvPtIdx.Branch(albr.AppendElement(1))[1], albr.AppendElement(3));
            #endregion  AL-BR

            //Curves BL-CR
            #region  BL-CR

            GH_Path blcr = new GH_Path(3);
            //L01
            SofistikCrvPtCoo.Add(SofistikCrvPtCoo.Branch(blcr.AppendElement(0))[0], blcr.AppendElement(2));
            SofistikCrvPtCoo.Add(SofistikCrvPtCoo.Branch(blcr.AppendElement(1))[0], blcr.AppendElement(2));

            SofistikCrvPtIdx.Add(SofistikCrvPtIdx.Branch(blcr.AppendElement(0))[0], blcr.AppendElement(2));
            SofistikCrvPtIdx.Add(SofistikCrvPtIdx.Branch(blcr.AppendElement(1))[0], blcr.AppendElement(2));
            //L02
            SofistikCrvPtCoo.Add(SofistikCrvPtCoo.Branch(blcr.AppendElement(0))[1], blcr.AppendElement(3));
            SofistikCrvPtCoo.Add(SofistikCrvPtCoo.Branch(blcr.AppendElement(1))[1], blcr.AppendElement(3));

            SofistikCrvPtIdx.Add(SofistikCrvPtIdx.Branch(blcr.AppendElement(0))[1], blcr.AppendElement(3));
            SofistikCrvPtIdx.Add(SofistikCrvPtIdx.Branch(blcr.AppendElement(1))[1], blcr.AppendElement(3));
            #endregion  BL-CR

            //Curves CL-AR
            #region  CL-AR

            GH_Path clar = new GH_Path(4);
            //L01
            SofistikCrvPtCoo.Add(SofistikCrvPtCoo.Branch(clar.AppendElement(0))[0], clar.AppendElement(2));
            SofistikCrvPtCoo.Add(SofistikCrvPtCoo.Branch(clar.AppendElement(1))[0], clar.AppendElement(2));

            SofistikCrvPtIdx.Add(SofistikCrvPtIdx.Branch(clar.AppendElement(0))[0], clar.AppendElement(2));
            SofistikCrvPtIdx.Add(SofistikCrvPtIdx.Branch(clar.AppendElement(1))[0], clar.AppendElement(2));
            //L02
            SofistikCrvPtCoo.Add(SofistikCrvPtCoo.Branch(clar.AppendElement(0))[1], clar.AppendElement(3));
            SofistikCrvPtCoo.Add(SofistikCrvPtCoo.Branch(clar.AppendElement(1))[1], clar.AppendElement(3));

            SofistikCrvPtIdx.Add(SofistikCrvPtIdx.Branch(clar.AppendElement(0))[1], clar.AppendElement(3));
            SofistikCrvPtIdx.Add(SofistikCrvPtIdx.Branch(clar.AppendElement(1))[1], clar.AppendElement(3));
            #endregion  CL-AR

            #endregion Sofistik Curves

        }

        public void SofistikCreateSurfaces()
        {

            List<Point3d> planarCrv0 = new List<Point3d>(AllSofistikPoints.Branch(0));
            planarCrv0.Add(planarCrv0[0]);

            List<Point3d> planarCrv1 = new List<Point3d>(AllSofistikPoints.Branch(1));
            planarCrv1.Add(planarCrv1[0]);

            SofistikSurfaces.Add(
                Brep.CreatePlanarBreps(PolylineCurve.CreateControlPointCurve(planarCrv0, 1)), new GH_Path(0)
                );
            SofistikSurfacesIdx.Add(0, new GH_Path(0));


            SofistikSurfaces.Add(
                 Brep.CreatePlanarBreps(PolylineCurve.CreateControlPointCurve(planarCrv1, 1)), new GH_Path(1)
                );
            SofistikSurfacesIdx.Add(1, new GH_Path(1));


            for (int i = 0; i < Curves.BranchCount; i++)
            {
                SofistikSurfaces.Add(Brep.CreateFromLoft(Curves.Branch(i), Point3d.Unset, Point3d.Unset, LoftType.Normal, false), new GH_Path(i + 2));
                SofistikSurfacesIdx.Add(i + 2, new GH_Path(i + 2));

            }

        }

    }
}
