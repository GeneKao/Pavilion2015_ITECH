using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;

using ICD;

namespace Pavillion2015
{
    public class GHC_PointDistributor2 : GH_Component
    {
        bool iReset = false;
        int iRandomSeed = -1;
        bool iPlay = false;
        double iRelaxationRate = double.NaN;
        bool iAddPointsAutomatically = false;
        bool iAddPointsManually = false;
        List<Circle> iApicalDiscs = new List<Circle>();
        Curve iBoundaryCurve = null;
        double iPointAdditionRate = double.NaN;
        int iMaxPointCount = 0;
        double iMinSeparation = double.NaN;
        double iMaxSeparation = double.NaN;
        double iGrowthRate = double.NaN;
        int iNumberColumnBasePlate = 1;
        int iGrowthOption = 1;

        string oInfo = string.Empty;
        List<Point3d> oDebugPoints1 = null;
        List<Point3d> oDebugPoints2 = null;
        List<Curve> oDebugCurves1 = null;
        List<Curve> oDebugCurves2 = null;
        List<Vector3d> oDebugVectors1 = null;
        List<Vector3d> oDebugVectors2 = null;
        List<double> oDebugNumbers1 = null;
        List<double> oDebugNumbers2 = null;

        List<PointAgent> pointAgents = null;
        int iteration = 0;

        //----------------------------------------------------------------------------------------

        public GHC_PointDistributor2()
            : base(
                "Point Distributor 2",
                "Point Distributor 2",
                "Point Distributor 2",
                "Pavillion 2015",
                "Pavillion 2015")
        {
        }


        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddBooleanParameter("Reset", "Reset - 00", "Reset", GH_ParamAccess.item, true);
            pManager.AddIntegerParameter("Random Seed", "Random Seed - 01", "Random Seed", GH_ParamAccess.item, -1);
            pManager.AddBooleanParameter("Play", "Play - 02", "Play", GH_ParamAccess.item, false);
            pManager.AddNumberParameter("Relaxation Rate", "Relaxation Rate -03", "Relaxation Rate", GH_ParamAccess.item, 0.4);
            pManager.AddBooleanParameter("Add Points Automatically", "Add Points Automatically - 04", "Add Points Automatically", GH_ParamAccess.item, true);
            pManager.AddBooleanParameter("Add Points Manually", "Add Points Manually - 05", "Add Points Automatically", GH_ParamAccess.item, true);
            pManager.AddCircleParameter("Apical Discs", "Apical Discs - 06", "Apical Discs", GH_ParamAccess.list);
            pManager.AddBooleanParameter("Fix Apical Discs", "Fix Apical Discs - 07", "Fix Apical Discs", GH_ParamAccess.item);
            pManager.AddCurveParameter("Boundary Curve", "Boundary Curve - 08", "Boundary Curve", GH_ParamAccess.item);
            pManager.AddNumberParameter("Point Addition Rate", "Point Addition Rate - 09", "Point Addition Rate", GH_ParamAccess.item, 1.0);
            pManager.AddIntegerParameter("Max. Point Count", "Max. Point Count - 10", "Max. Point Count", GH_ParamAccess.item, 300);
            pManager.AddNumberParameter("Min. Separation", "Min. Separation - 11", "Min. Separation", GH_ParamAccess.item, 0.1);
            pManager.AddNumberParameter("Max. Separation", "Max. Separation - 12", "Max. Separation", GH_ParamAccess.item, 1.0);
            pManager.AddNumberParameter("Growth Rate", "Growth Rate - 13", "Separation Increment Rate", GH_ParamAccess.item, 1.0);
            pManager.AddIntegerParameter("Number of Plates at Column Base", "Number of Plates at Column Base - 14", "Number of Plates at Column Base", GH_ParamAccess.item, 1);
            pManager.AddIntegerParameter("Growth Option", "Growth Option - 14", "1 = Half growth, 2 = Full growth", GH_ParamAccess.item, 1);
        }


        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddTextParameter("Info", "00 - Info", "Information", GH_ParamAccess.item);
            pManager.AddGenericParameter("Debug 1", "01 - Debug 1", "This output is reserved for debugging", GH_ParamAccess.list);
            pManager.AddGenericParameter("Debug 2", "02 - Debug 2", "This output is reserved for debugging", GH_ParamAccess.list);
            pManager.AddGenericParameter("Debug 3", "03 - Debug 3", "This output is reserved for debugging", GH_ParamAccess.list);
            pManager.AddGenericParameter("Debug 4", "04 - Debug 4", "This output is reserved for debugging", GH_ParamAccess.list);
            pManager.AddGenericParameter("Debug 5", "05 - Debug 5", "This output is reserved for debugging", GH_ParamAccess.list);
            pManager.AddPointParameter("Points", "06 - Points", "Points", GH_ParamAccess.list);
        }


        protected override void BeforeSolveInstance()
        {
            oInfo = string.Empty;
            oDebugPoints1 = new List<Point3d>();
            oDebugPoints2 = new List<Point3d>();
            oDebugCurves1 = new List<Curve>();
            oDebugCurves2 = new List<Curve>();
            oDebugVectors1 = new List<Vector3d>();
            oDebugVectors2 = new List<Vector3d>();
            oDebugNumbers1 = new List<double>();
            oDebugNumbers2 = new List<double>();
        }


        protected override void SolveInstance(IGH_DataAccess DA)
        {
            DA.GetData<bool>("Reset", ref iReset);
            DA.GetData<int>("Number of Plates at Column Base", ref iNumberColumnBasePlate);

            if (iReset)
            {
                DA.GetData<int>("Random Seed", ref iRandomSeed);
                iApicalDiscs.Clear();
                DA.GetDataList<Circle>("Apical Discs", iApicalDiscs);
                iteration = 0;
                Random random = iRandomSeed == 0 ? new Random() : new Random(iRandomSeed);

                pointAgents = new List<PointAgent>();
                foreach (Circle apicalDisc in iApicalDiscs)
                {

                    Curve apicalDiscCurve = apicalDisc.ToNurbsCurve();
                    Point3d[] divisionPoints = new Point3d[1];
                    pointAgents.Add(new PointAgent(apicalDisc));
                    apicalDiscCurve.DivideByCount(iNumberColumnBasePlate, true, out divisionPoints);
                    foreach (Point3d divisionPoint in divisionPoints)
                    {
                        pointAgents.Add(new PointAgent(divisionPoint, iMinSeparation));
                    }
                }

                goto Conclusion;
            }

            DA.GetData<bool>("Play", ref iPlay);
            if (!iPlay) goto Conclusion;

            iteration++;
            ExpireSolution(true);

            if (iApicalDiscs.Count == 0) goto Conclusion;

            ///////////////////////////////////////////////////         

            DA.GetData<Curve>("Boundary Curve", ref iBoundaryCurve);
            DA.GetData<double>("Min. Separation", ref iMinSeparation);
            DA.GetData<double>("Max. Separation", ref iMaxSeparation);
            DA.GetData<double>("Growth Rate", ref iGrowthRate);

            DA.GetData<int>("Growth Option", ref iGrowthOption);

            if (iGrowthOption != 1 && iGrowthOption != 2)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "The Growth Option must be set to 1 or 2");
                return;
            }

            Curve boundaryCurveXY = Curve.ProjectToPlane(iBoundaryCurve, Plane.WorldXY);


            // ==========================================================================================================
            // Generate new point agents
            // ==========================================================================================================

            DA.GetData<bool>("Add Points Automatically", ref iAddPointsAutomatically);

            if (iAddPointsAutomatically)
            {
                DA.GetData<double>("Point Addition Rate", ref iPointAdditionRate);
                DA.GetData<int>("Max. Point Count", ref iMaxPointCount);

                if (pointAgents.Count < iMaxPointCount)
                {
                    if (Utils.GetRandomDouble() < iPointAdditionRate)
                    {
                        Circle apicalDisc = iApicalDiscs[Utils.GetRandomInteger(0, iApicalDiscs.Count)];
                        double randomAngle = Utils.GetRandomDouble(0.0, Math.PI * 2.0);
                        Point3d center = apicalDisc.Center +
                            apicalDisc.Radius * new Point3d(Math.Cos(randomAngle), Math.Sin(randomAngle), 0.0);
                        pointAgents.Add(new PointAgent(center, iMinSeparation));
                    }
                }
            }


            // ==========================================================================================================
            // Compute the offset for each agent
            // (i.e. how much the agents push each other away when they are too close to each other)
            // ==========================================================================================================

            List<Vector3d> offsets = new List<Vector3d>();
            for (int i = 0; i < pointAgents.Count; i++)
                offsets.Add(Vector3d.Zero);

            DA.GetData<double>("Relaxation Rate", ref iRelaxationRate);

            for (int i = 0; i < pointAgents.Count; i++)
            {
                for (int j = i + 1; j < pointAgents.Count; j++)
                {
                    double d = Utils.Distance(pointAgents[i].Position, pointAgents[j].Position);
                    double R = pointAgents[i].Separation + pointAgents[j].Separation;
                    if (d < R)
                    {
                        Vector3d offset = (pointAgents[j].Position - pointAgents[i].Position) / d;
                        offset *= (R - d) * 0.5 * iRelaxationRate;
                        offsets[i] -= offset;
                        offsets[j] += offset;
                    }
                }
            }


            // ==========================================================================================================
            // Move the point agents to their new positions
            // ==========================================================================================================

            for (int i = iApicalDiscs.Count * iNumberColumnBasePlate + 1; i < pointAgents.Count; i++)
            {
                Point3d newPosition = pointAgents[i].Position + offsets[i];
                newPosition.Z = 0.0;

                // if the agent is outside the boundary then move it to the closet point on the boundary 
                if (boundaryCurveXY.Contains(newPosition, Plane.WorldXY) == PointContainment.Outside)
                {
                    double t;
                    boundaryCurveXY.ClosestPoint(newPosition, out t);

                    newPosition = boundaryCurveXY.PointAt(t);
                }

                newPosition.Z = 0.0;
                pointAgents[i].Position = newPosition;
            }


            // ==========================================================================================================
            // Update the separation strength for each point agent 
            // ==========================================================================================================

            for (int i = iApicalDiscs.Count; i < pointAgents.Count; i++)
            {
                PointAgent pointAgent = pointAgents[i];

                double minDist = 9999.0;

                foreach (Circle apicalDisc in iApicalDiscs)
                {
                    double d = Utils.Distance(apicalDisc.Center, pointAgent.Position) - apicalDisc.Radius;
                    if (d < minDist) minDist = d;
                }

                if (iGrowthOption == 2)
                {
                    double t;
                    boundaryCurveXY.ClosestPoint(pointAgent.Position, out t);
                    double d = Utils.Distance(pointAgent.Position, boundaryCurveXY.PointAt(t));
                    if (d < minDist) minDist = d;
                }

                pointAgent.Separation = iMinSeparation + (iMaxSeparation - iMinSeparation) *
                        (1 - Math.Pow(0.9, minDist * iGrowthRate));
            }


            // ==========================================================================================================
        // Conclusion
        // ==========================================================================================================

            Conclusion:

            oInfo += "Iteration: " + iteration.ToString();
            oInfo += "\nPoint Count: " + pointAgents.Count.ToString();

            List<Point3d> points = new List<Point3d>();
            List<double> separations = new List<double>();

            foreach (PointAgent pointAgent in pointAgents)
            {
                points.Add(pointAgent.Position);
                separations.Add(pointAgent.Separation);
            }

            DA.SetData(0, oInfo);
            DA.SetDataList(5, separations);
            DA.SetDataList(6, points);
        }

        protected override System.Drawing.Bitmap Icon { get { return null; } }
        public override Guid ComponentGuid { get { return new Guid("{f558159d-082f-4411-987f-bbfe40e8d279}"); } }
    }
}