using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;
using Rhino.Geometry.Intersect;

using ICD;

namespace Pavillion2015
{
    public class SpringMesh
    {
        public List<Vertex> Vertices = new List<Vertex>();
        public List<Edge> Edges = new List<Edge>();
        public List<Triangle> Triangles = new List<Triangle>();

        public int Integrator = 1; // 1 = Euler, 2 = Verlet
        public double DeltaT = 0.02;
        public double Damping = 0.2;
        public double GravityScale = 1;

        public double RestLengthScale = 1.0;
        public double RestLengthOffset = 1.0;
        public double Stiffness = 0.0;
        public double BendingStiffness = 0.0;

        public bool EnableBoundaryBendingStiffness = false;
        public double BoundaryBendingStiffness = 1.0;

        public bool EnableColumnBendingStiffness = false;
        public double ColumnnBendingStiffness = 1.0;

        public bool Is2D = false;

        //added by Gene
        public int EnableEquilateralTriangle = 0;
        public double EquilateralStrength = 1.0; 

        private int iteration = 0;


        public SpringMesh()
        {
        }


        public SpringMesh(SpringMesh springMesh)
        {
            Vertices = new List<Vertex>();
            foreach (Vertex vertex in springMesh.Vertices)
                Vertices.Add(new Vertex(vertex));

            Edges = new List<Edge>();
            foreach (Edge edge in springMesh.Edges)
                Edges.Add(new Edge(edge));

            Triangles = new List<Triangle>();
            foreach (Triangle triangle in springMesh.Triangles)
                Triangles.Add(new Triangle(triangle));

            iteration = 0;
        }


        public void Update()
        {
            // =======================================================================
            // Apply gravity and damping
            // =======================================================================

            foreach (Vertex vertex in Vertices)
                vertex.Force = 
                    new Vector3d(0.0, 0.0, vertex.Mass * -9.8 * GravityScale) -
                    vertex.Velocity * Damping;


            // =======================================================================
            // Apply spring force to the two vertices at the ends of each edge
            // =======================================================================

            foreach (Edge edge in Edges)
            {
                Vertex firstVertex = Vertices[edge.FirstVertexIndex];
                Vertex secondVertex = Vertices[edge.SecondVertexIndex];

                Vector3d springForce = firstVertex.Position - secondVertex.Position;
                double length = springForce.Length;
                springForce = springForce / length * (edge.Stiffness + Stiffness) * (length - (edge.RestLength * RestLengthScale + RestLengthOffset)); // Hooke's Law

                firstVertex.Force -= springForce;
                secondVertex.Force += springForce;
            }


            // =======================================================================
            // Apply bending resistance force
            // =======================================================================
            foreach (Edge edge in Edges)
            {
                if (edge.IsBoundaryEdge)
                {
                    Vertex v = Vertices[edge.FirstAdjacentVertexIndex];
                    if (v.IsFixed && edge.SecondAdjacentVertexIndex != -1)
                        v = Vertices[edge.SecondAdjacentVertexIndex];

                    Point3d M = Utils.ClosestPointOnLine(v.Position, Vertices[edge.FirstVertexIndex].Position, Vertices[edge.SecondVertexIndex].Position);

                    Vector3d force = M - v.Position;
                    force.Z = 0.0;

                    v.Force += force * BoundaryBendingStiffness;
                }
                else if (edge.IsColumnEdge)
                {
                    Vertex v = Vertices[edge.FirstAdjacentVertexIndex];
                    if (v.IsFixed && edge.SecondAdjacentVertexIndex != -1)
                        v = Vertices[edge.SecondAdjacentVertexIndex];

                    Point3d M = Utils.ClosestPointOnLine(v.Position, Vertices[edge.FirstVertexIndex].Position, Vertices[edge.SecondVertexIndex].Position);

                    Vector3d force = M - v.Position;
                    force.Z = 0.0;

                    v.Force += force * ColumnnBendingStiffness;
                }

                else
                {
                    Vertex A = Vertices[edge.FirstVertexIndex];
                    Vertex B = Vertices[edge.SecondVertexIndex];
                    Vertex M = Vertices[edge.FirstAdjacentVertexIndex];
                    Vertex N = Vertices[edge.SecondAdjacentVertexIndex];

                    Vector3d M_M = M.Position - Utils.ClosestPointOnLine(M.Position, A.Position, B.Position);
                    Vector3d N_N = N.Position - Utils.ClosestPointOnLine(N.Position, A.Position, B.Position);

                    M_M.Unitize();
                    N_N.Unitize();

                    //Vector3d F = (M_M + N_N) * (edge.BendingStiffness + BendingStiffness);

                    Vector3d F = (M_M + N_N);
                    F.Unitize();
                    F *= (edge.BendingStiffness + BendingStiffness)
                         * 2.0 * (Math.Cos(0.5 * Utils.AngleBetweenTwoTriangles(A.Position, B.Position, M.Position, N.Position)) - Math.Cos(0.5 * edge.RestAngle));

                    A.Force += F;
                    B.Force += F;
                    M.Force -= F;
                    N.Force -= F;
                }
            }

            if (Is2D)
                for (int i = 1; i < Vertices.Count - 1; i++)
                {
                    Point3d A = Vertices[i].Position;
                    Point3d B = Vertices[i - 1].Position;
                    Point3d C = Vertices[i + 1].Position;

                    Vector3d AB = (B - A); AB.Unitize();
                    Vector3d AC = (C - A); AC.Unitize();

                    Vector3d F = (((A + AB) + (A + AC)) - 2.0 * A);

                    Vertices[i].Force += 2 * F;
                    Vertices[i - 1].Force -= F;
                    Vertices[i + 1].Force -= F;
                }


            // ============================================================
            // Apply Equilateral Triangle  //// Added by Gene
            // ============================================================
            if (EnableEquilateralTriangle == 0)
                equilateralTriangleLength();
            else if (EnableEquilateralTriangle == 1)
                equilateralTriangleAngle();
            else if (EnableEquilateralTriangle == 2)
            {
                equilateralTriangleLength();
                equilateralTriangleAngle();
            }
            else
            {
            }

            // ============================================================
            // Update position and velocity of the vertices
            // ============================================================

            if (Integrator == 1) // Semi-Implicit Euler Integration
            {
                foreach (Vertex vertex in Vertices)
                {
                    if (vertex.IsFixed) continue;

                    vertex.Acceleration = vertex.Force / vertex.Mass;
                    vertex.Velocity += vertex.Acceleration * DeltaT;
                    vertex.Position += vertex.Velocity * DeltaT;
                }
            }
            else if (Integrator == 2) // Velocity Verlet Integration
            {
                foreach (Vertex vertex in Vertices)
                {
                    if (vertex.IsFixed) continue;

                    vertex.Acceleration = vertex.Force / vertex.Mass;

                    if (iteration > 0)
                        vertex.Velocity += 0.5 * vertex.Acceleration * DeltaT;
                    
                    vertex.Position += vertex.Velocity * DeltaT + 0.5 * vertex.Acceleration * DeltaT * DeltaT;

                    vertex.Velocity += 0.5 * DeltaT * vertex.Acceleration;
                }
            }

            iteration++;
        }

        //==================================================
        // added by Gene 
        // Equilateral Triangle
        //==================================================
        private void equilateralTriangleLength()
        {
            foreach (Triangle triangle in Triangles)
            {
                Vertex firstVertex = Vertices[triangle.FirstVertexIndex];
                Vertex secondVertex = Vertices[triangle.SecondVertexIndex];
                Vertex thirdVertex = Vertices[triangle.ThirdVertexIndex];

                double lengthfirstSecond = ICD.Utils.Distance(firstVertex.Position, secondVertex.Position);
                double lengthSecondThird = ICD.Utils.Distance(secondVertex.Position, thirdVertex.Position);
                double lengthThirdFirst = ICD.Utils.Distance(thirdVertex.Position, firstVertex.Position);

                double avg = 0.3333 * (lengthfirstSecond + lengthSecondThird + lengthThirdFirst);

                if (lengthfirstSecond - lengthSecondThird != 0 ||
                    lengthSecondThird - lengthThirdFirst != 0 || lengthThirdFirst - lengthfirstSecond != 0)
                {
                    Vector3d ForceFirst01 = (thirdVertex.Position - firstVertex.Position);
                    Vector3d ForceFirst02 = (secondVertex.Position - firstVertex.Position);
                    ForceFirst01.Unitize();
                    ForceFirst02.Unitize();
                    Vector3d ForceFirst = ForceFirst01 * (lengthThirdFirst - avg) + ForceFirst02 * (lengthfirstSecond - avg);

                    Vector3d ForceSecond01 = (thirdVertex.Position - secondVertex.Position);
                    Vector3d ForceSecond02 = (firstVertex.Position - secondVertex.Position);
                    ForceSecond01.Unitize();
                    ForceSecond02.Unitize();
                    Vector3d ForceSecond = ForceSecond01 * (lengthSecondThird - avg) + ForceSecond02 * (lengthfirstSecond - avg);

                    Vector3d ForceThird01 = (secondVertex.Position - thirdVertex.Position);
                    Vector3d ForceThird02 = (firstVertex.Position - thirdVertex.Position);
                    ForceThird01.Unitize();
                    ForceThird02.Unitize();

                    Vector3d ForceThird = ForceThird01 * (lengthSecondThird - avg) + ForceThird02 * (lengthThirdFirst - avg);

                    ForceFirst *= EquilateralStrength;
                    ForceSecond *= EquilateralStrength;
                    ForceThird *= EquilateralStrength;

                    firstVertex.Force += ForceFirst;
                    secondVertex.Force += ForceSecond;
                    thirdVertex.Force += ForceThird;
                }
            }
        }

        private void equilateralTriangleAngle()
        {
            for (int i = 0; i < Triangles.Count; i++)
            {
                Vertex firstVertex = Vertices[Triangles[i].FirstVertexIndex];
                Vertex secondVertex = Vertices[Triangles[i].SecondVertexIndex];
                Vertex thirdVertex = Vertices[Triangles[i].ThirdVertexIndex];

                double angle01 = ICD.Utils.AngleBetweenTwoVectors(
                    secondVertex.Position - firstVertex.Position, thirdVertex.Position - firstVertex.Position);
                double angle02 = ICD.Utils.AngleBetweenTwoVectors(
                    firstVertex.Position - secondVertex.Position, thirdVertex.Position - secondVertex.Position);
                double angle03 = ICD.Utils.AngleBetweenTwoVectors(
                    firstVertex.Position - thirdVertex.Position, secondVertex.Position - thirdVertex.Position);

                Vector3d triangleNormal = ComputeTriangleNormal(i);

                Vector3d F31 = Vector3d.CrossProduct(firstVertex.Position - thirdVertex.Position, triangleNormal);
                Vector3d F21 = Vector3d.CrossProduct(triangleNormal, firstVertex.Position - secondVertex.Position);
                Vector3d F12 = Vector3d.CrossProduct(secondVertex.Position - firstVertex.Position, triangleNormal);
                Vector3d F32 = Vector3d.CrossProduct(triangleNormal, secondVertex.Position - thirdVertex.Position);
                Vector3d F13 = Vector3d.CrossProduct(triangleNormal, thirdVertex.Position - firstVertex.Position);
                Vector3d F23 = Vector3d.CrossProduct(thirdVertex.Position - secondVertex.Position, triangleNormal);

                F31.Unitize();
                F21.Unitize();
                F12.Unitize();
                F32.Unitize();
                F13.Unitize();
                F23.Unitize();

                F31 *= (0.3333 * Math.PI - angle01) * EquilateralStrength;
                F21 *= (0.3333 * Math.PI - angle01) * EquilateralStrength;
                F12 *= (0.3333 * Math.PI - angle02) * EquilateralStrength;
                F32 *= (0.3333 * Math.PI - angle02) * EquilateralStrength;
                F13 *= (0.3333 * Math.PI - angle03) * EquilateralStrength;
                F23 *= (0.3333 * Math.PI - angle03) * EquilateralStrength;

                firstVertex.Force += F31 + F21;
                secondVertex.Force += F12 + F32;
                thirdVertex.Force += F13 + F23;
            }
        }
        //================================================== 

        public Point3d ComputeTriangleCentroid(int triangleIndex)
        {
            Triangle triangle = Triangles[triangleIndex];
            return 0.333333 * (
                Vertices[triangle.FirstVertexIndex].Position + 
                Vertices[triangle.SecondVertexIndex].Position + 
                Vertices[triangle.ThirdVertexIndex].Position
                );
        }


        public Point3d ComputeCircumscribedCircleCenter(int triangleIndex)
        {
            Triangle triangle = Triangles[triangleIndex];

            Circle circle = new Circle(
                Vertices[triangle.FirstVertexIndex].Position,
                Vertices[triangle.SecondVertexIndex].Position,
                Vertices[triangle.ThirdVertexIndex].Position
                );

            return circle.Center;
        }

        //--------------------------------------------------------------------------------------------------------------------------
        // Added by Gene
        public Circle ComputeDoubleLayerCircumscribedCircle(int triangleIndex, double scaleFactor, List<double> vertexThickness, bool flip, double flipDist)
        {
            List<Vector3d> normalList = ComputeVertexNormals();

            Triangle triangle = Triangles[triangleIndex];

            Point3d newVertexPos01;
            Point3d newVertexPos02;
            Point3d newVertexPos03;

            if (!flip)
            {
                newVertexPos01 = Vertices[triangle.FirstVertexIndex].Position +
                        normalList[triangle.FirstVertexIndex] * vertexThickness[triangle.FirstVertexIndex] * (1 - flipDist);
                newVertexPos02 = Vertices[triangle.SecondVertexIndex].Position +
                        normalList[triangle.SecondVertexIndex] * vertexThickness[triangle.SecondVertexIndex] * (1 - flipDist);
                newVertexPos03 = Vertices[triangle.ThirdVertexIndex].Position +
                        normalList[triangle.ThirdVertexIndex] * vertexThickness[triangle.ThirdVertexIndex] * (1 - flipDist);
            }
            else
            {
                newVertexPos01 = Vertices[triangle.FirstVertexIndex].Position -
                        normalList[triangle.FirstVertexIndex] * vertexThickness[triangle.FirstVertexIndex] * (1 + flipDist);
                newVertexPos02 = Vertices[triangle.SecondVertexIndex].Position -
                        normalList[triangle.SecondVertexIndex] * vertexThickness[triangle.SecondVertexIndex] * (1 + flipDist);
                newVertexPos03 = Vertices[triangle.ThirdVertexIndex].Position -
                        normalList[triangle.ThirdVertexIndex] * vertexThickness[triangle.ThirdVertexIndex] * (1 + flipDist);
            }
            Plane p = new Plane(newVertexPos01, newVertexPos02, newVertexPos03);

            Circle c = new Circle(
                newVertexPos01, newVertexPos02, newVertexPos03
                );

            return new Circle(p, c.Center, c.Radius * scaleFactor);
        }

        public Point3d ComputeTPI(int triangleIndex)
        {
            List<Vector3d> normalList = ComputeVertexNormals();

            Triangle triangle = Triangles[triangleIndex];

            Plane plane01 = new Plane(Vertices[triangle.FirstVertexIndex].Position, normalList[triangle.FirstVertexIndex]);
            Plane plane02 = new Plane(Vertices[triangle.SecondVertexIndex].Position, normalList[triangle.SecondVertexIndex]);
            Plane plane03 = new Plane(Vertices[triangle.ThirdVertexIndex].Position, normalList[triangle.ThirdVertexIndex]);

            Point3d intersectionPoint;
            Intersection.PlanePlanePlane(plane01, plane02, plane03, out intersectionPoint);

            return intersectionPoint;
        }
        //--------------------------------------------------------------------------------------------------------------------------
        // Added by Gene
        public Point3d ComputeDoubleLayerTPI(int triangleIndex, List<double> vertexThickness, bool flip, double flipDist)
        {
            List<Vector3d> normalList = ComputeVertexNormals();

            Triangle triangle = Triangles[triangleIndex];

            Point3d newVertexPos01;
            Point3d newVertexPos02;
            Point3d newVertexPos03;

            if (!flip)
            {
                newVertexPos01 = Vertices[triangle.FirstVertexIndex].Position +
                        normalList[triangle.FirstVertexIndex] * vertexThickness[triangle.FirstVertexIndex] * (1 - flipDist);
                newVertexPos02 = Vertices[triangle.SecondVertexIndex].Position +
                        normalList[triangle.SecondVertexIndex] * vertexThickness[triangle.SecondVertexIndex] * (1 - flipDist);
                newVertexPos03 = Vertices[triangle.ThirdVertexIndex].Position +
                        normalList[triangle.ThirdVertexIndex] * vertexThickness[triangle.ThirdVertexIndex] * (1 - flipDist);
            }
            else
            {
                newVertexPos01 = Vertices[triangle.FirstVertexIndex].Position -
                        normalList[triangle.FirstVertexIndex] * vertexThickness[triangle.FirstVertexIndex] * (1 + flipDist);
                newVertexPos02 = Vertices[triangle.SecondVertexIndex].Position -
                        normalList[triangle.SecondVertexIndex] * vertexThickness[triangle.SecondVertexIndex] * (1 + flipDist);
                newVertexPos03 = Vertices[triangle.ThirdVertexIndex].Position -
                        normalList[triangle.ThirdVertexIndex] * vertexThickness[triangle.ThirdVertexIndex] * (1 + flipDist);
            }
            Plane plane01 = new Plane(newVertexPos01, normalList[triangle.FirstVertexIndex]);
            Plane plane02 = new Plane(newVertexPos02, normalList[triangle.SecondVertexIndex]);
            Plane plane03 = new Plane(newVertexPos03, normalList[triangle.ThirdVertexIndex]);

            Point3d intersectionPoint;
            Intersection.PlanePlanePlane(plane01, plane02, plane03, out intersectionPoint);

            return intersectionPoint;
        }

        //--------------------------------------------------------------------------------------------------------------------------
        // Added by Gene
        public Circle ComputeIncircle(int triangleIndex, double scaleFactor)
        {
            Triangle triangle = Triangles[triangleIndex];

            Point3d v1 = Vertices[triangle.FirstVertexIndex].Position;
            Point3d v2 = Vertices[triangle.SecondVertexIndex].Position;
            Point3d v3 = Vertices[triangle.ThirdVertexIndex].Position;

            Plane p = new Plane(v1, v2, v3);
            Circle c = Circle.TryFitCircleTTT(new LineCurve(v2, v3), new LineCurve(v3, v1), new LineCurve(v1, v2), 0.001, 0.001, 0.001);

            return new Circle(p, c.Center, c.Radius * scaleFactor);
        }

        // Double layers
        public Circle ComputeDoubleLayerIncircle(int triangleIndex, double scaleFactor, List<double> vertexThickness, bool flip, double flipDist)
        {
            List<Vector3d> normalList = ComputeVertexNormals();

            Triangle triangle = Triangles[triangleIndex];

            Point3d newVertexPos01;
            Point3d newVertexPos02;
            Point3d newVertexPos03;

            if (!flip)
            {
                newVertexPos01 = Vertices[triangle.FirstVertexIndex].Position +
                        normalList[triangle.FirstVertexIndex] * vertexThickness[triangle.FirstVertexIndex] * (1 - flipDist);
                newVertexPos02 = Vertices[triangle.SecondVertexIndex].Position +
                        normalList[triangle.SecondVertexIndex] * vertexThickness[triangle.SecondVertexIndex] * (1 - flipDist);
                newVertexPos03 = Vertices[triangle.ThirdVertexIndex].Position +
                        normalList[triangle.ThirdVertexIndex] * vertexThickness[triangle.ThirdVertexIndex] * (1 - flipDist);
            }
            else
            {
                newVertexPos01 = Vertices[triangle.FirstVertexIndex].Position -
                        normalList[triangle.FirstVertexIndex] * vertexThickness[triangle.FirstVertexIndex] * (1 + flipDist);
                newVertexPos02 = Vertices[triangle.SecondVertexIndex].Position -
                        normalList[triangle.SecondVertexIndex] * vertexThickness[triangle.SecondVertexIndex] * (1 + flipDist);
                newVertexPos03 = Vertices[triangle.ThirdVertexIndex].Position -
                        normalList[triangle.ThirdVertexIndex] * vertexThickness[triangle.ThirdVertexIndex] * (1 + flipDist);
            }
            Plane p = new Plane(newVertexPos01, newVertexPos02, newVertexPos03);
            Circle c = Circle.TryFitCircleTTT(
                new LineCurve(newVertexPos02, newVertexPos03),
                new LineCurve(newVertexPos03, newVertexPos01),
                new LineCurve(newVertexPos01, newVertexPos02), 0.001, 0.001, 0.001);

            return new Circle(p, c.Center, c.Radius * scaleFactor);
        }

        public Plane ComputePlane(int triangleIndex)
        {
            Triangle triangle = Triangles[triangleIndex];

            Point3d v1 = Vertices[triangle.FirstVertexIndex].Position;
            Point3d v2 = Vertices[triangle.SecondVertexIndex].Position;
            Point3d v3 = Vertices[triangle.ThirdVertexIndex].Position;

            return new Plane(v1, v2, v3);
        }

        public Plane ComputeDoubleLayerPlane(int triangleIndex, List<double> vertexThickness, bool flip, double flipDist)
        {
            List<Vector3d> normalList = ComputeVertexNormals();

            Triangle triangle = Triangles[triangleIndex];

            Point3d newVertexPos01;
            Point3d newVertexPos02;
            Point3d newVertexPos03;

            if (!flip)
            {
                newVertexPos01 = Vertices[triangle.FirstVertexIndex].Position +
                        normalList[triangle.FirstVertexIndex] * vertexThickness[triangle.FirstVertexIndex] * (1 - flipDist);
                newVertexPos02 = Vertices[triangle.SecondVertexIndex].Position +
                        normalList[triangle.SecondVertexIndex] * vertexThickness[triangle.SecondVertexIndex] * (1 - flipDist);
                newVertexPos03 = Vertices[triangle.ThirdVertexIndex].Position +
                        normalList[triangle.ThirdVertexIndex] * vertexThickness[triangle.ThirdVertexIndex] * (1 - flipDist);
            }
            else
            {
                newVertexPos01 = Vertices[triangle.FirstVertexIndex].Position -
                        normalList[triangle.FirstVertexIndex] * vertexThickness[triangle.FirstVertexIndex] * (1 + flipDist);
                newVertexPos02 = Vertices[triangle.SecondVertexIndex].Position -
                        normalList[triangle.SecondVertexIndex] * vertexThickness[triangle.SecondVertexIndex] * (1 + flipDist);
                newVertexPos03 = Vertices[triangle.ThirdVertexIndex].Position -
                        normalList[triangle.ThirdVertexIndex] * vertexThickness[triangle.ThirdVertexIndex] * (1 + flipDist);
            }
            return new Plane(newVertexPos01, newVertexPos02, newVertexPos03);
        }

        //--------------------------------------------------------------------------------------------------------------------------


        public List<Point3d> ComputeTriangleThreePts(int triangleIndex, double factor)
        {

            Triangle triangle = Triangles[triangleIndex];

            Point3d centroid = ComputeTriangleCentroid(triangleIndex);

            Point3d A = centroid * factor + Vertices[triangle.FirstVertexIndex].Position * (1 - factor);
            Point3d B = centroid * factor + Vertices[triangle.SecondVertexIndex].Position * (1 - factor);
            Point3d C = centroid * factor + Vertices[triangle.ThirdVertexIndex].Position * (1 - factor);

            return new List<Point3d>() { A, B, C };
        }


        public Vector3d ComputeTriangleNormal(int triangleIndex)
        {
            Triangle triangle = Triangles[triangleIndex];
            return Utils.ComputePlaneNormal(
                Vertices[triangle.FirstVertexIndex].Position,
                Vertices[triangle.SecondVertexIndex].Position,
                Vertices[triangle.ThirdVertexIndex].Position
                );
        }


        public List<Vector3d> ComputeTriangleNormals()
        {
            List<Vector3d> triangleNormals = new List<Vector3d>();
            for (int i = 0; i < Triangles.Count; i++)
                triangleNormals.Add(ComputeTriangleNormal(i));
            return triangleNormals;
        }


        public List<Vector3d> ComputeVertexNormals()
        {
            List<Vector3d> vertexNormals = new List<Vector3d>();
            foreach (Vertex vertex in Vertices) vertexNormals.Add(Vector3d.Zero);

            for (int i = 0; i < Triangles.Count; i++)
            {
                Vector3d triangleNormal = ComputeTriangleNormal(i);
                vertexNormals[Triangles[i].FirstVertexIndex] += triangleNormal;
                vertexNormals[Triangles[i].SecondVertexIndex] += triangleNormal;
                vertexNormals[Triangles[i].ThirdVertexIndex] += triangleNormal;
            }

            for (int i = 0; i < Vertices.Count; i++)
            {
                Vector3d vertexNormal = vertexNormals[i];
                vertexNormal.Unitize();
                vertexNormals[i] = vertexNormal;
            }

            return vertexNormals;
        }


        public Mesh ConvertToRhinoMesh()
        {
            Mesh rhinoMesh = new Mesh();

            foreach (Vertex vertex in Vertices)
                rhinoMesh.Vertices.Add(vertex.Position);

            foreach (Triangle triangle in Triangles)
                rhinoMesh.Faces.AddFace(triangle.FirstVertexIndex, triangle.SecondVertexIndex, triangle.ThirdVertexIndex);

            return rhinoMesh;
        }
    
    }


    public class Vertex
    {
        public Point3d Position;
        public Vector3d Velocity;
        public List<int> NeighborVertexIndices;
        public double Mass;
        public bool IsFixed;
        public bool IsBoundaryVertex = false;
        public bool IsColumnVertex = false;

        public Vector3d Force;
        public Vector3d Acceleration;

        public double Factor;

        public Vertex()
        {
        }

        public Vertex(Vertex vertex)
        {
            Position = vertex.Position;
            Velocity = vertex.Velocity;
            NeighborVertexIndices = new List<int>(vertex.NeighborVertexIndices);
            Mass = vertex.Mass;
            IsFixed = vertex.IsFixed;
            IsBoundaryVertex = vertex.IsBoundaryVertex;
            IsColumnVertex = vertex.IsColumnVertex;  

            Force = vertex.Force;
            Acceleration = vertex.Acceleration;
        }


        public Vertex(
            Point3d position,
            Vector3d velocity,
            List<int> neighborVertexIndices,
            double mass = 1.0,
            bool isFixed = false,
            bool isBoundaryVertex = false,
            bool isColumnVertex = false)
        {
            Position = position;
            Velocity = velocity;
            NeighborVertexIndices = neighborVertexIndices;
            Mass = mass;
            IsFixed = isFixed;
            Acceleration = Vector3d.Zero;
            IsBoundaryVertex = isBoundaryVertex;
            IsColumnVertex = isColumnVertex;
        }
    }


    public class Edge
    {
        public int FirstVertexIndex = -1;
        public int SecondVertexIndex = -1;
        public int FirstTriangleIndex = -1;
        public int SecondTriangleIndex = -1;
        public int FirstAdjacentVertexIndex = -1;
        public int SecondAdjacentVertexIndex =- 1;
        public double RestLength = double.NaN;       
        public double Stiffness = double.NaN;
        public double BendingStiffness = double.NaN;
        public double RestAngle = double.NaN;
        public bool IsBoundaryEdge = false;
        public bool IsColumnEdge = false;

        public Edge()
        {
        }


        public Edge(Edge edge)
        {
            FirstVertexIndex = edge.FirstVertexIndex;
            SecondVertexIndex = edge.SecondVertexIndex;
            FirstTriangleIndex = edge.FirstTriangleIndex;
            SecondTriangleIndex = edge.SecondTriangleIndex;
            FirstAdjacentVertexIndex = edge.FirstAdjacentVertexIndex;
            SecondAdjacentVertexIndex = edge.SecondAdjacentVertexIndex;
            RestLength = edge.RestLength;           
            Stiffness = edge.Stiffness;
            RestAngle = edge.RestAngle;
            BendingStiffness = edge.BendingStiffness;
            IsBoundaryEdge = edge.IsBoundaryEdge;
            IsColumnEdge = edge.IsColumnEdge;
        }


        public Edge(int firstVertexIndex, int secondVertexIndex, double restLength, double stiffness = 1.0, double restAngle = Math.PI, double bendingStiffness = 0.0)
        {
            FirstVertexIndex = firstVertexIndex;
            SecondVertexIndex = secondVertexIndex;
            FirstTriangleIndex = -1;
            SecondTriangleIndex = -1;
            FirstAdjacentVertexIndex = -1;
            SecondAdjacentVertexIndex = -1;
            RestLength = restLength;
            Stiffness = stiffness;
            RestAngle = restAngle;
            BendingStiffness = bendingStiffness;
        }
    }


    public class Triangle
    {
        public int FirstVertexIndex = -1;
        public int SecondVertexIndex = -1;
        public int ThirdVertexIndex = -1;

        public int FirstEdgeIndex = -1;
        public int SecondEdgeIndex = -1;
        public int ThirdEdgeIndex = -1;

        // Gene Added =============================
        public int FirstAdjTriIndex = -1;
        public int SecondAdjTriIndex = -1;
        public int ThirdAdjTriIndex = -1;
        // ========================================

        public int FirstSecondAdjacentVertexIndex = -1;
        public int SecondThirdAdjacentVertexIndex = -1;
        public int ThirdFirstAdjacentVertexIndex = -1;

        public Triangle(Triangle triangle)
        {
            FirstVertexIndex = triangle.FirstVertexIndex;
            SecondVertexIndex = triangle.SecondVertexIndex;
            ThirdVertexIndex = triangle.ThirdVertexIndex;
            FirstEdgeIndex = triangle.FirstEdgeIndex;
            SecondEdgeIndex = triangle.SecondEdgeIndex;
            ThirdEdgeIndex = triangle.ThirdEdgeIndex;
            FirstSecondAdjacentVertexIndex = triangle.FirstSecondAdjacentVertexIndex;
            SecondThirdAdjacentVertexIndex = triangle.SecondThirdAdjacentVertexIndex;
            ThirdFirstAdjacentVertexIndex = triangle.ThirdFirstAdjacentVertexIndex;
            // Gene Added
            FirstAdjTriIndex = triangle.FirstAdjTriIndex;
            SecondAdjTriIndex = triangle.SecondAdjTriIndex;
            ThirdAdjTriIndex = triangle.ThirdAdjTriIndex;
        }

        public Triangle(
            int firstVertexIndex, 
            int secondVertexIndex, 
            int thirdVertexIndex
            )
        {
            FirstVertexIndex = firstVertexIndex;
            SecondVertexIndex = secondVertexIndex;
            ThirdVertexIndex = thirdVertexIndex;
        }
    }
    }
