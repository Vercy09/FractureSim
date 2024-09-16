using UnityEngine;
using System.Collections.Generic;
/// <summary>
/// Custom class for grouping necessary math functions.
/// </summary>
public static class MathUtils
{
    /// <summary>
    /// Returns true if the point is above or on the plane, false otherwise. Above is the half space in the direction of the plane normal.
    /// </summary>
    public static bool IsPointAbovePlane(Vector3 testPoint, Vector3 planeNormal, Vector3 planeOrigin)
    {
        // Px - test point, Po - plane origin, n - plane normal
        // (Px-Po) * n = 0 - vector plane equation (points that lie on the plane satisfy the equation)
        // v = Px-Po, the distance of Px from plane is the length of the projection of v on the n
        // proj_n(v) = (v*n)/||n||^2 * n,
        // if n is normalised: (v*n) * n,
        // where v*n is some scalar scaling the normalised n

        // (Px-Po)*n = 0: Point lies directly on the plane
        // (Px-Po)*n > 0: Point is above the plane
        // (Px-Po)*n < 0: Point is below the plane
        return (planeNormal.x * (testPoint.x - planeOrigin.x) + planeNormal.y * (testPoint.y - planeOrigin.y) + planeNormal.z * (testPoint.z - planeOrigin.z)) >= 0;
    }


    /// <summary>
    /// Finds the intersection between a line segment and a plane.
    /// Returns true plus the point of intersection and the parameter of the line segment if an intersection exists.
    /// Returns false if no intersection exists.
    /// </summary>
    public static bool LineSegPlaneIntersection(Vector3 p1,
                                                Vector3 p2,
                                                Vector3 planeNormal,
                                                Vector3 planeOrigin,
                                                out Vector3 intersectPoint,
                                                out float lineSegParam)
    {
        lineSegParam = 0;
        intersectPoint = Vector3.zero;

        if (p1 == p2 || planeNormal == Vector3.zero)
        {
            return false;
        }

        // d - lineSegParam, Po - plane origin, n - plane normal
        // d = (Po-p1)*n/(p2-p1)*n

        lineSegParam = Vector3.Dot(planeOrigin - p1, planeNormal) / Vector3.Dot(p2 - p1, planeNormal);

        // Intersection exists, find the point
        if (lineSegParam >= 0 && lineSegParam <= 1)
        {
            intersectPoint = p1 + (p2 - p1) * lineSegParam;
            return true;
        }
        else
        {
            return false;
        }
    }

    public static bool LineSegLineSegIntersection(Vector2 p1, Vector2 p2, Vector2 q1, Vector2 q2)
    {
        if (p1 == q1 || p1 == q2 || p2 == q1 || p2 == q2)
        {
            return false;
        }
        else if (p1 == p2 || q1 == q2)
        {
            return false;
        }
        else
        {
            Vector2 p12 = new Vector2(p2.x - p1.x, p2.y - p1.y);
            Vector2 q12 = new Vector2(q2.x - q1.x, q2.y - q1.y);

            float p1xq = (p1.x - q1.x) * q12.y - (p1.y - q1.y) * q12.x;
            float p2xq = (p2.x - q1.x) * q12.y - (p2.y - q1.y) * q12.x;
            float q1xp = (q1.x - p1.x) * p12.y - (q1.y - p1.y) * p12.x;
            float q2xp = (q2.x - p1.x) * p12.y - (q2.y - p1.y) * p12.x;

            return ((p1xq >= 0 && p2xq <= 0) || (p1xq <= 0 && p2xq >= 0)) &&
                   ((q1xp >= 0 && q2xp <= 0) || (q1xp <= 0 && q2xp >= 0));
        }
    }

    public static bool IsQuadConvex(Vector2 q1, Vector2 q2, Vector2 q3, Vector2 q4)
    {
        return LineSegLineSegIntersection(q1, q2, q3, q4);
    }

    /// <summary>
    /// Returns a normalized, interpolated 3d vector based on the provided normals and the interpolation factor.
    /// </summary>
    public static Vector3 InterpolateNormals(Vector3 normal1, Vector3 normal2, float s12)
    {
        return (normal1 + s12 * (normal2 - normal1)).normalized;
    }

    /// <summary>
    /// Returns an interpolated 2d vector based on the provided uv coords and the interpolation factor.
    /// </summary>
    public static Vector3 InterpolateUV(Vector3 uv1, Vector3 uv2, float s12)
    {
        return uv1 + s12 * (uv2 - uv1);
    }

    /// <summary>
    /// Returns true if the point p is on the right side or on the line segment p1,p2
    /// </summary>
    public static bool IsPointOnRightSideOfLineSeg(Vector2 p, Vector2 p1, Vector2 p2)
    {
        // Translate the points so that the p1 is on the origin and then
        // use the cross product between p2 and p to find if the angle is CW or CCW
        return ((p2.x - p1.x) * (p.y - p1.y) - (p2.y - p1.y) * (p.x - p1.x)) <= 0;
    }

    /// <summary>
    /// Returns true if the point p is on the left side or on the line segment p1,p2
    /// </summary>
    public static bool IsPointOnLeftSideOfLineSeg(Vector2 p, Vector2 p1, Vector2 p2)
    {
        return ((p2.x - p1.x) * (p.y - p1.y) - (p2.y - p1.y) * (p.x - p1.x)) >= 0;
    }

    public static bool LineLineIntersection(Vector2 p0, Vector2 v0, Vector2 p1, Vector2 v1, out float m0, out float m1)
    {
        var det = (v0.x * v1.y - v0.y * v1.x);

        if (Mathf.Abs(det) < 0.001f)
        {
            m0 = float.NaN;
            m1 = float.NaN;

            return false;
        }
        else
        {
            m0 = ((p0.y - p1.y) * v1.x - (p0.x - p1.x) * v1.y) / det;

            if (Mathf.Abs(v1.x) >= 0.001f)
            {
                m1 = (p0.x + m0 * v0.x - p1.x) / v1.x;
            }
            else
            {
                m1 = (p0.y + m0 * v0.y - p1.y) / v1.y;
            }

            return true;
        }
    }
    public static Vector2 LineLineIntersection(Vector2 p0, Vector2 v0, Vector2 p1, Vector2 v1)
    {
        float m0, m1;

        if (LineLineIntersection(p0, v0, p1, v1, out m0, out m1))
        {
            return p0 + m0 * v0;
        }
        else
        {
            return new Vector2(float.NaN, float.NaN);
        }
    }

    public static float GeomArea(IList<Vector2> polygon)
    {
        int numPoints = polygon.Count;
        if (numPoints < 3) return 0;

        float area = 0;

        for (int i = 0; i < numPoints; i++)
        {
            Vector2 current = polygon[i];
            Vector2 next = polygon[(i + 1) % numPoints];

            area += current.x * next.y;
            area -= current.y * next.x;
        }

        area = Mathf.Abs(area) / 2.0f;
        return area;
    }
}
