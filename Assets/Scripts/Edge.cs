using UnityEngine;

/// <summary>
/// Custom class for representing an edge constraint as vertex indices.
/// </summary>
public class Edge
{
    public int v1;
    public int v2;

    public int t1;
    public int t2;
    public int t1Edge;

    public Edge(int v1, int v2)
    {
        this.v1 = v1;
        this.v2 = v2;
        this.t1 = -1;
        this.t2 = -1;
    }

    public Edge(int v1, int v2, int triangle1, int triangle2, int edge1)
    {
        this.v1 = v1;
        this.v2 = v2;
        this.t1 = triangle1;
        this.t2 = triangle2;
        this.t1Edge = edge1;
    }

    public override bool Equals(object other)
    {
        if (other is Edge otherEdge)
        {
            return (v1 == otherEdge.v1 && v2 == otherEdge.v2) ||
                   (v1 == otherEdge.v2 && v2 == otherEdge.v1);
        }
        return false;
    }

    public override int GetHashCode()
    {
        return new { v1, v2 }.GetHashCode() + new { v2, v1 }.GetHashCode();
    }

    public static bool operator ==(Edge lhs, Edge rhs)
    {
        return lhs.Equals(rhs);
    }

    public static bool operator !=(Edge lhs, Edge rhs)
    {
        return !lhs.Equals(rhs);
    }


}
