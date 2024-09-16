using System.Collections;
using System.Collections.Generic;
using UnityEngine;

/// <summary>
/// Custom class for representing a point in the triangulation process.
/// </summary>
public class TriangulationPoint
{
    public int index;
    public Vector2 coordinates;
    public int bin;

    public TriangulationPoint(int index, Vector2 coordinates)
    {
        this.index = index;
        this.coordinates = coordinates;
        this.bin = -1;
    }
}
