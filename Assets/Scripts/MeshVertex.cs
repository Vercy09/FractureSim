using UnityEngine;

/// <summary>
/// Custom struct for storing the position, normal and UV data of a single vertex.
/// </summary>
public struct MeshVertex
{
    public Vector3 position;
    public Vector3 normal;
    public Vector2 uv;

    public MeshVertex(Vector3 position) 
    {
        this.position = position;
        this.normal = Vector3.zero;
        this.uv = Vector3.zero;
    }

    public MeshVertex(Vector3 position, Vector3 normal, Vector2 uv) 
    {
        this.position = position;
        this.normal = normal;
        this.uv = uv;
    }
}