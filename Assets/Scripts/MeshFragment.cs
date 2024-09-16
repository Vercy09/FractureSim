using System;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.Rendering;

/// <summary>
/// Custom class for storing the data of a single mesh fragment.
/// </summary>
public class MeshFragment
{

    // All vertices of the mesh fragment geometry
    private List<MeshVertex> vertices;

    // Vertices of the cut faces produced during the creation of this fragment
    private List<MeshVertex> cutFaceVertices;

    // Triangle indices of the non-cut faces
    private List<int> indices;

    // Triangle indices of the cut faces
    private List<int> cutFaceIndices;

    // Map of vertices in the source fragment to vertices in this fragment (sourceIndexMap[sourceInd] = newInd)
    // source.vertices[sourceInd] = this.vertices[newInd]
    private int[] sourceIndexMap;

    // List of fragment edges for trianglulation
    private List<Edge> edges;

    // AABB of the mesh fragment
    public Bounds bounds;


    public List<MeshVertex> GetVertices()
    {
        return vertices;
    }

    public List<MeshVertex> GetCutFaceVertices()
    {
        return cutFaceVertices;
    }

    public List<int> GetIndices()
    {
        return indices;
    }

    public List<int> GetCutFaceIndices()
    {
        return cutFaceIndices;
    }

    public List<Edge> GetEdges()
    {
        return edges;
    }

    public int VertexCount { get { return vertices.Count; } }

    public int CutFaceVertexCount { get { return cutFaceVertices.Count; } }

    public int TotalVertexCount { get { return vertices.Count + cutFaceVertices.Count; } }

    public int IndexCount { get { return indices.Count; } }

    public int CutFaceIndexCount { get { return cutFaceIndices.Count; } }

    public int TotalIndexCount { get { return indices.Count + cutFaceIndices.Count; } }

    public int[] IndexMap { get { return sourceIndexMap; } }


    /// <summary>
    /// Constructor for an empty fragment.
    /// Provide some initial vertex and index count to reduce list resizing.
    /// </summary>
    public MeshFragment(int vertexCount, int indexCount)
    {
        vertices = new List<MeshVertex>(vertexCount);
        cutFaceVertices = new List<MeshVertex>(vertexCount / 10);

        indices = new List<int>(indexCount);
        cutFaceIndices = new List<int>(indexCount / 10);

        sourceIndexMap = new int[vertexCount];
        edges = new List<Edge>();
    }

    public MeshFragment(Mesh mesh)
    {
        int vertexCount = mesh.vertices.Length;
        vertices = new List<MeshVertex>(mesh.vertexCount);
        cutFaceVertices = new List<MeshVertex>(mesh.vertexCount / 10);

        indices = new List<int>(mesh.GetTriangles(0));
        if(mesh.subMeshCount >= 2)
        {
            cutFaceIndices = new List<int>(mesh.GetTriangles(1));
        }
        else
        {
            cutFaceIndices = new List<int>(mesh.triangles.Length / 10);
        }
        

        edges = new List<Edge>();
        sourceIndexMap = new int[vertexCount];

        var meshVertices = mesh.vertices;
        var meshNormals = mesh.normals;
        var meshUV = mesh.uv;
        for (int i = 0; i < vertexCount; i++)
        {
            vertices.Add(new MeshVertex(meshVertices[i], meshNormals[i], meshUV[i]));
        }

        SetBounds();
    }

    /// <summary>
    /// Adds an existing vertex of the source fragment to this fragment.
    /// Performs index mapping.
    /// </summary>
    public void AddExistingVertex(MeshVertex vertex, int sourceIndex)
    {
        vertices.Add(vertex);
        sourceIndexMap[sourceIndex] = vertices.Count - 1;
    }


    /// <summary>
    /// Adds a new cut face vertex to this fragment.
    /// </summary>
    public void AddCutFaceVertex(MeshVertex vertex)
    {
        vertices.Add(vertex);
        cutFaceVertices.Add(vertex);
    }


    /// <summary>
    /// Adds a new triangle to this fragment based 
    /// on the indices relative to this fragment's vertices. 
    /// </summary>
    public void AddNewTriangle(int v1, int v2, int v3, bool isCutFace)
    {
        var subMesh = isCutFace ? cutFaceIndices : indices;
        subMesh.Add(v1);
        subMesh.Add(v2);
        subMesh.Add(v3);
    }


    /// <summary>
    /// Adds an existing triangle of the source fragment to this fragment.
    /// Do not use with triangles of cut faces;
    /// </summary>
    public void AddExistingTriangle(int v1, int v2, int v3)
    {
        indices.Add(sourceIndexMap[v1]);
        indices.Add(sourceIndexMap[v2]);
        indices.Add(sourceIndexMap[v3]);
    }

    /// <summary>
    /// Adds an existing cut face triangle of the source fragment to this fragment.
    /// Do not use with triangles of non cut faces;
    /// </summary>
    public void AddExistingCutFaceTriangle(int v1, int v2, int v3)
    {
        cutFaceIndices.Add(sourceIndexMap[v1]);
        cutFaceIndices.Add(sourceIndexMap[v2]);
        cutFaceIndices.Add(sourceIndexMap[v3]);
    }

    /// <summary>
    /// Adds an edge constraint to this fragment.
    /// </summary>
    public void AddEdge(Edge edge)
    {
        edges.Add(edge);
    }


    /// <summary>
    /// Removes dulpicate cut face vertices of the fragment.
    /// </summary>
    public void RemoveDuplicateCutFaceVerts()
    {
        List<MeshVertex> uniqueVertices = new List<MeshVertex>(CutFaceVertexCount);

        // Index map: cutFaceVert --> uniqueVert
        int[] indexMap = new int[CutFaceVertexCount];

        // Count of unique vertices found
        int uniqueCount = 0;

        for (int i = 0; i < CutFaceVertexCount; i++)
        {
            bool isDuplicate = false;
            for (int j = 0; j < uniqueVertices.Count; j++)
            {
                if (cutFaceVertices[i].position == uniqueVertices[j].position)
                {
                    indexMap[i] = j;
                    isDuplicate = true;
                    break;
                }
            }

            if (!isDuplicate)
            {
                uniqueVertices.Add(cutFaceVertices[i]);
                indexMap[i] = uniqueCount;
                uniqueCount++;
            }
        }

        // Map the edges to the indices of the new unique vertices
        for (int i = 0; i < edges.Count; i++)
        {
            var edge = edges[i];
            edge.v1 = indexMap[edge.v1];
            edge.v2 = indexMap[edge.v2];
        }

        uniqueVertices.TrimExcess();
        cutFaceVertices = new List<MeshVertex>(uniqueVertices);
    }

    public void SetBounds()
    {
        if (VertexCount < 1) return;

        Vector3 vPos = vertices[0].position;
        Vector3 min = new Vector3(vPos.x, vPos.y, vPos.z);
        Vector3 max = new Vector3(vPos.x, vPos.y, vPos.z);

        for (int i = 1; i < VertexCount; i++)
        {
            vPos = vertices[i].position;
            if (vPos.x < min.x) min.x = vPos.x;
            if (vPos.x > max.x) max.x = vPos.x;

            if (vPos.y < min.y) min.y = vPos.y;
            if (vPos.y > max.y) max.y = vPos.y;

            if (vPos.z < min.z) min.z = vPos.z;
            if (vPos.z > max.z) max.z = vPos.z;
        }

        bounds = new Bounds((max + min) / 2.0f, max - min);
    }

    public Mesh GetMesh()
    {
        Mesh mesh = new Mesh();
        
        var layout = new[]
        {
            new VertexAttributeDescriptor(VertexAttribute.Position, VertexAttributeFormat.Float32, 3),
            new VertexAttributeDescriptor(VertexAttribute.Normal, VertexAttributeFormat.Float32, 3),
            new VertexAttributeDescriptor(VertexAttribute.TexCoord0, VertexAttributeFormat.Float32, 2),
        };

        mesh.SetIndexBufferParams(TotalIndexCount, IndexFormat.UInt32);
        mesh.SetVertexBufferParams(TotalVertexCount, layout);

        mesh.SetVertexBufferData(vertices, 0, 0, VertexCount);

        mesh.SetVertexBufferData(cutFaceVertices, 0, VertexCount, CutFaceVertexCount);

        mesh.subMeshCount = 2;
        mesh.SetIndexBufferData(indices, 0, 0, IndexCount);
        mesh.SetSubMesh(0, new SubMeshDescriptor(0, IndexCount));

        mesh.SetIndexBufferData(cutFaceIndices, 0, IndexCount, CutFaceIndexCount);
        mesh.SetSubMesh(1, new SubMeshDescriptor(IndexCount, CutFaceIndexCount));

        mesh.RecalculateBounds();

        return mesh;
    }
}
