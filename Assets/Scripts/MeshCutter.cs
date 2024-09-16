using System.Threading.Tasks;
using UnityEngine;

/// <summary>
/// Class which handles cutting a mesh(defined by a MeshFragment) into two slices(MeshFragment).
/// </summary>
public class MeshCutter
{
    public static void Cut(MeshFragment sourceFragment,
                            Vector3 cutPlaneNormal,
                            Vector3 cutPlaneOrigin,
                            out MeshFragment topSlice,
                            out MeshFragment bottomSlice)
    {
        topSlice = new MeshFragment(sourceFragment.TotalVertexCount, sourceFragment.TotalIndexCount);
        bottomSlice = new MeshFragment(sourceFragment.TotalVertexCount, sourceFragment.TotalIndexCount);

        bool[] isAboveCut = new bool[sourceFragment.TotalVertexCount];

        // Assign the vertices of the sourceFragment to appropriate slices
        var vertices = sourceFragment.GetVertices();
        int vertexCount = sourceFragment.VertexCount;
        for (int i = 0; i < vertexCount; i++)
        {
            var vertex = vertices[i];
            isAboveCut[i] = MathUtils.IsPointAbovePlane(vertex.position, cutPlaneNormal, cutPlaneOrigin);
            if (isAboveCut[i]) topSlice.AddExistingVertex(vertex, i);
            else bottomSlice.AddExistingVertex(vertex, i);
        }

        // Do the same for any vertices that got produced in the previous cut
        var cutFaceVertices = sourceFragment.GetCutFaceVertices();
        int cutFaceVertexCount = sourceFragment.CutFaceVertexCount;
        for (int i = 0; i < cutFaceVertexCount; i++)
        {
            var vertex = cutFaceVertices[i];
            isAboveCut[i + vertexCount] = MathUtils.IsPointAbovePlane(vertex.position, cutPlaneNormal, cutPlaneOrigin);
            if (isAboveCut[i + vertexCount]) topSlice.AddExistingVertex(vertex, i + vertexCount);
            else bottomSlice.AddExistingVertex(vertex, i + vertexCount);
        }


        int[] indices = sourceFragment.GetIndices().ToArray();

        // Assign the triangles of the sourceFragment to appropriate slices
        int v1, v2, v3;
        for (int i = 0; i < sourceFragment.IndexCount; i += 3)
        {
            v1 = indices[i];
            v2 = indices[i + 1];
            v3 = indices[i + 2];

            if (isAboveCut[v1] && isAboveCut[v2] && isAboveCut[v3])
            {
                topSlice.AddExistingTriangle(v1, v2, v3);
            }
            else if (!isAboveCut[v1] && !isAboveCut[v2] && !isAboveCut[v3])
            {
                bottomSlice.AddExistingTriangle(v1, v2, v3);
            }
            // Triangle interesects the cutting plane, need to cut the triangle
            else
            {
                int inV1, inV2, inV3;
                bool isFlipped;
                // Two vertices below, one vertex above the cut plane
                if (!isAboveCut[v1] && !isAboveCut[v2] && isAboveCut[v3])
                {
                    inV1 = v1;
                    inV2 = v2;
                    inV3 = v3;
                    isFlipped = false;
                }
                else if (!isAboveCut[v2] && !isAboveCut[v3] && isAboveCut[v1])
                {
                    inV1 = v2;
                    inV2 = v3;
                    inV3 = v1;
                    isFlipped = false;
                }
                else if (!isAboveCut[v3] && !isAboveCut[v1] && isAboveCut[v2])
                {
                    inV1 = v3;
                    inV2 = v1;
                    inV3 = v2;
                    isFlipped = false;
                }


                // Two vertices above, one vertex below the cut plane
                else if (isAboveCut[v1] && isAboveCut[v2] && !isAboveCut[v3])
                {
                    inV1 = v1;
                    inV2 = v2;
                    inV3 = v3;
                    isFlipped = true;
                }
                else if (isAboveCut[v2] && isAboveCut[v3] && !isAboveCut[v1])
                {
                    inV1 = v2;
                    inV2 = v3;
                    inV3 = v1;
                    isFlipped = true;
                }
                else if (isAboveCut[v3] && isAboveCut[v1] && !isAboveCut[v2])
                {
                    inV1 = v3;
                    inV2 = v1;
                    inV3 = v2;
                    isFlipped = true;
                }


                else continue;

                CutTriangle(inV1, inV2, inV3, cutPlaneOrigin, cutPlaneNormal, sourceFragment, topSlice, bottomSlice, isFlipped);
            }
        }

        // Do the same for the cut face(s)
        int[] cutFaceIndices = sourceFragment.GetCutFaceIndices().ToArray();

        for (int i = 0; i < sourceFragment.CutFaceIndexCount; i += 3)
        {
            v1 = cutFaceIndices[i];
            v2 = cutFaceIndices[i + 1];
            v3 = cutFaceIndices[i + 2];

            if (isAboveCut[v1] && isAboveCut[v2] && isAboveCut[v3])
            {
                topSlice.AddExistingCutFaceTriangle(v1, v2, v3);
            }
            else if (!isAboveCut[v1] && !isAboveCut[v2] && !isAboveCut[v3])
            {
                bottomSlice.AddExistingCutFaceTriangle(v1, v2, v3);
            }
            // Triangle interesects the cutting plane, need to cut the triangle
            else
            {
                int inV1, inV2, inV3;
                bool isFlipped;
                // Two vertices below, one vertex above the cut plane
                if (!isAboveCut[v1] && !isAboveCut[v2] && isAboveCut[v3])
                {
                    inV1 = v1;
                    inV2 = v2;
                    inV3 = v3;
                    isFlipped = false;
                }
                else if (!isAboveCut[v2] && !isAboveCut[v3] && isAboveCut[v1])
                {
                    inV1 = v2;
                    inV2 = v3;
                    inV3 = v1;
                    isFlipped = false;
                }
                else if (!isAboveCut[v3] && !isAboveCut[v1] && isAboveCut[v2])
                {
                    inV1 = v3;
                    inV2 = v1;
                    inV3 = v2;
                    isFlipped = false;
                }

                // Two vertices above, one vertex below the cut plane
                else if (isAboveCut[v1] && isAboveCut[v2] && !isAboveCut[v3])
                {
                    inV1 = v1;
                    inV2 = v2;
                    inV3 = v3;
                    isFlipped = true;
                }
                else if (isAboveCut[v2] && isAboveCut[v3] && !isAboveCut[v1])
                {
                    inV1 = v2;
                    inV2 = v3;
                    inV3 = v1;
                    isFlipped = true;
                }
                else if (isAboveCut[v3] && isAboveCut[v1] && !isAboveCut[v2])
                {
                    inV1 = v3;
                    inV2 = v1;
                    inV3 = v2;
                    isFlipped = true;
                }
                else continue;

                CutTriangle(inV1, inV2, inV3, cutPlaneOrigin, cutPlaneNormal, sourceFragment, topSlice, bottomSlice, isFlipped, true);
            }
        }

        CoverCutFaces(topSlice, bottomSlice, -cutPlaneNormal);
    }


    private static void CutTriangle(int v1,
                                    int v2,
                                    int v3,
                                    Vector3 cutPlaneOrigin,
                                    Vector3 cutPlaneNormal,
                                    MeshFragment sourceFragment,
                                    MeshFragment topSlice,
                                    MeshFragment bottomSlice,
                                    bool isFlipped,
                                    bool isCutFace = false)
    {

        var vertices = sourceFragment.GetVertices();
        var cutFaceVertices = sourceFragment.GetCutFaceVertices();

        MeshVertex vert1;
        if (v1 >= sourceFragment.VertexCount)
            vert1 = cutFaceVertices[v1 - sourceFragment.VertexCount];
        else
            vert1 = vertices[v1];

        MeshVertex vert2;
        if (v2 >= sourceFragment.VertexCount)
            vert2 = cutFaceVertices[v2 - sourceFragment.VertexCount];
        else
            vert2 = vertices[v2];

        MeshVertex vert3;
        if (v3 >= sourceFragment.VertexCount)
            vert3 = cutFaceVertices[v3 - sourceFragment.VertexCount];
        else
            vert3 = vertices[v3];

        bool isEdge13intersected = MathUtils.LineSegPlaneIntersection(vert1.position, vert3.position, cutPlaneNormal, cutPlaneOrigin, out Vector3 p13, out float s13);
        bool isEdge23intersected = MathUtils.LineSegPlaneIntersection(vert2.position, vert3.position, cutPlaneNormal, cutPlaneOrigin, out Vector3 p23, out float s23);

        if (isEdge13intersected && isEdge23intersected)
        {
            // Calulate normals and uv coords for new points
            Vector3 normal13 = MathUtils.InterpolateNormals(vert1.normal, vert3.normal, s13);
            Vector3 normal23 = MathUtils.InterpolateNormals(vert2.normal, vert3.normal, s23);
            Vector2 uv13 = MathUtils.InterpolateUV(vert1.uv, vert3.uv, s13);
            Vector2 uv23 = MathUtils.InterpolateUV(vert2.uv, vert3.uv, s23);

            // The cut produces same vertices on both slices, so add to both
            MeshVertex vert13 = new MeshVertex(p13, normal13, uv13);
            MeshVertex vert23 = new MeshVertex(p23, normal23, uv23);
            topSlice.AddCutFaceVertex(vert13);
            topSlice.AddCutFaceVertex(vert23);
            bottomSlice.AddCutFaceVertex(vert13);
            bottomSlice.AddCutFaceVertex(vert23);

            // Indices for the new triangles to be added
            int v13_topSlice = topSlice.VertexCount - 2;
            int v23_topSlice = topSlice.VertexCount - 1;
            int v13_bottomSlice = bottomSlice.VertexCount - 2;
            int v23_bottomSlice = bottomSlice.VertexCount - 1;

            if (isFlipped)
            {
                // Quad above the cut, construct into two new triangles
                topSlice.AddNewTriangle(v23_topSlice, v13_topSlice, topSlice.IndexMap[v2], isCutFace);
                topSlice.AddNewTriangle(v13_topSlice, topSlice.IndexMap[v1], topSlice.IndexMap[v2], isCutFace);

                // Triangle below cut
                bottomSlice.AddNewTriangle(bottomSlice.IndexMap[v3], v13_bottomSlice, v23_bottomSlice, isCutFace);

                // Add edge constraints (CCW)
                topSlice.AddEdge(new Edge(topSlice.CutFaceVertexCount - 2, topSlice.CutFaceVertexCount - 1));
                bottomSlice.AddEdge(new Edge(bottomSlice.CutFaceVertexCount - 1, bottomSlice.CutFaceVertexCount - 2));
            }
            else
            {
                // Triangle above cut
                topSlice.AddNewTriangle(v13_topSlice, v23_topSlice, topSlice.IndexMap[v3], isCutFace);

                // Quad below the cut, construct into two new triangles
                bottomSlice.AddNewTriangle(bottomSlice.IndexMap[v1], bottomSlice.IndexMap[v2], v13_bottomSlice, isCutFace);
                bottomSlice.AddNewTriangle(bottomSlice.IndexMap[v2], v23_bottomSlice, v13_bottomSlice, isCutFace);

                // Add edge constraints (CCW)
                topSlice.AddEdge(new Edge(topSlice.CutFaceVertexCount - 1, topSlice.CutFaceVertexCount - 2));
                bottomSlice.AddEdge(new Edge(bottomSlice.CutFaceVertexCount - 2, bottomSlice.CutFaceVertexCount - 1));
            }

        }

    }

    private static void CoverCutFaces(MeshFragment topSlice, MeshFragment bottomSlice, Vector3 sliceNormal)
    {
        topSlice.RemoveDuplicateCutFaceVerts();

        if (topSlice.CutFaceVertexCount < 3) return;

        ConstrainedTriangulator triangulator = new ConstrainedTriangulator(topSlice.GetCutFaceVertices(), topSlice.GetEdges(), sliceNormal);

        int[] triangles = triangulator.ConstrainedTriangulate();

        var topCutFaceVertices = topSlice.GetCutFaceVertices();
        var bottomCutFaceVertices = bottomSlice.GetCutFaceVertices();
        var triPoints = triangulator.GetTriangulationPoints();

        for (int i = 0; i < topSlice.CutFaceVertexCount; i++)
        {
            MeshVertex vertex = topCutFaceVertices[i];
            TriangulationPoint point = triPoints[i];

            Vector2 uv = new Vector2(triangulator.NormalizationFactor * point.coordinates.x, triangulator.NormalizationFactor * point.coordinates.y);

            MeshVertex topVertex = vertex;
            topVertex.normal = sliceNormal;
            topVertex.uv = uv;
            topCutFaceVertices[i] = topVertex;

            MeshVertex bottomVertex = vertex;
            bottomVertex.normal = -sliceNormal;
            bottomVertex.uv = uv;            
            bottomCutFaceVertices[i] = bottomVertex;
        }

        int offsetTop = topSlice.VertexCount;
        int offsetBottom = bottomSlice.VertexCount;
        for (int i = 0; i < triangles.Length; i += 3)
        {
            topSlice.AddNewTriangle(
                offsetTop + triangles[i],
                offsetTop + triangles[i + 1],
                offsetTop + triangles[i + 2],
                true);

            bottomSlice.AddNewTriangle(
                offsetBottom + triangles[i],
                offsetBottom + triangles[i + 2],
                offsetBottom + triangles[i + 1],
                true);
        }
    }


}
