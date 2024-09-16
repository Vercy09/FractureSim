using System.Collections.Generic;
using UnityEngine;

public class ConstrainedTriangulator
{

    // List of edge constraints for the triangulation
    private List<Edge> constraints;

    // Map of vertices to the triangles that contain them
    // vertexTriangleMap[vertex] = triangle
    private int[] vertexTriangleMap;

    private bool[] visited;



    /**********************triangulator*******/
    /*================= Constants for indexing into the triangle and adjacency table ==================*/
    protected const int V1 = 0; // Vertex 1
    protected const int V2 = 1; // Vertex 2
    protected const int V3 = 2; // Vertex 3
    protected const int tE12 = 3; // Triangle adjacent to edge V1V2
    protected const int tE23 = 4; // Triangle adjacent to edge V2V3
    protected const int tE31 = 5; // Triangle adjacent to edge V3V1
    protected const int BOUNDARY_EDGE = -1; // Index for adjancency when the triangle is on the edge
    /*================================================================================================*/

    // An array of points in the triangulation process
    protected TriangulationPoint[] triPoints;

    // 2D array for storing triangle vertices, and triangle adjacencies
    // Row is the triangle index, Columns are V1,V2,V3 and tE12,tE23,tE31 
    protected int[,] triMatrix;

    protected int pointCount;
    protected int triangleCount;

    protected bool[] discardTriangle;

    protected float normalizationFactor = 1.0f;

    public TriangulationPoint[] GetTriangulationPoints()
    {
        return triPoints;
    }

    public float NormalizationFactor { get { return normalizationFactor; } }


    public int[] Triangulate(bool discardTriangles)
    {
        if (pointCount < 3) return new int[] { };

        // Add the points of the super triangle
        triPoints[pointCount] = new TriangulationPoint(pointCount, new Vector2(-100f, -100f));
        triPoints[pointCount + 1] = new TriangulationPoint(pointCount + 1, new Vector2(0f, 100f));
        triPoints[pointCount + 2] = new TriangulationPoint(pointCount + 2, new Vector2(100f, -100f));

        // Store the super triangle in the tri matrix as the first member
        triMatrix[0, V1] = pointCount;
        triMatrix[0, V2] = pointCount + 1;
        triMatrix[0, V3] = pointCount + 2;
        triMatrix[0, tE12] = BOUNDARY_EDGE;
        triMatrix[0, tE23] = BOUNDARY_EDGE;
        triMatrix[0, tE31] = BOUNDARY_EDGE;

        // Sort the points into bins based on their position
        // to help when searching for the triangle to insert to
        var sortedPoints = BinSortPoints();
        if (sortedPoints == null) Debug.LogWarning("Sorted points are null!");

        int tSearch = 0; // The current triangle being searched
        int tLast = 0; // The last triangle formed

        // Loop through each point and add it to the triangulation
        for (int i = 0; i < pointCount; i++)
        {
            TriangulationPoint point = sortedPoints[i];
            int numSearched = 0;

            while (numSearched <= tLast && tSearch != BOUNDARY_EDGE)
            {
                Vector2 p = point.coordinates;
                Vector2 p1 = triPoints[triMatrix[tSearch, V1]].coordinates;
                Vector2 p2 = triPoints[triMatrix[tSearch, V2]].coordinates;
                Vector2 p3 = triPoints[triMatrix[tSearch, V3]].coordinates;

                if (!MathUtils.IsPointOnRightSideOfLineSeg(p, p1, p2))
                {
                    tSearch = triMatrix[tSearch, tE12];
                }
                else if (!MathUtils.IsPointOnRightSideOfLineSeg(p, p2, p3))
                {
                    tSearch = triMatrix[tSearch, tE23];
                }
                else if (!MathUtils.IsPointOnRightSideOfLineSeg(p, p3, p1))
                {
                    tSearch = triMatrix[tSearch, tE31];
                }
                else
                {
                    InsertPointIntoTriangle(point, tSearch, tLast);
                    tLast += 2;
                    tSearch = tLast;
                    break;
                }
                numSearched++;
            }
        }

        if (discardTriangles)
        {
            DiscardTrianglesWithSuperTriangleVertices();
        }

        List<int> res = new List<int>(3 * triangleCount);
        for (int i = 0; i < triangleCount; i++)
        {
            // Add all triangles that don't contain a super-triangle vertex
            if (!discardTriangle[i])
            {
                res.Add(triMatrix[i, V1]);
                res.Add(triMatrix[i, V2]);
                res.Add(triMatrix[i, V3]);
            }
        }
        return res.ToArray();
    }

    protected TriangulationPoint[] BinSortPoints()
    {
        float xMin = triPoints[0].coordinates.x;
        float xMax = xMin;

        float yMin = triPoints[0].coordinates.y;
        float yMax = yMin;

        for (int i = 1; i < pointCount; i++)
        {
            TriangulationPoint point = triPoints[i];

            if (point.coordinates.x < xMin) xMin = point.coordinates.x;
            if (point.coordinates.x > xMax) xMax = point.coordinates.x;

            if (point.coordinates.y < yMin) yMin = point.coordinates.y;
            if (point.coordinates.y > yMax) yMax = point.coordinates.y;
        }

        normalizationFactor = Mathf.Max(xMax - xMin, yMax - yMin);

        for (int i = 0; i < pointCount; i++)
        {
            TriangulationPoint point = triPoints[i];
            triPoints[i].coordinates = new Vector2((point.coordinates.x - xMin) / normalizationFactor, (point.coordinates.y - yMin) / normalizationFactor);
        }


        // Compute the number of bins along each axis
        int n = Mathf.RoundToInt(Mathf.Pow((float)pointCount, 0.25f));

        // Total bin count
        int binCount = n * n;

        // Assign the bin numbers to each point
        for (int k = 0; k < pointCount; k++)
        {
            TriangulationPoint point = triPoints[k];
            int i = (int)(0.99f * n * point.coordinates.y);
            int j = (int)(0.99f * n * point.coordinates.x);
            point.bin = (i % 2 == 0) ? (i * n) + j : (i + 1) * n - j - 1;

        }

        int[] count = new int[binCount];
        TriangulationPoint[] res = new TriangulationPoint[triPoints.Length];

        if (binCount <= 1)
        {
            return triPoints;
        }

        // Ignore the super triangle when sorting
        for (int i = 0; i < pointCount; i++)
        {
            int j = triPoints[i].bin;
            count[j] += 1;
        }

        for (int i = 1; i < binCount; i++)
        {
            count[i] += count[i - 1];
        }

        for (int i = pointCount - 1; i >= 0; i--)
        {
            int j = triPoints[i].bin;
            count[j] -= 1;
            res[count[j]] = triPoints[i];
        }

        // Add the super triangle points
        for (int i = pointCount; i < res.Length; i++)
        {
            res[i] = triPoints[i];
        }

        return res;
    }

    protected void InsertPointIntoTriangle(TriangulationPoint point, int t, int currTriangleCount)
    {
        // When adding a new point an existing triangle is replaced by 3 new triangles.
        // Replace the existing triangle by one of the new ones, and add the rest to the back of the triMatrix.
        int t1 = t;
        int t2 = currTriangleCount + 1;
        int t3 = currTriangleCount + 2;

        // Add the vertices and adjacency information for the two new triangles
        triMatrix[t2, V1] = point.index;
        triMatrix[t2, V2] = triMatrix[t, V2];
        triMatrix[t2, V3] = triMatrix[t, V3];
        triMatrix[t2, tE12] = t3;
        triMatrix[t2, tE23] = triMatrix[t, tE23];
        triMatrix[t2, tE31] = t1;

        triMatrix[t3, V1] = point.index;
        triMatrix[t3, V2] = triMatrix[t, V1];
        triMatrix[t3, V3] = triMatrix[t, V2];
        triMatrix[t3, tE12] = t1;
        triMatrix[t3, tE23] = triMatrix[t, tE12];
        triMatrix[t3, tE31] = t2;

        // Update the adjacencies so that they reference the newly created triangles
        // No need to update for edge13 because the new triangle just replaced the existing one.
        UpdateAdjacency(triMatrix[t, tE12], t, t3);
        UpdateAdjacency(triMatrix[t, tE23], t, t2);

        // Replace the information of the existing triangle with the new one
        triMatrix[t1, V2] = triMatrix[t, V3];
        triMatrix[t1, V3] = triMatrix[t, V1];
        triMatrix[t1, V1] = point.index; // different order because of overwriting data
        triMatrix[t1, tE23] = triMatrix[t, tE31];
        triMatrix[t1, tE12] = t2;
        triMatrix[t1, tE31] = t3;

        // Check if the new triangles satisfy the Dealunay conditions
        RestoreDelauney(point, t1, t2, t3);
    }

    protected void UpdateAdjacency(int tUpdate, int tOld, int tNew)
    {
        // If the triangle to update is on the boundary return.
        if (tUpdate == BOUNDARY_EDGE)
        {
            return;
        }
        // Else we need to find which edge did the triangles share so that we need which one to update.
        else if (triMatrix[tUpdate, tE12] == tOld)
        {
            triMatrix[tUpdate, tE12] = tNew;
        }
        else if (triMatrix[tUpdate, tE23] == tOld)
        {
            triMatrix[tUpdate, tE23] = tNew;
        }
        else if (triMatrix[tUpdate, tE31] == tOld)
        {
            triMatrix[tUpdate, tE31] = tNew;
        }
    }

    protected void RestoreDelauney(TriangulationPoint p, int t1, int t2, int t3)
    {
        // Stack of quadrilaterals represented as triangle pairs
        Stack<(int, int)> stack = new Stack<(int, int)>();

        // We need to check all the neighboring triangles of the triangle
        // where the point was inserted.
        // Since the new point has been added as vertex 1 all 
        // the oposing edges are on the tE23 index.
        stack.Push((t1, triMatrix[t1, tE23]));
        stack.Push((t2, triMatrix[t2, tE23]));
        stack.Push((t3, triMatrix[t3, tE23]));

        while (stack.Count > 0)
        {
            // Pop a quadrilateral off the stack.
            // Reuse variable t1 as the triangle containing the point and 
            // and t2 as the neighboring triangle.
            (t1, t2) = stack.Pop();

            // If t2 is the boundary edge ignore the quadrilateral
            if (t2 == BOUNDARY_EDGE)
            {
                continue;
            }
            // Else check if the quadrilateral diagonal needs to be flipped and flip it if neccessary
            else if (IsQuadDiagonalFlipped(p.index, t1, t2, out t3, out int t4))
            {
                // If the flip was performed add the new triangles so that we can 
                // check if they violate the Dealunay condition.
                stack.Push((t1, t3));
                stack.Push((t2, t4));
            }
        }
    }

    protected bool IsQuadDiagonalFlipped(int p, int t1, int t2, out int t3, out int t4)
    {
        // We need to find the points of the quadrilateral, one of the points is p and the rest can 
        // be found in t2.
        // t2's points may be defined in any order so we must check which edge t2 shares with t1 to
        // correctly identify the vertices and t2's neighboring triangles.

        int v1, v2, v3;
        if (triMatrix[t2, tE12] == t1)
        {
            v1 = triMatrix[t2, V2];
            v2 = triMatrix[t2, V1];
            v3 = triMatrix[t2, V3];

            t3 = triMatrix[t2, tE23];
            t4 = triMatrix[t2, tE31];
        }
        else if (triMatrix[t2, tE23] == t1)
        {
            v1 = triMatrix[t2, V3];
            v2 = triMatrix[t2, V2];
            v3 = triMatrix[t2, V1];

            t3 = triMatrix[t2, tE31];
            t4 = triMatrix[t2, tE12];
        }
        else if (triMatrix[t2, tE31] == t1)
        {
            v1 = triMatrix[t2, V1];
            v2 = triMatrix[t2, V3];
            v3 = triMatrix[t2, V2];

            t3 = triMatrix[t2, tE12];
            t4 = triMatrix[t2, tE23];
        }
        else
        {
            Debug.Log("Error in IsQuadDiagonalFlipped() function");
            t3 = -1; t4 = -1;
            return false;
        }

        // Check the Sloan's conditions to see if the circle around triangle v1,v2,v3 contains p
        if (IsSwapped(triPoints[v1].coordinates, triPoints[v2].coordinates, triPoints[v3].coordinates, triPoints[p].coordinates))
        {
            // After the swap t3 and tE31 will switch neighbours from t2 to t1 and vice versa
            UpdateAdjacency(t3, t2, t1);
            UpdateAdjacency(triMatrix[t1, tE31], t1, t2);

            // Construct the new t1 and t2 by re-assigning vertices
            triMatrix[t1, V1] = p;
            triMatrix[t1, V2] = v1;
            triMatrix[t1, V3] = v3;

            triMatrix[t2, V1] = p;
            triMatrix[t2, V2] = v3;
            triMatrix[t2, V3] = v2;

            // Update the other adjancencies (ORDER IS IMPORTANT!!!)
            triMatrix[t2, tE12] = t1;
            triMatrix[t2, tE23] = t4;
            triMatrix[t2, tE31] = triMatrix[t1, tE31];

            triMatrix[t1, tE23] = t3;
            triMatrix[t1, tE31] = t2;
            // triMatrix[t1, tE12] - no change

            return true;
        }
        else
        {
            return false;
        }
    }


    protected bool IsSwapped(Vector2 v1, Vector2 v2, Vector2 v3, Vector2 v4)
    {
        // S.W. Sloan's swap conditions
        float x13 = v1.x - v3.x;
        float y13 = v1.y - v3.y;

        float x23 = v2.x - v3.x;
        float y23 = v2.y - v3.y;

        float x14 = v1.x - v4.x;
        float y14 = v1.y - v4.y;

        float x24 = v2.x - v4.x;
        float y24 = v2.y - v4.y;

        float cosA = x13 * x23 + y13 * y23;
        float cosB = x24 * x14 + y24 * y14;

        if (cosA >= 0 && cosB >= 0)
        {
            return false;
        }
        else if (cosA < 0 && cosB < 0)
        {
            return true;
        }
        else
        {
            float sinA = x13 * y23 - x23 * y13;
            float sinB = x24 * y14 - x14 * y24;
            float sinAB = sinA * cosB + sinB * cosA;
            return sinAB < 0;
        }
    }

    protected void DiscardTrianglesWithSuperTriangleVertices()
    {
        for (int i = 0; i < triangleCount; i++)
        {
            // Find all the triangles that contain a super triangle vertex and mark them to be discarded
            if (TriangleContainsVertex(i, pointCount) ||
                TriangleContainsVertex(i, pointCount + 1) ||
                TriangleContainsVertex(i, pointCount + 2))
            {
                discardTriangle[i] = true;
            }
        }
    }

    protected bool TriangleContainsVertex(int t, int v)
    {
        return triMatrix[t, V1] == v || triMatrix[t, V2] == v || triMatrix[t, V3] == v;
    }
    /**************************************************************************/



    public ConstrainedTriangulator(List<MeshVertex> vertices, List<Edge> edges, Vector3 normal)
    {
        if (vertices == null || vertices.Count < 3)
        {
            return;
        }

        pointCount = vertices.Count;
        triangleCount = 2 * pointCount + 1;
        triMatrix = new int[triangleCount, 6];
        discardTriangle = new bool[triangleCount];
        triPoints = new TriangulationPoint[pointCount + 3]; //3 extra points for super triangle

        // Find basis vectors of the triangulation plane
        Vector3 b1 = (vertices[0].position - vertices[1].position).normalized;
        Vector3 b2 = Vector3.Cross(b1, normal.normalized).normalized;

        // Project the vertices onto the triangulation plane
        for (int i = 0; i < pointCount; i++)
        {
            var coords = new Vector2(Vector3.Dot(vertices[i].position, b1), Vector3.Dot(vertices[i].position, b2));
            triPoints[i] = new TriangulationPoint(i, coords);
        }


        constraints = edges;
    }

    public int[] ConstrainedTriangulate()
    {
        Triangulate(false);

        if (constraints.Count > 0)
        {
            ApplyConstraints();
            DiscardTrianglesViolatingConstraints();
        }

        DiscardTrianglesWithSuperTriangleVertices();

        List<int> triangles = new List<int>(3 * triangleCount);
        for (int i = 0; i < triangleCount; i++)
        {
            if (!discardTriangle[i])
            {
                triangles.Add(triMatrix[i, V1]);
                triangles.Add(triMatrix[i, V2]);
                triangles.Add(triMatrix[i, V3]);
            }
        }

        return triangles.ToArray();
    }

    private void ApplyConstraints()
    {
        // Map each vertex to a triangle that contains it
        visited = new bool[triMatrix.GetLength(0)];
        vertexTriangleMap = new int[pointCount + 3];
        for (int i = 0; i < triMatrix.GetLength(0); i++)
        {
            vertexTriangleMap[triMatrix[i, V1]] = i;
            vertexTriangleMap[triMatrix[i, V2]] = i;
            vertexTriangleMap[triMatrix[i, V3]] = i;
        }

        // Loop through each edge constraint
        foreach (Edge constraint in constraints)
        {
            // Ignore degenerate edges
            if (constraint.v1 == constraint.v2) continue;

            RemoveIntersectingEdges(constraint);
        }
    }

    private void RemoveIntersectingEdges(Edge constraint)
    {
        Queue<Edge> intersectingEdges = FindIntersectingEdges(constraint);
        List<Edge> newEdges = new List<Edge>();
        Edge edge, newEdge;

        int numVisited = 0;

        while (intersectingEdges.Count > 0 && numVisited <= intersectingEdges.Count)
        {
            edge = intersectingEdges.Dequeue();

            if (QuadFromEdge(edge.t1, edge.t1Edge, out Quad quad))
            {
                // If the quad is convex we can swap the diagonal
                
                /*if (MathUtils.LineSegLineSegIntersection(triPoints[quad.q4].coordinates,
                        triPoints[quad.q3].coordinates,
                        triPoints[quad.q1].coordinates,
                        triPoints[quad.q2].coordinates))*/
                
                if(MathUtils.IsQuadConvex(triPoints[quad.q1].coordinates, 
                                          triPoints[quad.q2].coordinates,
                                          triPoints[quad.q3].coordinates,
                                          triPoints[quad.q4].coordinates))
                {
                    SwapDiagonal(quad);
                    UpdateEdges(intersectingEdges, quad);
                    UpdateEdges(newEdges, quad);
                    UpdateEdges(constraints, quad);

                    newEdge = new Edge(quad.q3, quad.q4, quad.t1, quad.t2, tE31);

                    // Check if the flipped diagonal intersects the constraint
                    if (MathUtils.LineSegLineSegIntersection(triPoints[constraint.v1].coordinates,
                            triPoints[constraint.v2].coordinates,
                            triPoints[quad.q3].coordinates,
                            triPoints[quad.q4].coordinates))
                    {
                        intersectingEdges.Enqueue(newEdge);
                    }
                    else
                    {
                        numVisited = 0;
                        newEdges.Add(newEdge);
                    }
                }
                else
                {
                    intersectingEdges.Enqueue(edge);
                }
            }

            numVisited++;
        }

        RestoreDelauneyTriangulation(constraint, newEdges);
    }

    private Queue<Edge> FindIntersectingEdges(Edge constraint)
    {
        if (!FindStartEdge(constraint, out Edge startEdge))
            return new Queue<Edge>();

        Queue<Edge> res = new Queue<Edge>();
        res.Enqueue(startEdge);

        int tPrevious = startEdge.t1;
        int tSearch = triMatrix[tPrevious, startEdge.t1Edge];
        while (true)
        {

            Vector2 c1 = triPoints[constraint.v1].coordinates;
            Vector2 c2 = triPoints[constraint.v2].coordinates;
            Vector2 v1 = triPoints[triMatrix[tSearch, V1]].coordinates;
            Vector2 v2 = triPoints[triMatrix[tSearch, V2]].coordinates;
            Vector2 v3 = triPoints[triMatrix[tSearch, V3]].coordinates;

            if (TriangleContainsVertex(tSearch, constraint.v2))
            {
                break;
            }
            else if (triMatrix[tSearch, tE12] != tPrevious && MathUtils.LineSegLineSegIntersection(c1, c2, v1, v2))
            {
                res.Enqueue(new Edge(triMatrix[tSearch, V1], triMatrix[tSearch, V2], tSearch, triMatrix[tSearch, tE12], tE12));
                tPrevious = tSearch;
                tSearch = triMatrix[tSearch, tE12];
            }
            else if (triMatrix[tSearch, tE23] != tPrevious && MathUtils.LineSegLineSegIntersection(c1, c2, v2, v3))
            {
                res.Enqueue(new Edge(triMatrix[tSearch, V2], triMatrix[tSearch, V3], tSearch, triMatrix[tSearch, tE23], tE23));
                tPrevious = tSearch;
                tSearch = triMatrix[tSearch, tE23];
            }
            else if (triMatrix[tSearch, tE31] != tPrevious && MathUtils.LineSegLineSegIntersection(c1, c2, v3, v1))
            {
                res.Enqueue(new Edge(triMatrix[tSearch, V3], triMatrix[tSearch, V1], tSearch, triMatrix[tSearch, tE31], tE31));
                tPrevious = tSearch;
                tSearch = triMatrix[tSearch, tE31];
            }
            else
            {
                break;
            }
        }

        return res;
    }

    private bool FindStartEdge(Edge constraint, out Edge startEdge)
    {
        //bool[] visited = new bool[triangleCount];

        for (int i = 0; i < visited.Length; i++)
        {
            visited[i] = false;
        }

        int tSearch = vertexTriangleMap[constraint.v1];

        while (true)
        {
            visited[tSearch] = true;

            // Check if the edge constraint is the triangle's side
            if ((triMatrix[tSearch, V1] == constraint.v1 || triMatrix[tSearch, V2] == constraint.v1 || triMatrix[tSearch, V3] == constraint.v1) &&
               (triMatrix[tSearch, V1] == constraint.v2 || triMatrix[tSearch, V2] == constraint.v2 || triMatrix[tSearch, V3] == constraint.v2))
            {
                startEdge = new Edge(-1, -1);
                return false;
            }
            // Check if any of the triangle's sides intersect the constraint edge
            else
            {
                Vector2 c1 = triPoints[constraint.v1].coordinates;
                Vector2 c2 = triPoints[constraint.v2].coordinates;
                Vector2 v1 = triPoints[triMatrix[tSearch, V1]].coordinates;
                Vector2 v2 = triPoints[triMatrix[tSearch, V2]].coordinates;
                Vector2 v3 = triPoints[triMatrix[tSearch, V3]].coordinates;

                if (MathUtils.LineSegLineSegIntersection(c1, c2, v1, v2))
                {
                    int e1 = triMatrix[tSearch, V1];
                    int e2 = triMatrix[tSearch, V2];
                    int triangle2 = triMatrix[tSearch, tE12];

                    startEdge = new Edge(e1, e2, tSearch, triangle2, tE12);
                    return true;
                }
                else if (MathUtils.LineSegLineSegIntersection(c1, c2, v2, v3))
                {
                    int e1 = triMatrix[tSearch, V2];
                    int e2 = triMatrix[tSearch, V3];
                    int triangle2 = triMatrix[tSearch, tE23];

                    startEdge = new Edge(e1, e2, tSearch, triangle2, tE23);
                    return true;
                }
                else if (MathUtils.LineSegLineSegIntersection(c1, c2, v3, v1))
                {
                    int e1 = triMatrix[tSearch, V3];
                    int e2 = triMatrix[tSearch, V1];
                    int triangle2 = triMatrix[tSearch, tE31];

                    startEdge = new Edge(e1, e2, tSearch, triangle2, tE31);
                    return true;
                }
                // No intersection found for current triangle
                // Check the neighbouring triangles
                else
                {
                    
                    /*
                    int[] tNext = { triMatrix[tSearch, tE12], triMatrix[tSearch, tE23], triMatrix[tSearch, tE31] };
                    bool found = false;

                    foreach (int t in tNext)
                    {
                        if (!visited[t] && t != BOUNDARY_EDGE && TriangleContainsVertex(t, constraint.v1))
                        {
                            tSearch = t;
                            found = true;
                            break;
                        }
                    }
                    
                    // Did not find any appropriate neighbour
                    if (!found)
                    {
                        startEdge = new Edge(-1, -1);
                        return false;
                    }*/

                    
                    int triE12 = triMatrix[tSearch, tE12];
                    int triE23 = triMatrix[tSearch, tE23];
                    int triE31 = triMatrix[tSearch, tE31];
                    if (triE12 != BOUNDARY_EDGE && !visited[triE12] && TriangleContainsVertex(triE12, constraint.v1))
                    {
                        tSearch = triE12;
                    }
                    else if (triE23 != BOUNDARY_EDGE && !visited[triE23] && TriangleContainsVertex(triE23, constraint.v1))
                    {
                        tSearch = triE23;
                    }
                    else if (triE31 != BOUNDARY_EDGE && !visited[triE31] && TriangleContainsVertex(triE31, constraint.v1))
                    {
                        tSearch = triE31;
                    }
                    else
                    {
                        startEdge = new Edge(-1, -1);
                        return false;
                    }

                }

            }
        }
    }

    private bool QuadFromEdge(int t1, int t1SharedEdge, out Quad quad)
    {

        int q1, q2, q3, q4;
        int t1L, t1R, t2L, t2R;

        int t2 = triMatrix[t1, t1SharedEdge];
        int t2SharedEdge;

        if (t2 == BOUNDARY_EDGE)
        {
            quad = new Quad();
            return false;
        }
        else if (triMatrix[t2, tE12] == t1)
        {
            t2SharedEdge = tE12;
        }
        else if (triMatrix[t2, tE23] == t1)
        {
            t2SharedEdge = tE23;
        }
        else if (triMatrix[t2, tE31] == t1)
        {
            t2SharedEdge = tE31;
        }
        else
        {
            quad = new Quad();
            return false;
        }

        if (t1SharedEdge == tE12)
        {
            q4 = triMatrix[t1, V3];
            t1L = triMatrix[t1, tE31];
            t1R = triMatrix[t1, tE23];
        }
        else if (t1SharedEdge == tE23)
        {
            q4 = triMatrix[t1, V1];
            t1L = triMatrix[t1, tE12];
            t1R = triMatrix[t1, tE31];
        }
        else
        {
            q4 = triMatrix[t1, V2];
            t1L = triMatrix[t1, tE23];
            t1R = triMatrix[t1, tE12];
        }

        if (t2SharedEdge == tE12)
        {
            q2 = triMatrix[t2, V1];
            q1 = triMatrix[t2, V2];
            q3 = triMatrix[t2, V3];
            t2L = triMatrix[t2, tE23];
            t2R = triMatrix[t2, tE31];
        }
        else if (t2SharedEdge == tE23)
        {
            q2 = triMatrix[t2, V2];
            q1 = triMatrix[t2, V3];
            q3 = triMatrix[t2, V1];
            t2L = triMatrix[t2, tE31];
            t2R = triMatrix[t2, tE12];
        }
        else
        {
            q2 = triMatrix[t2, V3];
            q1 = triMatrix[t2, V1];
            q3 = triMatrix[t2, V2];
            t2L = triMatrix[t2, tE12];
            t2R = triMatrix[t2, tE23];
        }

        quad = new Quad(q1, q2, q3, q4, t1, t2, t1L, t1R, t2L, t2R);
        return true;
    }

    private void SwapDiagonal(Quad quad)
    {
        triMatrix[quad.t1, V1] = quad.q4;
        triMatrix[quad.t1, V2] = quad.q1;
        triMatrix[quad.t1, V3] = quad.q3;
        triMatrix[quad.t1, tE12] = quad.t1L;
        triMatrix[quad.t1, tE23] = quad.t2L;
        triMatrix[quad.t1, tE31] = quad.t2;

        triMatrix[quad.t2, V1] = quad.q4;
        triMatrix[quad.t2, V2] = quad.q3;
        triMatrix[quad.t2, V3] = quad.q2;
        triMatrix[quad.t2, tE12] = quad.t1;
        triMatrix[quad.t2, tE23] = quad.t2R;
        triMatrix[quad.t2, tE31] = quad.t1R;

        UpdateAdjacency(quad.t2L, quad.t2, quad.t1);
        UpdateAdjacency(quad.t1R, quad.t1, quad.t2);

        vertexTriangleMap[quad.q1] = quad.t1;
        vertexTriangleMap[quad.q2] = quad.t2;
    }

    internal void UpdateEdges(IEnumerable<Edge> edges, Quad quad)
    {
        foreach (Edge edge in edges)
        {
            if (edge.t1 == quad.t1)
            {
                if (edge.t2 == quad.t1R)
                {
                    edge.t1 = quad.t2;
                    edge.t2 = quad.t1R;
                    edge.t1Edge = tE31;
                }
                else if (edge.t2 == quad.t1L)
                {
                    edge.t1Edge = tE12;
                }
            }
            else if (edge.t1 == quad.t2)
            {
                if (edge.t2 == quad.t2L)
                {
                    edge.t1 = quad.t1;
                    edge.t2 = quad.t2L;
                    edge.t1Edge = tE23;
                }
                else if (edge.t2 == quad.t2R)
                {
                    edge.t1Edge = tE23;
                }
            }
            else if (edge.t2 == quad.t1)
            {
                if (edge.t1 == quad.t1R)
                {
                    edge.t2 = quad.t2;
                }
                else if (edge.t1 == quad.t2L)
                {
                    edge.t2 = quad.t1;
                }
            }
        }
    }

    private void RestoreDelauneyTriangulation(Edge constraint, List<Edge> newEdges)
    {
        bool didSwap;
        do
        {
            didSwap = false;

            for (int i = 0; i < newEdges.Count; i++)
            {
                Edge edge = newEdges[i];

                if (edge == constraint)
                    continue;

                if (QuadFromEdge(edge.t1, edge.t1Edge, out Quad quad))
                {
                    if (IsSwapped(triPoints[quad.q1].coordinates, triPoints[quad.q2].coordinates, triPoints[quad.q3].coordinates, triPoints[quad.q4].coordinates))
                    {
                        SwapDiagonal(quad);
                        UpdateEdges(newEdges, quad);
                        UpdateEdges(constraints, quad);

                        newEdges[i] = new Edge(quad.q3, quad.q4, quad.t1, quad.t2, tE31);
                        didSwap = true;
                    }
                }
            }
        } while (didSwap);
    }

    private void DiscardTrianglesViolatingConstraints()
    {
        //bool[] visited = new bool[triangleCount];
        for (int i = 0; i < visited.Length; i++)
        {
            visited[i] = false;
        }

        for (int i = 0; i < triangleCount; i++)
        {
            discardTriangle[i] = true;
        }

        HashSet<(int, int)> boundaries = new HashSet<(int, int)>();
        foreach (Edge constraint in constraints)
        {
            boundaries.Add((constraint.v1, constraint.v2));
        }

        Queue<int> search = new Queue<int>();

        int v1, v2, v3;
        bool isE12Boundary, isE23Boundary, isE31Boundary;
        for (int i = 0; i < triangleCount; i++)
        {
            if (visited[i]) continue;
            //visited[i] = true;

            v1 = triMatrix[i, V1];
            v2 = triMatrix[i, V2];
            v3 = triMatrix[i, V3];
            isE12Boundary = boundaries.Contains((v1, v2));
            isE23Boundary = boundaries.Contains((v2, v3));
            isE31Boundary = boundaries.Contains((v3, v1));

            if (isE12Boundary || isE23Boundary || isE31Boundary)
            {
                discardTriangle[i] = false;
                //visited[i] = true;

                search.Clear();
                if (!isE12Boundary)
                {
                    search.Enqueue(triMatrix[i, tE12]);
                }
                if (!isE23Boundary)
                {
                    search.Enqueue(triMatrix[i, tE23]);
                }
                if (!isE31Boundary)
                {
                    search.Enqueue(triMatrix[i, tE31]);
                }

                while (search.Count > 0)
                {
                    int tKeep = search.Dequeue();

                    if (tKeep == BOUNDARY_EDGE || visited[tKeep]) continue;

                    visited[tKeep] = true;
                    discardTriangle[tKeep] = false;


                    v1 = triMatrix[tKeep, V1];
                    v2 = triMatrix[tKeep, V2];
                    v3 = triMatrix[tKeep, V3];

                    if (!boundaries.Contains((v1, v2)))
                    {
                        search.Enqueue(triMatrix[tKeep, tE12]);
                    }
                    if (!boundaries.Contains((v2, v3)))
                    {
                        search.Enqueue(triMatrix[tKeep, tE23]);
                    }
                    if (!boundaries.Contains((v3, v1)))
                    {
                        search.Enqueue(triMatrix[tKeep, tE31]);
                    }
                }
            }
        }
    }
}

public class Quad
{
    public int q1, q2, q3, q4;
    public int t1, t2;
    public int t1L, t1R, t2L, t2R;

    public Quad(int q1, int q2, int q3, int q4, int t1, int t2, int t1L, int t1R, int t2L, int t2R)
    {
        this.q1 = q1;
        this.q2 = q2;
        this.q3 = q3;
        this.q4 = q4;
        this.t1 = t1;
        this.t2 = t2;
        this.t1L = t1L;
        this.t1R = t1R;
        this.t2L = t2L;
        this.t2R = t2R;
    }

    public Quad() { }

}
