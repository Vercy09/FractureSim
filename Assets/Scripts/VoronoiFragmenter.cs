using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;


public class VoronoiFragmenter
{

    private VoronoiDiagram voronoiDiag;
    private List<Vector2> pointsIn = new List<Vector2>();
    private List<Vector2> pointsOut = new List<Vector2>();

    public float SiteCount
    {
        get
        {
            return voronoiDiag.Sites.Count;
        }
    }


    public VoronoiFragmenter(string diagramFileName)
    {
        TextAsset jsonFile = Resources.Load<TextAsset>(diagramFileName);

        if (jsonFile != null)
        {
            // Parse JSON content
            voronoiDiag = JsonUtility.FromJson<VoronoiDiagram>(jsonFile.text);
        }
        else
        {
            Debug.LogError("Failed to load JSON file.");
        }
    }

    public void ClipSite(IList<Vector2> polygon, int site, ref List<Vector2> clipped, Vector2 offset)
    {
        pointsIn.Clear();
        pointsIn.AddRange(polygon);

        int firstEdge, lastEdge;

        if (site == voronoiDiag.Sites.Count - 1)
        {
            firstEdge = voronoiDiag.FirstEdgeBySite[site];
            lastEdge = voronoiDiag.Edges.Count - 1;
        }
        else
        {
            firstEdge = voronoiDiag.FirstEdgeBySite[site];
            lastEdge = voronoiDiag.FirstEdgeBySite[site + 1] - 1;
        }

        for (int ei = firstEdge; ei <= lastEdge; ei++)
        {
            pointsOut.Clear();

            var edge = voronoiDiag.Edges[ei];

            Vector2 lp, ld;

            if (edge.Type == VoronoiDiagram.EdgeType.RayCCW || edge.Type == VoronoiDiagram.EdgeType.RayCW)
            {
                lp = voronoiDiag.Vertices[edge.Vert0];
                ld = edge.Direction;

                if (edge.Type == VoronoiDiagram.EdgeType.RayCW)
                {
                    ld *= -1;
                }
            }
            else if (edge.Type == VoronoiDiagram.EdgeType.Segment)
            {
                Vector2 lp0 = voronoiDiag.Vertices[edge.Vert0];
                Vector2 lp1 = voronoiDiag.Vertices[edge.Vert1];

                lp = lp0;
                ld = lp1 - lp0;
            }
            else
            {
                Debug.Assert(false);
                return;
            }
            lp = lp + offset;
            for (int pi0 = 0; pi0 < pointsIn.Count; pi0++)
            {
                var pi1 = pi0 == pointsIn.Count - 1 ? 0 : pi0 + 1;

                Vector2 p0 = pointsIn[pi0];
                Vector2 p1 = pointsIn[pi1];

                var p0Inside = MathUtils.IsPointOnLeftSideOfLineSeg(p0, lp, lp + ld);
                var p1Inside = MathUtils.IsPointOnLeftSideOfLineSeg(p1, lp, lp + ld);

                if (p0Inside && p1Inside)
                {
                    pointsOut.Add(p1);
                }
                else if (!p0Inside && !p1Inside)
                {
                    // Do nothing, both are outside
                }
                else
                {
                    Vector2 intersection = MathUtils.LineLineIntersection(lp, ld.normalized, p0, (p1 - p0).normalized);

                    if (p0Inside)
                    {
                        pointsOut.Add(intersection);
                    }
                    else if (p1Inside)
                    {
                        pointsOut.Add(intersection);
                        pointsOut.Add(p1);
                    }
                    else
                    {
                        Debug.Assert(false);
                    }
                }
            }

            List<Vector2> tmp = pointsIn;
            pointsIn = pointsOut;
            pointsOut = tmp;
        }

        if (clipped == null)
        {
            clipped = new List<Vector2>();
        }
        else
        {
            clipped.Clear();
        }

        clipped.AddRange(pointsIn);
    }
}


[System.Serializable]
public class VoronoiDiagram
{
    public List<Vector2> Sites;
    public List<Vector2> Vertices;
    public List<Edge> Edges;
    public List<int> FirstEdgeBySite;

   [System.Serializable]
    public enum EdgeType
    {
        Line,
        RayCCW,
        RayCW,
        Segment
    }

    [System.Serializable]
    public struct Edge
    {
        public EdgeType Type;
        public int Site;
        public int Vert0;
        public int Vert1;
        public Vector2 Direction;

        public Edge(EdgeType type, int site, int vert0, int vert1, Vector2 direction)
        {
            this.Type = type;
            this.Site = site;
            this.Vert0 = vert0;
            this.Vert1 = vert1;
            this.Direction = direction;
        }
    }
}
