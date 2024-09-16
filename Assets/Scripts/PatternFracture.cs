
using System;
using System.Collections.Generic;
using System.IO;
using UnityEngine;
using UnityEngine.TestTools;
using Random = UnityEngine.Random;

[RequireComponent(typeof(MeshFilter))]
[RequireComponent(typeof(MeshRenderer))]
[RequireComponent(typeof(Rigidbody))]
[RequireComponent(typeof(MeshCollider))]
public class PatternFracture : MonoBehaviour
{
    private GameObject rootFragmentObj;
    public MeshFilter Filter { get; private set; }
    public MeshRenderer Renderer { get; private set; }
    public MeshCollider Collider { get; private set; }
    public Rigidbody Rigidbody { get; private set; }

    [SerializeField]
    public List<string> fractureCausingTags;

    [SerializeField]
    float forceToFracture = 0.0f;

    public List<Vector2> Polygon;

    [HideInInspector]
    public float thickness = 1.0f;
    public float minBreakArea = 1.0f;
    public float _Area = -1.0f;


    [SerializeField]
    private string fileName = "GlassVoronoiDiag";

    public float Area
    {
        get
        {
            if (_Area < 0.0f)
            {
                _Area = MathUtils.GeomArea(Polygon);
            }

            return _Area;
        }
    }

    void Start()
    {
        if (Filter == null) Filter = GetComponent<MeshFilter>();
        if (Renderer == null) Renderer = GetComponent<MeshRenderer>();
        if (Collider == null) Collider = GetComponent<MeshCollider>();
        if (Rigidbody == null) Rigidbody = GetComponent<Rigidbody>();

        if (Polygon.Count == 0)
        {
            Vector3 scale = 0.5f * transform.localScale;

            Polygon.Add(new Vector2(-scale.x, -scale.y));
            Polygon.Add(new Vector2(scale.x, -scale.y));
            Polygon.Add(new Vector2(scale.x, scale.y));
            Polygon.Add(new Vector2(-scale.x, scale.y));

            thickness = 2.0f * scale.z;
            transform.localScale = Vector3.one;
        }

        Mesh mesh = PolygonToMesh(Polygon, thickness);

        Filter.sharedMesh = mesh;
        Collider.sharedMesh = mesh;
    }

    void OnCollisionEnter(Collision collision)
    {
        if (collision.contactCount > 0)
        {
            float collisionForce = collision.impulse.magnitude / Time.fixedDeltaTime;
            if (collisionForce > forceToFracture && fractureCausingTags.Contains(collision.GetContact(0).otherCollider.gameObject.tag))
            {
                Vector3 impactPoint = collision.GetContact(0).point;
                Vector2 fracturePoint = (Vector2)transform.InverseTransformPoint(impactPoint);
                Fragment(fracturePoint);
            }
        }
    }

    public void Fragment(Vector2 fracturePoint)
    {
        float area = Area;
        if (area < 0.0f) Debug.Log("negative area");
        if (area > minBreakArea)
        {
            VoronoiFragmenter clip = new VoronoiFragmenter(fileName);
            List<Vector2> clipped = new List<Vector2>();

            for (int i = 0; i < clip.SiteCount; i++)
            {
                clip.ClipSite(Polygon, i, ref clipped, fracturePoint);

                if (clipped.Count > 0)
                {
                    GameObject newFragment = Instantiate(gameObject, transform.parent);
                    newFragment.transform.SetLocalPositionAndRotation(transform.localPosition, transform.localRotation);

                    if (!rootFragmentObj)
                    {
                        rootFragmentObj = new GameObject($"{name}'s fragments");
                        rootFragmentObj.transform.SetParent(transform.parent);
                        rootFragmentObj.transform.SetPositionAndRotation(transform.position, transform.rotation);
                        rootFragmentObj.transform.localScale = new Vector3(1.0f, 1.0f, 1.0f);
                    }
                    newFragment.transform.SetParent(rootFragmentObj.transform);

                    var fragmentPatternFracture = newFragment.GetComponent<PatternFracture>();
                    fragmentPatternFracture.thickness = thickness;
                    fragmentPatternFracture.Polygon.Clear();
                    fragmentPatternFracture.Polygon.AddRange(clipped);

                    float childArea = fragmentPatternFracture.Area;

                    var fragmentRigidBody = fragmentPatternFracture.GetComponent<Rigidbody>();
                    fragmentRigidBody.mass = Rigidbody.mass * (childArea / area);
                }
            }

            gameObject.SetActive(false);
            Destroy(gameObject);
        }
    }

    static Mesh PolygonToMesh(List<Vector2> polygon, float thickness)
    {
        int n = polygon.Count;

        Vector3[] vertices = new Vector3[6 * n];
        // #Triangles:
        // TOP: n-2
        // BOT: n-2
        // SIDE: 2
        // TOTAL: (n-2) + (n-2) + n * 2 = 4n - 4
        int[] indices = new int[3 * (4 * n - 4)];

        Vector3[] normals = new Vector3[6 * n];
        // @todo: add uv

        int vertInd = 0;
        int normalInd = 0;
        int indInd = 0;

        float extrusion = thickness / 2.0f;

        // TOP: n vertices
        for (int i = 0; i < n; i++)
        {
            vertices[vertInd++] = new Vector3(polygon[i].x, polygon[i].y, extrusion);
            normals[normalInd++] = Vector3.forward;
        }

        // BOT: n vertices
        for (int i = 0; i < n; i++)
        {
            vertices[vertInd++] = new Vector3(polygon[i].x, polygon[i].y, -extrusion);
            normals[normalInd++] = Vector3.back;
        }

        // SIDES: n sides, 4 vertices per side
        for (int i = 0; i < n; i++)
        {
            int iNext = (i + 1) % n;

            vertices[vertInd++] = new Vector3(polygon[i].x, polygon[i].y, extrusion);
            vertices[vertInd++] = new Vector3(polygon[i].x, polygon[i].y, -extrusion);
            vertices[vertInd++] = new Vector3(polygon[iNext].x, polygon[iNext].y, -extrusion);
            vertices[vertInd++] = new Vector3(polygon[iNext].x, polygon[iNext].y, extrusion);


            Vector3 normal = Vector3.Cross(polygon[iNext] - polygon[i], Vector3.forward).normalized;
            normals[normalInd++] = normal;
            normals[normalInd++] = normal;
            normals[normalInd++] = normal;
            normals[normalInd++] = normal;
        }


        for (int vert = 2; vert < n; vert++)
        {
            indices[indInd++] = 0;
            indices[indInd++] = vert - 1;
            indices[indInd++] = vert;
        }

        for (int vert = 2; vert < n; vert++)
        {
            indices[indInd++] = n;
            indices[indInd++] = n + vert;
            indices[indInd++] = n + vert - 1;
        }

        for (int vert = 0; vert < n; vert++)
        {
            int offset = 2 * n + 4 * vert;

            indices[indInd++] = offset;
            indices[indInd++] = offset + 1;
            indices[indInd++] = offset + 2;

            indices[indInd++] = offset;
            indices[indInd++] = offset + 2;
            indices[indInd++] = offset + 3;
        }

        Debug.Assert(indInd == indices.Length);
        Debug.Assert(vertInd == vertices.Length);

        Mesh mesh = new Mesh();
        mesh.vertices = vertices;
        mesh.triangles = indices;
        mesh.normals = normals;

        return mesh;
    }
}

