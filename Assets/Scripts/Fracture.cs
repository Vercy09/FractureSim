using System;
using System.Collections.Generic;
using Unity.Collections;
using UnityEngine;
using UnityEngine.Rendering;
using Random = UnityEngine.Random;

[RequireComponent(typeof(MeshFilter))]
[RequireComponent(typeof(MeshRenderer))]
[RequireComponent(typeof(Rigidbody))]
[RequireComponent(typeof(Collider))]
public class Fracture : MonoBehaviour
{
    private GameObject rootFragmentObj;

    [SerializeField]
    public List<string> fractureCausingObjectTags;

    [SerializeField]
    protected float forceToFracture = 0.0f;

    [SerializeField]
    protected bool isLocalized = false;

    [SerializeField]
    protected bool separateIslands = false;

    [SerializeField]
    protected float fracturePointSampleOffset = 0.1f;

    [SerializeField]
    protected int fragmentCount = 10;

    [SerializeField]
    protected bool logVertexTriangleCount = false;

    [SerializeField]
    public int fractureDepth = 0;

    [HideInInspector]
    public int currentFractureDepth = 0;

    [SerializeField]
    private Material innerMaterial;

    protected void OnCollisionEnter(Collision collision)
    {
        if (collision.contactCount > 0)
        {
            float collisionForce = collision.impulse.magnitude / Time.fixedDeltaTime;
            if (collisionForce > forceToFracture && fractureCausingObjectTags.Contains(collision.GetContact(0).otherCollider.gameObject.tag))
            {
                var mesh = GetComponent<MeshFilter>().sharedMesh;
                if (logVertexTriangleCount)
                {
                    Debug.Log("Vertices before fracture: " + mesh.vertexCount);
                    Debug.Log("Triangles before fracture: " + mesh.triangles.Length / 3);
                }
                Vector3 impactPoint = collision.GetContact(0).point;
                impactPoint = transform.InverseTransformPoint(impactPoint);
                if (mesh)
                {
                    if (!rootFragmentObj)
                    {
                        rootFragmentObj = new GameObject($"{name}'s fragments");
                        rootFragmentObj.transform.SetParent(transform.parent);
                        rootFragmentObj.transform.SetPositionAndRotation(transform.position, transform.rotation);
                        rootFragmentObj.transform.localScale = new Vector3(1.0f, 1.0f, 1.0f);
                    }
                    GameObject newFragmentModel = new GameObject();
                    AddComponents(ref newFragmentModel);

                    Fragment(newFragmentModel, rootFragmentObj.transform, impactPoint);

                    Destroy(newFragmentModel);
                    gameObject.SetActive(false);
                }
            }
        }

    }

    private void AddComponents(ref GameObject newFragmentModel)
    {
        newFragmentModel.name = name + "'s fragment";
        newFragmentModel.tag = tag;

        // Mesh filter
        newFragmentModel.AddComponent<MeshFilter>();
        var meshRenderer = newFragmentModel.AddComponent<MeshRenderer>();
        meshRenderer.sharedMaterials = new Material[2] {
            GetComponent<MeshRenderer>().sharedMaterial,
            innerMaterial
        };

        // Collider
        var collider = GetComponent<Collider>();
        var fragmentCollider = newFragmentModel.AddComponent<MeshCollider>();
        fragmentCollider.convex = true;
        fragmentCollider.isTrigger = collider.isTrigger;
        fragmentCollider.providesContacts = collider.providesContacts;
        fragmentCollider.sharedMaterial = collider.sharedMaterial;


        // Rigid Body
        var rigidBody = GetComponent<Rigidbody>();
        var fragmentRigidBody = newFragmentModel.AddComponent<Rigidbody>();
        fragmentRigidBody.velocity = rigidBody.velocity;
        fragmentRigidBody.angularVelocity = rigidBody.angularVelocity;
        fragmentRigidBody.drag = rigidBody.drag;
        fragmentRigidBody.angularDrag = rigidBody.angularDrag;
        fragmentRigidBody.useGravity = rigidBody.useGravity;

        // Fracture
        if (currentFractureDepth < fractureDepth)
        {
            var fractureComponent = newFragmentModel.AddComponent<Fracture>();
            fractureComponent.rootFragmentObj = rootFragmentObj;
            fractureComponent.fractureCausingObjectTags = fractureCausingObjectTags;
            fractureComponent.forceToFracture = forceToFracture;
            fractureComponent.fragmentCount = fragmentCount;
            fractureComponent.fractureDepth = fractureDepth;
            fractureComponent.currentFractureDepth = currentFractureDepth + 1;
            fractureComponent.innerMaterial = innerMaterial;
            fractureComponent.isLocalized = isLocalized;
            fractureComponent.fracturePointSampleOffset = fracturePointSampleOffset;
            fractureComponent.separateIslands = separateIslands;
            //fractureComponent.logVertexTriangleCount = logVertexTriangleCount;
        }
    }

    public void Fragment(GameObject fragmentTemplate, Transform parent, Vector3 impactPoint)
    {
        var fragmentsToCut = new Queue<MeshFragment>();
        var initialFragment = new MeshFragment(gameObject.GetComponent<MeshFilter>().sharedMesh);
        fragmentsToCut.Enqueue(initialFragment);


        while (fragmentsToCut.Count < fragmentCount)
        {
            MeshFragment fragmentToCut = fragmentsToCut.Dequeue();
            fragmentToCut.SetBounds();


            Vector3 samplePointOffset = new Vector3(Random.Range(-fracturePointSampleOffset, fracturePointSampleOffset),
                                                    Random.Range(-fracturePointSampleOffset, fracturePointSampleOffset),
                                                    Random.Range(-fracturePointSampleOffset, fracturePointSampleOffset));

            Vector3 fracturePoint = isLocalized ? impactPoint : fragmentToCut.bounds.center;
            fracturePoint += samplePointOffset;


            Vector3 normal = new Vector3(Random.Range(-1.0f, 1.0f), Random.Range(-1.0f, 1.0f), Random.Range(-1.0f, 1.0f));
            MeshCutter.Cut(fragmentToCut,
                           normal,
                           fracturePoint,
                           out MeshFragment topSlice,
                           out MeshFragment bottomSlice);

            fragmentsToCut.Enqueue(topSlice);
            fragmentsToCut.Enqueue(bottomSlice);
        }

        int vertexCount = 0;
        int indexCount = 0;
        foreach (var fragment in fragmentsToCut)
        {
            vertexCount += fragment.TotalVertexCount;
            indexCount += fragment.TotalIndexCount;
            CreateFragment(fragment, fragmentTemplate, parent);
        }
        if (logVertexTriangleCount)
        {
            Debug.Log("Vertices after fracture: " + vertexCount);
            Debug.Log("Triangles after fracture: " + indexCount / 3);
        }
    }

    private void CreateFragment(MeshFragment fragment, GameObject fragmentTemplate, Transform parent)
    {
        if (fragment.TotalVertexCount <= 0) return;


        Mesh[] meshes;
        Mesh fragmentMesh = fragment.GetMesh();

        if (separateIslands)
        {
            meshes = SeparateMeshIslands(fragmentMesh);
            //meshes = new Mesh[] { };
        }
        else
        {
            meshes = new Mesh[] { fragmentMesh };
        }

        var parentSize = gameObject.GetComponent<MeshFilter>().sharedMesh.bounds.size;
        float parentVolume = parentSize.x * parentSize.y * parentSize.z;
        var parentMass = gameObject.GetComponent<Rigidbody>().mass;

        int vertexCount = 0;
        foreach (var mesh in meshes)
        {
            if (mesh.vertexCount < 4) continue;
            GameObject fragmentObj = GameObject.Instantiate(fragmentTemplate, parent);
            fragmentObj.name = "Fragment";
            fragmentObj.transform.SetLocalPositionAndRotation(Vector3.zero, Quaternion.identity);
            fragmentObj.transform.localScale = transform.localScale;

            mesh.name = Guid.NewGuid().ToString();
            fragmentObj.GetComponent<MeshFilter>().sharedMesh = mesh;

            var collider = fragmentObj.GetComponent<MeshCollider>();
            collider.sharedMesh = mesh;
            collider.convex = true;
            collider.sharedMaterial = fragmentObj.GetComponent<Collider>().sharedMaterial;

            var rigidBody = fragmentObj.GetComponent<Rigidbody>();
            var size = fragmentMesh.bounds.size;
            float volume = size.x * size.y * size.z;
            float density = parentVolume / parentMass;
            rigidBody.mass = volume / density;

            vertexCount += mesh.vertexCount;
        }
    }

    public static Mesh[] SeparateMeshIslands(Mesh mesh)
    {
        List<Mesh> meshIslands = new List<Mesh>();

        var vertices = mesh.vertices;
        var triangles = mesh.triangles;
        var normals = mesh.normals;
        var uvs = mesh.uv;

        int[] triangleSubMesh = new int[triangles.Length / 3];
        int subMeshIndex = 0;
        int subMeshSize = mesh.GetTriangles(subMeshIndex).Length / 3;
        for (int i = 0; i < triangles.Length / 3; i++)
        {
            if (i >= subMeshSize)
            {
                subMeshIndex++;
                subMeshSize += mesh.GetTriangles(subMeshIndex).Length / 3;
            }
            triangleSubMesh[i] = subMeshIndex;
        }

        List<int>[] coincidentVertices = new List<int>[vertices.Length];
        for (int i = 0; i < vertices.Length; i++)
        {
            coincidentVertices[i] = new List<int>();
        }
        for (int i = 0; i < vertices.Length; i++)
        {
            for (int k = i + 1; k < vertices.Length; k++)
            {
                if (vertices[i] == vertices[k])
                {
                    coincidentVertices[k].Add(i);
                    coincidentVertices[i].Add(k);
                }
            }
        }

        List<int>[] vertexTriangles = new List<int>[vertices.Length];

        for (int i = 0; i < vertices.Length; i++)
        {
            vertexTriangles[i] = new List<int>();
        }

        int v1, v2, v3;
        for (int i = 0; i < triangles.Length; i += 3)
        {
            int t = i / 3;

            v1 = triangles[i];
            v2 = triangles[i + 1];
            v3 = triangles[i + 2];

            vertexTriangles[v1].Add(t);
            vertexTriangles[v2].Add(t);
            vertexTriangles[v3].Add(t);
        }


        bool[] visitedVertices = new bool[vertices.Length];
        bool[] visitedTriangles = new bool[triangles.Length];
        Queue<int> search = new Queue<int>();

        NativeArray<MeshVertex> islandVertices = new NativeArray<MeshVertex>(vertices.Length, Allocator.Temp, NativeArrayOptions.UninitializedMemory);

        int[][] islandTriangles = new int[mesh.subMeshCount][];
        for (int i = 0; i < mesh.subMeshCount; i++)
        {
            islandTriangles[i] = new int[triangles.Length];
        }

        int vertexCount = 0;
        int totalIndexCount = 0;
        int[] subMeshIndexCounts = new int[mesh.subMeshCount];

        for (int i = 0; i < vertices.Length; i++)
        {
            if (visitedVertices[i]) continue;

            vertexCount = 0;
            totalIndexCount = 0;
            for (int j = 0; j < mesh.subMeshCount; j++)
            {
                subMeshIndexCounts[j] = 0;
            }

            search.Enqueue(i);

            int[] vertexMap = new int[vertices.Length];
            for (int j = 0; j < vertices.Length; j++)
            {
                vertexMap[j] = -1;
            }

            while (search.Count > 0)
            {
                int k = search.Dequeue();

                if (visitedVertices[k])
                {
                    continue;
                }
                else
                {
                    visitedVertices[k] = true;
                }

                vertexMap[k] = vertexCount;
                islandVertices[vertexCount++] = new MeshVertex(vertices[k], normals[k], uvs[k]);

                foreach (int t in vertexTriangles[k])
                {
                    if (!visitedTriangles[t])
                    {
                        visitedTriangles[t] = true;

                        for (int m = t * 3; m < t * 3 + 3; m++)
                        {
                            int v = triangles[m];
                            subMeshIndex = triangleSubMesh[t];
                            islandTriangles[subMeshIndex][subMeshIndexCounts[subMeshIndex]++] = v;
                            totalIndexCount++;

                            search.Enqueue(v);

                            foreach (int cv in coincidentVertices[v])
                            {
                                search.Enqueue(cv);
                            }
                        }
                    }
                }
            }

            if (vertexCount > 0)
            {
                Mesh meshIsland = new Mesh();

                meshIsland.SetIndexBufferParams(totalIndexCount, IndexFormat.UInt32);
                meshIsland.SetVertexBufferParams(vertexCount, layout);
                meshIsland.SetVertexBufferData(islandVertices, 0, 0, vertexCount);


                meshIsland.subMeshCount = mesh.subMeshCount;
                int indexStart = 0;
                for (subMeshIndex = 0; subMeshIndex < mesh.subMeshCount; subMeshIndex++)
                {
                    var subMeshIndexBuffer = islandTriangles[subMeshIndex];
                    var subMeshIndexCount = subMeshIndexCounts[subMeshIndex];


                    for (int k = 0; k < subMeshIndexCount; k++)
                    {
                        int originalIndex = subMeshIndexBuffer[k];
                        subMeshIndexBuffer[k] = vertexMap[originalIndex];
                    }

                    meshIsland.SetIndexBufferData(subMeshIndexBuffer, 0, indexStart, (int)subMeshIndexCount);
                    meshIsland.SetSubMesh(subMeshIndex, new SubMeshDescriptor(indexStart, subMeshIndexCount));

                    indexStart += subMeshIndexCount;
                }

                meshIsland.RecalculateBounds();
                meshIslands.Add(meshIsland);
            }
        }

        return meshIslands.ToArray();
    }

    private static VertexAttributeDescriptor[] layout = new[]
    {
        new VertexAttributeDescriptor(VertexAttribute.Position, VertexAttributeFormat.Float32, 3),
        new VertexAttributeDescriptor(VertexAttribute.Normal, VertexAttributeFormat.Float32, 3),
        new VertexAttributeDescriptor(VertexAttribute.TexCoord0, VertexAttributeFormat.Float32, 2),
    };
}
