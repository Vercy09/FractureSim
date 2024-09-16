using System;
using System.Collections;
using System.Collections.Generic;
using Unity.VisualScripting;
using UnityEngine;

public class Test : MonoBehaviour
{
    // Start is called before the first frame update
    void Start()
    {
        /*
        List<MeshVertex> vertices = new List<MeshVertex>
        {
            new MeshVertex(new Vector3(100, -100, 40)),
            new MeshVertex(new Vector3(200, -100, 30)),
            new MeshVertex(new Vector3(15, -10, 70)),
            new MeshVertex(new Vector3(-50, 100, 40)),
            new MeshVertex(new Vector3(400, 50, 120)),
            new MeshVertex(new Vector3(30, -30, -40)),
            new MeshVertex(new Vector3(330, -230, 100)),
            new MeshVertex(new Vector3(30, -30, 120)),
            new MeshVertex(new Vector3(303, 210, -2)),
            new MeshVertex(new Vector3(12, 320, -12)),
        };

        List<Edge> constraints = new List<Edge>
        {
            new Edge(0, 1),
            new Edge(1, 2),
            new Edge(4, 5),
            new Edge(2, 3),
        };

        Vector3 n = new Vector3(100, 20, 30);

        //Triangulator t = new Triangulator(vertices, n);
        ConstrainedTriangulator t = new ConstrainedTriangulator(vertices, constraints, n);
        var res = t.Triangulate(true);


        for (int i = 0; i < res.Length - 3; i += 3)
        {
            string str = "{" + res[i] + "," + res[i + 1] + "," + res[i + 2] + "}";
            Debug.Log(str);
        }*/

    }

    void Awake()
    {
        QualitySettings.vSyncCount = 0;
        Application.targetFrameRate = 300;
    }
}
