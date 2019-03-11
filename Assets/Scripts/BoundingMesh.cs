using System.Collections;
using System.Collections.Generic;
using UnityEngine;

#if UNITY_EDITOR
using UnityEditor;
#endif

public class BoundingMesh : MonoBehaviour {

    public MeshRenderer R;
    public MeshCollider C;

    public int Count = 0;

	// Use this for initialization
	void Start () {
		
	}
	
	// Update is called once per frame
	void Update () {
		
	}


}


#if UNITY_EDITOR

[CustomEditor(typeof(BoundingMesh))]
class BoundingMeshEditor : Editor
{
    public override void OnInspectorGUI()
    {
        base.OnInspectorGUI();
        if (GUILayout.Button("构建模型"))
        {
            var self = this.target as BoundingMesh;
            var filter = self.R.GetComponent<MeshFilter>();
            var oldMesh = filter.sharedMesh;
            var mesh = new Mesh();

            mesh.vertices = oldMesh.vertices;
            mesh.uv = oldMesh.uv;
            mesh.triangles = oldMesh.triangles;

            mesh.RecalculateNormals();
            self.C.sharedMesh = mesh;
        }
        if (GUILayout.Button("裁剪碰撞"))
        {
            var self = this.target as BoundingMesh;
            var filter = self.R.GetComponent<MeshFilter>();
            var oldMesh = filter.sharedMesh;

            boundingmesh.diagnosis.Assert(oldMesh.triangles.Length % 3 == 0);

            var mesh = new boundingmesh.Mesh();
            foreach (var v in oldMesh.vertices)
            {
                mesh.AddVertex(v);
            }
            for(int i=0; i<oldMesh.triangles.Length; i+=3)
            {
                mesh.addTriangle(oldMesh.triangles[i], oldMesh.triangles[i+1], oldMesh.triangles[i+2]);
            }

            var decimator = new boundingmesh.Decimate();
            decimator.SetMesh(mesh);
            decimator.setTargetVertices(oldMesh.triangles.Length / 3 - self.Count);
            var deciMesh = decimator.Compute();

            UnityEngine.Debug.Log("模型顶点个数:" + deciMesh.nVertices() + ", 旧的：" + oldMesh.vertices.Length);
            UnityEngine.Debug.Log("模型三角形个数:" + deciMesh.nTriangles() + ", 旧的：" + (oldMesh.triangles.Length / 3));

            var newMesh = new Mesh();
            var verts = new List<Vector3>();
            var trs = new List<int>();
            foreach (var v in deciMesh.vertices_)
            {
                verts.Add(v.position());
            }
            foreach (var v in deciMesh.triangles_)
            {
                trs.AddRange(v.vertices);
            }
            newMesh.vertices = verts.ToArray();
            newMesh.triangles = trs.ToArray();
            newMesh.RecalculateNormals();
            self.C.sharedMesh = newMesh;
        }
    }
}


#endif