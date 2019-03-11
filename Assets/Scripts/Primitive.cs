using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using boundingmesh;
namespace boundingmesh
{
	using Real = System.Single; //System.Double
	using Index = System.Int32;

    public static class Utils
    {
        public static Matrix4x4 Add(this Matrix4x4 self, Matrix4x4 right)
        {
            var ret = self;
            ret.m00 += right.m00; ret.m01 += right.m01; ret.m02 += right.m02; ret.m03 += right.m03;
            ret.m10 += right.m10; ret.m11 += right.m11; ret.m12 += right.m12; ret.m13 += right.m13;
            ret.m20 += right.m20; ret.m21 += right.m21; ret.m22 += right.m22; ret.m23 += right.m23;
            ret.m30 += right.m30; ret.m31 += right.m31; ret.m32 += right.m32; ret.m33 += right.m33;

            return ret;
        }
        public static Matrix4x4 Scale(this Matrix4x4 self, float scalar)
        {
            var ret = self;

            ret.m00 *= scalar; ret.m01 *= scalar; ret.m02 *= scalar; ret.m03 *= scalar;
            ret.m10 *= scalar; ret.m11 *= scalar; ret.m12 *= scalar; ret.m13 *= scalar;
            ret.m20 *= scalar; ret.m21 *= scalar; ret.m22 *= scalar; ret.m23 *= scalar;
            ret.m30 *= scalar; ret.m31 *= scalar; ret.m32 *= scalar; ret.m33 *= scalar;

            return ret;
        }
        public static void CopyArray(ref Index[] dst, Index[] src)
        {
            for (int i = 0; i < src.Length; i++)
            {
                dst[i] = src[i];
            }
        }
        public static float Norm(Matrix4x4 mat)
        {
            float sum = 0;

            sum += Mathf.Pow(mat.m00, 2); sum += Mathf.Pow(mat.m01, 2); sum += Mathf.Pow(mat.m02, 2);
            sum += Mathf.Pow(mat.m10, 2); sum += Mathf.Pow(mat.m11, 2); sum += Mathf.Pow(mat.m12, 2);
            sum += Mathf.Pow(mat.m20, 2); sum += Mathf.Pow(mat.m21, 2); sum += Mathf.Pow(mat.m22, 2);


            return Mathf.Sqrt(sum);
        }
        public static Matrix4x4 VvMatrix(Vector3 v1, Vector3 v2)
        {
            var m1 = new Matrix4x4();
            m1.SetColumn(0, v1);
            var m2 = new Matrix4x4();
            m2.SetRow(0, v2);

            return m1 * m2;
        }
    }


	public class Plane
	{
	    public Vector3 normal;
	    public Real d;

        public Plane Clone()
        {
            Plane ret = new Plane();

            ret.normal = normal;
            ret.d = d;

            return ret;
        }
	
	    public Plane()
	    {
	        normal = Vector3.zero;
	        d = 0;
	    }
	    public Plane(Vector3 n, Real dd)
	    {
	        normal = n.normalized;
	        d = dd;
	    }
	    public Plane(Vector3 n, Vector3 point)
	    {
	        normal = n.normalized;
	        d = -1 * Vector3.Dot(normal, point);
	    }

        public Vector3 p1;
        public Vector3 p2;
        public Vector3 p3;
	    public Plane(Vector3 point1, Vector3 point2, Vector3 point3)
	    {
            p1 = point1;
            p2 = point2;
            p3 = point3;

	        normal = Vector3.Cross((point2 - point1), (point3 - point1));
	        normal.Normalize();
	        d = -1 * Vector3.Dot(normal, point1);
	    }
	
	
	    public Real distance(Vector3 point)
	    {
	        return Vector3.Dot(normal, point);
	    }
	    public Vector4 vector()
	    {
	        return new Vector4(normal.x, normal.y, normal.z, d);
	    }
	    public Matrix4x4 distanceMatrix()
	    {
            //var m1 = new Matrix4x4(vector(), new Vector4(0, 1, 0, 0), new Vector4(0, 0, 1, 0), new Vector4(0, 0, 0, 1));
            //var m2 = new Matrix4x4(vector(), new Vector4(0, 1, 0, 0), new Vector4(0, 0, 1, 0), new Vector4(0, 0, 0, 1));
            var v = vector();
            var m1 = new Matrix4x4();
            var m2 = new Matrix4x4();
            m1.SetColumn(0, v);
            m2.SetRow(0, v);

            var m = m1 * m2;

            return m;
	    }
	
	}
	
	public class Vertex
	{
	    Vector3 position_ = Vector3.zero;
	    public List<Index> triangles_ = new List<Index>();
	    public List<Index> edges_ = new List<Index>();

        public Vertex Clone()
        {
            Vertex ret = new Vertex();
            ret.position_ = position_;
            ret.triangles_.AddRange(triangles_);
            ret.edges_.AddRange(edges_);

            return ret;
        }

	    public Vertex()
	    {
	        
	    }
	    public Vertex(Vector3 pos)
	    {
	        position_ = pos;
	    }
	    public Vertex(Real x, Real y, Real z)
	    {
	        position_ = new Vector3(x, y, z);
	    }
	
	    public Vector3 position()
	    {
	        return position_;
	    }
	    public int nTriangles()
	    {
	        return triangles_.Count;
	    }
	    public int nEdges()
	    {
	        return edges_.Count;
	    }
	    public Index triangle(int i)
	    {
	        return triangles_[i];
	    }
	    public Index edge(int i)
	    {
	        return edges_[i];
	    }
	
	}
	
	public class Edge
	{
	    public Index[] vertices_ = { 0, 0 };
	    public List<Index> triangles_ = new List<Index>();

        public Edge Clone()
        {
            Edge ret = new Edge();

            Utils.CopyArray(ref ret.vertices_, vertices_);
            ret.triangles_.AddRange(this.triangles_);            

            return ret;
        }
        public Edge()
        { }
        public Edge(Index i1, Index i2)
	    {
	        vertices_[0] = i1;
	        vertices_[1] = i2;
	    }
	
	    public int nTriangles()
	    {
	        return triangles_.Count;
	    }
	
	    public bool border()
	    {
	        return triangles_.Count == 1;
	    }
	
	    public Index vertex(int i)
	    {
	        return vertices_[i];
	    }
	
	    public Index triangle(int i)
	    {
	        return triangles_[i];
	    }
	}
	
	public class Triangle
	{
        public Triangle Clone()
        {
            var ret = new Triangle();

            Utils.CopyArray(ref ret.edges_, edges_);
            Utils.CopyArray(ref ret.vertices, vertices);
            ret.plane_ = plane_.Clone();

            return ret;
        }
        public Triangle()
        { }
        public Triangle(Index i1, Index i2, Index i3)
	    {
	        vertices[0] = i1;
	        vertices[1] = i2;
	        vertices[2] = i3;
	    }
	    public Index vertex(int i)
	    {
	        return vertices[i];
	    }
	    public Index edge(int i)
	    {
	        return edges_[i];
	    }
	    public Plane plane()
	    {
	        return plane_;
	    }
	
	    public Index[] vertices = { 0, 0, 0 };
        public Index[] edges_ = { 0, 0, 0 };
        public Plane plane_ = new Plane();
	}
}
