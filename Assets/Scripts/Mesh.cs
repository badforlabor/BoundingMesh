using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using boundingmesh;

namespace boundingmesh
{
    using Real = System.Single; //System.Double
    using Index = System.Int32;

    static class diagnosis
    {
        public static void Assert(bool b)
        {
            if (!b)
            {
                throw new Exception("assert failed!");
            }
        }
    }


    class VertexPosition
    {
        public const Real vertex_epsilon = 0.000001f;

        public VertexPosition(Vector3 position)
        {
            position_ = position;
        }
        public VertexPosition(Index idx, Vector3 pos)
        {
            index_ = idx;
            position_ = pos;
        }
        public Index index()
        {
            return index_;
        }

        public static bool operator<(VertexPosition left, VertexPosition right)
        {
            if (Math.Abs(left.position_.x - right.position_.x) > vertex_epsilon)
            {
                return left.position_.x < right.position_.x;
            }
            else if (Math.Abs(left.position_.y - right.position_.y) > vertex_epsilon)
            {
                return left.position_.y < right.position_.y;
            }
            else
            {
                return left.position_.z < right.position_.z;
            }
        }
        public static bool operator >(VertexPosition left, VertexPosition right)
        {
            return right < left;
        }

        Index index_;
        Vector3 position_;
        bool searching;
    }

    class VertexPositionSet
    {
        public VertexPositionSet(Mesh mesh)
        {
            mesh_ = mesh;
        }
        public Index addVertex(Vector3 position)
        {
            throw new System.Exception("todo. undefine");
        }
        Mesh mesh_;
        HashSet<VertexPosition> set_ = new HashSet<VertexPosition>();
    }

    class Deque<T> : List<T>
    {
        public T Top()
        {
            return this[0];
        }
        public void Pop()
        {
            this.RemoveAt(0);
        }
        public void PushBack(T obj)
        {
            this.Add(obj);
        }

        public T this[uint idx]
        {
            get
            {
                return this[(int)idx];
            }
            set
            {
                this[(int)idx] = value;
            }
        }
    }
    class Stack<T> : Deque<T>
    {

    }

    class Mesh
    {
        public Mesh Clone()
        {
            Mesh ret = new Mesh();

            ret.deleted_edges.AddRange(deleted_edges);
            ret.deleted_triangles_.AddRange(deleted_triangles_);
            ret.deleted_vertices.AddRange(deleted_vertices);
            ret.vertices_.AddRange(vertices_);
            ret.edges_.AddRange(edges_);
            ret.triangles_.AddRange(triangles_);

            ret.n_valid_edges_ = n_valid_edges_;
            ret.n_valid_triangles_ = n_valid_triangles_;
            ret.n_valid_vertices_ = n_valid_vertices_;

            return ret;
        }

        class AutoDeleteIndex
        {
            public AutoDeleteIndex(Stack<Index> deleted_vertices)
            {
                deleted_vertices_sorted = new List<Index>(deleted_vertices.Count);
                deleted_vertices_sorted.AddRange(deleted_vertices);
                deleted_vertices_sorted.Sort();
            }

            int deleted_i = 0;
            int next_index = 0;
            List<Index> deleted_vertices_sorted = null;

            public int NextIndex()
            {
                while (deleted_i < deleted_vertices_sorted.Count 
                    && next_index == deleted_vertices_sorted[deleted_i])
                {
                    deleted_i++;
                    next_index++;
                }
                var ret = next_index;
                next_index++;

                return ret;
            }
        }

        public void cleanAndRenumber()
        {
            AutoDeleteIndex autoVertex = new AutoDeleteIndex(deleted_vertices);
            var new_vertices = new Deque<Vertex>();
            for (int i = 0; i < n_valid_vertices_; i++)
            {
                var next_index = autoVertex.NextIndex();
                var new_index = new_vertices.Count;

                var vt = vertex(next_index);

                foreach (var ie in vt.edges_)
                {
                    var e = edge(ie);
                    for (int k = 0; k < 2; k++)
                    {
                        if (e.vertex(k) == next_index)
                        {
                            e.vertices_[k] = new_index;
                            break;
                        }
                    }
                }

                foreach (var it in vt.triangles_)
                {
                    var t = triangles_[it];
                    for (int k = 0; k < 3; k++)
                    {
                        if (t.vertex(k) == next_index)
                        {
                            t.vertices[k] = new_index;
                            break;
                        }
                    }
                }

                new_vertices.Add(vt);
            }
            vertices_ = new_vertices;

            AutoDeleteIndex autoTriangle = new AutoDeleteIndex(deleted_triangles_);
            var new_triangles = new Deque<Triangle>();
            for (int i = 0; i < n_valid_triangles_; i++)
            {
                var next_index = autoTriangle.NextIndex();
                var new_index = new_triangles.Count;
                var tri = triangles_[next_index];

                foreach(var iv in tri.vertices)
                {
                    var vt = vertices_[iv];
                    for(int k=0; k<vt.triangles_.Count; k++)
                    {
                        if (vt.triangles_[k] == next_index)
                        {                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
                            vt.triangles_[k] = new_index;
                            break;
                        }
                    }
                }

                foreach(var ie in tri.edges_)
                {
                    var e = edge(ie);
                    for (int k = 0; k < e.nTriangles(); k++)
                    {
                        if (e.triangles_[k] == next_index)
                        {
                            e.triangles_[k] = new_index;
                            break;
                        }
                    }
                }

                new_triangles.Add(tri);
            }
            triangles_ = new_triangles;

            var autoEdge = new AutoDeleteIndex(deleted_edges);
            var new_edges = new Deque<Edge>();
            for (int i = 0; i < n_valid_edges_; i++)
            {
                var next_index = autoEdge.NextIndex();
                var new_index = new_edges.Count;
                var e = edge(next_index);

                for (int k = 0; k < 2; k++)
                {
                    var ver = vertices_[e.vertices_[k]];
                    for (int q = 0; q < ver.edges_.Count; q++)
                    {
                        if (ver.edges_[q] == next_index)
                        {
                            ver.edges_[q] = new_index;
                            break;
                        }
                    }
                }
                for (int k = 0; k < e.triangles_.Count; k++)
                {
                    var tri = triangles_[e.triangles_[k]];
                    for (int q = 0; q < tri.edges_.Length; q++)
                    {
                        if (tri.edges_[q] == next_index)
                        {
                            tri.edges_[q] = new_index;
                            break;
                        }
                    }
                }

                new_edges.Add(e);
            }
            edges_ = new_edges;
            
        }
        public Index AddVertex(Vector3 position)
        {
            Index new_index = 0;
            if (deleted_vertices.Count > 0)
            {
                new_index = deleted_vertices.Top();
                deleted_vertices.Pop();
                vertices_[new_index] = new Vertex(position);
            }
            else
            {
                new_index = vertices_.Count;
                vertices_.PushBack(new Vertex(position));
            }
            n_valid_vertices_++;

            return new_index;
        }
        public Index addTriangle(Index vertex1, Index vertex2, Index vertex3)
        {
            //if (vertex1 == vertex2 || vertex1 == vertex3 || vertex2 == vertex3)
            //{
            //    return 0;
            //}

            diagnosis.Assert(vertex1 < vertices_.Count);
            diagnosis.Assert(vertex2 < vertices_.Count);
            diagnosis.Assert(vertex3 < vertices_.Count);
            diagnosis.Assert(vertex1 != vertex2);
            diagnosis.Assert(vertex1 != vertex3);
            diagnosis.Assert(vertex2 != vertex3);
            Index new_index = 0;
            if (deleted_triangles_.Count > 0)
            {
                new_index = deleted_triangles_.Top();
                deleted_triangles_.Pop();
                triangles_[new_index] = new Triangle();
            }
            else
            {
                new_index = triangles_.Count;
                triangles_.PushBack(new Triangle());
            }

            Triangle triangle = new Triangle(vertex1, vertex2, vertex3);
            n_valid_triangles_++;

            triangle.plane_ = new Plane(vertices_[vertex1].position(),
                vertices_[vertex2].position(),
                vertices_[vertex3].position());

            triangle.edges_[0] = registerEdge(vertex1, vertex2, new_index);
            triangle.edges_[1] = registerEdge(vertex2, vertex3, new_index);
            triangle.edges_[2] = registerEdge(vertex3, vertex1, new_index);
            triangles_[new_index] = triangle;

            vertices_[vertex1].triangles_.Add(new_index);
            vertices_[vertex2].triangles_.Add(new_index);
            vertices_[vertex3].triangles_.Add(new_index);

            return new_index;
        }
        public Index registerEdge(Index vertex1, Index vertex2, Index triangle)
        {
            diagnosis.Assert(vertex1 != vertex2);
            diagnosis.Assert(vertex1 < vertices_.Count);
            diagnosis.Assert(vertex2 < vertices_.Count);
            diagnosis.Assert(triangle < triangles_.Count);

            Index first_vertex = vertex1;
            Index second_vertex = vertex2;

            // 规定，edge的两个点，依次是由小到大的。
            if (vertex1 > vertex2)
            {
                first_vertex = vertex2;
                second_vertex = vertex1;
            }

            bool already_exists = false;
            Index existing_edge_index = 0;
            for (int i = 0; i < vertices_[first_vertex].nEdges(); i++)
            {
                var edge_index = vertices_[first_vertex].edge(i);
                if (edges_[edge_index].vertex(1) == second_vertex)
                {
                    already_exists = true;
                    existing_edge_index = edge_index;
                    break;
                }
            }

            Index new_index = 0;
            if (!already_exists)
            {
                var edge = new Edge(first_vertex, second_vertex);
                edge.triangles_.Add(triangle);
                if (deleted_edges.Count > 0)
                {
                    new_index = deleted_edges.Top();
                    deleted_edges.Pop();
                    edges_[new_index] = edge;
                }
                else
                {
                    new_index = edges_.Count;
                    edges_.PushBack(edge);
                }

                n_valid_edges_++;
                vertices_[vertex1].edges_.Add(new_index);
                vertices_[vertex2].edges_.Add(new_index);
            }
            else
            {
                new_index = existing_edge_index;
                diagnosis.Assert(!edges_[existing_edge_index].triangles_.Contains(triangle));
                edges_[existing_edge_index].triangles_.Add(triangle);
            }

            return new_index;
        }
        public void RemoveVertex(Index index)
        {
            if (this.vertex(index) == null)
            {
                return;
            }

            while (this.vertex(index).triangles_.Count > 0)
            {
                RemoveTriangle(this.vertex(index).triangles_[0]);
            }

            deleted_vertices.Add(index);
            vertices_[index] = new Vertex();
            n_valid_vertices_--;
            dirty_ = true;
        }
        public void RemoveTriangle(Index index)
        {
            var t = triangle(index);
            if (t == null)
            {
                return;
            }
            for (int i = 0; i < 3; i++)
            {
                var v = this.vertex(t.vertex(i));
                v.triangles_.Remove(index);
            }

            for (int i = 0; i < 3; i++)
            {
                var edge_index = t.edge(i);
                var e = this.edge(edge_index);

                if (e.border())
                {
                    var v1 = this.vertex(e.vertex(0));
                    v1.edges_.Remove(edge_index);

                    var v2 = this.vertex(e.vertex(1));
                    v2.edges_.Remove(edge_index);

                    deleted_edges.Add(edge_index);
                    this.edges_[edge_index] = new Edge();
                    n_valid_edges_--;
                }
                else
                {
                    e.triangles_.Remove(index);
                }
            }

            this.triangles_[index] = new Triangle();
            deleted_triangles_.Add(index);
            n_valid_triangles_--;
            dirty_ = true;
        }
        public int nVertices()
        {
            return n_valid_vertices_;
        }
        public Vertex vertex(Index i)
        {
            return vertices_[i];
        }
        public int nTriangles()
        {
            return n_valid_triangles_;
        }
        public Triangle triangle(Index i)
        {
            return triangles_[i];
        }
        public int nEdges()
        {
            return n_valid_edges_;
        }
        public Edge edge(Index i)
        {
            return edges_[i];
        }


        // 这两个是一对。用来构建可复用的数组。
        public Deque<Vertex> vertices_ = new Deque<Vertex>();
        public Stack<Index> deleted_vertices = new Stack<Index>();

        Deque<Edge> edges_ = new Deque<Edge>();
        Stack<Index> deleted_edges = new Stack<Index>();
        public Deque<Triangle> triangles_ = new Deque<Triangle>();
        Stack<Index> deleted_triangles_ = new Stack<Index>();
        int n_valid_triangles_ = 0;
        int n_valid_vertices_ = 0;
        int n_valid_edges_ = 0;
        bool dirty_ = false;
    }
}
