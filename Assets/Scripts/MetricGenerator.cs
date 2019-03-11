using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using boundingmesh;

namespace boundingmesh
{
    using Real = System.Single; //System.Double
    using Index = System.Int32;

    enum Metric
    {
        ClassicQEM,
        ModifiedQEM,
        MinimizedConstant,
        Diagonalization,
        Average
    };

    enum Initialization
    {
        DistancePrimitives,
        Midpoint
    };

    class MatrixArray : Deque<Matrix4x4>
    { }

    class MetricGenerator
    {
        public MetricGenerator()
        {

        }

        public void SetMesh(Mesh mesh)
        {
            mesh_ = mesh;
            Initialize();
        }
        public void contractEdge(Index edge_index)
        {
            var edge = mesh_.edge(edge_index);
            switch (metric_)
            {
                case Metric.ClassicQEM:
                    qems_[edge.vertex(1)] = getErrorMetric(edge_index);
                    break;

                case Metric.Average:
                    qems_merge_[edge.vertex(1)] = getErrorMetric(edge_index);
                    break;

                default:
                    // todo.
                    throw new System.Exception("undefine");
            }
        }
        public void cleanAndRenumber()
        {
            switch (metric_)
            {
                case Metric.ClassicQEM:
                    shrinkIndexedArray(ref qems_, mesh_.deleted_vertices);
                    break;

                case Metric.Average:
                    shrinkIndexedArray(ref qems_merge_, mesh_.deleted_vertices);
                    break;

                default:
                    throw new System.Exception("not defined");
            }
        }
        void shrinkIndexedArray(ref MatrixArray array, Stack<Index> deleted_indices)
        {
            MatrixArray newArray = new MatrixArray();
            List<Index> deleted_indices_sorted = new List<Index>(deleted_indices.Count);

            deleted_indices_sorted.AddRange(deleted_indices);
            deleted_indices_sorted.Sort();

            if (array.Count <= deleted_indices_sorted.Count)
            {
                return;
            }

            var next_index = 0;
            var deleted_i = 0;
            var valid_num = array.Count - deleted_indices_sorted.Count;
            for (int i = 0; i < valid_num; i++)
            {
                while (deleted_i < deleted_indices_sorted.Count && next_index == deleted_indices_sorted[deleted_i])
                {
                    next_index++;
                    deleted_i++;
                }

                newArray.Add(array[next_index]);
                next_index++;
            }

            array = newArray;
        }

        void Initialize()
        {
            if (mesh_ == null)
            {
                return;
            }

            switch (metric_)
            {
                case Metric.ClassicQEM:
                    qems_.Clear();
                    for (int i = 0; i < mesh_.nVertices(); i++)
                    {
                        qems_.Add(computeQEM(i));
                    }
                    break;

                case Metric.Average:
                    qems_merge_.Clear();
                    foreach (var v in mesh_.vertices_)
                    {
                        qems_merge_.Add(computeInitialMergeMetric(v));
                    }
                    break;
                default:
                    // todo.
                    throw new System.Exception("unsupported.");
            }
        }
        Matrix4x4 computeInitialMergeMetric(Vertex vertex)
        {
            var qem = Matrix4x4.zero;

            if (initialization_ == Initialization.Midpoint)
            {
                foreach (var t in vertex.triangles_)
                {
                    var tri = mesh_.triangle(t);
                    var pt = Vector3.zero;
                    foreach (var v in tri.vertices)
                    {
                        var vt = mesh_.vertex(v);
                        pt = pt + vt.position();
                    }
                    pt = pt / tri.vertices.Length;
                    var distance_point = Matrix4x4.identity;
                    distance_point.SetColumn(3, new Vector4(-pt.x, -pt.y, -pt.z, 0));
                    distance_point = distance_point.transpose * distance_point;
                    qem = mergeMatrices(qem, distance_point);
                }

                {
                    var pt = vertex.position();
                    var distance_point = Matrix4x4.identity;
                    distance_point.SetColumn(3, new Vector4(-pt.x, -pt.y, -pt.z, 0));
                    distance_point = distance_point.transpose * distance_point;
                    qem = mergeMatrices(qem, distance_point);
                }

            }
            else
            {
                throw new System.Exception("undefine");
            }

            return qem;
        }
        Matrix4x4 mergeMatrices(Matrix4x4 a, Matrix4x4 b)
        {
            switch (metric_)
            {
                case Metric.Average:
                    return mergeAverage(a, b);
                default:
                    throw new System.Exception("undefine");
            }
        }
        Matrix4x4 mergeAverage(Matrix4x4 a, Matrix4x4 b)
        {
            var avg = a.Add(b).Scale(0.5f);
            return avg;
        }
        public Matrix4x4 getErrorMetric(Index edge_index)
        {
            var qem = Matrix4x4.zero;
            var edge = mesh_.edge(edge_index);
            switch (metric_)
            {
                case Metric.ClassicQEM:
                    qem = qems_[edge.vertex(0)].Add(qems_[edge.vertex(1)]);
                    break;

                case Metric.Average:
                    qem = mergeMatrices(qems_merge_[edge.vertex(0)], qems_merge_[edge.vertex(1)]);
                    break;

                default:
                    throw new System.Exception("undefine");
            }
            
            return qem;
        }

        Matrix4x4 computeQEM(Index vertex_index)
        {
            Matrix4x4 qem = Matrix4x4.zero;
            var vertex = mesh_.vertex(vertex_index);
            for (int i = 0; i < vertex.nTriangles(); i++)
            {
                qem = qem.Add(mesh_.triangle(i).plane().distanceMatrix());
            }
            return qem;
        }


        Metric metric_ = Metric.Average;
        Initialization initialization_ = Initialization.Midpoint;
        Mesh mesh_ = null;

        MatrixArray qems_ = new MatrixArray();
        MatrixArray qems_merge_ = new MatrixArray();
    }
}
