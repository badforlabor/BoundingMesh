using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using boundingmesh;
using MathNet.Numerics.LinearAlgebra;

namespace boundingmesh
{
    using Real = System.Single; //System.Double
    using Index = System.Int32;
    
    enum DecimationDirection
    {
        Outward,
        Inward,
        Any
    };

    class Decimate
    {
        public const int default_target_vertices = 1000;
        public const Real default_maximum_error = 1.0f;

        public const Real pi = (Real)Math.PI;
        public const Real phi = (2 * pi) * (30 / 360);

        public void SetMesh(Mesh mesh)
        {
            result_mesh_ = mesh.Clone();
            metric_generator_.SetMesh(result_mesh_);

            recomputeQueue();
        }

        public void setTargetVertices(int cnt)
        {
            target_vertices = cnt;
            target_vertices_used_ = true;
        }

        public Mesh Compute()
        {
            if (queue_.Size() == 0)
            {
                result_mesh_ = new Mesh();
                return result_mesh_;
            }

            if (!target_vertices_used_ && !maximum_error_used_)
            {
                return result_mesh_;
            }
            int cnt = result_mesh_.nVertices();
            while (cnt >= 0 && (queue_.Size() > 0) && (!target_vertices_used_ || result_mesh_.nTriangles() > target_vertices)
                && (!maximum_error_used_ || queue_.First().cost() < maximum_error_))
            {
                cnt--;
                var c = queue_.First();
                var changed = executeEdgeContraction(c);
                if (!changed)
                {
                    break;
                }
            }

            cleanAndRenumber();

            return result_mesh_;
        }

        void cleanAndRenumber()
        {
            ContractionQueue new_queue = new ContractionQueue();
            Index next_index = 0;
            foreach (var ci in queue_.indices_)
            {
                if (ci == null)
                {
                    continue;
                }
                new_queue.Insert(new EdgeContraction(next_index, ci.new_point(), ci.cost(), ci.qem()));
                next_index++;
            }

            queue_ = new_queue;

            metric_generator_.cleanAndRenumber();
            result_mesh_.cleanAndRenumber();
        }

        bool executeEdgeContraction(EdgeContraction c)
        {
            var changed = false;
            current_error_ = c.cost();
            var edge = result_mesh_.edge(c.edge());
            if (edge.nTriangles() == 0)
            {
                queue_.Remove(c.edge());
                return true;
            }
            metric_generator_.contractEdge(c.edge());

            HashSet<Index> edge_to_remove = new HashSet<Index>();
            List<List<Index>> hole_border = new List<List<Index>>();
            collectRemovalData(edge.vertex(0), edge.vertex(1), edge_to_remove, hole_border);
            collectRemovalData(edge.vertex(1), edge.vertex(0), edge_to_remove, hole_border);

            foreach (var it in edge_to_remove)
            {
                queue_.Remove(it);
                changed = true;
            }

            result_mesh_.RemoveVertex(edge.vertex(0));
            result_mesh_.RemoveVertex(edge.vertex(1));

            var new_vertex_index = result_mesh_.AddVertex(c.new_point());
            for (int i = 0; i < hole_border.Count; i++)
            {
                result_mesh_.addTriangle(new_vertex_index, hole_border[i][0], hole_border[i][1]);
            }
            hole_border.Clear();

            HashSet<Index> edges_to_add = new HashSet<Index>();
            foreach (var ti in result_mesh_.vertex(new_vertex_index).triangles_)
            {
                var t = result_mesh_.triangle(ti);
                foreach (var e in t.edges_)
                {
                    edges_to_add.Add(e);
                }
            }

            foreach (var e in edges_to_add)
            {
                var cc = computeEdgeContraction(e);
                if (this.result_mesh_.edge(cc.edge()).nTriangles() == 0)
                {
                    diagnosis.Assert(false);
                }
                queue_.Insert(cc);
                changed = true;
            }

            return changed;
        }

        void collectRemovalData(Index index, Index other_index, HashSet<Index> edge_to_remove, List<List<Index>> hole_border)
        {
            var vertex = result_mesh_.vertex(index);
            for (int i = 0; vertex != null && i < vertex.triangles_.Count; i++)
            {
                var triangle = result_mesh_.triangle(vertex.triangle(i));
                bool is_shared_triangle = false;
                for (int j = 0; j < 3; j++)
                {
                    edge_to_remove.Add(triangle.edge(j));
                    if (triangle.vertex(j) == other_index)
                    {
                        is_shared_triangle = true;
                    }
                }
                if (!is_shared_triangle)
                {
                    var border = new List<Index>();
                    if (triangle.vertex(0) == index)
                    {
                        border.Add(triangle.vertex(1));
                        border.Add(triangle.vertex(2));
                    }
                    else if (triangle.vertex(1) == index)
                    {
                        border.Add(triangle.vertex(2));
                        border.Add(triangle.vertex(0));
                    }
                    else if (triangle.vertex(2) == index)
                    {
                        border.Add(triangle.vertex(0));
                        border.Add(triangle.vertex(1));
                    }
                    hole_border.Add(border);
                }
            }
        }

        void recomputeQueue()
        {
            queue_ = new ContractionQueue();
            for (int i = 0; i < result_mesh_.nEdges(); i++)
            {
                var c = computeEdgeContraction(i);
                if (this.result_mesh_.edge(c.edge()).nTriangles() == 0)
                {
                    diagnosis.Assert(false);
                }
                queue_.Insert(c);
            }
        }
        EdgeContraction computeEdgeContraction(Index edge_index)
        {
            var edge = result_mesh_.edge(edge_index);
            var vertex1 = result_mesh_.vertex(edge.vertex(0));
            var vertex2 = result_mesh_.vertex(edge.vertex(1));

            var qem = metric_generator_.getErrorMetric(edge_index);
            var new_point = Vector3.zero;

            var border_vertices = 0;
            for (int i = 0; i < vertex1.nEdges(); i++)
            {
                if (result_mesh_.edge(i).border())
                {
                    new_point = vertex1.position();
                    border_vertices++;
                    break;
                }
            }
            for (int i = 0; i < vertex2.nEdges(); i++)
            {
                if (result_mesh_.edge(i).border())
                {
                    new_point = vertex2.position();
                    border_vertices++;
                    break;
                }
            }
            if (border_vertices == 2 && !edge.border())
            {
                return new EdgeContraction(edge_index, Vector3.zero, Real.MaxValue, Matrix4x4.zero);
            }

            if (direction_ == DecimationDirection.Any)
            {
                if (border_vertices == 2 || border_vertices == 0)
                {
                    var found_valid = solveConstrainedMinimization(qem, new List<Plane>(), new List<Index>(), direction_, ref new_point);
                    diagnosis.Assert(found_valid);
                }

                var new_point_homegeneous = new Vector4(new_point.x, new_point.y, new_point.z, 1);
                var cost = MathUtil.VtmV(new_point_homegeneous, qem);
                return new EdgeContraction(edge_index, new_point, cost, qem);
            }

            var constraints = new List<Plane>();
            for (int i = 0; i < vertex1.triangles_.Count; i++)
            {
                var t = result_mesh_.triangle(vertex1.triangles_[i]);
                constraints.Add(t.plane());
            }

            for (int i = 0; i < vertex2.triangles_.Count; i++)
            {
                var t = result_mesh_.triangle(vertex2.triangles_[i]);
                bool shared = false;
                foreach (var vt in t.vertices)
                {
                    if (vt == edge.vertex(0))
                    {
                        shared = true;
                        break;
                    }
                }
                if (!shared)
                {
                    constraints.Add(t.plane());
                }
            }

            if (border_vertices == 1)
            {
                var valid_solution = true;
                var new_point_homogeneous = new Vector4(new_point.x, new_point.y, new_point.z, 1);
                for (int i = 0; i < constraints.Count; i++)
                {
                    if (direction_ == DecimationDirection.Outward)
                    {
                        var cost = Vector4.Dot(constraints[i].vector(), new_point_homogeneous);
                        if (cost + Mathf.Epsilon < 0)
                        {
                            valid_solution = false;
                            break;
                        }
                    }
                    else if (direction_ == DecimationDirection.Inward)
                    {
                        var cost = Vector4.Dot(constraints[i].vector(), new_point_homogeneous);
                        if (cost - Mathf.Epsilon > 0)
                        {
                            valid_solution = false;
                            break;
                        }
                    }
                }
                if (valid_solution)
                {
                    var new_point_homegeneous = new Vector4(new_point.x, new_point.y, new_point.z, 1);
                    var cost = MathUtil.VtmV(new_point_homegeneous, qem);
                    return new EdgeContraction(edge_index, new_point, cost, qem);
                }
                else
                {
                    return new EdgeContraction(edge_index, Vector3.zero, Real.MaxValue, Matrix4x4.zero);
                }
            }

            Real minimal_cost = 0;
            var found = solveConstrainedMinimizationInequalities(qem, constraints, direction_, ref new_point, ref minimal_cost);
            if (found)
            {
                foreach (var c in constraints)
                {
                    if (minimal_cost < Single.MaxValue && new_point == Vector3.zero)
                    {
                        if (Vector3.Distance(c.p1, new_point) < minimal_cost)
                        {
                            Debug.Log("cost error. 1");
                        }
                        if (Vector3.Distance(c.p2, new_point) < minimal_cost)
                        {
                            Debug.Log("cost error. 2");
                        }
                        if (Vector3.Distance(c.p3, new_point) < minimal_cost)
                        {
                            Debug.Log("cost error. 3");
                        }
                    }
                }


                return new EdgeContraction(edge_index, new_point, minimal_cost, qem);
            }
            else
            {
                return new EdgeContraction(edge_index, Vector3.zero, Real.MaxValue, Matrix4x4.zero);
            }
        }
        bool solveConstrainedMinimizationInequalities(Matrix4x4 qem, List<Plane> constraints, DecimationDirection direction,
            ref Vector3 result, ref Real result_cost)
        {
            var found = false;
            var mini_cost = Real.MaxValue;
            int used_m = -1;
            var new_point = Vector3.zero;

            for (int m = 3; m >= 0; m--)
            {
                int n_subsets = nSubsets(m, constraints.Count);
                List<int> subset = new List<int>(m);
                for (int i = 0; i < m; i++)
                {
                    subset.Add(i);
                }
                for (int i = 0; i < n_subsets; i++)
                {
                    Vector3 mini_result = Vector3.zero;
                    bool found_result = solveConstrainedMinimization(qem, constraints, subset, direction_, ref mini_result);
                    var result_homogeneous = new Vector4(mini_result.x, mini_result.y, mini_result.z, 1);
                    if (found_result)
                    {
                        found = true;
                        var new_cost = MathUtil.VtmV(result_homogeneous, qem);
                        if (new_cost < mini_cost)
                        {
                            mini_cost = new_cost;
                            new_point = mini_result;
                            used_m = m;
                        }
                    }
                    nextSubset(ref subset, constraints.Count);
                }
                if (found)
                {
                    break;
                }
            }
            if (mini_cost < -Mathf.Epsilon)
            {
                mini_cost = Real.MaxValue;
            }

            result_cost = mini_cost;
            result = new_point;

            return found;
        }
        void nextSubset(ref List<int> indices_subset, int total_size)
        {
            for (int i = indices_subset.Count - 1; i >= 0; i--)
            {
                var max_next_index = indices_subset[i] + 1 + (indices_subset.Count - 1 - i);
                if (max_next_index < total_size)
                {
                    indices_subset[i] = indices_subset[i] + 1;
                    for (int j = i + 1; j < indices_subset.Count; j++)
                    {
                        indices_subset[j] = indices_subset[j - 1] + 1;
                    }
                    break;
                }
            }
        }
        int nSubsets(int subset_size, int total_size)
        {
            int result = 1;
            for (int i = 0; i < subset_size; i++)
            {
                result *= (total_size - i) / (i + 1);
            }
            return result;
        }
        bool solveConstrainedMinimization(Matrix4x4 qem, List<Plane> constrains, List<int> subset,
            DecimationDirection direction, ref Vector3 result)
        {
            Vector3 result_position = Vector3.zero;
            switch (subset.Count)
            {
                case 0:
                    result_position = minimizeSubspace(qem);
                    break;
                case 1:
                    result_position = minimizeSubspace(qem, constrains[subset[0]]);
                    break;
                case 2:
                    result_position = minimizeSubspace(qem, constrains[subset[0]], constrains[subset[1]]);
                    break;
                case 3:
                    result_position = minimizeSubspace(qem, constrains[subset[0]], constrains[subset[1]], constrains[subset[2]]);
                    break;
                default:
                    throw new System.Exception("undefined.");                    
            }
            //if (result_position == Vector3.zero)
            //{
            //    return false;
            //}

            Vector4 result_homogeneous = new Vector4(result_position.x, result_position.y, result_position.z, 1);
            bool valid_resultion = true;
            foreach (var c in constrains)
            {
                if (direction == DecimationDirection.Outward)
                {
                    if (Vector4.Dot(c.vector(), result_homogeneous) + Mathf.Epsilon < 0)
                    {
                        valid_resultion = false;
                        break;
                    }
                }
                else if (direction == DecimationDirection.Inward)
                {
                    if (Vector4.Dot(c.vector(), result_homogeneous) - Mathf.Epsilon > 0)
                    {
                        valid_resultion = false;
                        break;
                    }
                }
            }
            result = result_position;

            return valid_resultion;
        }
        Vector3 minimizeSubspace(Matrix4x4 quadratic_cost)
        {
            var E = Matrix4x4.identity;
            var f = new Vector4(0, 0, 0, 1);
            var A = E.transpose * quadratic_cost * E;
            var b = -1 * (E.transpose * quadratic_cost * f);
            var result = Cholesky(A).inverse * b;

            return new Vector3(result.x, result.y, result.z);
        }
        Matrix4x4 Cholesky(Matrix4x4 m)
        {

            float[,] m1 = 
                {
                {m.m00, m.m01, m.m02, m.m03 },
                {m.m10, m.m11, m.m12, m.m13 },
                {m.m20, m.m21, m.m22, m.m23 },
                {m.m30, m.m31, m.m32, m.m33 }
                };
            var m2 = Cholesky(m1);

            var ret = new Matrix4x4();
            ret.m00 = m2[0, 0];
            ret.m01 = m2[0, 1];
            ret.m02 = m2[0, 2];
            ret.m03 = m2[0, 3];
            ret.m10 = m2[1, 0];
            ret.m11 = m2[1, 1];
            ret.m12 = m2[1, 2];
            ret.m13 = m2[1, 3];
            ret.m20 = m2[2, 0];
            ret.m21 = m2[2, 1];
            ret.m22 = m2[2, 2];
            ret.m23 = m2[2, 3];
            ret.m30 = m2[3, 0];
            ret.m31 = m2[3, 1];
            ret.m32 = m2[3, 2];
            ret.m33 = m2[3, 3];

            return ret;
        }
        public static float[,] Cholesky(float[,] a)
        {
            int n = (int)Math.Sqrt(a.Length);

            float[,] ret = new float[n, n];
            for (int r = 0; r < n; r++)
                for (int c = 0; c <= r; c++)
                {
                    if (c == r)
                    {
                        float sum = 0;
                        for (int j = 0; j < c; j++)
                        {
                            sum += ret[c, j] * ret[c, j];
                        }
                        ret[c, c] = Mathf.Sqrt(a[c, c] - sum);
                    }
                    else
                    {
                        float sum = 0;
                        for (int j = 0; j < c; j++)
                            sum += ret[r, j] * ret[c, j];
                        ret[r, c] = 1.0f / ret[c, c] * (a[r, c] - sum);
                    }
                }

            return ret;
        }
        Vector3 minimizeSubspace(Matrix4x4 quadratic_cost, Plane plane)
        {
            var some_unit = new Vector3(1, 0, 0);
            var some_other_unit = new Vector3(0, 1, 0);
            var plane_direction_1 = Vector3.zero;
            var plane_direction_2 = Vector3.zero;

            if (Vector3.Dot(plane.normal, some_unit) < Vector3.Dot(plane.normal, some_other_unit))            
            //if(Utils.Norm(Utils.VvMatrix(plane.normal,  some_unit)) < Utils.Norm(Utils.VvMatrix(plane.normal, some_other_unit)))
            {
                plane_direction_1 = Vector3.Cross(plane.normal, some_unit);
            }
            else
            {
                plane_direction_1 = Vector3.Cross(plane.normal, some_other_unit);
            }

            plane_direction_1.Normalize();
            plane_direction_2 = Vector3.Cross(plane.normal, plane_direction_1);
            plane_direction_2.Normalize();

            var E = new Matrix4x4(plane_direction_1, plane_direction_2, new Vector4(0, 0, 0, 0), new Vector4(0, 0, 0, 0));
            var v3 = -1 * plane.d * plane.normal;
            var f = new Vector4(v3.x, v3.y, v3.z, 1);

            var m = (E.transpose * quadratic_cost * E).inverse * (-1 * (E.transpose * quadratic_cost * f));
            m.z = 0; m.w = 0;
            var minimizer = E * m + f;

            return new Vector3(minimizer.x, minimizer.y, minimizer.z);
        }
        Vector3 minimizeSubspace(Matrix4x4 quadratic_cost, Plane plane1, Plane plane2)
        {
            var point_on_p = plane1.normal * (-plane1.d);
            var along_p_to_edge = Vector3.Cross(plane1.normal, Vector3.Cross(plane1.normal, plane2.normal));

            var mu = (-plane2.d - Vector3.Dot(plane2.normal, point_on_p)) / (Vector3.Dot(plane2.normal, along_p_to_edge));
            var point_on_edge = point_on_p + mu * along_p_to_edge;

            var a1 = Vector3.Dot(point_on_edge, plane1.normal);
            var a2 = Vector3.Dot(point_on_edge, plane2.normal);

            if (Mathf.Abs(a1 + plane1.d) > Mathf.Epsilon
                || Mathf.Abs(a2 + plane2.d) > Mathf.Epsilon)
            {
                return Vector3.zero;
            }

            var ep1 = Vector3.Cross(plane1.normal, plane2.normal);
            var E = new Vector4(ep1.x, ep1.y, ep1.z, 0);
            var f = new Vector4(point_on_edge.x, point_on_edge.y, point_on_edge.z, 1);

            Real m = -1 * (MathUtil.VtmV(E, quadratic_cost, f)) / (MathUtil.VtmV(E, quadratic_cost, E));
            var minimizer = E * m + f;

            return new Vector3(minimizer.x, minimizer.y, minimizer.z);
        }

        // 三个平面相交时，取交点
        Vector3 minimizeSubspace(Matrix4x4 quadratic_cost, Plane plane1, Plane plane2, Plane plane3)
        {
#if true
            Matrix<Real> m1 = MathUtil.NewMatrix3x4(plane1.normal, plane2.normal, plane3.normal);
            Vector<Real> b1 = Vector<Real>.Build.Dense(new Real[] { -plane1.d, -plane2.d, -plane3.d });
            m1 = m1.Inverse();
            var n = m1.FrobeniusNorm();
            var b2 = m1.Multiply(b1);
            if (n > 1e3)
            {
                return Vector3.zero;
            }
            return MathUtil.ToVector3(b2);
#else

            var m = new Matrix4x4(plane1.normal, plane2.normal, plane3.normal, new Vector4(0,0,0,1));
            var b = new Vector3(-plane1.d, -plane2.d, -plane3.d);

            m = m.inverse;

            var n = Utils.Norm(m);
            if (n > 1e3)
            {
                return Vector3.zero;
            }

            return m * b;
#endif
        }



        Mesh result_mesh_ = null;
        MetricGenerator metric_generator_ = new MetricGenerator();
        ContractionQueue queue_ = new ContractionQueue();
        DecimationDirection direction_ = DecimationDirection.Outward;
        bool target_vertices_used_ = false;
        bool maximum_error_used_ = false;
        int target_vertices = 0;
        Real maximum_error_ = 0;
        Real current_error_ = 0;
    }
}
