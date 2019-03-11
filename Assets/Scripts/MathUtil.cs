using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using boundingmesh;
using MathNet.Numerics.LinearAlgebra;

namespace boundingmesh
{
    using Real = System.Single; //System.Double

    public static class MathUtil
    {
        public static float VtmV(Vector4 v, Matrix4x4 m)
        {
            // vt * m * v;
            return VtmV(v, m, v);
        }
        public static float VtmV(Vector4 v, Matrix4x4 m, Vector4 v2)
        {
            // vt * m * v;
            return Vector4.Dot(m.transpose * v, v2);
        }
        public static Vector3 ToVector3(Vector<Real> v)
        {
            return new Vector3(v.At(0), v.At(1), v.At(2));
        }
        public static Matrix<Real> NewMatrix3x4(Vector3 v1, Vector3 v2, Vector3 v3)
        {
            Matrix < Real > m = Matrix<Real>.Build.Dense(3, 3);
            m[0, 0] = v1.x; m[0, 1] = v1.y; m[0, 2] = v1.z;
            m[1, 0] = v2.x; m[1, 1] = v2.y; m[1, 2] = v2.z;
            m[2, 0] = v3.x; m[2, 1] = v3.y; m[2, 2] = v3.z;

            return m;
        }
    }
}
