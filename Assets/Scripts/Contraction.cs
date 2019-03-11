using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using boundingmesh;

namespace boundingmesh
{
    using Real = System.Single; //System.Double
    using Index = System.Int32;

    class EdgeContractionCompare : IComparer<EdgeContraction>
    {
        public int Compare(EdgeContraction x, EdgeContraction y)
        {
            if (x.cost() < y.cost())
            {
                return -1;
            }
            else if (x.cost() > y.cost())
            {
                return 1;
            }
            else
            {
                return 0;
            }
        }
    }
    class EdgeContraction
    {
        public EdgeContraction()
        { }

        public EdgeContraction(Index edge, Vector3 new_point, Real cost, Matrix4x4 qem)
        {
            edge_ = edge;
            new_point_ = new_point;
            cost_ = cost;
            qem_ = qem;
        }

        public static bool operator<(EdgeContraction left, EdgeContraction right)
        {
            return left.cost_ < right.cost_;
        }
        public static bool operator>(EdgeContraction left, EdgeContraction right)
        {
            return right < left;
        }

        public Index edge()
        {
            return edge_;
        }
        public Vector3 new_point()
        {
            return new_point_;
        }
        public Real cost()
        {
            return cost_;
        }
        public Matrix4x4 qem()
        {
            return qem_;
        }

        Index edge_;
        Vector3 new_point_;
        Real cost_;
        Matrix4x4 qem_;
    }

    class ContractionIndex
    {
        public ContractionIndex(Index index)
        {
            index_ = index;
            searching_ = true;
            iterator = null;
        }
        public ContractionIndex(Index index, EdgeContraction c)
        {
            index_ = index;
            searching_ = false;
            iterator = c;
        }
        public Index index_;
        public EdgeContraction iterator;
        bool searching_ = false;
    }

    class ContractionIndexCompare : IComparer<ContractionIndex>
    {
        public int Compare(ContractionIndex x, ContractionIndex y)
        {
            if (x.index_ < y.index_)
            {
                return -1;
            }
            else if (x.index_ > y.index_)
            {
                return 1;
            }
            else
            {
                return 0;
            }
        }
    }

    class ContractionQueue
    {
        public int Size()
        {
            return contractions_.Count;
        }
        public EdgeContraction First()
        {
            return contractions_[0];
        }
        public void Insert(EdgeContraction c)
        {
            while (c.edge() >= indices_.Count)
            {
                indices_.Add(null);
            }

            diagnosis.Assert(indices_[c.edge()] == null);
            indices_[c.edge()] = c;

            contractions_.Add(c);

            {
                var comp = new EdgeContractionCompare();
                contractions_.Sort(comp);
            }
        }
        public void Remove(Index index)
        {
            EdgeContraction ci = null;
            {
                ci = indices_[index];
                if (ci == null)
                {
                    return;
                }
                indices_[index] = null;
            }

            {
                var comp = new EdgeContractionCompare();
                //var cj = contractions_.BinarySearch(ci, comp);
                //if (cj >= 0)
                {
                    //contractions_.RemoveAt(cj);
                    contractions_.Remove(ci);
                    contractions_.Sort(comp);
                }
            }
        }

        public List<EdgeContraction> contractions_ = new List<EdgeContraction>();
        public List<EdgeContraction> indices_ = new List<EdgeContraction>();
    }


}
