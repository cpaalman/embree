// ======================================================================== //
// Copyright 2009-2013 Intel Corporation                                    //
//                                                                          //
// Licensed under the Apache License, Version 2.0 (the "License");          //
// you may not use this file except in compliance with the License.         //
// You may obtain a copy of the License at                                  //
//                                                                          //
//     http://www.apache.org/licenses/LICENSE-2.0                           //
//                                                                          //
// Unless required by applicable law or agreed to in writing, software      //
// distributed under the License is distributed on an "AS IS" BASIS,        //
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. //
// See the License for the specific language governing permissions and      //
// limitations under the License.                                           //
// ======================================================================== //

#pragma once

#include "builders/primrefalloc.h"
#include "builders/primrefblock.h"
#include "builders/primrefgen.h"
#include "builders/splitter.h"
#include "builders/splitter_parallel.h"
#include "geometry/primitive.h"

#define logBlockSize 2

namespace embree
{
  /* BVH builder. The builder is multi-threaded and implements 3
   * different build strategies: 1) Small tasks are finished in a
   * single thread (BuildTask) 2) Medium sized tasks are split into
   * two tasks using a single thread (SplitTask) and 3) Large tasks are
   * split using multiple threads on one processor. */

  class BVH4Builder2 : public Builder
  {
    ALIGNED_CLASS;
  public:

    /*! Type shortcuts */
    typedef typename BVH4::Node    Node;
    typedef typename BVH4::NodeRef NodeRef;

    /*! Split type of the split heuristic. */
    //typedef typename Heuristic::Split Split;
    //typedef typename Heuristic::PrimInfo PrimInfo;
    //typedef Splitter<Heuristic> SplitterNormal;

    class ObjectSplitBinner
  {
    /*! Maximal number of bins. */
    static const size_t maxBins = 32;

  public:

    static const std::string name() { return "objectsplit"; }
    
    /*! stores bounding information for a set of primitives */
    class PrimInfo
    {
    public:
      __forceinline PrimInfo () 
        : num(0), geomBounds(empty), centBounds(empty) {}

      /*! returns the number of primitives */
      __forceinline size_t size() const { 
        return num; 
      }

      /*! return the surface area heuristic when creating a leaf */
      __forceinline float sah () const { 
        return halfArea(geomBounds)*blocks(num); 
      }

    public:
      size_t num;          //!< number of primitives
      BBox3fa geomBounds;   //!< geometry bounds of primitives
      BBox3fa centBounds;   //!< centroid bounds of primitives
    };
    
    /*! mapping from bounding boxes to bins */
    class Mapping
    {
    public:
      
      /*! default constructor */
      __forceinline Mapping () {}
      
      /*! construct from primitive info */
      __forceinline Mapping (const PrimInfo& pinfo)
      {
        num   = min(maxBins,size_t(4.0f + 0.05f*pinfo.num));
        ofs   = pinfo.centBounds.lower;
        scale = rcp(pinfo.centBounds.upper - pinfo.centBounds.lower) * Vec3fa(float(num));
      }
      
      /*! returns number of bins */
      __forceinline size_t size() const { return num; }
      
      /*! Computes the bin numbers for each dimension for a box. */
      __forceinline Vec3ia bin_unsafe(const BBox3fa& box) const {
        return Vec3ia((center2(box) - ofs)*scale-Vec3fa(0.5f));
      }
      
      /*! Computes the bin numbers for each dimension for a box. */
      __forceinline Vec3ia bin(const BBox3fa& box) const {
        return clamp(bin_unsafe(box),Vec3ia(0),Vec3ia(int(num-1)));
      }
      
    private:
      size_t num;       //!< number of bins to use
      Vec3fa ofs;        //!< offset to compute bin
      Vec3fa scale;      //!< scaling factor to compute bin
    };
    
    /*! Stores information about an object split. */
    class Split
    {
    public:
      
      /*! create an invalid split by default */
      __forceinline Split () : dim(0), pos(0), cost(inf) {}
      
      /*! return SAH cost of performing the split */
      __forceinline float sah() const { return cost; } 

      /*! splitting of the primitive if required */
      __forceinline void split(const PrimRef& prim, PrimRef& lprim_o, PrimRef& rprim_o) const {
        new (&lprim_o) PrimRef(empty,-1,-1);
        new (&rprim_o) PrimRef(empty,-1,-1);
        assert(false);
      }
      
      void split(size_t threadIndex, PrimRefAlloc* alloc, atomic_set<PrimRefBlock>& prims, atomic_set<PrimRefBlock>& lprims, atomic_set<PrimRefBlock>& rprims);

    public:
      PrimInfo pinfo;
      Mapping mapping;    //!< Mapping to bins
      int dim;            //!< Best object split dimension
      int pos;            //!< Best object split position
      float cost;         //!< SAH cost of performing best object split
    };
    
    /*! default constructor */
    ObjectSplitBinner (atomic_set<PrimRefBlock>& prims);
  
  private:

    void add(atomic_set<PrimRefBlock>& prims);
    
    void bin(atomic_set<PrimRefBlock>& prims);

    /*! calculate the best possible split */
    void best(Split& split_o);

    void info(const PrimRef* prims, size_t num);

    /*! bin an array of primitives */
    void bin(const PrimRef* prim, size_t num);
    
    
    /*! Compute the number of blocks occupied for each dimension. */
    __forceinline static Vec3ia blocks(const Vec3ia& a) { return (a+Vec3ia((1 << logBlockSize)-1)) >> logBlockSize; }
    
    /*! Compute the number of blocks occupied in one dimension. */
    __forceinline static int  blocks(size_t a) { return (int)((a+((1LL << logBlockSize)-1)) >> logBlockSize); }
  
  public:
    PrimInfo pinfo;                //!< bounding information of geometry
    Mapping mapping;               //!< mapping from geometry to the bins
    
    /* initialize binning counter and bounds */
    Vec3ia   counts    [maxBins];    //< number of primitives mapped to bin
    BBox3fa geomBounds[maxBins][4]; //< bounds for every bin in every dimension

    Split split;
  };

  
  public:

    /*! builder entry point */
    void build(size_t threadIndex, size_t threadCount);

    /*! Constructor. */
    BVH4Builder2 (BVH4* bvh, BuildSource* source, void* geometry, const size_t minLeafSize = 1, const size_t maxLeafSize = inf);

    /*! creates a leaf node */
    NodeRef createLeaf(size_t threadIndex, atomic_set<PrimRefBlock>& prims, const ObjectSplitBinner::Split& split);

    NodeRef recurse(size_t threadIndex, size_t depth, atomic_set<PrimRefBlock>& prims, const ObjectSplitBinner::Split& split);

  private:
    BuildSource* source;      //!< build source interface
    void* geometry;           //!< input geometry
    
  public:
    const PrimitiveType& primTy;          //!< triangle type stored in BVH4
    size_t minLeafSize;                 //!< minimal size of a leaf
    size_t maxLeafSize;                 //!< maximal size of a leaf
    PrimRefAlloc alloc;                 //!< Allocator for primitive blocks
    TaskScheduler::QUEUE taskQueue;     //!< Task queue to use

  public:
    BVH4* bvh;                      //!< Output BVH4
  };
}
