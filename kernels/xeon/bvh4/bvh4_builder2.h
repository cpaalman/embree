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
#include "geometry/bezier1.h"

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

    typedef PrimRefBlockT<PrimRef> TriRefBlock;
    typedef atomic_set<TriRefBlock> TriRefList;

    typedef PrimRefBlockT<Bezier1> BezierRefBlock;
    typedef atomic_set<BezierRefBlock> BezierRefList;

    /*! stores bounding information for a set of primitives */
    class PrimInfo
    {
    public:
      __forceinline PrimInfo () 
        : num(0), numTriangles(0), numBeziers(0), geomBounds(empty), centBounds(empty) {}

      __forceinline PrimInfo (size_t numTriangles, size_t numBeziers, const BBox3fa& geomBounds, const BBox3fa& centBounds) 
        : num(numTriangles+numBeziers), numTriangles(numTriangles), numBeziers(numBeziers), geomBounds(geomBounds), centBounds(centBounds) {}

      /*! returns the number of primitives */
      __forceinline size_t size() const { 
        return num; 
      }

      __forceinline float leafSAH(float triCost, float bezierCost) const { 
        return halfArea(geomBounds)*(triCost*blocks(numTriangles) + bezierCost*numBeziers); 
      }

    public:
      size_t num;          //!< number of primitives
      size_t numTriangles, numBeziers;
      BBox3fa geomBounds;   //!< geometry bounds of primitives
      BBox3fa centBounds;   //!< centroid bounds of primitives
    };

    class ObjectSplitBinner
  {
    /*! Maximal number of bins. */
    static const size_t maxBins = 32;

  public:
    
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
      __forceinline float splitSAH() const { return cost; } 

      void split(size_t threadIndex, PrimRefBlockAlloc<PrimRef>* alloc, TriRefList& prims, TriRefList& lprims, TriRefList& rprims);

      void split(size_t threadIndex, PrimRefBlockAlloc<Bezier1>* alloc, BezierRefList& prims, BezierRefList& lprims, BezierRefList& rprims);

    public:
      Mapping mapping;    //!< Mapping to bins
      int dim;            //!< Best object split dimension
      int pos;            //!< Best object split position
      float cost;         //!< SAH cost of performing best object split
      LinearSpace3fa space;
    };
    
    /*! default constructor */
    ObjectSplitBinner (TriRefList& prims, float triCost);
    ObjectSplitBinner (TriRefList& prims, float triCost, BezierRefList& beziers, float bezierCost);
  
  private:

    void add(TriRefList& prims);
    void add(BezierRefList& prims);

    void setup_binning();
    
    void bin(TriRefList& prims);

    void bin(BezierRefList& prims);

    /*! calculate the best possible split */
    void best();

    void info(const PrimRef* prims, size_t num);

    void info(const Bezier1* prims, size_t num);

    /*! bin an array of primitives */
    void bin(const PrimRef* prim, size_t num);

    /*! bin an array of beziers */
    void bin(const Bezier1* prim, size_t num);
      
  public:
    PrimInfo pinfo;                //!< bounding information of geometry
    Mapping mapping;               //!< mapping from geometry to the bins
    
    /* initialize binning counter and bounds */
    Vec3ia  triCounts[maxBins];    
    Vec3ia  bezierCounts[maxBins];    
    BBox3fa geomBounds[maxBins][4]; 

    Split split;
    float triCost;
    float bezierCost;
  };

    class ObjectSplitBinnerUnaligned
  {
    /*! Maximal number of bins. */
    static const size_t maxBins = 32;

  public:
    
    typedef ObjectSplitBinner::Mapping Mapping;

    /*! Stores information about an object split. */
    class Split
    {
    public:
      
      /*! create an invalid split by default */
      __forceinline Split () : dim(0), pos(0), cost(inf) {}
      
      /*! return SAH cost of performing the split */
      __forceinline float splitSAH() const { return cost; } 

      void split(size_t threadIndex, PrimRefBlockAlloc<PrimRef>* alloc, TriRefList& prims, TriRefList& lprims, TriRefList& rprims);

      void split(size_t threadIndex, PrimRefBlockAlloc<Bezier1>* alloc, BezierRefList& prims, BezierRefList& lprims, BezierRefList& rprims);

    public:
      Mapping mapping;    //!< Mapping to bins
      int dim;            //!< Best object split dimension
      int pos;            //!< Best object split position
      float cost;         //!< SAH cost of performing best object split
      LinearSpace3fa space;
    };
    
    /*! default constructor */
    ObjectSplitBinnerUnaligned (const LinearSpace3fa& space, TriRefList& prims, float triCost);
    ObjectSplitBinnerUnaligned (const LinearSpace3fa& space, TriRefList& prims, float triCost, BezierRefList& beziers, float bezierCost);
  
  private:

    void add(TriRefList& prims);
    void add(BezierRefList& prims);

    void setup_binning();
    
    void bin(TriRefList& prims);

    void bin(BezierRefList& prims);

    /*! calculate the best possible split */
    void best();

    void info(const PrimRef* prims, size_t num);

    void info(const Bezier1* prims, size_t num);

    /*! bin an array of primitives */
    void bin(const PrimRef* prim, size_t num);

    /*! bin an array of beziers */
    void bin(const Bezier1* prim, size_t num);
      
  public:
    PrimInfo pinfo;                //!< bounding information of geometry
    Mapping mapping;               //!< mapping from geometry to the bins
    LinearSpace3fa space;
    
    /* initialize binning counter and bounds */
    Vec3ia  triCounts[maxBins];    
    Vec3ia  bezierCounts[maxBins];    
    BBox3fa geomBounds[maxBins][4]; 

    Split split;
    float triCost;
    float bezierCost;
  };

    struct SpatialSplit
    {
      /*! Maximal number of bins. */
      static const size_t BINS = 16;

    public:

      SpatialSplit (TriRefList& tris, float triCost, BezierRefList& beziers, float bezierCost);

      class Split
      {
      public:
        
        /*! create an invalid split by default */
        __forceinline Split () : dim(0), pos(0), ipos(0), cost(inf), numSpatialSplits(0) {}
        
        /*! return SAH cost of performing the split */
        __forceinline float splitSAH() const { return cost; } 
        
        void split(size_t threadIndex, PrimRefBlockAlloc<PrimRef>* alloc, TriRefList& prims, TriRefList& lprims, TriRefList& rprims) const;
        
        void split(size_t threadIndex, PrimRefBlockAlloc<Bezier1>* alloc, BezierRefList& prims, BezierRefList& lprims, BezierRefList& rprims) const;
      
      public:
        //LinearSpace3fa space;
        float pos;
        int ipos;
        int dim;
        float cost;
        ssef ofs,scale;
        size_t numSpatialSplits;
      };

      void bin(TriRefList& tris);
      void bin(BezierRefList& beziers);
      void best();

      Split split;

    BBox3fa geomBounds;
      ssef ofs;
      ssef diag;
      ssef scale;

    BBox3fa bounds[BINS][4];
    Vec3ia    numTriBegin[BINS];
    Vec3ia    numTriEnd[BINS];
    Vec3ia    numBezierBegin[BINS];
    Vec3ia    numBezierEnd[BINS];

      float triCost;
      float bezierCost;
    };

  class ObjectTypePartitioning
  {
  public:
    ObjectTypePartitioning (TriRefList& prims, float triCost, BezierRefList& beziers, float bezierCost);

  public:
    class Split 
    { 
    public:
      /*! return SAH cost of performing the split */
      __forceinline float splitSAH() const { return cost; } 

      void split(size_t threadIndex, PrimRefBlockAlloc<PrimRef>* alloc, TriRefList& prims, TriRefList& lprims, TriRefList& rprims);
      void split(size_t threadIndex, PrimRefBlockAlloc<Bezier1>* alloc, BezierRefList& prims, BezierRefList& lprims, BezierRefList& rprims);

    public:
      float cost;
    };

  public:
    Split split;
    float triCost;
    float bezierCost;
  };
    
    class __aligned(16) GeneralSplit
    {
      enum Type { OBJECT_SPLIT, OBJECT_SPLIT_UNALIGNED, SPATIAL_SPLIT, TYPE_SPLIT, FALLBACK_SPLIT };

    public:

      __forceinline GeneralSplit () {}

      __forceinline GeneralSplit (size_t N) 
        : type(FALLBACK_SPLIT), aligned(true), split_sah(inf), num(N) {}

      __forceinline GeneralSplit(const ObjectSplitBinner::Split& split, bool aligned_in) {
        type = OBJECT_SPLIT; aligned = aligned_in;
        split_sah = split.splitSAH();
        new (data) ObjectSplitBinner::Split(split);
      }

      __forceinline GeneralSplit(const ObjectSplitBinnerUnaligned::Split& split) {
        type = OBJECT_SPLIT_UNALIGNED; aligned = false;
        split_sah = split.splitSAH();
        new (data) ObjectSplitBinnerUnaligned::Split(split);
      }

      __forceinline GeneralSplit(const SpatialSplit::Split& split, bool aligned_in) {
        type = SPATIAL_SPLIT; aligned = aligned_in;
        split_sah = split.splitSAH();
        new (data) SpatialSplit::Split(split);
      }

      __forceinline GeneralSplit(const ObjectTypePartitioning::Split& split) {
        type = TYPE_SPLIT; aligned = true;
        split_sah = split.splitSAH();
        new (data) ObjectTypePartitioning::Split(split);
      }

      __forceinline void split(size_t threadIndex, PrimRefBlockAlloc<PrimRef>* alloc, TriRefList& prims, TriRefList& lprims, TriRefList& rprims) 
      {
        switch (type) {
        case OBJECT_SPLIT : ((ObjectSplitBinner::Split* )data)->split(threadIndex,alloc,prims,lprims,rprims); break;          
        case OBJECT_SPLIT_UNALIGNED : ((ObjectSplitBinnerUnaligned::Split* )data)->split(threadIndex,alloc,prims,lprims,rprims); break;          
        case SPATIAL_SPLIT: ((SpatialSplit::Split*)data)->split(threadIndex,alloc,prims,lprims,rprims); break;          
        case TYPE_SPLIT   : ((ObjectTypePartitioning::Split*   )data)->split(threadIndex,alloc,prims,lprims,rprims); break;         
        default           : split_fallback(threadIndex,alloc,prims,lprims,rprims); break;
        }
      }

      void split(size_t threadIndex, PrimRefBlockAlloc<Bezier1>* alloc, BezierRefList& prims, BezierRefList& lprims, BezierRefList& rprims)
      {
        switch (type) {
        case OBJECT_SPLIT : ((ObjectSplitBinner::Split* )data)->split(threadIndex,alloc,prims,lprims,rprims); break;          
        case OBJECT_SPLIT_UNALIGNED : ((ObjectSplitBinnerUnaligned::Split* )data)->split(threadIndex,alloc,prims,lprims,rprims); break;          
        case SPATIAL_SPLIT: ((SpatialSplit::Split*)data)->split(threadIndex,alloc,prims,lprims,rprims); break;          
        case TYPE_SPLIT   : ((ObjectTypePartitioning::Split*   )data)->split(threadIndex,alloc,prims,lprims,rprims); break;         
        default           : split_fallback(threadIndex,alloc,prims,lprims,rprims); break;
        }
      }

      __forceinline float splitSAH() const { return split_sah; }
      
      bool aligned;
    private:
      Type type;
      float split_sah;
      size_t num;
      __aligned(16) char data[256];
    };
  
    /*! stores all info to build a subtree */
    struct BuildTask
    {
      __forceinline BuildTask () {}

      __forceinline BuildTask (BVH4::NodeRef* dst, size_t depth, TriRefList& tris, BezierRefList& beziers, const PrimInfo& pinfo, const GeneralSplit& split, const NAABBox3fa& bounds)
        : dst(dst), depth(depth), tris(tris), beziers(beziers), pinfo(pinfo), split(split), nodeBounds(bounds) {}

    public:
      __forceinline friend bool operator< (const BuildTask& a, const BuildTask& b) {
        return area(a.pinfo.geomBounds) < area(b.pinfo.geomBounds);
      //  return a.pinfo.size() < b.pinfo.size();
        }

      /*! execute single task and create subtasks */
      void process(size_t threadIndex, BVH4Builder2* builder, BuildTask task_o[BVH4::N], size_t& N);

      /*! recursive build function for aligned and non-aligned bounds */
      void recurse(size_t threadIndex, BVH4Builder2* builder);
      
    public:
      BVH4::NodeRef* dst;
      size_t depth;
      TriRefList tris;
      BezierRefList beziers;
      PrimInfo pinfo;
      GeneralSplit split;
      NAABBox3fa nodeBounds;
    };

  public:

    static const PrimInfo computePrimInfo(BezierRefList& beziers);
    static const PrimInfo computePrimInfo(TriRefList& tris, BezierRefList& beziers);

    /*! calculate bounds for range of primitives */
    static const BBox3fa computeAlignedBounds(TriRefList& tris);
    static const BBox3fa computeAlignedBounds(BezierRefList& beziers);
    static const BBox3fa computeAlignedBounds(TriRefList& tris, BezierRefList& beziers);
    
    /*! calculate bounds for range of primitives */
    static const NAABBox3fa computeAlignedBounds(TriRefList& tris, const LinearSpace3fa& space);
    static const NAABBox3fa computeAlignedBounds(BezierRefList& beziers, const LinearSpace3fa& space);
    static const NAABBox3fa computeAlignedBounds(TriRefList& tris, BezierRefList& beziers, const LinearSpace3fa& space);

    static const LinearSpace3fa computeHairSpace(BezierRefList& prims);

    static void split_fallback(size_t threadIndex, PrimRefBlockAlloc<PrimRef>* alloc, TriRefList& prims, TriRefList& lprims, TriRefList& rprims);
    static void split_fallback(size_t threadIndex, PrimRefBlockAlloc<Bezier1>* alloc, BezierRefList& prims, BezierRefList& lprims, BezierRefList& rprims);

    /*! builder entry point */
    void build(size_t threadIndex, size_t threadCount);

    /*! Constructor. */
    BVH4Builder2 (BVH4* bvh, BuildSource* source, void* geometry, const size_t minLeafSize = 1, const size_t maxLeafSize = inf);

    /*! creates a leaf node */
    NodeRef leaf   (size_t threadIndex, size_t depth, TriRefList& prims   , const PrimInfo& pinfo);
    NodeRef leaf   (size_t threadIndex, size_t depth, BezierRefList& prims, const PrimInfo& pinfo);
    NodeRef leaf   (size_t threadIndex, size_t depth, TriRefList& tris, BezierRefList& beziers, const PrimInfo& pinfo);

    void heuristic(TriRefList& tris, BezierRefList& beziers, GeneralSplit& split, const NAABBox3fa& nodeBounds);

    TASK_RUN_FUNCTION(BVH4Builder2,task_build_parallel);

  private:
    BuildSource* source;      //!< build source interface
    void* geometry;           //!< input geometry
    

  public:
    size_t minLeafSize;                 //!< minimal size of a leaf
    size_t maxLeafSize;                 //!< maximal size of a leaf
    PrimRefBlockAlloc<Bezier1> allocBezierRefs;                 //!< Allocator for primitive blocks
    PrimRefBlockAlloc<PrimRef> allocTriRefs;                 //!< Allocator for primitive blocks

    /*! Compute the number of blocks occupied for each dimension. */
    __forceinline static Vec3ia blocks(const Vec3ia& a) { return (a+Vec3ia((1 << logBlockSize)-1)) >> logBlockSize; }
    
    /*! Compute the number of blocks occupied in one dimension. */
    __forceinline static int  blocks(size_t a) { return (int)((a+((1LL << logBlockSize)-1)) >> logBlockSize); }

  public:
    BVH4* bvh;                      //!< Output BVH4

    MutexSys taskMutex;
    volatile atomic_t numActiveTasks;
    volatile atomic_t remainingSpatialSplits;
    std::vector<BuildTask> tasks;

  };
}
