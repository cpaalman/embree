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

    /*! Block of build primitives */
    template<typename PrimRef>
    class PrimRefBlock
    {
    public:

      typedef PrimRef T;

      /*! Number of primitive references inside a block */
      static const size_t blockSize = 511;

      /*! default constructor */
      PrimRefBlock () : num(0) {}

      /*! frees the block */
      __forceinline void clear(size_t n = 0) { num = n; }
      
      /*! return base pointer */
      __forceinline PrimRef* base() { return ptr; }

      /*! returns number of elements */
      __forceinline size_t size() const { return num; }

      /*! inserts a primitive reference */
      __forceinline bool insert(const PrimRef& ref) {
        if (unlikely(num >= blockSize)) return false;
        ptr[num++] = ref;
        return true;
      }

      /*! access the i-th primitive reference */
      __forceinline       PrimRef& operator[] (size_t i)       { return ptr[i]; }
      __forceinline const PrimRef& operator[] (size_t i) const { return ptr[i]; }
    
      /*! access the i-th primitive reference */
      __forceinline       PrimRef& at (size_t i)       { return ptr[i]; }
      __forceinline const PrimRef& at (size_t i) const { return ptr[i]; }
    
    private:
      PrimRef ptr[blockSize];   //!< Block with primitive references
      size_t num;               //!< Number of primitive references in block
    };

    typedef PrimRefBlock<PrimRef> TriRefBlock;
    typedef atomic_set<TriRefBlock> TriRefList;

    typedef PrimRefBlock<Bezier1> BezierRefBlock;
    typedef atomic_set<BezierRefBlock> BezierRefList;

    template<typename PrimRef>
      class PrimRefBlockAlloc : public AllocatorBase
    {
      ALIGNED_CLASS;
    public:
   
      struct __aligned(4096) ThreadPrimBlockAllocator 
      {
        ALIGNED_CLASS_(4096);
      public:
      
        __forceinline typename atomic_set<PrimRefBlock<PrimRef> >::item* malloc(size_t thread, AllocatorBase* alloc) 
        {
          /* try to take a block from local list */
          atomic_set<PrimRefBlock<PrimRef> >::item* ptr = local_free_blocks.take_unsafe();
          if (ptr) return new (ptr) typename atomic_set<PrimRefBlock<PrimRef> >::item();
          
          /* if this failed again we have to allocate more memory */
          ptr = (typename atomic_set<PrimRefBlock<PrimRef> >::item*) alloc->malloc(sizeof(typename atomic_set<PrimRefBlock<PrimRef> >::item));
          
          /* return first block */
          return new (ptr) typename atomic_set<PrimRefBlock<PrimRef> >::item();
        }
      
        __forceinline void free(typename atomic_set<PrimRefBlock<PrimRef> >::item* ptr) {
          local_free_blocks.insert_unsafe(ptr);
        }
        
      public:
        atomic_set<PrimRefBlock<PrimRef> > local_free_blocks; //!< only accessed from one thread
      };
      
    public:
      
      /*! Allocator default construction. */
      PrimRefBlockAlloc () {
        threadPrimBlockAllocator = new ThreadPrimBlockAllocator[getNumberOfLogicalThreads()];
      }
      
      /*! Allocator destructor. */
      virtual ~PrimRefBlockAlloc() {
        delete[] threadPrimBlockAllocator; threadPrimBlockAllocator = NULL;
      }
      
      /*! Allocate a primitive block */
      __forceinline typename atomic_set<PrimRefBlock<PrimRef> >::item* malloc(size_t thread) {
        return threadPrimBlockAllocator[thread].malloc(thread,this);
      }
      
      /*! Frees a primitive block */
      __forceinline void free(size_t thread, typename atomic_set<PrimRefBlock<PrimRef> >::item* block) {
        return threadPrimBlockAllocator[thread].free(block);
      }
      
    private:
      ThreadPrimBlockAllocator* threadPrimBlockAllocator;  //!< Thread local allocator
    };

    class ObjectSplitBinner
  {
    /*! Maximal number of bins. */
    static const size_t maxBins = 32;

  public:

    /*! stores bounding information for a set of primitives */
    class PrimInfo
    {
    public:
      __forceinline PrimInfo () 
        : num(0), numTriangles(0), numBeziers(0), geomBounds(empty), centBounds(empty) {}

      /*! returns the number of primitives */
      __forceinline size_t size() const { 
        return num; 
      }

    public:
      size_t num;          //!< number of primitives
      size_t numTriangles, numBeziers;
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
      
      __forceinline float leafSAH() const { return halfArea(pinfo.geomBounds)*(blocks(pinfo.numTriangles) + pinfo.numBeziers); }

      /*! return SAH cost of performing the split */
      __forceinline float splitSAH() const { return cost; } 

      /*! splitting of the primitive if required */
      __forceinline void split(const PrimRef& prim, PrimRef& lprim_o, PrimRef& rprim_o) const {
        new (&lprim_o) PrimRef(empty,-1,-1);
        new (&rprim_o) PrimRef(empty,-1,-1);
        assert(false);
      }
      
      void split(size_t threadIndex, PrimRefBlockAlloc<PrimRef>* alloc, TriRefList& prims, TriRefList& lprims, TriRefList& rprims);

      void split(size_t threadIndex, PrimRefBlockAlloc<Bezier1>* alloc, BezierRefList& prims, BezierRefList& lprims, BezierRefList& rprims);

    public:
      PrimInfo pinfo;
      Mapping mapping;    //!< Mapping to bins
      int dim;            //!< Best object split dimension
      int pos;            //!< Best object split position
      float cost;         //!< SAH cost of performing best object split
    };
    
    /*! default constructor */
    ObjectSplitBinner (TriRefList& prims);
    ObjectSplitBinner (TriRefList& prims, BezierRefList& beziers);
  
  private:

    void add(TriRefList& prims);

    void add(BezierRefList& prims);

    void setup_binning();
    
    void bin(TriRefList& prims);

    void bin(BezierRefList& prims);

    /*! calculate the best possible split */
    void best(Split& split_o);

    void info(const PrimRef* prims, size_t num);

    void info(const Bezier1* prims, size_t num);

    /*! bin an array of primitives */
    void bin(const PrimRef* prim, size_t num);

    /*! bin an array of beziers */
    void bin(const Bezier1* prim, size_t num);
    
    /*! Compute the number of blocks occupied for each dimension. */
    __forceinline static Vec3ia blocks(const Vec3ia& a) { return (a+Vec3ia((1 << logBlockSize)-1)) >> logBlockSize; }
    
    /*! Compute the number of blocks occupied in one dimension. */
    __forceinline static int  blocks(size_t a) { return (int)((a+((1LL << logBlockSize)-1)) >> logBlockSize); }
  
  public:
    PrimInfo pinfo;                //!< bounding information of geometry
    Mapping mapping;               //!< mapping from geometry to the bins
    
    /* initialize binning counter and bounds */
    Vec3ia  triCounts[maxBins];    
    Vec3ia  bezierCounts[maxBins];    
    BBox3fa geomBounds[maxBins][4]; 

    Split split;
  };

  
  public:

    /*! calculate bounds for range of primitives */
    static const BBox3fa computeAlignedBounds(TriRefList& tris);
    static const BBox3fa computeAlignedBounds(BezierRefList& beziers);
    static const BBox3fa computeAlignedBounds(TriRefList& tris, BezierRefList& beziers);
    
    /*! calculate bounds for range of primitives */
    static const NAABBox3fa computeAlignedBounds(TriRefList& tris, const LinearSpace3fa& space);
    static const NAABBox3fa computeAlignedBounds(BezierRefList& beziers, const LinearSpace3fa& space);
    static const NAABBox3fa computeAlignedBounds(TriRefList& tris, BezierRefList& beziers, const LinearSpace3fa& space);


    /*! builder entry point */
    void build(size_t threadIndex, size_t threadCount);

    /*! Constructor. */
    BVH4Builder2 (BVH4* bvh, BuildSource* source, void* geometry, const size_t minLeafSize = 1, const size_t maxLeafSize = inf);

    /*! creates a leaf node */
    NodeRef leaf(size_t threadIndex, size_t depth, TriRefList& prims, const ObjectSplitBinner::Split& split);
    NodeRef leaf(size_t threadIndex, size_t depth, BezierRefList& prims, const ObjectSplitBinner::Split& split);
    NodeRef leaf(size_t threadIndex, size_t depth, TriRefList& tris, BezierRefList& beziers, const ObjectSplitBinner::Split& split);

    NodeRef recurse(size_t threadIndex, size_t depth, TriRefList& tris, BezierRefList& beziers, const ObjectSplitBinner::Split& split);

  private:
    BuildSource* source;      //!< build source interface
    void* geometry;           //!< input geometry
    
  public:
    const PrimitiveType& primTy;          //!< triangle type stored in BVH4
    size_t minLeafSize;                 //!< minimal size of a leaf
    size_t maxLeafSize;                 //!< maximal size of a leaf
    PrimRefBlockAlloc<Bezier1> allocBezierRefs;                 //!< Allocator for primitive blocks
    PrimRefBlockAlloc<PrimRef> allocTriRefs;                 //!< Allocator for primitive blocks
    //TaskScheduler::QUEUE taskQueue;     //!< Task queue to use

  public:
    BVH4* bvh;                      //!< Output BVH4
  };
}
