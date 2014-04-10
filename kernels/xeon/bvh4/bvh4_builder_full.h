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

#include "geometry/primitive.h"
#include "geometry/bezier1.h"

namespace embree
{
  class BVH4BuilderFull : public Builder
  {
    ALIGNED_CLASS;
  public:

    /*! builder entry point */
    void build(size_t threadIndex, size_t threadCount);

    /*! Constructor. */
    BVH4BuilderFull (BVH4* bvh, Scene* scene);

  private:

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

    typedef PrimRefBlock<Triangle1v> TriRefBlock;
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
      
        __forceinline atomic_set<PrimRefBlock<PrimRef> >::item* malloc(size_t thread, AllocatorBase* alloc) 
        {
          /* try to take a block from local list */
          atomic_set<PrimRefBlock<PrimRef> >::item* ptr = local_free_blocks.take_unsafe();
          if (ptr) return new (ptr) atomic_set<PrimRefBlock<PrimRef> >::item();
          
          /* if this failed again we have to allocate more memory */
          ptr = (atomic_set<PrimRefBlock<PrimRef> >::item*) alloc->malloc(sizeof(atomic_set<PrimRefBlock<PrimRef> >::item));
          
          /* return first block */
          return new (ptr) atomic_set<PrimRefBlock<PrimRef> >::item();
        }
      
        __forceinline void free(atomic_set<PrimRefBlock<PrimRef> >::item* ptr) {
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
      __forceinline atomic_set<PrimRefBlock<PrimRef> >::item* malloc(size_t thread) {
        return threadPrimBlockAllocator[thread].malloc(thread,this);
      }
      
      /*! Frees a primitive block */
      __forceinline void free(size_t thread, atomic_set<PrimRefBlock<PrimRef> >::item* block) {
        return threadPrimBlockAllocator[thread].free(block);
      }
      
    private:
      ThreadPrimBlockAllocator* threadPrimBlockAllocator;  //!< Thread local allocator
    };

    /*! stores all info to build a subtree */
    struct BuildTask
    {
      __forceinline BuildTask () {}

      __forceinline BuildTask (BVH4::NodeRef* dst, size_t depth, size_t size, bool makeleaf, TriRefList& tris, BezierRefList& beziers, const NAABBox3fa& bounds)
        : dst(dst), depth(depth), size(size), makeleaf(makeleaf), prims(prims), bounds(bounds) {}

    public:
      __forceinline friend bool operator< (const BuildTask& a, const BuildTask& b) {
        return area(a.bounds.bounds) < area(b.bounds.bounds);
        //return area(a.bounds.bounds)/double(a.size) < area(b.bounds.bounds)/double(b.size);
      }

    public:
      BVH4::NodeRef* dst;
      size_t depth;
      size_t size;
      bool makeleaf;
      TriRefList tris;
      BezierRefList beziers;
      NAABBox3fa bounds;
    };

    /*! Performs standard object binning */
    struct ObjectSplit
    {
      /*! number of bins */
      static const size_t BINS = 16;

    public:

      /*! default constructor */
      __forceinline ObjectSplit ()
        : dim(-1), pos(0), cost(inf), num0(0), num1(0), bounds0(inf), bounds1(inf) {}
      
      /*! calculates standard surface area heuristic for the split */
      __forceinline float standardSAH() const {
        return BVH4::intCost*float(num0)*embree::area(bounds0.bounds) + BVH4::intCost*float(num1)*embree::area(bounds1.bounds);
      }

      __forceinline bool operator() (const Triangle1v& prim) const
      {
        const Vec3fa center = prim.center(space);
        //const ssei bin = clamp(floori((ssef(center) - ofs)*scale),ssei(0),ssei(BINS-1));
        const ssei bin = floori((ssef(center)-ofs)*scale);
        return bin[dim] < pos;
      }

      __forceinline bool operator() (const Bezier1& prim) const
      {
        const Vec3fa center = prim.center(space);
        //const ssei bin = clamp(floori((ssef(center) - ofs)*scale),ssei(0),ssei(BINS-1));
        const ssei bin = floori((ssef(center)-ofs)*scale);
        return bin[dim] < pos;
      }

      /*! performs object binning to the the best partitioning */
      static ObjectSplit find(size_t threadIndex, size_t depth, BVH4BuilderFull* parent, TriRefList& tris, BezierRefList& beziers, const LinearSpace3fa& space);

      /*! splits hairs into two sets */
      void split(size_t threadIndex, BVH4BuilderFull* parent, 
                 TriRefList& tris,    BezierRefList& beziers, 
                 TriRefList& ltris_o, BezierRefList& lbeziers_o,
                 TriRefList& rtris_o, BezierRefList& rbeziers_o) const;
      
    public:
      LinearSpace3fa space;
      NAABBox3fa bounds0, bounds1;
      int dim;
      int pos;
      float cost;
      size_t num0,num1;
      ssef ofs,scale;
    };

    /*! Performs fallback splits */
    struct FallBackSplit
    {
      __forceinline FallBackSplit (const NAABBox3fa& bounds0, size_t num0, const NAABBox3fa& bounds1, size_t num1)
        : bounds0(bounds0), num0(num0), bounds1(bounds1), num1(num1) {}

      /*! finds some partitioning */
      static FallBackSplit find(size_t threadIndex, BVH4BuilderFull* parent, 
                                TriRefList& tris,    BezierRefList& beziers, 
                                TriRefList& ltris_o, BezierRefList& lbeziers_o,
                                TriRefList& rtris_o, BezierRefList& rbeziers_o);

    public:
      size_t num0, num1;
      NAABBox3fa bounds0, bounds1;
    };

  private:

    template<typename List>
      void insert(size_t threadIndex, List& prims_i, List& prims_o);

    template<typename List, typename Left>
      void split(size_t threadIndex, List& prims, const Left& left, 
                 List& lprims_o, size_t& lnum_o, List& rprims_o, size_t& rnum_o);

    /*! creates a leaf node */
    BVH4::NodeRef leaf(size_t threadIndex, size_t depth, atomic_set<PrimRefBlock>& prims, const NAABBox3fa& bounds);

    bool split(size_t threadIndex, size_t depth, 
               atomic_set<PrimRefBlock>& prims, const NAABBox3fa& bounds, size_t size,
               atomic_set<PrimRefBlock>& lprims, size_t& lsize,
               atomic_set<PrimRefBlock>& rprims, size_t& rsize,
               bool& isAligned);

    /*! execute single task and create subtasks */
    void processTask(size_t threadIndex, BuildTask& task, BuildTask task_o[BVH4::N], size_t& N);

    /*! recursive build function for aligned and non-aligned bounds */
    void recurseTask(size_t threadIndex, BuildTask& task);

    TASK_RUN_FUNCTION(BVH4BuilderFull,task_build_parallel);

  public:
    Scene* scene;          //!< source
    size_t minLeafSize;    //!< minimal size of a leaf
    size_t maxLeafSize;    //!< maximal size of a leaf

    size_t numAlignedObjectSplits;
    size_t numAlignedSpatialSplits;
    size_t numUnalignedObjectSplits;
    size_t numUnalignedSpatialSplits;
    size_t numStrandSplits;
    size_t numFallbackSplits;

    bool enableAlignedObjectSplits;
    bool enableAlignedSpatialSplits;
    bool enableUnalignedObjectSplits;
    bool enableUnalignedSpatialSplits;
    bool enableStrandSplits;
    int enablePreSubdivision;

    BVH4* bvh;         //!< output
    PrimRefBlockAlloc alloc;                 //!< Allocator for primitive blocks

    MutexSys taskMutex;
    volatile atomic_t numActiveTasks;
    volatile atomic_t numGeneratedPrims;
    volatile atomic_t remainingReplications;
    std::vector<BuildTask> tasks;
  };
}
