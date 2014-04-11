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

namespace embree
{
  /* BVH builder. The builder is multi-threaded and implements 3
   * different build strategies: 1) Small tasks are finished in a
   * single thread (BuildTask) 2) Medium sized tasks are split into
   * two tasks using a single thread (SplitTask) and 3) Large tasks are
   * split using multiple threads on one processor. */

  template<typename Heuristic>
    class BVH4Builder2 : public Builder
  {
    ALIGNED_CLASS;
  public:

    /*! Type shortcuts */
    typedef typename BVH4::Node    Node;
    typedef typename BVH4::NodeRef NodeRef;

    /*! Split type of the split heuristic. */
    typedef typename Heuristic::Split Split;
    typedef typename Heuristic::PrimInfo PrimInfo;
    
  public:

    /*! builder entry point */
    void build(size_t threadIndex, size_t threadCount);

    /*! Constructor. */
    BVH4Builder2 (BVH4* bvh, BuildSource* source, void* geometry, const size_t minLeafSize = 1, const size_t maxLeafSize = inf);

    /*! creates a leaf node */
    NodeRef createLeaf(size_t threadIndex, atomic_set<PrimRefBlock>& prims, const Split& split);

    NodeRef recurse(size_t threadIndex, size_t depth, atomic_set<PrimRefBlock>& prims, const Split& split);

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
