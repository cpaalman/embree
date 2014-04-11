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
// distributed under the License is dist ributed on an "AS IS" BASIS,        //
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. //
// See the License for the specific language governing permissions and      //
// limitations under the License.                                           //
// ======================================================================== //

#include "bvh4.h"
#include "bvh4_builder2.h"
#include "bvh4_statistics.h"

#include "builders/heuristic_binning.h"
#include "builders/heuristic_binning2.h"
#include "builders/splitter_fallback.h"

#include "common/scene_triangle_mesh.h"

namespace embree
{
  template<typename Heuristic>
  BVH4Builder2<Heuristic>::BVH4Builder2 (BVH4* bvh, BuildSource* source, void* geometry, const size_t minLeafSize, const size_t maxLeafSize)
    : source(source), geometry(geometry), primTy(*bvh->primTy[0]), minLeafSize(minLeafSize), maxLeafSize(maxLeafSize), bvh(bvh)
  {
    size_t maxLeafPrims = BVH4::maxLeafBlocks*primTy.blockSize;
    if (maxLeafPrims < this->maxLeafSize) 
      this->maxLeafSize = maxLeafPrims;
  }

  template<typename Heuristic>
  void BVH4Builder2<Heuristic>::build(size_t threadIndex, size_t threadCount) 
  {
    size_t numPrimitives = source->size();
    bvh->init(numPrimitives);
    if (source->isEmpty()) 
      return;

    if (g_verbose >= 2) 
      std::cout << "building " + bvh->name() + " with " << Heuristic::name() << " SAH builder ... " << std::flush;

    double t0 = 0.0, t1 = 0.0f;
    if (g_verbose >= 2 || g_benchmark)
      t0 = getSeconds();
    
    /* first generate primrefs */
    size_t numTriangles = 0;
    size_t numVertices = 0;
    atomic_set<PrimRefBlock> tris;
    BBox3fa geomBounds = empty;
    BBox3fa centBounds = empty;
    Scene* scene = (Scene*)geometry;
    for (size_t i=0; i<scene->size(); i++) 
    {
      Geometry* geom = scene->get(i);
      if (!geom->isEnabled()) continue;

      if (geom->type == TRIANGLE_MESH) 
      {
        TriangleMesh* set = (TriangleMesh*) geom;
        if (set->numTimeSteps != 1) continue;
        numTriangles += set->numTriangles;
        //numVertices  += set->numVertices;
        for (size_t j=0; j<set->numTriangles; j++) {
          const BBox3fa bounds = set->bounds(j);
          geomBounds.extend(bounds);
          centBounds.extend(center2(bounds));
          const PrimRef prim = PrimRef(bounds,i,j);
           atomic_set<PrimRefBlock>::item* block = tris.head();
          if (block == NULL || !block->insert(prim)) {
            block = tris.insert(alloc.malloc(threadIndex));
            block->insert(prim);
          }
        }
      }
    }
    
    /* perform binning */
    PrimInfo pinfo(numTriangles,geomBounds,centBounds);
    Heuristic heuristic(pinfo);
    atomic_set<PrimRefBlock>::iterator i=tris;
    while (PrimRefBlock* block = i.next())
      heuristic.bin(block->base(),block->size());
    
    Split split; heuristic.best(split);

    /* perform binning */
    bvh->numPrimitives = numTriangles;
    bvh->bounds = geomBounds;
    if (primTy.needVertices) bvh->numVertices = numVertices;
    else                     bvh->numVertices = 0;

    bvh->root = recurse(threadIndex,1,tris,split);

    /* free all temporary blocks */
    Alloc::global.clear();

    if (g_verbose >= 2 || g_benchmark) 
      t1 = getSeconds();
    
    if (g_verbose >= 2) {
      std::cout << "[DONE]" << std::endl;
      std::cout << "  dt = " << 1000.0f*(t1-t0) << "ms, perf = " << 1E-6*double(source->size())/(t1-t0) << " Mprim/s" << std::endl;
      std::cout << BVH4Statistics(bvh).str();
    }
  }
  
  template<typename Heuristic>
  typename BVH4Builder2<Heuristic>::NodeRef BVH4Builder2<Heuristic>::createLeaf(size_t threadIndex, atomic_set<PrimRefBlock>& prims, const Split& split)
  {
    /* allocate leaf node */
    size_t blocks = primTy.blocks(split.pinfo.size());
    char* leaf = bvh->allocPrimitiveBlocks(threadIndex,0,blocks);
    assert(blocks <= (size_t)BVH4::maxLeafBlocks);

    /* insert all triangles */
    atomic_set<PrimRefBlock>::block_iterator_unsafe iter(prims);
    for (size_t i=0; i<blocks; i++) {
      primTy.pack(leaf+i*primTy.bytes,iter,geometry);
    }
    assert(!iter);
    
    /* free all primitive blocks */
    while (atomic_set<PrimRefBlock>::item* block = prims.take())
      alloc.free(threadIndex,block);
        
    return bvh->encodeLeaf(leaf,blocks,0);
  }

  template<typename Heuristic>
  typename BVH4Builder2<Heuristic>::NodeRef BVH4Builder2<Heuristic>::recurse(size_t threadIndex, size_t depth, atomic_set<PrimRefBlock>& prims, const Split& split)
  {
    /*! compute leaf and split cost */
    const float leafSAH  = primTy.intCost*split.pinfo.sah();
    const float splitSAH = BVH4::travCost*halfArea(split.pinfo.geomBounds)+primTy.intCost*split.sah();
    assert(atomic_set<PrimRefBlock>::block_iterator_unsafe(prims).size() == split.pinfo.size());
    assert(split.pinfo.size() == 0 || leafSAH >= 0 && splitSAH >= 0);
    
    /*! create a leaf node when threshold reached or SAH tells us to stop */
    if (split.pinfo.size() <= minLeafSize || depth > BVH4::maxBuildDepth || (split.pinfo.size() <= maxLeafSize && leafSAH <= splitSAH)) {
      return createLeaf(threadIndex,prims,split);
    }
    
    /*! initialize child list */
    atomic_set<PrimRefBlock> cprims[BVH4::N]; cprims[0] = prims;
    Split                    csplit[BVH4::N]; csplit[0] = split;
    size_t numChildren = 1;
    
    /*! split until node is full or SAH tells us to stop */
    do {
      
      /*! find best child to split */
      float bestSAH = 0; 
      ssize_t bestChild = -1;
      for (size_t i=0; i<numChildren; i++) 
      {
        float dSAH = csplit[i].sah()-csplit[i].pinfo.sah();
        if (csplit[i].pinfo.size() <= minLeafSize) continue; 
        if (csplit[i].pinfo.size() > maxLeafSize) dSAH = min(0.0f,dSAH); //< force split for large jobs
        if (dSAH <= bestSAH) { bestChild = i; bestSAH = dSAH; }
      }
      if (bestChild == -1) break;
      
      /*! perform best found split and find new splits */
      Split lsplit,rsplit;
      PrimInfo linfo,rinfo;
      atomic_set<PrimRefBlock> lprims,rprims;
      csplit[bestChild].split(threadIndex,&alloc,
                              cprims[bestChild],
                              lprims,lsplit,
                              rprims,rsplit);

      Heuristic lheuristic;
      lheuristic.add(lprims);
      Split lsplit1;
      lheuristic.best(lsplit1);

      Heuristic rheuristic;
      rheuristic.add(rprims);
      Split rsplit1;
      rheuristic.best(rsplit1);
                              
      cprims[bestChild  ] = lprims; csplit[bestChild  ] = lsplit1;
      cprims[numChildren] = rprims; csplit[numChildren] = rsplit1;
      numChildren++;
      
    } while (numChildren < BVH4::N);
    
    /*! create an inner node */
    BVH4::UANode* node = bvh->allocUANode(threadIndex);
    for (size_t i=0; i<numChildren; i++) node->set(i,csplit[i].pinfo.geomBounds,recurse(threadIndex,depth+1,cprims[i],csplit[i]));
    return bvh->encodeNode(node);
  }
  
  Builder* BVH4Builder2ObjectSplit4 (void* accel, BuildSource* source, void* geometry, const size_t minLeafSize, const size_t maxLeafSize) {
    return new BVH4Builder2<HeuristicBinning2<2> >((BVH4*)accel,source,geometry,minLeafSize,maxLeafSize);
  }
}
