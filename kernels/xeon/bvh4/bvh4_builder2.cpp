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

//#include "builders/heuristic_binning.h"
//#include "builders/heuristic_binning2.h"
#include "builders/splitter_fallback.h"

#include "common/scene_triangle_mesh.h"

namespace embree
{
  const size_t BVH4Builder2::ObjectSplitBinner::maxBins;

  BVH4Builder2::ObjectSplitBinner::ObjectSplitBinner(atomic_set<PrimRefBlock>& prims) 
  {
    add(prims);
    bin(prims);
    best(split);
  }

  void BVH4Builder2::ObjectSplitBinner::add(atomic_set<PrimRefBlock>& prims)
  {
    atomic_set<PrimRefBlock>::iterator i=prims;
    while (PrimRefBlock* block = i.next()) {
      info(block->base(),block->size());
    }
  }

  void BVH4Builder2::ObjectSplitBinner::bin(atomic_set<PrimRefBlock>& prims)
  {
    new (&mapping) Mapping(pinfo);
    for (size_t i=0; i<mapping.size(); i++) {
      counts[i] = 0;
      geomBounds[i][0] = geomBounds[i][1] = geomBounds[i][2] = empty;
    }

    atomic_set<PrimRefBlock>::iterator j=prims;
    while (PrimRefBlock* block = j.next())
      bin(block->base(),block->size());
  }

  void BVH4Builder2::ObjectSplitBinner::info(const PrimRef* prims, size_t num)
  {
    BBox3fa geomBounds = empty;
    BBox3fa centBounds = empty;
    for (size_t i=0; i<num; i++)
    {
      const BBox3fa bounds = prims[i].bounds(); 
      geomBounds.extend(bounds);
      centBounds.extend(center2(bounds));
    }
    pinfo.num += num;
    pinfo.geomBounds.extend(geomBounds);
    pinfo.centBounds.extend(centBounds);
  }

  void BVH4Builder2::ObjectSplitBinner::bin(const PrimRef* prims, size_t num)
  {
    if (num == 0) return;
    
    size_t i; for (i=0; i<num-1; i+=2)
    {
      /*! map even and odd primitive to bin */
      const BBox3fa prim0 = prims[i+0].bounds(); const Vec3ia bin0 = mapping.bin(prim0); const Vec3fa center0 = Vec3fa(center2(prim0));
      const BBox3fa prim1 = prims[i+1].bounds(); const Vec3ia bin1 = mapping.bin(prim1); const Vec3fa center1 = Vec3fa(center2(prim1));
      
      /*! increase bounds for bins for even primitive */
      const int b00 = bin0.x; counts[b00][0]++; geomBounds[b00][0].extend(prim0); 
      const int b01 = bin0.y; counts[b01][1]++; geomBounds[b01][1].extend(prim0); 
      const int b02 = bin0.z; counts[b02][2]++; geomBounds[b02][2].extend(prim0); 
      
      /*! increase bounds of bins for odd primitive */
      const int b10 = bin1.x; counts[b10][0]++; geomBounds[b10][0].extend(prim1); 
      const int b11 = bin1.y; counts[b11][1]++; geomBounds[b11][1].extend(prim1); 
      const int b12 = bin1.z; counts[b12][2]++; geomBounds[b12][2].extend(prim1); 
    }
    
    /*! for uneven number of primitives */
    if (i < num)
    {
      /*! map primitive to bin */
      const BBox3fa prim0 = prims[i].bounds(); const Vec3ia bin0 = mapping.bin(prim0); const Vec3fa center0 = Vec3fa(center2(prim0));
      
      /*! increase bounds of bins */
      const int b00 = bin0.x; counts[b00][0]++; geomBounds[b00][0].extend(prim0); 
      const int b01 = bin0.y; counts[b01][1]++; geomBounds[b01][1].extend(prim0); 
      const int b02 = bin0.z; counts[b02][2]++; geomBounds[b02][2].extend(prim0); 
    }
  }
  
  void BVH4Builder2::ObjectSplitBinner::best(Split& split)
  {
    Vec3fa rAreas [maxBins];      //!< area of bounds of primitives on the right
    Vec3fa rCounts[maxBins];      //!< blocks of primitives on the right
    
    /* sweep from right to left and compute parallel prefix of merged bounds */
    assert(mapping.size() > 0);
    Vec3ia count = 0; BBox3fa bx = empty; BBox3fa by = empty; BBox3fa bz = empty;
    for (size_t i=mapping.size()-1; i>0; i--)
    {
      count += counts[i];
      rCounts[i] = Vec3fa(blocks(count));
      bx = merge(bx,geomBounds[i][0]); rAreas[i][0] = halfArea(bx);
      by = merge(by,geomBounds[i][1]); rAreas[i][1] = halfArea(by);
      bz = merge(bz,geomBounds[i][2]); rAreas[i][2] = halfArea(bz);
    }
    
    /* sweep from left to right and compute SAH */
    Vec3ia ii = 1; Vec3fa bestSAH = pos_inf; Vec3ia bestPos = 0;
    count = 0; bx = empty; by = empty; bz = empty;
    for (size_t i=1; i<mapping.size(); i++, ii+=1)
    {
      count += counts[i-1];
      bx = merge(bx,geomBounds[i-1][0]); float Ax = halfArea(bx);
      by = merge(by,geomBounds[i-1][1]); float Ay = halfArea(by);
      bz = merge(bz,geomBounds[i-1][2]); float Az = halfArea(bz);
      const Vec3fa lCount = Vec3fa(blocks(count));
      const Vec3fa lArea = Vec3fa(Ax,Ay,Az);
      const Vec3fa sah = lArea*lCount + rAreas[i]*rCounts[i];
      bestPos = select(lt_mask(sah,bestSAH),ii ,bestPos);
      bestSAH = select(lt_mask(sah,bestSAH),sah,bestSAH);
    }

    /* find best dimension */
    for (int i=0; i<3; i++) 
    {
      if (unlikely(pinfo.centBounds.lower[i] >= pinfo.centBounds.upper[i])) 
        continue;
      
      if (bestSAH[i] < split.cost && bestPos[i] != 0) {
        split.dim = i;
        split.pos = bestPos[i];
        split.cost = bestSAH[i];
      }
    }

    split.mapping = mapping;
    split.pinfo = pinfo;
  }
  
  void BVH4Builder2::ObjectSplitBinner::Split::split(size_t thread, PrimRefAlloc* alloc, 
                                                     atomic_set<PrimRefBlock>& prims, 
                                                     atomic_set<PrimRefBlock>& lprims, 
                                                     atomic_set<PrimRefBlock>& rprims)
  {
    atomic_set<PrimRefBlock>::item* lblock = lprims.insert(alloc->malloc(thread));
    atomic_set<PrimRefBlock>::item* rblock = rprims.insert(alloc->malloc(thread));
    
    while (atomic_set<PrimRefBlock>::item* block = prims.take()) 
    {
      for (size_t i=0; i<block->size(); i++) 
      {
        const PrimRef& prim = block->at(i); 
                
        if (mapping.bin_unsafe(prim.bounds())[dim] < pos) 
        {
          if (likely(lblock->insert(prim))) continue; 
          lblock = lprims.insert(alloc->malloc(thread));
          lblock->insert(prim);
        } 
        else 
        {
          if (likely(rblock->insert(prim))) continue;
          rblock = rprims.insert(alloc->malloc(thread));
          rblock->insert(prim);
        }
      }
      alloc->free(thread,block);
    }
  }

  BVH4Builder2::BVH4Builder2 (BVH4* bvh, BuildSource* source, void* geometry, const size_t minLeafSize, const size_t maxLeafSize)
    : source(source), geometry(geometry), primTy(*bvh->primTy[0]), minLeafSize(minLeafSize), maxLeafSize(maxLeafSize), bvh(bvh)
  {
    size_t maxLeafPrims = BVH4::maxLeafBlocks*primTy.blockSize;
    if (maxLeafPrims < this->maxLeafSize) 
      this->maxLeafSize = maxLeafPrims;
  }

  void BVH4Builder2::build(size_t threadIndex, size_t threadCount) 
  {
    size_t numPrimitives = source->size();
    bvh->init(numPrimitives);
    if (source->isEmpty()) 
      return;

    if (g_verbose >= 2) 
      std::cout << "building " + bvh->name() + " with SAH builder ... " << std::flush;

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
    ObjectSplitBinner heuristic(tris); 

#if 0
    float bestSAH = inf;
    ObjectBinning object_binning_aligned(tris,beziers);
    bestSAH = min(bestSAH,object_binning_aligned.sah());

    SpatialBinning spatial_binning_aligned(tris,beziers);
    bestSAH = min(bestSAH,spatial_binning_aligned.sah());
    
    Space hairspace = computeSpace(tris,beziers);
    ObjectBinning object_binning_unaligned(tris,beziers,hairspace);
    bestSAH = min(bestSAH,object_binning_unaligned.sah());
    
    SpatialBinning spatial_binning_unaligned(tris,beziers,hairspace);
    bestSAH = min(bestSAH,spatial_binning_unaligned.sah());

    Split split;
    if (bestSAH == object_binning_aligned.sah()) 
      new (&split) Split(object_binning_aligned.split);
    else if (bestSAH == spatial_binning_aligned.sah()) 
      new (&split) Split(spatial_binning_aligned.split);
    else if (bestSAH == object_binning_unaligned.sah()) 
      new (&split) Split(object_binning_unaligned.split);
    else if (bestSAH == spatial_binning_unaligned.sah()) 
      new (&split) Split(spatial_binning_unaligned.split);
    
#endif

    /* perform binning */
    bvh->numPrimitives = numTriangles;
    bvh->bounds = geomBounds;
    if (primTy.needVertices) bvh->numVertices = numVertices;
    else                     bvh->numVertices = 0;

    bvh->root = recurse(threadIndex,1,tris,heuristic.split);

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
  
  typename BVH4Builder2::NodeRef BVH4Builder2::createLeaf(size_t threadIndex, atomic_set<PrimRefBlock>& prims, const ObjectSplitBinner::Split& split)
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

  typename BVH4Builder2::NodeRef BVH4Builder2::recurse(size_t threadIndex, size_t depth, atomic_set<PrimRefBlock>& prims, const ObjectSplitBinner::Split& split)
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
    ObjectSplitBinner::Split                    csplit[BVH4::N]; csplit[0] = split;
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
      atomic_set<PrimRefBlock> lprims,rprims;
      csplit[bestChild].split(threadIndex,&alloc,cprims[bestChild],lprims,rprims);
      ObjectSplitBinner lheuristic(lprims); 
      ObjectSplitBinner rheuristic(rprims);
      cprims[bestChild  ] = lprims; csplit[bestChild  ] = lheuristic.split;
      cprims[numChildren] = rprims; csplit[numChildren] = rheuristic.split;
      numChildren++;
      
    } while (numChildren < BVH4::N);
    
    /*! create an inner node */
    BVH4::UANode* node = bvh->allocUANode(threadIndex);
    for (size_t i=0; i<numChildren; i++) node->set(i,csplit[i].pinfo.geomBounds,recurse(threadIndex,depth+1,cprims[i],csplit[i]));
    return bvh->encodeNode(node);
  }
  
  Builder* BVH4Builder2ObjectSplit4 (void* accel, BuildSource* source, void* geometry, const size_t minLeafSize, const size_t maxLeafSize) {
    return new BVH4Builder2((BVH4*)accel,source,geometry,minLeafSize,maxLeafSize);
  }
}
