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
#include "common/scene_bezier_curves.h"
#include "common/scene_triangle_mesh.h"
#include "geometry/bezier1.h"
#include "geometry/triangle4.h"

//#include "builders/heuristic_binning.h"
//#include "builders/heuristic_binning2.h"
#include "builders/splitter_fallback.h"

#include "common/scene_triangle_mesh.h"

namespace embree
{
  const size_t BVH4Builder2::ObjectSplitBinner::maxBins;

  Scene* g_scene = NULL; // FIXME: remove me
  BVH4Builder2* g_builder2 = NULL;

  BBox3fa primBounds(const LinearSpace3fa& space, const PrimRef& prim)
  {
    TriangleMesh* mesh = (TriangleMesh*) g_scene->get(prim.geomID());
    TriangleMesh::Triangle tri = mesh->triangle(prim.primID());
    BBox3fa bounds = empty;
    bounds.extend(xfmPoint(space,mesh->vertex(tri.v[0])));
    bounds.extend(xfmPoint(space,mesh->vertex(tri.v[1])));
    bounds.extend(xfmPoint(space,mesh->vertex(tri.v[2])));
    return bounds;
  }

  BVH4Builder2::ObjectSplitBinner::ObjectSplitBinner(const LinearSpace3fa& space, TriRefList& triangles) 
    : space(space)
  {
    add(triangles);
    setup_binning();
    bin(triangles);
    best(split);
  }

  BVH4Builder2::ObjectSplitBinner::ObjectSplitBinner(const LinearSpace3fa& space, TriRefList& triangles, BezierRefList& beziers) 
    : space(space)
  {
    add(triangles);
    add(beziers);
    setup_binning();
    bin(triangles);
    bin(beziers);
    best(split);
  }

  void BVH4Builder2::ObjectSplitBinner::add(TriRefList& prims)
  {
    TriRefList::iterator i=prims;
    while (TriRefBlock* block = i.next()) {
      info(block->base(),block->size());
    }
  }

  void BVH4Builder2::ObjectSplitBinner::add(BezierRefList& prims)
  {
    BezierRefList::iterator i=prims;
    while (BezierRefBlock* block = i.next()) {
      info(block->base(),block->size());
    }
  }

  void BVH4Builder2::ObjectSplitBinner::setup_binning()
  {
    new (&mapping) Mapping(pinfo);
    for (size_t i=0; i<mapping.size(); i++) {
      triCounts[i] = bezierCounts[i] = 0;
      geomBounds[i][0] = geomBounds[i][1] = geomBounds[i][2] = empty;
    }
  }

  void BVH4Builder2::ObjectSplitBinner::bin(TriRefList& prims)
  {
    TriRefList::iterator j=prims;
    while (TriRefBlock* block = j.next())
      bin(block->base(),block->size());
  }

  void BVH4Builder2::ObjectSplitBinner::bin(BezierRefList& prims)
  {
    BezierRefList::iterator j=prims;
    while (BezierRefBlock* block = j.next())
      bin(block->base(),block->size());
  }

  void BVH4Builder2::ObjectSplitBinner::info(const PrimRef* prims, size_t num)
  {
    BBox3fa geomBounds = empty;
    BBox3fa centBounds = empty;
    for (size_t i=0; i<num; i++)
    {
      //const BBox3fa bounds = prims[i].bounds(space); 
      const BBox3fa bounds = primBounds(space,prims[i]); 
      geomBounds.extend(bounds);
      centBounds.extend(center2(bounds));
    }
    pinfo.num += num;
    pinfo.numTriangles += num;
    pinfo.geomBounds.extend(geomBounds);
    pinfo.centBounds.extend(centBounds);
  }

  void BVH4Builder2::ObjectSplitBinner::info(const Bezier1* prims, size_t num)
  {
    BBox3fa geomBounds = empty;
    BBox3fa centBounds = empty;
    for (size_t i=0; i<num; i++)
    {
      const BBox3fa bounds = prims[i].bounds(space); 
      geomBounds.extend(bounds);
      centBounds.extend(center2(bounds));
    }
    pinfo.num += num;
    pinfo.numBeziers += num;
    pinfo.geomBounds.extend(geomBounds);
    pinfo.centBounds.extend(centBounds);
  }

  void BVH4Builder2::ObjectSplitBinner::bin(const PrimRef* prims, size_t num)
  {
    if (num == 0) return;
    
    size_t i; for (i=0; i<num-1; i+=2)
    {
      /*! map even and odd primitive to bin */
      //const BBox3fa prim0 = prims[i+0].bounds(space); const Vec3ia bin0 = mapping.bin(prim0); const Vec3fa center0 = Vec3fa(center2(prim0));
      //const BBox3fa prim1 = prims[i+1].bounds(space); const Vec3ia bin1 = mapping.bin(prim1); const Vec3fa center1 = Vec3fa(center2(prim1));
      const BBox3fa prim0 = primBounds(space,prims[i+0]); const Vec3ia bin0 = mapping.bin(prim0); const Vec3fa center0 = Vec3fa(center2(prim0));
      const BBox3fa prim1 = primBounds(space,prims[i+1]); const Vec3ia bin1 = mapping.bin(prim1); const Vec3fa center1 = Vec3fa(center2(prim1));
      
      /*! increase bounds for bins for even primitive */
      const int b00 = bin0.x; triCounts[b00][0]++; geomBounds[b00][0].extend(prim0); 
      const int b01 = bin0.y; triCounts[b01][1]++; geomBounds[b01][1].extend(prim0); 
      const int b02 = bin0.z; triCounts[b02][2]++; geomBounds[b02][2].extend(prim0); 
      
      /*! increase bounds of bins for odd primitive */
      const int b10 = bin1.x; triCounts[b10][0]++; geomBounds[b10][0].extend(prim1); 
      const int b11 = bin1.y; triCounts[b11][1]++; geomBounds[b11][1].extend(prim1); 
      const int b12 = bin1.z; triCounts[b12][2]++; geomBounds[b12][2].extend(prim1); 
    }
    
    /*! for uneven number of primitives */
    if (i < num)
    {
      /*! map primitive to bin */
      //const BBox3fa prim0 = prims[i].bounds(space); const Vec3ia bin0 = mapping.bin(prim0); const Vec3fa center0 = Vec3fa(center2(prim0));
      const BBox3fa prim0 = primBounds(space,prims[i]); const Vec3ia bin0 = mapping.bin(prim0); const Vec3fa center0 = Vec3fa(center2(prim0));
      
      /*! increase bounds of bins */
      const int b00 = bin0.x; triCounts[b00][0]++; geomBounds[b00][0].extend(prim0); 
      const int b01 = bin0.y; triCounts[b01][1]++; geomBounds[b01][1].extend(prim0); 
      const int b02 = bin0.z; triCounts[b02][2]++; geomBounds[b02][2].extend(prim0); 
    }
  }

  void BVH4Builder2::ObjectSplitBinner::bin(const Bezier1* prims, size_t num)
  {
    for (size_t i=0; i<num; i++)
    {
      /*! map even and odd primitive to bin */
      const BBox3fa prim0 = prims[i+0].bounds(space); const Vec3ia bin0 = mapping.bin(prim0); const Vec3fa center0 = Vec3fa(center2(prim0));
      
      /*! increase bounds for bins for even primitive */
      const int b00 = bin0.x; bezierCounts[b00][0]++; geomBounds[b00][0].extend(prim0); 
      const int b01 = bin0.y; bezierCounts[b01][1]++; geomBounds[b01][1].extend(prim0); 
      const int b02 = bin0.z; bezierCounts[b02][2]++; geomBounds[b02][2].extend(prim0); 
    }
  }
  
  void BVH4Builder2::ObjectSplitBinner::best(Split& split)
  {
    Vec3fa rAreas [maxBins];      //!< area of bounds of primitives on the right
    Vec3fa rCounts[maxBins];      //!< blocks of primitives on the right
    Vec3fa rBezierCounts[maxBins];      //!< blocks of primitives on the right
    
    /* sweep from right to left and compute parallel prefix of merged bounds */
    assert(mapping.size() > 0);
    Vec3ia triCount = 0, bezierCount = 0; 
    BBox3fa bx = empty; BBox3fa by = empty; BBox3fa bz = empty;
    for (size_t i=mapping.size()-1; i>0; i--)
    {
      triCount += triCounts[i];
      bezierCount += bezierCounts[i];
      rCounts[i] = Vec3fa(blocks(triCount) + bezierCount);
      bx = merge(bx,geomBounds[i][0]); rAreas[i][0] = halfArea(bx);
      by = merge(by,geomBounds[i][1]); rAreas[i][1] = halfArea(by);
      bz = merge(bz,geomBounds[i][2]); rAreas[i][2] = halfArea(bz);
    }
    
    /* sweep from left to right and compute SAH */
    Vec3ia ii = 1; Vec3fa bestSAH = pos_inf; Vec3ia bestPos = 0;
    triCount = 0; bezierCount = 0; bx = empty; by = empty; bz = empty;
    for (size_t i=1; i<mapping.size(); i++, ii+=1)
    {
      triCount += triCounts[i-1];
      bezierCount += bezierCounts[i-1];
      bx = merge(bx,geomBounds[i-1][0]); float Ax = halfArea(bx);
      by = merge(by,geomBounds[i-1][1]); float Ay = halfArea(by);
      bz = merge(bz,geomBounds[i-1][2]); float Az = halfArea(bz);
      const Vec3fa lCount = Vec3fa(blocks(triCount) + bezierCount);
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
    split.space = space;
  }

  void BVH4Builder2::split_fallback(size_t threadIndex, PrimRefBlockAlloc<PrimRef>* alloc, TriRefList& prims, TriRefList& lprims, TriRefList& rprims)
  {  
    size_t num = 0;
    atomic_set<TriRefBlock>::item* lblock = lprims.insert(alloc->malloc(threadIndex));
    atomic_set<TriRefBlock>::item* rblock = rprims.insert(alloc->malloc(threadIndex));
    
    while (atomic_set<TriRefBlock>::item* block = prims.take()) 
    {
      for (size_t i=0; i<block->size(); i++) 
      {
        const PrimRef& prim = block->at(i); 
        
        if (num&1) 
        {
          num++;
          if (likely(lblock->insert(prim))) continue; 
          lblock = lprims.insert(alloc->malloc(threadIndex));
          lblock->insert(prim);
        } else {
          num++;
          if (likely(rblock->insert(prim))) continue;
          rblock = rprims.insert(alloc->malloc(threadIndex));
          rblock->insert(prim);
        }
      }
      alloc->free(threadIndex,block);
    }
  }

  void BVH4Builder2::split_fallback(size_t threadIndex, PrimRefBlockAlloc<Bezier1>* alloc, BezierRefList& prims, BezierRefList& lprims, BezierRefList& rprims)
  {  
    size_t num = 0;
    atomic_set<BezierRefBlock>::item* lblock = lprims.insert(alloc->malloc(threadIndex));
    atomic_set<BezierRefBlock>::item* rblock = rprims.insert(alloc->malloc(threadIndex));
    
    while (atomic_set<BezierRefBlock>::item* block = prims.take()) 
    {
      for (size_t i=0; i<block->size(); i++) 
      {
        const Bezier1& prim = block->at(i); 
        
        if (num&1) 
        {
          num++;
          if (likely(lblock->insert(prim))) continue; 
          lblock = lprims.insert(alloc->malloc(threadIndex));
          lblock->insert(prim);
        } else {
          num++;
          if (likely(rblock->insert(prim))) continue;
          rblock = rprims.insert(alloc->malloc(threadIndex));
          rblock->insert(prim);
        }
      }
      alloc->free(threadIndex,block);
    }
  }

  void BVH4Builder2::ObjectSplitBinner::Split::split(size_t thread, PrimRefBlockAlloc<PrimRef>* alloc, TriRefList& prims, TriRefList& lprims, TriRefList& rprims)
  {
    if (cost == float(inf)) {
      g_builder2->split_fallback(thread,alloc,prims,lprims,rprims);
      return;
    }

    TriRefList::item* lblock = lprims.insert(alloc->malloc(thread));
    TriRefList::item* rblock = rprims.insert(alloc->malloc(thread));
    
    while (TriRefList::item* block = prims.take()) 
    {
      for (size_t i=0; i<block->size(); i++) 
      {
        const PrimRef& prim = block->at(i); 
                
        if (mapping.bin_unsafe(primBounds(space,prim))[dim] < pos) 
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

  void BVH4Builder2::ObjectSplitBinner::Split::split(size_t thread, PrimRefBlockAlloc<Bezier1>* alloc, BezierRefList& prims, BezierRefList& lprims, BezierRefList& rprims)
  {
    if (cost == float(inf)) {
      g_builder2->split_fallback(thread,alloc,prims,lprims,rprims);
      return;
    }

    BezierRefList::item* lblock = lprims.insert(alloc->malloc(thread));
    BezierRefList::item* rblock = rprims.insert(alloc->malloc(thread));
    
    while (BezierRefList::item* block = prims.take()) 
    {
      for (size_t i=0; i<block->size(); i++) 
      {
        const Bezier1& prim = block->at(i); 
                
        if (mapping.bin_unsafe(prim.bounds(space))[dim] < pos) 
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
    : source(source), geometry(geometry), /*primTy(*bvh->primTy[0]),*/ minLeafSize(minLeafSize), maxLeafSize(maxLeafSize), bvh(bvh)
  {
    size_t maxLeafTris    = BVH4::maxLeafBlocks*bvh->primTy[0]->blockSize;
    size_t maxLeafBeziers = BVH4::maxLeafBlocks*bvh->primTy[1]->blockSize;
    if (maxLeafTris < this->maxLeafSize) this->maxLeafSize = maxLeafTris;
    if (maxLeafBeziers < this->maxLeafSize) this->maxLeafSize = maxLeafBeziers; // FIXME: keep separate for tris and beziers
  }

  const BBox3fa BVH4Builder2::computeAlignedBounds(TriRefList& tris)
  {
    BBox3fa bounds = empty;
    for (TriRefList::block_iterator_unsafe i=tris; i; i++)
      bounds.extend(i->bounds());
    return bounds;
  }

  const BBox3fa BVH4Builder2::computeAlignedBounds(BezierRefList& beziers)
  {
    BBox3fa bounds = empty;
    for (BezierRefList::block_iterator_unsafe i=beziers; i; i++)
      bounds.extend(i->bounds());
    return bounds;
  }

  const BBox3fa BVH4Builder2::computeAlignedBounds(TriRefList& tris, BezierRefList& beziers) {
    return merge(computeAlignedBounds(tris),computeAlignedBounds(beziers));
  }

  const NAABBox3fa BVH4Builder2::computeAlignedBounds(TriRefList& tris, const LinearSpace3fa& space)
  {
    BBox3fa bounds = empty;
    for (TriRefList::block_iterator_unsafe i=tris; i; i++)
      bounds.extend(primBounds(space,*i));
    return NAABBox3fa(space,bounds);
  }

  const NAABBox3fa BVH4Builder2::computeAlignedBounds(BezierRefList& beziers, const LinearSpace3fa& space)
  {
    BBox3fa bounds = empty;
    for (BezierRefList::block_iterator_unsafe i=beziers; i; i++)
      bounds.extend(i->bounds(space));
    return NAABBox3fa(space,bounds);
  }

  const NAABBox3fa BVH4Builder2::computeAlignedBounds(TriRefList& tris, BezierRefList& beziers, const LinearSpace3fa& space) {
    return NAABBox3fa(space,merge(computeAlignedBounds(tris,space).bounds,computeAlignedBounds(beziers,space).bounds));
  }

  const LinearSpace3fa BVH4Builder2::computeHairSpace(BezierRefList& prims)
  {
    size_t N = BezierRefList::block_iterator_unsafe(prims).size();
    if (N == 0)
      return one; // FIXME: can cause problems with compression

    float bestArea = inf;
    LinearSpace3fa bestSpace = one;
    BBox3fa bestBounds = empty;

    size_t k=0;
    for (BezierRefList::block_iterator_unsafe i = prims; i; i++)
    {
      if ((k++) % ((N+3)/4)) continue;
      //size_t k = begin + rand() % (end-begin);
      const Vec3fa axis = normalize(i->p3 - i->p0);
      if (length(i->p3 - i->p0) < 1E-9) continue;
      const LinearSpace3fa space0 = LinearSpace3fa::rotate(Vec3fa(0,0,1),2.0f*float(pi)*drand48())*frame(axis).transposed();
      const LinearSpace3fa space = clamp(space0);
      BBox3fa bounds = empty;
      float area = 0.0f;
      for (BezierRefList::block_iterator_unsafe j = prims; j; j++) {
        const BBox3fa cbounds = j->bounds(space);
        area += embree::area(cbounds);
        bounds.extend(cbounds);
      }

      if (area <= bestArea) {
        bestBounds = bounds;
        bestSpace = space;
        bestArea = area;
      }
    }
    //assert(bestArea != (float)inf); // FIXME: can get raised if all selected curves are points
#ifdef DEBUG
    if (bestArea == (float)inf)
      {
        std::cout << "WARNING: bestArea == (float)inf" << std::endl; 
      }
#endif

    return bestSpace;
  }

  void BVH4Builder2::build(size_t threadIndex, size_t threadCount) 
  {
    size_t numBeziers = 0;
    size_t numTriangles = 0;
    Scene* scene = (Scene*) geometry;
    g_scene = scene; // FIXME: hack
    g_builder2 = this; // FIXME: hack
    for (size_t i=0; i<scene->size(); i++) 
    {
      Geometry* geom = scene->get(i);
      if (!geom->isEnabled()) continue;

      if (geom->type == BEZIER_CURVES) 
      {
        BezierCurves* set = (BezierCurves*) geom;
        numBeziers  += set->numCurves;
      }

      if (geom->type == TRIANGLE_MESH) 
      {
        TriangleMesh* set = (TriangleMesh*) geom;
        if (set->numTimeSteps != 1) continue;
        numTriangles += set->numTriangles;
      }
    }
    size_t numPrimitives = numBeziers + numTriangles;

    bvh->init(numPrimitives); //+(size_t)(g_hair_builder_replication_factor*numPrimitives));
    if (numPrimitives == 0) return;

    if (g_verbose >= 2) 
      std::cout << "building " + bvh->name() + " with SAH builder ... " << std::flush;

    double t0 = 0.0, t1 = 0.0f;
    if (g_verbose >= 2 || g_benchmark)
      t0 = getSeconds();
    
    /* first generate primrefs */
    //size_t numTriangles = 0;
    size_t numVertices = 0;
    TriRefList tris;
    BezierRefList beziers;
    BBox3fa geomBounds = empty;
    //BBox3fa centBounds = empty;

    for (size_t i=0; i<scene->size(); i++) 
    {
      Geometry* geom = scene->get(i);
      if (!geom->isEnabled()) continue;

      if (geom->type == TRIANGLE_MESH) 
      {
        TriangleMesh* set = (TriangleMesh*) geom;
        if (set->numTimeSteps != 1) continue;
        //numTriangles += set->numTriangles;
        //numVertices  += set->numVertices;
        for (size_t j=0; j<set->numTriangles; j++) {
          const BBox3fa bounds = set->bounds(j);
          geomBounds.extend(bounds);
          //centBounds.extend(center2(bounds));
          const PrimRef prim = PrimRef(bounds,i,j);
           TriRefList::item* block = tris.head();
          if (block == NULL || !block->insert(prim)) {
            block = tris.insert(allocTriRefs.malloc(threadIndex));
            block->insert(prim);
          }
        }
      }

      if (geom->type == BEZIER_CURVES) 
      {
        BezierCurves* set = (BezierCurves*) geom;
        //numBeziers  += set->numCurves;
        //numVertices += set->numVertices;
        for (size_t j=0; j<set->numCurves; j++) {
          const int ofs = set->curve(j);
          const Vec3fa& p0 = set->vertex(ofs+0);
          const Vec3fa& p1 = set->vertex(ofs+1);
          const Vec3fa& p2 = set->vertex(ofs+2);
          const Vec3fa& p3 = set->vertex(ofs+3);
          const Bezier1 bezier(p0,p1,p2,p3,0,1,i,j);
          //bounds.extend(subdivideAndAdd(threadIndex,prims,bezier,enablePreSubdivision));
          geomBounds.extend(bezier.bounds());
          BezierRefList::item* block = beziers.head();
          if (block == NULL || !block->insert(bezier)) {
            block = beziers.insert(allocBezierRefs.malloc(threadIndex));
            block->insert(bezier);
          }
        }
      }
    }
    
    /* perform binning */
    ObjectSplitBinner heuristic(one,tris,beziers); 

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
    bvh->numPrimitives = numPrimitives;
    bvh->bounds = geomBounds;
    //if (primTy.needVertices) bvh->numVertices = numVertices; // FIXME
    //else                     bvh->numVertices = 0;

    bvh->root = recurse(threadIndex,1,tris,beziers,heuristic.split);

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

  void BVH4Builder2::heuristic(TriRefList& tris, BezierRefList& beziers, ObjectSplitBinner::Split& split)
  {
    /* perform binning */
#if 0
    ObjectSplitBinner heuristic(one,tris,beziers); 
    split = heuristic.split;
    split.aligned = true;

#else

    float bestSAH = inf;
    ObjectSplitBinner object_binning_aligned(one,tris,beziers);
    bestSAH = min(bestSAH,object_binning_aligned.split.splitSAH());

    //SpatialBinning spatial_binning_aligned(tris,beziers);
    //bestSAH = min(bestSAH,spatial_binning_aligned.split.splitSAH());
    
    const LinearSpace3fa hairspace = computeHairSpace(beziers);
    ObjectSplitBinner object_binning_unaligned(hairspace,tris,beziers);
    bestSAH = min(bestSAH,object_binning_unaligned.split.splitSAH());
    
    //SpatialBinning spatial_binning_unaligned(tris,beziers,hairspace);
    //bestSAH = min(bestSAH,spatial_binning_unaligned.split.splitSAH());

    if (bestSAH == object_binning_aligned.split.splitSAH()) {
      new (&split) ObjectSplitBinner::Split(object_binning_aligned.split);
      split.aligned = true;
    }
    //else if (bestSAH == spatial_binning_aligned.split.splitSAH()) {
    //new (&split) Split(spatial_binning_aligned.split);
    //split.aligned = true;
    //}
    else if (bestSAH == object_binning_unaligned.split.splitSAH()) {
      new (&split) ObjectSplitBinner::Split(object_binning_unaligned.split);
      split.aligned = false;
    }
    //else if (bestSAH == spatial_binning_unaligned.split.splitSAH()) {
    //new (&split) Split(spatial_binning_unaligned.split);
    //split.aligned = true;
    //}
    
#endif
  }

  BVH4::NodeRef BVH4Builder2::leaf(size_t threadIndex, size_t depth, TriRefList& prims, const ObjectSplitBinner::Split& split)
  {
    size_t N = bvh->primTy[0]->blocks(TriRefList::block_iterator_unsafe(prims).size());

    if (N > (size_t)BVH4::maxLeafBlocks) {
      std::cout << "WARNING: Loosing " << N-BVH4::maxLeafBlocks << " primitives during build!" << std::endl;
      //std::cout << "!" << std::flush;
      N = (size_t)BVH4::maxLeafBlocks;
    }
    Scene* scene = (Scene*) geometry;

    if (bvh->primTy[0] == &SceneTriangle4::type) 
    {
      char* leaf = bvh->allocPrimitiveBlocks(threadIndex,0,N);

      TriRefList::block_iterator_unsafe iter(prims);
      for (size_t j=0; j<N; j++) 
      {
        void* This = leaf+j*bvh->primTy[0]->bytes;
        ssei geomID = -1, primID = -1, mask = -1;
        sse3f v0 = zero, v1 = zero, v2 = zero;
    
        for (size_t i=0; i<4 && iter; i++, iter++)
        {
          const PrimRef& prim = *iter;
          const TriangleMesh* mesh = scene->getTriangleMesh(prim.geomID());
          const TriangleMesh::Triangle& tri = mesh->triangle(prim.primID());
          const Vec3fa& p0 = mesh->vertex(tri.v[0]);
          const Vec3fa& p1 = mesh->vertex(tri.v[1]);
          const Vec3fa& p2 = mesh->vertex(tri.v[2]);
          geomID [i] = prim.geomID();
          primID [i] = prim.primID();
          mask   [i] = mesh->mask;
          v0.x[i] = p0.x; v0.y[i] = p0.y; v0.z[i] = p0.z;
          v1.x[i] = p1.x; v1.y[i] = p1.y; v1.z[i] = p1.z;
          v2.x[i] = p2.x; v2.y[i] = p2.y; v2.z[i] = p2.z;
        }
        new (This) Triangle4(v0,v1,v2,geomID,primID,mask);
      }
      assert(!iter);
    
      /* free all primitive blocks */
      while (TriRefList::item* block = prims.take())
        allocTriRefs.free(threadIndex,block);
        
      return bvh->encodeLeaf(leaf,N,0);
    } 
    else 
      throw std::runtime_error("unknown primitive type");
  }

  BVH4::NodeRef BVH4Builder2::leaf(size_t threadIndex, size_t depth, BezierRefList& prims, const ObjectSplitBinner::Split& split)
  {
    size_t N = BezierRefList::block_iterator_unsafe(prims).size();

    if (N > (size_t)BVH4::maxLeafBlocks) {
      std::cout << "WARNING: Loosing " << N-BVH4::maxLeafBlocks << " primitives during build!" << std::endl;
      //std::cout << "!" << std::flush;
      N = (size_t)BVH4::maxLeafBlocks;
    }
    Scene* scene = (Scene*) geometry;

    if (bvh->primTy[1] == &Bezier1Type::type) 
    { 
      Bezier1* leaf = (Bezier1*) bvh->allocPrimitiveBlocks(threadIndex,1,N);
      BezierRefList::block_iterator_unsafe iter(prims);
      for (size_t i=0; i<N; i++) { leaf[i] = *iter; iter++; }
      assert(!iter);

      /* free all primitive blocks */
      while (BezierRefList::item* block = prims.take())
        allocBezierRefs.free(threadIndex,block);

      return bvh->encodeLeaf((char*)leaf,N,1);
    } 
    else if (bvh->primTy[1] == &SceneBezier1i::type) 
    {
      Bezier1i* leaf = (Bezier1i*) bvh->allocPrimitiveBlocks(threadIndex,1,N);
      BezierRefList::block_iterator_unsafe iter(prims);
      for (size_t i=0; i<N; i++) {
        const Bezier1& curve = *iter; iter++;
        const BezierCurves* in = (BezierCurves*) scene->get(curve.geomID);
        const Vec3fa& p0 = in->vertex(in->curve(curve.primID));
        leaf[i] = Bezier1i(&p0,curve.geomID,curve.primID,-1); // FIXME: support mask
      }

      /* free all primitive blocks */
      while (BezierRefList::item* block = prims.take())
        allocBezierRefs.free(threadIndex,block);

      return bvh->encodeLeaf((char*)leaf,N,1);
    }
    else 
      throw std::runtime_error("unknown primitive type");
  }

  BVH4::NodeRef BVH4Builder2::leaf(size_t threadIndex, size_t depth, TriRefList& tris, BezierRefList& beziers, const ObjectSplitBinner::Split& split)
  {
    size_t Ntris    = TriRefList   ::block_iterator_unsafe(tris   ).size();
    size_t Nbeziers = BezierRefList::block_iterator_unsafe(beziers).size();
    if (Ntris == 0) 
      return leaf(threadIndex,depth,beziers,split);
    else if (Nbeziers == 0)
      return leaf(threadIndex,depth,tris,split);
    
    BVH4::UANode* node = bvh->allocUANode(threadIndex);
    node->set(0,computeAlignedBounds(tris   ),leaf(threadIndex,depth+1,tris   ,split));
    node->set(1,computeAlignedBounds(beziers),leaf(threadIndex,depth+1,beziers,split));
    return bvh->encodeNode(node);
  }

  typename BVH4Builder2::NodeRef BVH4Builder2::recurse(size_t threadIndex, size_t depth, TriRefList& tris, BezierRefList& beziers, const ObjectSplitBinner::Split& split)
  {
    //PRINT(depth);

    /*! compute leaf and split cost */
    const float leafSAH  = /*primTy.intCost*/split.leafSAH ();
    const float splitSAH = /*primTy.intCost*/split.splitSAH() + BVH4::travCost*halfArea(split.pinfo.geomBounds);
    //PRINT(leafSAH);
    //PRINT(splitSAH);
    //PRINT(split.pinfo.size());
    //PRINT(maxLeafSize);
    assert(TriRefList::block_iterator_unsafe(prims).size() == split.pinfo.size());
    assert(split.pinfo.size() == 0 || leafSAH >= 0 && splitSAH >= 0);
    
    /*! create a leaf node when threshold reached or SAH tells us to stop */
    if (split.pinfo.size() <= minLeafSize || depth > BVH4::maxBuildDepth || (split.pinfo.size() <= maxLeafSize && leafSAH <= splitSAH)) {
      return leaf(threadIndex,depth,tris,beziers,split);
    }
    
    /*! initialize child list */
    TriRefList    ctris   [BVH4::N]; ctris   [0] = tris;
    BezierRefList cbeziers[BVH4::N]; cbeziers[0] = beziers;
    ObjectSplitBinner::Split csplit[BVH4::N]; csplit[0] = split;
    size_t numChildren = 1;
    bool aligned = true;

    /*! split until node is full or SAH tells us to stop */
    do {
      
      /*! find best child to split */
      float bestSAH = 0; 
      ssize_t bestChild = -1;
      for (size_t i=0; i<numChildren; i++) 
      {
        float dSAH = csplit[i].splitSAH()-csplit[i].leafSAH();
        if (csplit[i].pinfo.size() <= minLeafSize) continue; 
        if (csplit[i].pinfo.size() > maxLeafSize) dSAH = min(0.0f,dSAH); //< force split for large jobs
        if (dSAH <= bestSAH) { bestChild = i; bestSAH = dSAH; }
      }
      if (bestChild == -1) break;
      
      /*! perform best found split and find new splits */
      aligned &= csplit[bestChild].aligned;
      TriRefList    ltris,   rtris;    csplit[bestChild].split(threadIndex,&allocTriRefs,   ctris[bestChild],   ltris,   rtris);
      BezierRefList lbeziers,rbeziers; csplit[bestChild].split(threadIndex,&allocBezierRefs,cbeziers[bestChild],lbeziers,rbeziers);
      ObjectSplitBinner::Split lsplit; heuristic(ltris,lbeziers,lsplit);
      ObjectSplitBinner::Split rsplit; heuristic(rtris,rbeziers,rsplit);
      ctris[bestChild  ] = ltris; cbeziers[bestChild  ] = lbeziers; csplit[bestChild  ] = lsplit;
      ctris[numChildren] = rtris; cbeziers[numChildren] = rbeziers; csplit[numChildren] = rsplit;
      numChildren++;
      
    } while (numChildren < BVH4::N);
    
    /*! create an aligned node */
    if (aligned)
    {
      BVH4::UANode* node = bvh->allocUANode(threadIndex);
      for (size_t i=0; i<numChildren; i++) {
        const BBox3fa bounds = computeAlignedBounds(ctris[i],cbeziers[i]);
        node->set(i,bounds,recurse(threadIndex,depth+1,ctris[i],cbeziers[i],csplit[i]));
      }
      return bvh->encodeNode(node);
    } 

    /*! create an unaligned node */
    else
    {
      BVH4::UUNode* node = bvh->allocUUNode(threadIndex);
      for (size_t i=0; i<numChildren; i++) {
        const LinearSpace3fa hairspace = computeHairSpace(cbeziers[i]);
        const NAABBox3fa bounds = computeAlignedBounds(ctris[i],cbeziers[i],hairspace);
        node->set(i,bounds,recurse(threadIndex,depth+1,ctris[i],cbeziers[i],csplit[i]));
      }
      return bvh->encodeNode(node);
    }
  }
  
  Builder* BVH4Builder2ObjectSplit4 (void* accel, BuildSource* source, void* geometry, const size_t minLeafSize, const size_t maxLeafSize) {
    return new BVH4Builder2((BVH4*)accel,source,geometry,minLeafSize,maxLeafSize);
  }
}
