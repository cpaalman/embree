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

#define BVH4HAIR_COMPRESS_ALIGNED_NODES 0
#define BVH4HAIR_COMPRESS_UNALIGNED_NODES 0

#include "bvh4.h"
#include "bvh4_builder_full.h"
#include "bvh4_statistics.h"
#include "common/scene_bezier_curves.h"
#include "geometry/triangle4.h"

namespace embree
{
  extern double g_hair_builder_replication_factor;
  
  BVH4BuilderFull::BVH4BuilderFull (BVH4* bvh, Scene* scene)
    : scene(scene), minLeafSize(1), maxLeafSize(inf), bvh(bvh), remainingReplications(0)
  {
    if (BVH4::maxLeafBlocks < this->maxLeafSize) 
      this->maxLeafSize = BVH4::maxLeafBlocks;

    enableAlignedObjectSplits = false;
    enableAlignedSpatialSplits = false;
    enableUnalignedObjectSplits = false;
    enableUnalignedSpatialSplits = false;
    enableStrandSplits = false;
    enablePreSubdivision = 0;
    
    for (size_t i=0; i<g_hair_accel_mode.size();)
    {
      if      (g_hair_accel_mode.substr(i,2) == "P0" ) { enablePreSubdivision = 0; i+=2; } 
      else if (g_hair_accel_mode.substr(i,2) == "P1" ) { enablePreSubdivision = 1; i+=2; } 
      else if (g_hair_accel_mode.substr(i,2) == "P2" ) { enablePreSubdivision = 2; i+=2; } 
      else if (g_hair_accel_mode.substr(i,2) == "P3" ) { enablePreSubdivision = 3; i+=2; } 
      else if (g_hair_accel_mode.substr(i,2) == "P4" ) { enablePreSubdivision = 4; i+=2; } 
      else if (g_hair_accel_mode.substr(i,2) == "aO" ) { enableAlignedObjectSplits = true; i+=2; } 
      else if (g_hair_accel_mode.substr(i,2) == "uO" ) { enableUnalignedObjectSplits = true; i+=2; } 
      else if (g_hair_accel_mode.substr(i,3) == "auO" ) { enableAlignedObjectSplits = enableUnalignedObjectSplits = true; i+=3; } 
      else if (g_hair_accel_mode.substr(i,3) == "uST") { enableStrandSplits = true; i+=3; } 
      else if (g_hair_accel_mode.substr(i,3) == "aSP") { enableAlignedSpatialSplits = true; i+=3; } 
      else if (g_hair_accel_mode.substr(i,3) == "uSP") { enableUnalignedSpatialSplits = true; i+=3; } 
      else if (g_hair_accel_mode.substr(i,4) == "auSP") { enableAlignedSpatialSplits = enableUnalignedSpatialSplits = true; i+=4; } 
      else throw std::runtime_error("invalid hair accel mode");
    }
  }

  void BVH4BuilderFull::build(size_t threadIndex, size_t threadCount) 
  {
    /* fast path for empty BVH */
    size_t numPrimitives = scene->numCurves << enablePreSubdivision;
    bvh->init(numPrimitives+(size_t)(g_hair_builder_replication_factor*numPrimitives));
    if (numPrimitives == 0) return;
    numGeneratedPrims = 0;
    numAlignedObjectSplits = 0;
    numAlignedSpatialSplits = 0;
    numUnalignedObjectSplits = 0;
    numUnalignedSpatialSplits = 0;
    numStrandSplits = 0;
    numFallbackSplits = 0;

    double t0 = 0.0;
    if (g_verbose >= 2) 
    {
      PRINT(enableAlignedObjectSplits);
      PRINT(enableAlignedSpatialSplits);
      PRINT(enableUnalignedObjectSplits);
      PRINT(enableUnalignedSpatialSplits);
      PRINT(enableStrandSplits);
      PRINT(enablePreSubdivision);

      std::cout << "building " << bvh->name() + " using BVH4BuilderFull ..." << std::flush;
      t0 = getSeconds();
    }

    size_t N = 0;
    float r = 0;

    /* create initial curve list */
    BBox3fa bounds = empty;
    size_t numTriangles = 0;
    size_t numVertices = 0;
    TriRefList tris;
    BezierRefList beziers;
    for (size_t i=0; i<scene->size(); i++) 
    {
      Geometry* geom = scene->get(i);
      if (!geom->isEnabled()) continue;
      if (geom->type == BEZIER_CURVES) {
        BezierCurves* set = (BezierCurves*) geom;
        numVertices += set->numVertices;
        for (size_t j=0; j<set->numCurves; j++) {
          const int ofs = set->curve(j);
          const Vec3fa& p0 = set->vertex(ofs+0);
          const Vec3fa& p1 = set->vertex(ofs+1);
          const Vec3fa& p2 = set->vertex(ofs+2);
          const Vec3fa& p3 = set->vertex(ofs+3);
          const Bezier1 bezier(p0,p1,p2,p3,0,1,i,j);
          //bounds.extend(subdivideAndAdd(threadIndex,prims,bezier,enablePreSubdivision));
          bounds.extend(bezier.bounds());
          BezierRefList::item* block = beziers.head();
          if (block == NULL || !block->insert(bezier)) {
            block = beziers.insert(allocBezierRefs.malloc(threadIndex));
            block->insert(bezier);
          }
        }
      }
      if (geom->type == TRIANGLE_MESH) 
      {
        TriangleMesh* set = (TriangleMesh*) geom;
        if (set->numTimeSteps != 1) continue;
        numTriangles += set->numTriangles;
        numVertices += set->numVertices;
        for (size_t j=0; j<set->numTriangles; j++) {
          const TriangleMesh::Triangle tri = set->triangle(j);
          const Vec3fa& p0 = set->vertex(tri.v[0]);
          const Vec3fa& p1 = set->vertex(tri.v[1]);
          const Vec3fa& p2 = set->vertex(tri.v[2]);
          const Triangle1v triangle(p0,p1,p2,i,j,0);
          //bounds.extend(subdivideAndAdd(threadIndex,prims,bezier,enablePreSubdivision));
          bounds.extend(triangle.bounds());
          TriRefList::item* block = tris.head();
          if (block == NULL || !block->insert(triangle)) {
            block = tris.insert(allocTriRefs.malloc(threadIndex));
            block->insert(triangle);
          }
        }
      }
    }

    bvh->numPrimitives = scene->numCurves + numTriangles;
    bvh->numVertices = 0;
    if (bvh->primTy[0] == &SceneBezier1i::type) bvh->numVertices = numVertices; // FIXME

    /* start recursive build */
    remainingReplications = g_hair_builder_replication_factor*numPrimitives;
    //const NAABBox3fa ubounds = computeUnalignedBounds(prims);
    const BBox3fa abounds = computeAlignedBounds(tris,beziers);
    BuildTask task(&bvh->root,0,numPrimitives,false,tris,beziers,abounds);
    bvh->bounds = bounds;

#if 1
    recurseTask(threadIndex,task);
#else
    numActiveTasks = 1;
    tasks.push_back(task);
    push_heap(tasks.begin(),tasks.end());
    TaskScheduler::executeTask(threadIndex,threadCount,_task_build_parallel,this,threadCount,"BVH4Builder2::build_parallel");
#endif
    
    if (g_verbose >= 2) {
      double t1 = getSeconds();
      std::cout << " [DONE]" << std::endl;
      std::cout << "  dt = " << 1000.0f*(t1-t0) << "ms, perf = " << 1E-6*double(numPrimitives)/(t1-t0) << " Mprim/s" << std::endl;
      PRINT(numAlignedObjectSplits);
      PRINT(numAlignedSpatialSplits);
      PRINT(numUnalignedObjectSplits);
      PRINT(numUnalignedSpatialSplits);
      PRINT(numStrandSplits);
      PRINT(numFallbackSplits);
      std::cout << BVH4Statistics(bvh).str();
    }
  }

#if 0
  void BVH4BuilderFull::task_build_parallel(size_t threadIndex, size_t threadCount, size_t taskIndex, size_t taskCount, TaskScheduler::Event* event) 
  {
    while (numActiveTasks) 
    {
      taskMutex.lock();
      if (tasks.size() == 0) {
        taskMutex.unlock();
        continue;
      }

      /* take next task from heap */
      BuildTask task = tasks.front();
      pop_heap(tasks.begin(),tasks.end());
      tasks.pop_back();
      taskMutex.unlock();

      /* recursively finish task */
      if (task.size < 512) {
        atomic_add(&numActiveTasks,-1);
        recurseTask(threadIndex,task);
      }
      
      /* execute task and add child tasks */
      else 
      {
        size_t numChildren;
        BuildTask ctasks[BVH4::N];
        processTask(threadIndex,task,ctasks,numChildren);
        taskMutex.lock();
        for (size_t i=0; i<numChildren; i++) {
          atomic_add(&numActiveTasks,+1);
          tasks.push_back(ctasks[i]);
          push_heap(tasks.begin(),tasks.end());
        }
        atomic_add(&numActiveTasks,-1);
        taskMutex.unlock();
      }
    }
  }
#endif

  void BVH4BuilderFull::insert(size_t threadIndex, TriRefList& prims_i, TriRefList& prims_o)
  {
    while (TriRefList::item* block = prims_i.take()) {
      if (block->size()) prims_o.insert(block);
      else allocTriRefs.free(threadIndex,block);
    }
  }

  void BVH4BuilderFull::insert(size_t threadIndex, BezierRefList& prims_i, BezierRefList& prims_o)
  {
    while (BezierRefList::item* block = prims_i.take()) {
      if (block->size()) prims_o.insert(block);
      else allocBezierRefs.free(threadIndex,block);
    }
  }

  template<typename Left>
  void BVH4BuilderFull::split(size_t threadIndex, 
                              TriRefList& prims, const Left& left, 
                              TriRefList& lprims_o, size_t& lnum_o, 
                              TriRefList& rprims_o, size_t& rnum_o)
  {
    //lnum_o = rnum_o = 0;
    TriRefList::item* lblock = lprims_o.insert(allocTriRefs.malloc(threadIndex));
    TriRefList::item* rblock = rprims_o.insert(allocTriRefs.malloc(threadIndex));
    
    while (TriRefList::item* block = prims.take()) 
    {
      for (size_t i=0; i<block->size(); i++) 
      {
        const Triangle1v& prim = block->at(i); 
        if (left(prim)) 
        {
          lnum_o++;
          if (likely(lblock->insert(prim))) continue; 
          lblock = lprims_o.insert(allocTriRefs.malloc(threadIndex));
          lblock->insert(prim);
        } 
        else 
        {
          rnum_o++;
          if (likely(rblock->insert(prim))) continue;
          rblock = rprims_o.insert(allocTriRefs.malloc(threadIndex));
          rblock->insert(prim);
        }
      }
      allocTriRefs.free(threadIndex,block);
    }
  }

  template<typename Left>
  void BVH4BuilderFull::split(size_t threadIndex, 
                              BezierRefList& prims, const Left& left, 
                              BezierRefList& lprims_o, size_t& lnum_o, 
                              BezierRefList& rprims_o, size_t& rnum_o)
  {
    //lnum_o = rnum_o = 0;
    BezierRefList::item* lblock = lprims_o.insert(allocBezierRefs.malloc(threadIndex));
    BezierRefList::item* rblock = rprims_o.insert(allocBezierRefs.malloc(threadIndex));
    
    while (BezierRefList::item* block = prims.take()) 
    {
      for (size_t i=0; i<block->size(); i++) 
      {
        const Bezier1& prim = block->at(i); 
        if (left(prim)) 
        {
          lnum_o++;
          if (likely(lblock->insert(prim))) continue; 
          lblock = lprims_o.insert(allocBezierRefs.malloc(threadIndex));
          lblock->insert(prim);
        } 
        else 
        {
          rnum_o++;
          if (likely(rblock->insert(prim))) continue;
          rblock = rprims_o.insert(allocBezierRefs.malloc(threadIndex));
          rblock->insert(prim);
        }
      }
      allocBezierRefs.free(threadIndex,block);
    }
  }

  const BBox3fa BVH4BuilderFull::computeAlignedBounds(TriRefList& tris)
  {
    BBox3fa bounds = empty;
    for (TriRefList::block_iterator_unsafe i=tris; i; i++)
      bounds.extend(i->bounds());
    return bounds;
  }

  const BBox3fa BVH4BuilderFull::computeAlignedBounds(BezierRefList& beziers)
  {
    BBox3fa bounds = empty;
    for (BezierRefList::block_iterator_unsafe i=beziers; i; i++)
      bounds.extend(i->bounds());
    return bounds;
  }

  const BBox3fa BVH4BuilderFull::computeAlignedBounds(TriRefList& tris, BezierRefList& beziers) {
    return merge(computeAlignedBounds(tris),computeAlignedBounds(beziers));
  }

  const NAABBox3fa BVH4BuilderFull::computeAlignedBounds(TriRefList& tris, const LinearSpace3fa& space)
  {
    BBox3fa bounds = empty;
    for (TriRefList::block_iterator_unsafe i=tris; i; i++)
      bounds.extend(i->bounds(space));
    return NAABBox3fa(space,bounds);
  }

  const NAABBox3fa BVH4BuilderFull::computeAlignedBounds(BezierRefList& beziers, const LinearSpace3fa& space)
  {
    BBox3fa bounds = empty;
    for (BezierRefList::block_iterator_unsafe i=beziers; i; i++)
      bounds.extend(i->bounds(space));
    return NAABBox3fa(space,bounds);
  }

  const NAABBox3fa BVH4BuilderFull::computeAlignedBounds(TriRefList& tris, BezierRefList& beziers, const LinearSpace3fa& space) {
    return NAABBox3fa(space,merge(computeAlignedBounds(tris,space).bounds,computeAlignedBounds(beziers,space).bounds));
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////

  __forceinline BVH4BuilderFull::ObjectSplit BVH4BuilderFull::ObjectSplit::find(size_t threadIndex, size_t depth, BVH4BuilderFull* parent, 
                                                                                TriRefList& tris, BezierRefList& beziers, const LinearSpace3fa& space)
  {
    /* calculate geometry and centroid bounds */
    BBox3fa centBounds = empty;
    BBox3fa geomBounds = empty;

    for (TriRefList::block_iterator_unsafe i=tris; i; i++) {
      geomBounds.extend(i->bounds(space));
      centBounds.extend(i->center(space));
    }

    for (BezierRefList::block_iterator_unsafe i=beziers; i; i++) {
      geomBounds.extend(i->bounds(space));
      centBounds.extend(i->center(space));
    }
    
    /* calculate binning function */
    const ssef ofs  = (ssef) centBounds.lower;
    const ssef diag = (ssef) centBounds.size();
    const ssef scale = select(diag != 0.0f,rcp(diag) * ssef(BINS * 0.99f),ssef(0.0f));

    /* initialize bins */
    BBox3fa bounds[BINS][4];
    ssei    counts[BINS];
    for (size_t i=0; i<BINS; i++) {
      bounds[i][0] = bounds[i][1] = bounds[i][2] = bounds[i][3] = empty;
      counts[i] = 0;
    }
 
    /* perform binning of triangles */
    for (TriRefList::block_iterator_unsafe i=tris; i; i++)
    {
      const BBox3fa cbounds = i->bounds(space);
      const Vec3fa  center  = i->center(space);
      const ssei bin = clamp(floori((ssef(center) - ofs)*scale),ssei(0),ssei(BINS-1));
      //const ssei bin = floori((ssef(center) - ofs)*scale);
      assert(bin[0] >=0 && bin[0] < BINS);
      assert(bin[1] >=0 && bin[1] < BINS);
      assert(bin[2] >=0 && bin[2] < BINS);
      const int b0 = bin[0]; counts[b0][0]++; bounds[b0][0].extend(cbounds);
      const int b1 = bin[1]; counts[b1][1]++; bounds[b1][1].extend(cbounds);
      const int b2 = bin[2]; counts[b2][2]++; bounds[b2][2].extend(cbounds);
    }

    /* perform binning of bezier curves */
    for (BezierRefList::block_iterator_unsafe i=beziers; i; i++)
    {
      const BBox3fa cbounds = i->bounds(space);
      const Vec3fa  center  = i->center(space);
      const ssei bin = clamp(floori((ssef(center) - ofs)*scale),ssei(0),ssei(BINS-1));
      //const ssei bin = floori((ssef(center) - ofs)*scale);
      assert(bin[0] >=0 && bin[0] < BINS);
      assert(bin[1] >=0 && bin[1] < BINS);
      assert(bin[2] >=0 && bin[2] < BINS);
      const int b0 = bin[0]; counts[b0][0]++; bounds[b0][0].extend(cbounds);
      const int b1 = bin[1]; counts[b1][1]++; bounds[b1][1].extend(cbounds);
      const int b2 = bin[2]; counts[b2][2]++; bounds[b2][2].extend(cbounds);
    }
    
    /* sweep from right to left and compute parallel prefix of merged bounds */
    ssef rAreas[BINS];
    ssei rCounts[BINS];
    ssei count = 0; BBox3fa bx = empty; BBox3fa by = empty; BBox3fa bz = empty;
    for (size_t i=BINS-1; i>0; i--)
    {
      count += counts[i];
      rCounts[i] = count;
      bx.extend(bounds[i][0]); rAreas[i][0] = area(bx);
      by.extend(bounds[i][1]); rAreas[i][1] = area(by);
      bz.extend(bounds[i][2]); rAreas[i][2] = area(bz);
    }
    
    /* sweep from left to right and compute SAH */
    ssei ii = 1; ssef bestSAH = pos_inf; ssei bestPos = 0; ssei bestLeft = 0; ssei bestRight = 0;
    count = 0; bx = empty; by = empty; bz = empty;
    for (size_t i=1; i<BINS; i++, ii+=1)
    {
      count += counts[i-1];
      bx.extend(bounds[i-1][0]); float Ax = area(bx);
      by.extend(bounds[i-1][1]); float Ay = area(by);
      bz.extend(bounds[i-1][2]); float Az = area(bz);
      const ssef lArea = ssef(Ax,Ay,Az,Az);
      const ssef rArea = rAreas[i];
      const ssei lCount = (count     +ssei(3)) >> 2;
      const ssei rCount = (rCounts[i]+ssei(3)) >> 2;
      const ssef sah = lArea*ssef(lCount) + rArea*ssef(rCount);
      bestPos = select(sah < bestSAH,ii ,bestPos);
      bestLeft= select(sah < bestSAH,count,bestLeft);
      bestRight=select(sah < bestSAH,rCounts[i],bestRight);
      bestSAH = select(sah < bestSAH,sah,bestSAH);
    }
    
    /* find best dimension */
    ObjectSplit split;
    split.space = space;
    split.ofs = ofs;
    split.scale = scale;
    split.cost = inf;

    for (size_t dim=0; dim<3; dim++) 
    {
      /* ignore zero sized dimensions */
      if (unlikely(scale[dim] == 0.0f)) 
        continue;
      
      /* test if this is a better dimension */
      if (bestSAH[dim] < split.cost && bestPos[dim] != 0) {
        split.dim = dim;
        split.pos = bestPos[dim];
        split.cost = bestSAH[dim];
        split.num0 = bestLeft[dim];
        split.num1 = bestRight[dim];
      }
    }

    if (split.dim == -1) {
      split.num0 = split.num1 = 1;
      split.bounds0 = split.bounds1 = BBox3fa(inf);
      return split;
    }
    
    size_t lnum = 0, rnum = 0;
    TriRefList ltris, rtris; 
    BezierRefList lbeziers, rbeziers; 
    parent->split(threadIndex,tris,split,ltris,lnum,rtris,rnum); // FIXME: compute bounds from binning info
    parent->split(threadIndex,beziers,split,lbeziers,lnum,rbeziers,rnum); // FIXME: compute bounds from binning info
    split.bounds0 = computeAlignedBounds(ltris,lbeziers,space);
    split.bounds1 = computeAlignedBounds(rtris,rbeziers,space);
    parent->insert(threadIndex,ltris,tris);
    parent->insert(threadIndex,rtris,tris);
    parent->insert(threadIndex,lbeziers,beziers);
    parent->insert(threadIndex,rbeziers,beziers);
    return split;
  }

  __forceinline void BVH4BuilderFull::ObjectSplit::split(size_t threadIndex, BVH4BuilderFull* parent, 
                                                         TriRefList& tris, BezierRefList& beziers, 
                                                         TriRefList& ltris_o, BezierRefList& lbeziers_o,
                                                         TriRefList& rtris_o, BezierRefList& rbeziers_o) const
  {
    size_t lnum = 0, rnum = 0;
    parent->split(threadIndex,tris   ,*this,ltris_o   ,lnum,rtris_o   ,rnum);
    parent->split(threadIndex,beziers,*this,lbeziers_o,lnum,rbeziers_o,rnum);
    assert(lnum == num0);
    assert(rnum == num1);
  }
  

  __forceinline BVH4BuilderFull::FallBackSplit BVH4BuilderFull::FallBackSplit::find(size_t threadIndex, BVH4BuilderFull* parent, 
                                                                                    TriRefList& tris,    BezierRefList& beziers, 
                                                                                    TriRefList& ltris_o, BezierRefList& lbeziers_o,
                                                                                    TriRefList& rtris_o, BezierRefList& rbeziers_o)
  {
    size_t num = 0, lnum = 0, rnum = 0;
    {
      TriRefList::item* lblock = ltris_o.insert(parent->allocTriRefs.malloc(threadIndex));
      TriRefList::item* rblock = rtris_o.insert(parent->allocTriRefs.malloc(threadIndex));
    
      while (TriRefList::item* block = tris.take()) 
      {
        for (size_t i=0; i<block->size(); i++) 
        {
          const Triangle1v& prim = block->at(i); 
          if ((num++)%2) 
          {
            lnum++;
            if (likely(lblock->insert(prim))) continue; 
            lblock = ltris_o.insert(parent->allocTriRefs.malloc(threadIndex));
            lblock->insert(prim);
          } 
          else 
          {
            rnum++;
            if (likely(rblock->insert(prim))) continue;
            rblock = rtris_o.insert(parent->allocTriRefs.malloc(threadIndex));
            rblock->insert(prim);
          }
        }
        parent->allocTriRefs.free(threadIndex,block);
      }
    }

    {
      BezierRefList::item* lblock = lbeziers_o.insert(parent->allocBezierRefs.malloc(threadIndex));
      BezierRefList::item* rblock = rbeziers_o.insert(parent->allocBezierRefs.malloc(threadIndex));
    
      while (BezierRefList::item* block = beziers.take()) 
      {
        for (size_t i=0; i<block->size(); i++) 
        {
          const Bezier1& prim = block->at(i); 
          if ((num++)%2) 
          {
            lnum++;
            if (likely(lblock->insert(prim))) continue; 
            lblock = lbeziers_o.insert(parent->allocBezierRefs.malloc(threadIndex));
            lblock->insert(prim);
          } 
          else 
          {
            rnum++;
            if (likely(rblock->insert(prim))) continue;
            rblock = rbeziers_o.insert(parent->allocBezierRefs.malloc(threadIndex));
            rblock->insert(prim);
          }
        }
        parent->allocBezierRefs.free(threadIndex,block);
      }
    }
    const NAABBox3fa bounds0 = computeAlignedBounds(ltris_o,lbeziers_o);
    const NAABBox3fa bounds1 = computeAlignedBounds(rtris_o,rbeziers_o);
    return FallBackSplit(bounds0,lnum,bounds1,rnum);
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////

  BVH4::NodeRef BVH4BuilderFull::leaf(size_t threadIndex, size_t depth, TriRefList& prims, const NAABBox3fa& bounds)
  {
    //size_t N = end-begin;
    size_t N = bvh->primTy[0]->blocks(TriRefList::block_iterator_unsafe(prims).size());

    if (N > (size_t)BVH4::maxLeafBlocks) {
      //std::cout << "WARNING: Loosing " << N-BVH4::maxLeafBlocks << " primitives during build!" << std::endl;
      std::cout << "!" << std::flush;
      N = (size_t)BVH4::maxLeafBlocks;
    }
    size_t numGeneratedPrimsOld = atomic_add(&numGeneratedPrims,N); 
    if (numGeneratedPrimsOld%10000 > (numGeneratedPrimsOld+N)%10000) std::cout << "." << std::flush; 
    //assert(N <= (size_t)BVH4::maxLeafBlocks);

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
          const Triangle1v& prim = *iter;
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

  BVH4::NodeRef BVH4BuilderFull::leaf(size_t threadIndex, size_t depth, BezierRefList& prims, const NAABBox3fa& bounds)
  {
    //size_t N = end-begin;
    size_t N = BezierRefList::block_iterator_unsafe(prims).size();

    if (N > (size_t)BVH4::maxLeafBlocks) {
      //std::cout << "WARNING: Loosing " << N-BVH4::maxLeafBlocks << " primitives during build!" << std::endl;
      std::cout << "!" << std::flush;
      N = (size_t)BVH4::maxLeafBlocks;
    }
    size_t numGeneratedPrimsOld = atomic_add(&numGeneratedPrims,N); 
    if (numGeneratedPrimsOld%10000 > (numGeneratedPrimsOld+N)%10000) std::cout << "." << std::flush; 
    //assert(N <= (size_t)BVH4::maxLeafBlocks);
    if (bvh->primTy[1] == &Bezier1Type::type) { 
      Bezier1* leaf = (Bezier1*) bvh->allocPrimitiveBlocks(threadIndex,0,N);
      BezierRefList::block_iterator_unsafe iter(prims);
      for (size_t i=0; i<N; i++) { leaf[i] = *iter; iter++; }
      assert(!iter);

      /* free all primitive blocks */
      while (BezierRefList::item* block = prims.take())
        allocBezierRefs.free(threadIndex,block);

      return bvh->encodeLeaf((char*)leaf,N,1);
    } 
    else if (bvh->primTy[1] == &SceneBezier1i::type) {
      Bezier1i* leaf = (Bezier1i*) bvh->allocPrimitiveBlocks(threadIndex,0,N);
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

  BVH4::NodeRef BVH4BuilderFull::leaf(size_t threadIndex, size_t depth, TriRefList& tris, BezierRefList& beziers, const NAABBox3fa& bounds)
  {
    size_t Ntris    = TriRefList   ::block_iterator_unsafe(tris   ).size();
    size_t Nbeziers = BezierRefList::block_iterator_unsafe(beziers).size();
    if (Ntris == 0) 
      return leaf(threadIndex,depth,beziers,bounds);
    else if (Nbeziers == 0)
      return leaf(threadIndex,depth,tris,bounds);
    
    BVH4::UANode* node = bvh->allocUANode(threadIndex);
    node->set(0,computeAlignedBounds(tris   ),leaf(threadIndex,depth+1,tris   ,bounds));
    node->set(1,computeAlignedBounds(beziers),leaf(threadIndex,depth+1,beziers,bounds));
    return bvh->encodeNode(node);
  }

  bool BVH4BuilderFull::split(size_t threadIndex, size_t depth, 
                              TriRefList& tris,    BezierRefList& beziers, const NAABBox3fa& bounds, size_t size,
                              TriRefList& ltris_o, BezierRefList& lbeziers_o, size_t& lsize,
                              TriRefList& rtris_o, BezierRefList& rbeziers_o, size_t& rsize,
                              bool& isAligned)
  {
    /* variable to track the SAH of the best splitting approach */
    bool enableSpatialSplits = remainingReplications > 0;
    const int travCostAligned = isAligned ? BVH4::travCostAligned : BVH4::travCostUnaligned;
    const float leafSAH = BVH4::intCost*float(size)*embree::area(bounds.bounds);
    float bestSAH = leafSAH;
    
    /* perform standard binning in aligned space */
    ObjectSplit alignedObjectSplit = ObjectSplit::find(threadIndex,depth,this,tris,beziers,one);
    float alignedObjectSAH = travCostAligned*embree::area(bounds.bounds) + alignedObjectSplit.standardSAH();
    bestSAH = min(bestSAH,alignedObjectSAH);

    /* perform fallback split */
    if (bestSAH == float(inf)) {
      if (size <= maxLeafSize) return false;
      numFallbackSplits++;
      const FallBackSplit fallbackSplit = FallBackSplit::find(threadIndex,this,tris,beziers,ltris_o,lbeziers_o,rtris_o,rbeziers_o);
      lsize = fallbackSplit.num0;
      rsize = fallbackSplit.num1;
      return true;
    }

    /* perform aligned object split */
    else if (bestSAH == alignedObjectSAH) {
      numAlignedObjectSplits++;
      alignedObjectSplit.split(threadIndex,this,tris,beziers,ltris_o,lbeziers_o,rtris_o,rbeziers_o);
      lsize = alignedObjectSplit.num0;
      rsize = alignedObjectSplit.num1;
      return true;
    }
 
    else {
      throw std::runtime_error("bvh4hair_builder: internal error");
      return true;
    }
  }

  void BVH4BuilderFull::processTask(size_t threadIndex, BuildTask& task, BuildTask task_o[BVH4::N], size_t& numTasks_o)
  {
    /* create enforced leaf */
    if (task.size <= minLeafSize || task.depth >= BVH4::maxBuildDepth || task.makeleaf) {
      *task.dst = leaf(threadIndex,task.depth,task.tris,task.beziers,task.bounds);
      numTasks_o = 0;
      return;
    }

    /*! initialize child list */
    bool isAligned = true;
    NAABBox3fa cbounds[BVH4::N];
    TriRefList ctris[BVH4::N];
    BezierRefList cbeziers[BVH4::N];
    size_t csize[BVH4::N];
    bool isleaf[BVH4::N];
    ctris[0] = task.tris;
    cbeziers[0] = task.beziers;
    cbounds[0] = task.bounds;
    csize[0] = task.size;
    isleaf[0] = false;
    size_t numChildren = 1;
    
    /*! split until node is full or SAH tells us to stop */
    do {
      
      /*! find best child to split */
      float bestArea = neg_inf; 
      ssize_t bestChild = -1;
      for (size_t i=0; i<numChildren; i++) 
      {
        size_t Ntris    = TriRefList::block_iterator_unsafe(ctris[i]).size(); // FIXME: slow
        size_t Nbeziers = BezierRefList::block_iterator_unsafe(cbeziers[i]).size(); // FIXME: slow
        size_t N = Ntris + Nbeziers;
        float A = embree::area(cbounds[i].bounds);
        if (N <= minLeafSize) continue;  
        if (isleaf[i]) continue;
        if (A > bestArea) { bestChild = i; bestArea = A; }
      }
      if (bestChild == -1) break;

      /*! split selected child */
      size_t lsize, rsize;
      TriRefList ltris, rtris;
      BezierRefList lbeziers, rbeziers;
      bool done = split(threadIndex,task.depth,ctris[bestChild],cbeziers[bestChild],cbounds[bestChild],csize[bestChild],
                        ltris,lbeziers,lsize,rtris,rbeziers,rsize,isAligned);
      if (!done) { isleaf[bestChild] = true; continue; }
      ctris[numChildren] = rtris; cbeziers[numChildren] = rbeziers; isleaf[numChildren] = false; csize[numChildren] = rsize;
      ctris[bestChild  ] = ltris; cbeziers[bestChild  ] = lbeziers; isleaf[bestChild  ] = false; csize[bestChild  ] = lsize;

      //cbounds[numChildren] = computeUnalignedBounds(cprims[numChildren]);
      //cbounds[bestChild  ] = computeUnalignedBounds(cprims[bestChild  ]);
      cbounds[numChildren] = computeAlignedBounds(ctris[numChildren],cbeziers[numChildren]);
      cbounds[bestChild  ] = computeAlignedBounds(ctris[bestChild  ],cbeziers[bestChild  ]);
      numChildren++;
      
    } while (numChildren < BVH4::N);

    /* create aligned node */
    //if (isAligned) 
    {
      BVH4::UANode* node = bvh->allocUANode(threadIndex);

      BBox3fa bounds = empty;
      NAABBox3fa abounds[BVH4::N];
      for (size_t i=0; i<numChildren; i++) {
        abounds[i] = computeAlignedBounds(ctris[i],cbeziers[i]);
        bounds.extend(abounds[i].bounds);
      }

      for (ssize_t i=0; i<numChildren; i++) {
        node->set(i,abounds[i].bounds);
	//const NAABBox3fa ubounds = computeUnalignedBounds(cprims[i]);
        const NAABBox3fa ubounds = computeAlignedBounds(ctris[i],cbeziers[i]);
        new (&task_o[i]) BuildTask(&node->child(i),task.depth+1,csize[i],isleaf[i],ctris[i],cbeziers[i],ubounds);
      }
      numTasks_o = numChildren;
      *task.dst = bvh->encodeNode(node);
    }
    
    /* create unaligned node */
    /*else 
    {
      BVH4::UUNode* node = bvh->allocUUNode(threadIndex);
      for (ssize_t i=numChildren-1; i>=0; i--) {
        node->set(i,cbounds[i]);
        new (&task_o[i]) BuildTask(&node->child(i),task.depth+1,csize[i],isleaf[i],cprims[i],cbounds[i]);
      }
      numTasks_o = numChildren;
      *task.dst = bvh->encodeNode(node);
      }*/
  }

  void BVH4BuilderFull::recurseTask(size_t threadIndex, BuildTask& task)
  {
    size_t numChildren;
    BuildTask tasks[BVH4::N];
    processTask(threadIndex,task,tasks,numChildren);
    for (size_t i=0; i<numChildren; i++) 
      recurseTask(threadIndex,tasks[i]);
  }

  Builder* BVH4BuilderFull_ (BVH4* accel, Scene* scene) {
    return new BVH4BuilderFull(accel,scene);
  }
}
