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
#include "bvh4_builder_hair.h"
#include "bvh4_statistics.h"
#include "common/scene_bezier_curves.h"

namespace embree
{
  extern double g_hair_builder_replication_factor;
  
  /*! scales orthonormal transformation into the range -127 to +127 */
  __forceinline const LinearSpace3fa compressTransform(const LinearSpace3fa& xfm)
  {
#if BVH4HAIR_COMPRESS_UNALIGNED_NODES
    assert(xfm.vx.x >= -1.0f && xfm.vx.x <= 1.0f);
    assert(xfm.vx.y >= -1.0f && xfm.vx.y <= 1.0f);
    assert(xfm.vx.z >= -1.0f && xfm.vx.z <= 1.0f);
    assert(xfm.vy.x >= -1.0f && xfm.vy.x <= 1.0f);
    assert(xfm.vy.y >= -1.0f && xfm.vy.y <= 1.0f);
    assert(xfm.vy.z >= -1.0f && xfm.vy.z <= 1.0f);
    assert(xfm.vz.x >= -1.0f && xfm.vz.x <= 1.0f);
    assert(xfm.vz.y >= -1.0f && xfm.vz.y <= 1.0f);
    assert(xfm.vz.z >= -1.0f && xfm.vz.z <= 1.0f);
    return LinearSpace3fa (clamp(trunc(127.0f*xfm.vx),Vec3fa(-127.0f),Vec3fa(+127.0f))/127.0f,
                           clamp(trunc(127.0f*xfm.vy),Vec3fa(-127.0f),Vec3fa(+127.0f))/127.0f,
                           clamp(trunc(127.0f*xfm.vz),Vec3fa(-127.0f),Vec3fa(+127.0f))/127.0f);
#else
    return xfm;
#endif
  }

  BVH4BuilderHair::BVH4BuilderHair (BVH4* bvh, Scene* scene)
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

  void BVH4BuilderHair::build(size_t threadIndex, size_t threadCount) 
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

      std::cout << "building " << bvh->name() + " using BVH4BuilderHair ..." << std::flush;
      t0 = getSeconds();
    }

    size_t N = 0;
    float r = 0;

    /* create initial curve list */
    BBox3fa bounds = empty;
    size_t numVertices = 0;
    atomic_set<PrimRefBlock> prims;
    for (size_t i=0; i<scene->size(); i++) 
    {
      Geometry* geom = scene->get(i);
      if (geom->type != BEZIER_CURVES) continue;
      if (!geom->isEnabled()) continue;
      BezierCurves* set = (BezierCurves*) geom;
      numVertices += set->numVertices;
      for (size_t j=0; j<set->numCurves; j++) {
        const int ofs = set->curve(j);
        const Vec3fa& p0 = set->vertex(ofs+0);
        const Vec3fa& p1 = set->vertex(ofs+1);
        const Vec3fa& p2 = set->vertex(ofs+2);
        const Vec3fa& p3 = set->vertex(ofs+3);
        const Bezier1 bezier(p0,p1,p2,p3,0,1,i,j);
        bounds.extend(subdivideAndAdd(threadIndex,prims,bezier,enablePreSubdivision));
      }
    }

    bvh->numPrimitives = scene->numCurves;
    bvh->numVertices = 0;
    if (bvh->primTy[0] == &SceneBezier1i::type) bvh->numVertices = numVertices;

    /* start recursive build */
    PRINT(g_hair_builder_replication_factor);
    remainingReplications = g_hair_builder_replication_factor*numPrimitives;
    const NAABBox3fa ubounds = computeUnalignedBounds(prims);
    BuildTask task(&bvh->root,0,numPrimitives,false,prims,ubounds);
    bvh->bounds = bounds;

#if 0
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

  const BBox3fa BVH4BuilderHair::subdivideAndAdd(size_t threadIndex, atomic_set<PrimRefBlock>& prims, const Bezier1& bezier, size_t depth)
  {
    if (depth == 0) {
      atomic_set<PrimRefBlock>::item* block = prims.head();
      if (block == NULL || !block->insert(bezier)) {
        block = prims.insert(alloc.malloc(threadIndex));
        block->insert(bezier);
      }
      return bezier.bounds();
    }

    Bezier1 bezier0,bezier1;
    bezier.subdivide(bezier0,bezier1);
    const BBox3fa bounds0 = subdivideAndAdd(threadIndex,prims,bezier0,depth-1);
    const BBox3fa bounds1 = subdivideAndAdd(threadIndex,prims,bezier1,depth-1);
    return merge(bounds0,bounds1);
  }

  void BVH4BuilderHair::task_build_parallel(size_t threadIndex, size_t threadCount, size_t taskIndex, size_t taskCount, TaskScheduler::Event* event) 
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

  void BVH4BuilderHair::insert(size_t threadIndex, atomic_set<PrimRefBlock>& prims_i, atomic_set<PrimRefBlock>& prims_o)
  {
    while (atomic_set<PrimRefBlock>::item* block = prims_i.take()) {
      if (block->size()) prims_o.insert(block);
      else alloc.free(threadIndex,block);
    }
  }

  template<typename Left>
  void BVH4BuilderHair::split(size_t threadIndex, atomic_set<PrimRefBlock>& prims, const Left& left, 
                               atomic_set<PrimRefBlock>& lprims_o, size_t& lnum_o, 
                               atomic_set<PrimRefBlock>& rprims_o, size_t& rnum_o)
  {
    lnum_o = rnum_o = 0;
    atomic_set<PrimRefBlock>::item* lblock = lprims_o.insert(alloc.malloc(threadIndex));
    atomic_set<PrimRefBlock>::item* rblock = rprims_o.insert(alloc.malloc(threadIndex));
    
    while (atomic_set<PrimRefBlock>::item* block = prims.take()) 
    {
      for (size_t i=0; i<block->size(); i++) 
      {
        const PrimRef& prim = block->at(i); 
        if (left(prim)) 
        {
          lnum_o++;
          if (likely(lblock->insert(prim))) continue; 
          lblock = lprims_o.insert(alloc.malloc(threadIndex));
          lblock->insert(prim);
        } 
        else 
        {
          rnum_o++;
          if (likely(rblock->insert(prim))) continue;
          rblock = rprims_o.insert(alloc.malloc(threadIndex));
          rblock->insert(prim);
        }
      }
      alloc.free(threadIndex,block);
    }
  }

  const BBox3fa BVH4BuilderHair::computeAlignedBounds(atomic_set<PrimRefBlock>& prims)
  {
    float area = 0.0f;
    BBox3fa bounds = empty;
    for (atomic_set<PrimRefBlock>::block_iterator_unsafe i = prims; i; i++)
    {
      const BBox3fa cbounds = i->bounds();
      area += embree::area(cbounds);
      bounds.extend(cbounds);
    }
    bounds.upper.w = area;
    return bounds;
  }

  const NAABBox3fa BVH4BuilderHair::computeAlignedBounds(atomic_set<PrimRefBlock>& prims, const LinearSpace3fa& space)
  {
    float area = 0.0f;
    BBox3fa bounds = empty;
    for (atomic_set<PrimRefBlock>::block_iterator_unsafe i = prims; i; i++)
    {
      const BBox3fa cbounds = i->bounds(space);
      area += embree::area(cbounds);
      bounds.extend(cbounds);
    }
    bounds.upper.w = area;
    return NAABBox3fa(space,bounds);
  }

  const NAABBox3fa BVH4BuilderHair::computeUnalignedBounds(atomic_set<PrimRefBlock>& prims)
  {
    size_t N = atomic_set<PrimRefBlock>::block_iterator_unsafe(prims).size();
    if (N == 0)
      return NAABBox3fa(empty); // FIXME: can cause problems with compression

    float bestArea = inf;
    LinearSpace3fa bestSpace = one;
    BBox3fa bestBounds = empty;

    size_t k=0;
    for (atomic_set<PrimRefBlock>::block_iterator_unsafe i = prims; i; i++)
    {
      if ((k++) % ((N+3)/4)) continue;
      //size_t k = begin + rand() % (end-begin);
      const Vec3fa axis = normalize(i->p3 - i->p0);
      if (length(i->p3 - i->p0) < 1E-9) continue;
      const LinearSpace3fa space0 = LinearSpace3fa::rotate(Vec3fa(0,0,1),2.0f*float(pi)*drand48())*frame(axis).transposed();
      const LinearSpace3fa space = compressTransform(clamp(space0));
      BBox3fa bounds = empty;
      float area = 0.0f;
      for (atomic_set<PrimRefBlock>::block_iterator_unsafe j = prims; j; j++) {
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
/*#ifdef DEBUG
    if (bestArea == (float)inf)
      {
        std::cout << "WARNING: bestArea == (float)inf" << std::endl; 
      }
      #endif*/

    bestBounds.upper.w = bestArea;
    return NAABBox3fa(bestSpace,bestBounds);
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////

  struct StrandSplitFunction
  {
    __forceinline StrandSplitFunction (const Vec3fa& axis0, const Vec3fa& axis1)
      : axis0(axis0), axis1(axis1) {}

    __forceinline bool operator() (const Bezier1& prim) const 
    {
      const Vec3fa axisi = normalize(prim.p3-prim.p0);
      const float cos0 = abs(dot(axisi,axis0));
      const float cos1 = abs(dot(axisi,axis1));
      return cos0 > cos1;
    }

    const Vec3fa axis0;
    const Vec3fa axis1;
  };

  __forceinline BVH4BuilderHair::StrandSplit::StrandSplit (const NAABBox3fa& bounds0, const Vec3fa& axis0, const size_t num0,
                                                           const NAABBox3fa& bounds1, const Vec3fa& axis1, const size_t num1)
    : bounds0(bounds0), bounds1(bounds1), axis0(axis0), axis1(axis1), num0(num0), num1(num1) {}
  
  __forceinline const BVH4BuilderHair::StrandSplit BVH4BuilderHair::StrandSplit::find(size_t threadIndex, BVH4BuilderHair* parent, atomic_set<PrimRefBlock>& prims)
  {
    /* first try to split two hair strands */
    atomic_set<PrimRefBlock>::block_iterator_unsafe i = prims;
    Vec3fa axis0 = normalize(i->p3 - i->p0);
    float bestCos = 1.0f;
    Bezier1 bestI = *i;

    for (i; i; i++) {
      Vec3fa axisi = i->p3 - i->p0;
      float leni = length(axisi);
      if (leni == 0.0f) continue;
      axisi /= leni;
      float cos = abs(dot(axisi,axis0));
      if (cos < bestCos) { bestCos = cos; bestI = *i; }
    }
    Vec3fa axis1 = normalize(bestI.p3-bestI.p0);

    /* partition the two strands */
    size_t num0, num1;
    atomic_set<PrimRefBlock> lprims, rprims; 
    parent->split(threadIndex,prims,StrandSplitFunction(axis0,axis1),lprims,num0,rprims,num1);

    NAABBox3fa naabb0(one,inf);
    NAABBox3fa naabb1(one,inf);
    if (num0 == 0 || num1 == 0) {
      num0 = num1 = 1;
    } 
    else {
      naabb0 = computeUnalignedBounds(lprims);
      naabb1 = computeUnalignedBounds(rprims);
    }

    /* merge lists again */
    parent->insert(threadIndex,lprims,prims);
    parent->insert(threadIndex,rprims,prims);

    return StrandSplit(naabb0,axis0,num0,naabb1,axis1,num1);
  }

  __forceinline void BVH4BuilderHair::StrandSplit::split(size_t threadIndex, BVH4BuilderHair* parent, 
                                                          atomic_set<PrimRefBlock>& prims, atomic_set<PrimRefBlock>& lprims_o, atomic_set<PrimRefBlock>& rprims_o) const 
  {
    size_t lnum,rnum;
    parent->split(threadIndex,prims,StrandSplitFunction(axis0,axis1),lprims_o,lnum,rprims_o,rnum);
    assert(lnum == num0);
    assert(rnum == num1);
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////

  __forceinline BVH4BuilderHair::ObjectSplit BVH4BuilderHair::ObjectSplit::find(size_t threadIndex, size_t depth, BVH4BuilderHair* parent, 
                                                                                  atomic_set<PrimRefBlock>& prims, const LinearSpace3fa& space)
  {
    /* calculate geometry and centroid bounds */
    BBox3fa centBounds = empty;
    BBox3fa geomBounds = empty;
    for (atomic_set<PrimRefBlock>::block_iterator_unsafe i = prims; i; i++) {
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
 
    /* perform binning of curves */
    for (atomic_set<PrimRefBlock>::block_iterator_unsafe i = prims; i; i++)
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
#if BVH4HAIR_WIDTH == 8
      const ssei lCount = (count     +ssei(7)) >> 3;
      const ssei rCount = (rCounts[i]+ssei(7)) >> 3;
#else
      const ssei lCount = (count     +ssei(3)) >> 2;
      const ssei rCount = (rCounts[i]+ssei(3)) >> 2;
#endif
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
    
    size_t lnum, rnum;
    atomic_set<PrimRefBlock> lprims, rprims; 
    parent->split(threadIndex,prims,split,lprims,lnum,rprims,rnum);
    split.bounds0 = computeAlignedBounds(lprims,space);
    split.bounds1 = computeAlignedBounds(rprims,space);
    parent->insert(threadIndex,lprims,prims);
    parent->insert(threadIndex,rprims,prims);
    return split;
  }

  __forceinline void BVH4BuilderHair::ObjectSplit::split(size_t threadIndex, BVH4BuilderHair* parent, 
                                                          atomic_set<PrimRefBlock>& prims, atomic_set<PrimRefBlock>& lprims_o, atomic_set<PrimRefBlock>& rprims_o) const
  {
    size_t lnum,rnum;
    parent->split(threadIndex,prims,*this,lprims_o,lnum,rprims_o,rnum);
    assert(lnum == num0);
    assert(rnum == num1);
  }
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////

  const BVH4BuilderHair::SpatialSplit BVH4BuilderHair::SpatialSplit::find(size_t threadIndex, size_t depth, size_t size, BVH4BuilderHair* parent, atomic_set<PrimRefBlock>& prims, const LinearSpace3fa& space)
  {
    /* calculate geometry and centroid bounds */
    BBox3fa geomBounds = empty;
    for (atomic_set<PrimRefBlock>::block_iterator_unsafe i = prims; i; i++)
      geomBounds.extend(i->bounds(space));

    /* calculate binning function */
    const ssef ofs  = (ssef) geomBounds.lower;
    const ssef diag = (ssef) geomBounds.size();
    const ssef scale = select(diag != 0.0f,rcp(diag) * ssef(BINS * 0.99f),ssef(0.0f));

    /* initialize bins */
    BBox3fa bounds[BINS][4];
    float   areas [BINS][4];
    ssei    numBegin[BINS];
    ssei    numEnd[BINS];
    for (size_t i=0; i<BINS; i++) {
      bounds[i][0] = bounds[i][1] = bounds[i][2] = bounds[i][3] = empty;
      areas [i][0] = areas [i][1] = areas [i][2] = areas [i][3] = 0.0f;
      numBegin[i] = numEnd[i] = 0;
    }
 
    /* perform binning of curves */
    for (atomic_set<PrimRefBlock>::block_iterator_unsafe i = prims; i; i++)
    {
      const Vec3fa v0 = xfmPoint(space,i->p0);
      const ssei bin0 = clamp(floori((ssef(v0)-ofs)*scale),ssei(0),ssei(BINS-1));
      //const ssei bin0 = floori((ssef(v0)-ofs)*scale);
      assert(bin0[0] >=0 && bin0[0] < BINS);
      assert(bin0[1] >=0 && bin0[1] < BINS);
      assert(bin0[2] >=0 && bin0[2] < BINS);
      const Vec3fa v1 = xfmPoint(space,i->p3);
      const ssei bin1 = clamp(floori((ssef(v1)-ofs)*scale),ssei(0),ssei(BINS-1));
      //const ssei bin1 = floori((ssef(v1)-ofs)*scale);
      assert(bin1[0] >=0 && bin1[0] < BINS);
      assert(bin1[1] >=0 && bin1[1] < BINS);
      assert(bin1[2] >=0 && bin1[2] < BINS);
      const ssei startbin = min(bin0,bin1);
      const ssei endbin   = max(bin0,bin1);

      for (size_t dim=0; dim<3; dim++) 
      {
        size_t bin;
        Bezier1 curve = *i;
        for (bin=startbin[dim]; bin<endbin[dim]; bin++) // FIXME: one can prevent many transformations in this loop here !!!
        {
          const float pos = float(bin+1)/scale[dim]+ofs[dim];
          const Vec3fa plane(space.vx[dim],space.vy[dim],space.vz[dim],-pos);
          Bezier1 bincurve,restcurve; 
          if (curve.split(plane,bincurve,restcurve)) {
            const BBox3fa cbounds = bincurve.bounds(space);
            bounds[bin][dim].extend(cbounds);
            areas [bin][dim] += embree::area(cbounds); // FIXME: not correct
            curve = restcurve;
          }
        }
        numBegin[startbin[dim]][dim]++;
        numEnd  [endbin  [dim]][dim]++;
        const BBox3fa cbounds = curve.bounds(space);
        bounds[bin][dim].extend(cbounds);
        areas [bin][dim] += embree::area(cbounds);  // FIXME: not correct
      }
    }
    
    /* sweep from right to left and compute parallel prefix of merged bounds */
    ssef rAreas[BINS];
    ssei rCounts[BINS];
    ssei count = 0; BBox3fa bx = empty; BBox3fa by = empty; BBox3fa bz = empty;
    for (size_t i=BINS-1; i>0; i--)
    {
      count += numEnd[i];
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
      count += numBegin[i-1];
      bx.extend(bounds[i-1][0]); float Ax = area(bx);
      by.extend(bounds[i-1][1]); float Ay = area(by);
      bz.extend(bounds[i-1][2]); float Az = area(bz);
      const ssef lArea = ssef(Ax,Ay,Az,Az);
      const ssef rArea = rAreas[i];
#if BVH4HAIR_WIDTH == 8
      const ssei lCount = (count     +ssei(7)) >> 3;
      const ssei rCount = (rCounts[i]+ssei(7)) >> 3;
#else
      const ssei lCount = (count     +ssei(3)) >> 2;
      const ssei rCount = (rCounts[i]+ssei(3)) >> 2;
#endif
      const ssef sah = lArea*ssef(lCount) + rArea*ssef(rCount);
      bestPos  = select(sah < bestSAH,ii ,bestPos);
      bestLeft = select(sah < bestSAH,count,bestLeft);
      bestRight= select(sah < bestSAH,rCounts[i],bestRight);
      bestSAH  = select(sah < bestSAH,sah,bestSAH);
    }
    
    /* find best dimension */
    SpatialSplit split;
    split.space = space;
    split.ofs = ofs;
    split.scale = scale;
    split.cost = inf;
    split.dim = -1;
    split.pos = 0.0f;
    split.num0 = split.num1 = 1;
    split.bounds0 = split.bounds1 = BBox3fa(inf);

    float bestCost = inf;
    for (size_t dim=0; dim<3; dim++) 
    {
      /* ignore zero sized dimensions */
      if (unlikely(scale[dim] == 0.0f)) 
        continue;
      
      /* test if this is a better dimension */
      if (bestSAH[dim] < bestCost && bestPos[dim] != 0) {
        split.dim = dim;
        split.pos = bestPos[dim]/scale[dim]+ofs[dim];
        split.cost = bestSAH[dim];
        split.num0 = bestLeft[dim];
        split.num1 = bestRight[dim];
        bestCost = bestSAH[dim];
      }
    }

    /* compute bounds of left and right side */
    if (split.dim == -1)
      return split;

    /* calculate bounding box of left and right side */
    BBox3fa lbounds = empty, rbounds = empty;
    float   larea   = 0.0f,  rarea   = 0.0f;
    size_t lnum = 0, rnum = 0;

    for (atomic_set<PrimRefBlock>::block_iterator_unsafe i = prims; i; i++)
    {
      const Vec3fa plane(space.vx[split.dim],space.vy[split.dim],space.vz[split.dim],-split.pos);
      
      const float p0p = dot(i->p0,plane)+plane.w;
      const float p3p = dot(i->p3,plane)+plane.w;

      /* sort to the left side */
      if (p0p <= 0.0f && p3p <= 0.0f) {
        const BBox3fa bounds = i->bounds(space);
        lbounds.extend(bounds);
        larea += embree::area(bounds);
        lnum++;
        continue;
      }
      
      /* sort to the right side */
      if (p0p >= 0.0f && p3p >= 0.0f) {
        const BBox3fa bounds = i->bounds(space);
        rbounds.extend(bounds);
        rarea += embree::area(bounds);
        rnum++;
        continue;
      }

      Bezier1 left,right; 
      if (i->split(plane,left,right)) {
        const BBox3fa lcbounds = left.bounds(space);
        const BBox3fa rcbounds = right.bounds(space);
        lbounds.extend(lcbounds); larea += embree::area(lcbounds); lnum++;
        rbounds.extend(rcbounds); rarea += embree::area(rcbounds); rnum++;
        continue;
      }
      
      const BBox3fa bounds = i->bounds(space);
      lbounds.extend(bounds); larea += embree::area(bounds); lnum++;
    }
    lbounds.upper.w = larea;
    rbounds.upper.w = rarea;
    split.bounds0 = NAABBox3fa(space,lbounds);
    split.bounds1 = NAABBox3fa(space,rbounds);
    //assert(lnum == split.num0);
    //assert(rnum == split.num1);
    split.num0 = lnum;
    split.num1 = rnum;
    split.numReplications = split.num0 + split.num1 - size;
    assert(split.numReplications >= 0);
    //assert(lnum > 0);
    //assert(rnum > 0);

    if (split.num0 == 0 || split.num1 == 0) {
      split.cost = inf;
      split.num0 = split.num1 = 1;
      split.bounds0 = split.bounds1 = BBox3fa(inf);
    }

    return split;
  }
      
  void BVH4BuilderHair::SpatialSplit::split(size_t threadIndex, BVH4BuilderHair* parent, atomic_set<PrimRefBlock>& prims, atomic_set<PrimRefBlock>& lprims_o, atomic_set<PrimRefBlock>& rprims_o) const
  {
    /* calculate splitting plane */
    const Vec3fa plane(space.vx[dim],space.vy[dim],space.vz[dim],-pos);
    
    /* sort each curve to left, right, or left and right */
    atomic_set<PrimRefBlock>::item* lblock = lprims_o.insert(parent->alloc.malloc(threadIndex));
    atomic_set<PrimRefBlock>::item* rblock = rprims_o.insert(parent->alloc.malloc(threadIndex));
    
    while (atomic_set<PrimRefBlock>::item* block = prims.take()) 
    {
      for (size_t i=0; i<block->size(); i++) 
      {
        const PrimRef& prim = block->at(i); 
        const float p0p = dot(prim.p0,plane)+plane.w;
        const float p3p = dot(prim.p3,plane)+plane.w;

        /* sort to the left side */
        if (p0p <= 0.0f && p3p <= 0.0f)
        {
          if (likely(lblock->insert(prim))) continue; 
          lblock = lprims_o.insert(parent->alloc.malloc(threadIndex));
          lblock->insert(prim);
          continue;
        }

        /* sort to the right side */
        if (p0p >= 0.0f && p3p >= 0.0f)
        {
          if (likely(rblock->insert(prim))) continue;
          rblock = rprims_o.insert(parent->alloc.malloc(threadIndex));
          rblock->insert(prim);
          continue;
        }

        /* split and sort to left and right */
        Bezier1 left,right;
        if (prim.split(plane,left,right)) 
        {
          if (!lblock->insert(left)) {
            lblock = lprims_o.insert(parent->alloc.malloc(threadIndex));
            lblock->insert(left);
          }
          if (!rblock->insert(right)) {
            rblock = rprims_o.insert(parent->alloc.malloc(threadIndex));
            rblock->insert(right);
          }
          continue;
        }

        /* insert to left side as fallback */
        if (!lblock->insert(prim)) {
          lblock = lprims_o.insert(parent->alloc.malloc(threadIndex));
          lblock->insert(prim);
        }
      }
      parent->alloc.free(threadIndex,block);
    }
    //assert(atomic_set<PrimRefBlock>::block_iterator_unsafe(lprims).size() == num0);
    //assert(atomic_set<PrimRefBlock>::block_iterator_unsafe(rprims).size() == num1);
  }

  __forceinline BVH4BuilderHair::FallBackSplit BVH4BuilderHair::FallBackSplit::find(size_t threadIndex, BVH4BuilderHair* parent, atomic_set<PrimRefBlock>& prims, 
                                                                                      atomic_set<PrimRefBlock>& lprims_o, atomic_set<PrimRefBlock>& rprims_o)
  {
    size_t num = 0, lnum = 0, rnum = 0;
    atomic_set<PrimRefBlock>::item* lblock = lprims_o.insert(parent->alloc.malloc(threadIndex));
    atomic_set<PrimRefBlock>::item* rblock = rprims_o.insert(parent->alloc.malloc(threadIndex));
    
    while (atomic_set<PrimRefBlock>::item* block = prims.take()) 
    {
      for (size_t i=0; i<block->size(); i++) 
      {
        const PrimRef& prim = block->at(i); 
        if ((num++)%2) 
        {
          lnum++;
          if (likely(lblock->insert(prim))) continue; 
          lblock = lprims_o.insert(parent->alloc.malloc(threadIndex));
          lblock->insert(prim);
        } 
        else 
        {
          rnum++;
          if (likely(rblock->insert(prim))) continue;
          rblock = rprims_o.insert(parent->alloc.malloc(threadIndex));
          rblock->insert(prim);
        }
      }
      parent->alloc.free(threadIndex,block);
    }
    const NAABBox3fa bounds0 = computeAlignedBounds(lprims_o);
    const NAABBox3fa bounds1 = computeAlignedBounds(rprims_o);
    return FallBackSplit(bounds0,lnum,bounds1,rnum);
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////

  BVH4::NodeRef BVH4BuilderHair::leaf(size_t threadIndex, size_t depth, atomic_set<PrimRefBlock>& prims, const NAABBox3fa& bounds)
  {
    //size_t N = end-begin;
    size_t N = atomic_set<PrimRefBlock>::block_iterator_unsafe(prims).size();

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
      atomic_set<PrimRefBlock>::block_iterator_unsafe iter(prims);
      for (size_t i=0; i<N; i++) { leaf[i] = *iter; iter++; }
      assert(!iter);

      /* free all primitive blocks */
      while (atomic_set<PrimRefBlock>::item* block = prims.take())
        alloc.free(threadIndex,block);

      return bvh->encodeLeaf((char*)leaf,N,1);
    } 
    else if (bvh->primTy[1] == &SceneBezier1i::type) {
      Bezier1i* leaf = (Bezier1i*) bvh->allocPrimitiveBlocks(threadIndex,0,N);
      atomic_set<PrimRefBlock>::block_iterator_unsafe iter(prims);
      for (size_t i=0; i<N; i++) {
        const Bezier1& curve = *iter; iter++;
        const BezierCurves* in = (BezierCurves*) scene->get(curve.geomID);
        const Vec3fa& p0 = in->vertex(in->curve(curve.primID));
        leaf[i] = Bezier1i(&p0,curve.geomID,curve.primID,-1); // FIXME: support mask
      }

      /* free all primitive blocks */
      while (atomic_set<PrimRefBlock>::item* block = prims.take())
        alloc.free(threadIndex,block);

      return bvh->encodeLeaf((char*)leaf,N,1);
    }
    else 
      throw std::runtime_error("unknown primitive type");
  }

  bool BVH4BuilderHair::split(size_t threadIndex, size_t depth, 
                               atomic_set<PrimRefBlock>& prims, const NAABBox3fa& bounds, size_t size,
                               atomic_set<PrimRefBlock>& lprims_o, size_t& lsize,
                               atomic_set<PrimRefBlock>& rprims_o, size_t& rsize,
                               bool& isAligned)
  {
    /* variable to track the SAH of the best splitting approach */
    float bestSAH = inf;
    bool enableSpatialSplits = remainingReplications > 0;
    const int travCostAligned = isAligned ? BVH4::travCostAligned : BVH4::travCostUnaligned;
    const float leafSAH = BVH4::intCost*countfunc(size)*embree::area(bounds.bounds);
    
    /* perform standard binning in aligned space */
    ObjectSplit alignedObjectSplit;
    float alignedObjectSAH = neg_inf;
    if (enableAlignedObjectSplits) {
      alignedObjectSplit = ObjectSplit::find(threadIndex,depth,this,prims,one);
      alignedObjectSAH = travCostAligned*embree::area(bounds.bounds) + alignedObjectSplit.standardSAH();
      bestSAH = min(bestSAH,alignedObjectSAH);
    }

    /* perform spatial split in aligned space */
    SpatialSplit alignedSpatialSplit;
    float alignedSpatialSAH = neg_inf;
    if (enableSpatialSplits && enableAlignedSpatialSplits) {
      alignedSpatialSplit = SpatialSplit::find(threadIndex,depth,size,this,prims,one);
      alignedSpatialSAH = travCostAligned*embree::area(bounds.bounds) + alignedSpatialSplit.standardSAH();
      bestSAH = min(bestSAH,alignedSpatialSAH);
    }

    /* perform standard binning in unaligned space */
    ObjectSplit unalignedObjectSplit;
    float unalignedObjectSAH = neg_inf;
    if (enableUnalignedObjectSplits) {
      unalignedObjectSplit = ObjectSplit::find(threadIndex,depth,this,prims,bounds.space);
      unalignedObjectSAH = BVH4::travCostUnaligned*embree::area(bounds.bounds) + unalignedObjectSplit.standardSAH();
      bestSAH = min(bestSAH,unalignedObjectSAH);
    }

    /* perform spatial split in unaligned space */
    SpatialSplit unalignedSpatialSplit;
    float unalignedSpatialSAH = neg_inf;
    if (enableSpatialSplits && enableUnalignedSpatialSplits) {
      unalignedSpatialSplit = SpatialSplit::find(threadIndex,depth,size,this,prims,bounds.space);
      unalignedSpatialSAH = BVH4::travCostUnaligned*embree::area(bounds.bounds) + unalignedSpatialSplit.standardSAH();
      bestSAH = min(bestSAH,unalignedSpatialSAH);
    }

    /* perform splitting into two strands */
    StrandSplit strandSplit;
    float strandSAH = neg_inf;
    if (enableStrandSplits) {
      strandSplit = StrandSplit::find(threadIndex,this,prims);
      strandSAH = BVH4::travCostUnaligned*embree::area(bounds.bounds) + strandSplit.standardSAH();
      bestSAH = min(bestSAH,strandSAH);
    }

    /* perform fallback split */
    if (bestSAH == float(inf)) {
      //if (N <= maxLeafSize) return false;
      numFallbackSplits++;
      const FallBackSplit fallbackSplit = FallBackSplit::find(threadIndex,this,prims,lprims_o,rprims_o);
      assert(atomic_set<PrimRefBlock>::block_iterator_unsafe(lprims_o).size());
      assert(atomic_set<PrimRefBlock>::block_iterator_unsafe(rprims_o).size());
      lsize = fallbackSplit.num0;
      rsize = fallbackSplit.num1;
      return true;
    }

    /* perform aligned object split */
    else if (bestSAH == alignedObjectSAH && enableAlignedObjectSplits) {
      numAlignedObjectSplits++;
      alignedObjectSplit.split(threadIndex,this,prims,lprims_o,rprims_o);
      assert(atomic_set<PrimRefBlock>::block_iterator_unsafe(lprims_o).size());
      assert(atomic_set<PrimRefBlock>::block_iterator_unsafe(rprims_o).size());
      lsize = alignedObjectSplit.num0;
      rsize = alignedObjectSplit.num1;
      return true;
    }

    /* perform aligned spatial split */
    else if (bestSAH == alignedSpatialSAH && enableSpatialSplits && enableAlignedSpatialSplits) {
      numAlignedSpatialSplits++;
      alignedSpatialSplit.split(threadIndex,this,prims,lprims_o,rprims_o);
      assert(atomic_set<PrimRefBlock>::block_iterator_unsafe(lprims_o).size());
      assert(atomic_set<PrimRefBlock>::block_iterator_unsafe(rprims_o).size());
      lsize = alignedSpatialSplit.num0;
      rsize = alignedSpatialSplit.num1;
      atomic_add(&remainingReplications,-alignedSpatialSplit.numReplications);
      return true;
    }

    /* perform unaligned object split */
    else if (bestSAH == unalignedObjectSAH && enableUnalignedObjectSplits) {
      numUnalignedObjectSplits++;
      unalignedObjectSplit.split(threadIndex,this,prims,lprims_o,rprims_o);
      assert(atomic_set<PrimRefBlock>::block_iterator_unsafe(lprims_o).size());
      assert(atomic_set<PrimRefBlock>::block_iterator_unsafe(rprims_o).size());
      lsize = unalignedObjectSplit.num0;
      rsize = unalignedObjectSplit.num1;
      isAligned = false;
      return true;
    }

    /* perform unaligned spatial split */
    else if (bestSAH == unalignedSpatialSAH && enableSpatialSplits && enableUnalignedSpatialSplits) {
      numUnalignedSpatialSplits++;
      unalignedSpatialSplit.split(threadIndex,this,prims,lprims_o,rprims_o);
      assert(atomic_set<PrimRefBlock>::block_iterator_unsafe(lprims_o).size());
      assert(atomic_set<PrimRefBlock>::block_iterator_unsafe(rprims_o).size());
      lsize = unalignedSpatialSplit.num0;
      rsize = unalignedSpatialSplit.num1;
      atomic_add(&remainingReplications,-unalignedSpatialSplit.numReplications);
      isAligned = false;
      return true;
    }

    /* perform strand split */
    else if (bestSAH == strandSAH && enableStrandSplits) {
      numStrandSplits++;
      strandSplit.split(threadIndex,this,prims,lprims_o,rprims_o);
      assert(atomic_set<PrimRefBlock>::block_iterator_unsafe(lprims_o).size());
      assert(atomic_set<PrimRefBlock>::block_iterator_unsafe(rprims_o).size());
      lsize = strandSplit.num0;
      rsize = strandSplit.num1;
      isAligned = false;
      return true;
    }
 
    else {
      throw std::runtime_error("bvh4hair_builder: internal error");
      return true;
    }
  }

  void BVH4BuilderHair::processTask(size_t threadIndex, BuildTask& task, BuildTask task_o[BVH4::N], size_t& numTasks_o)
  {
    /* create enforced leaf */
    if (task.size <= minLeafSize || task.depth >= BVH4::maxBuildDepth || task.makeleaf) {
      *task.dst = leaf(threadIndex,task.depth,task.prims,task.bounds);
      numTasks_o = 0;
      return;
    }

    /*! initialize child list */
    bool isAligned = true;
    NAABBox3fa cbounds[BVH4::N];
    atomic_set<PrimRefBlock> cprims[BVH4::N];
    size_t csize[BVH4::N];
    bool isleaf[BVH4::N];
    cprims[0] = task.prims;
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
        size_t N = atomic_set<PrimRefBlock>::block_iterator_unsafe(cprims[i]).size(); // FIXME: slow
        float A = embree::area(cbounds[i].bounds);
        if (N <= minLeafSize) continue;  
        if (isleaf[i]) continue;
        if (A > bestArea) { bestChild = i; bestArea = A; }
      }
      if (bestChild == -1) break;

      /*! split selected child */
      size_t lsize, rsize;
      atomic_set<PrimRefBlock> lprims, rprims;
      bool done = split(threadIndex,task.depth,cprims[bestChild],cbounds[bestChild],csize[bestChild],lprims,lsize,rprims,rsize,isAligned);
      if (!done) { isleaf[bestChild] = true; continue; }
      cprims[numChildren] = rprims; isleaf[numChildren] = false; csize[numChildren] = rsize;
      cprims[bestChild  ] = lprims; isleaf[bestChild  ] = false; csize[bestChild  ] = lsize;

#if BVH4HAIR_COMPRESS_UNALIGNED_NODES
      cbounds[numChildren] = computeAlignedBounds(cprims[numChildren],task.bounds.space);
      cbounds[bestChild  ] = computeAlignedBounds(cprims[bestChild  ],task.bounds.space);
#else
      cbounds[numChildren] = computeUnalignedBounds(cprims[numChildren]);
      cbounds[bestChild  ] = computeUnalignedBounds(cprims[bestChild  ]);
#endif
      numChildren++;
      
    } while (numChildren < BVH4::N);

    /* create aligned node */
    if (isAligned) 
    {
      BVH4::UANode* node = bvh->allocUANode(threadIndex);

      BBox3fa bounds = empty;
      NAABBox3fa abounds[BVH4::N];
      for (size_t i=0; i<numChildren; i++) {
        abounds[i] = computeAlignedBounds(cprims[i]);
        bounds.extend(abounds[i].bounds);
      }

#if BVH4HAIR_COMPRESS_ALIGNED_NODES
      node->set(bounds);
#endif
      for (ssize_t i=0; i<numChildren; i++) {
        node->set(i,abounds[i].bounds);
	const NAABBox3fa ubounds = computeUnalignedBounds(cprims[i]);
        new (&task_o[i]) BuildTask(&node->child(i),task.depth+1,csize[i],isleaf[i],cprims[i],ubounds);
      }
      numTasks_o = numChildren;
      *task.dst = bvh->encodeNode(node);
    }
    
    /* create unaligned node */
    else {
      BVH4::UUNode* node = bvh->allocUUNode(threadIndex);
#if BVH4HAIR_COMPRESS_UNALIGNED_NODES
      node->set(task.bounds);
      for (ssize_t i=0; i<numChildren; i++) {
        node->set(i,cbounds[i].bounds);
        const NAABBox3fa ubounds = computeUnalignedBounds(cprims[i]);
        new (&task_o[i]) BuildTask(&node->child(i),task.depth+1,csize[i],isleaf[i],cprims[i],ubounds);
      }
#else
      for (ssize_t i=numChildren-1; i>=0; i--) {
        node->set(i,cbounds[i]);
        new (&task_o[i]) BuildTask(&node->child(i),task.depth+1,csize[i],isleaf[i],cprims[i],cbounds[i]);
      }
#endif
      numTasks_o = numChildren;
      *task.dst = bvh->encodeNode(node);
    }
  }

  void BVH4BuilderHair::recurseTask(size_t threadIndex, BuildTask& task)
  {
    size_t numChildren;
    BuildTask tasks[BVH4::N];
    processTask(threadIndex,task,tasks,numChildren);
    for (size_t i=0; i<numChildren; i++) 
      recurseTask(threadIndex,tasks[i]);
  }

  Builder* BVH4BuilderHair_ (BVH4* accel, Scene* scene) {
    return new BVH4BuilderHair(accel,scene);
  }
}
