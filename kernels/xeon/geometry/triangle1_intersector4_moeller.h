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

#include "triangle1.h"
#include "common/ray4.h"
#include "geometry/filter.h"

namespace embree
{
  /*! Intersector for individual precomputed triangles with 4
   *  rays. This intersector implements a modified version of the
   *  Moeller Trumbore intersector from the paper "Fast, Minimum
   *  Storage Ray-Triangle Intersection". In contrast to the paper we
   *  precalculate some factors and factor the calculations
   *  differently to allow precalculating the cross product e1 x
   *  e2. */
  struct Triangle1Intersector4MoellerTrumbore
  {
    typedef Triangle1 Primitive;
    
    static __forceinline void intersect(const sseb& valid_i, Ray4& ray, const Triangle1* __restrict__ tris, size_t num, const void* geom)
    {
      for (size_t i=0; i<num; i++) 
      {
        STAT3(normal.trav_prims,1,popcnt(valid_i),4);

        sseb valid = valid_i;
        const sse3f org = ray.org;
        const sse3f dir = ray.dir;
        const ssef zero = ssef(0.0f);
        const Triangle1& tri = tris[i];
        
        /* load vertices and calculate edges */
        const ssef v0 = load4f(&tri.v0.x);
        const ssef v1 = load4f(&tri.v1.x);
        const ssef v2 = load4f(&tri.v2.x);
        const ssef e1 = v0-v1;
        const ssef e2 = v2-v0;
        
        /* calculate denominator */
        const sse3f _v0 = sse3f(shuffle<0>(v0),shuffle<1>(v0),shuffle<2>(v0));
        const sse3f C =  _v0 - org;
        const sse3f Ng = sse3f(tri.Ng);
        const ssef den = dot(dir,Ng);
        const ssef sgnDen = signmsk(den);
        const ssef absDen = abs(den);
#if defined(__BACKFACE_CULLING__)
        valid &= den > ssef(zero);
#else
        valid &= den != ssef(zero);
#endif

        /* perform edge tests */
        const sse3f R = cross(dir,C);
        const sse3f _e1(shuffle<0>(e1),shuffle<1>(e1),shuffle<2>(e1));
        const ssef V = dot(R,_e1)^sgnDen;
        const sse3f _e2(shuffle<0>(e2),shuffle<1>(e2),shuffle<2>(e2));
        const ssef U = dot(R,_e2)^sgnDen;
        valid &= V >= zero & U >= zero & U+V <= absDen;
        if (unlikely(none(valid))) continue;
      
        /* perform depth test */
        const ssef T = dot(C,Ng) ^ sgnDen;
        valid &= (T >= absDen*ray.tnear) & (absDen*ray.tfar >= T);
        if (unlikely(none(valid))) continue;
        
        /* ray masking test */
#if defined(__USE_RAY_MASK__)
        valid &= (tri.mask() & ray.mask) != 0;
        if (unlikely(none(valid))) continue;
#endif

        /* calculate hit information */
        const ssef rcpAbsDen = rcp(absDen);
        const ssef u = U*rcpAbsDen;
        const ssef v = V*rcpAbsDen;
        const ssef t = T*rcpAbsDen;
        const int geomID = tri.geomID();
        const int primID = tri.primID();

        /* intersection filter test */
#if defined(__INTERSECTION_FILTER__)
        Geometry* geometry = ((Scene*)geom)->get(geomID);
        if (unlikely(geometry->hasIntersectionFilter4())) {
          runIntersectionFilter4(valid,geometry,ray,u,v,t,Ng,geomID,primID);
          continue;
        }
#endif

        /* update hit information */
        store4f(valid,&ray.u,u);
        store4f(valid,&ray.v,v);
        store4f(valid,&ray.tfar,t);
        store4i(valid,&ray.geomID,geomID);
        store4i(valid,&ray.primID,primID);
        store4f(valid,&ray.Ng.x,Ng.x);
        store4f(valid,&ray.Ng.y,Ng.y);
        store4f(valid,&ray.Ng.z,Ng.z);
      }
    }

    static __forceinline sseb occluded(const sseb& valid_i, Ray4& ray, const Triangle1* __restrict__ tris, size_t num, const void* geom)
    {
      sseb valid0 = valid_i;

      for (size_t i=0; i<num; i++) 
      {
        STAT3(shadow.trav_prims,1,popcnt(valid0),4);

        sseb valid = valid0;
        const sse3f org = ray.org;
        const sse3f dir = ray.dir;
        const ssef zero = ssef(0.0f);
        const Triangle1& tri = tris[i];
        
        /* load vertices and calculate edges */
        const ssef v0 = load4f(&tri.v0.x);
        const ssef v1 = load4f(&tri.v1.x);
        const ssef v2 = load4f(&tri.v2.x);
        const ssef e1 = v0-v1;
        const ssef e2 = v2-v0;
        
        /* calculate denominator */
        const sse3f _v0 = sse3f(shuffle<0>(v0),shuffle<1>(v0),shuffle<2>(v0));
        const sse3f C =  _v0 - org;
        const sse3f Ng = sse3f(tri.Ng);
        const ssef den = dot(dir,Ng);
        const ssef sgnDen = signmsk(den);
        const ssef absDen = abs(den);
#if defined(__BACKFACE_CULLING__)
        valid &= den > ssef(zero);
#else
        valid &= den != ssef(zero);
#endif
        
        /* perform edge tests */
        const sse3f R = cross(dir,C);
        const sse3f _e1(shuffle<0>(e1),shuffle<1>(e1),shuffle<2>(e1));
        const ssef V = dot(R,_e1)^sgnDen;
        const sse3f _e2(shuffle<0>(e2),shuffle<1>(e2),shuffle<2>(e2));
        const ssef U = dot(R,_e2)^sgnDen;
        valid &= V >= zero & U >= zero & U+V <= absDen;
        if (unlikely(none(valid))) continue;
      
        /* perform depth test */
        const ssef T = dot(C,Ng) ^ sgnDen;
        valid &= (T >= absDen*ray.tnear) & (absDen*ray.tfar >= T);
        if (unlikely(none(valid))) continue;

        /* ray masking test */
#if defined(__USE_RAY_MASK__)
        valid &= (tri.mask() & ray.mask) != 0;
        if (unlikely(none(valid))) continue;
#endif

        /* intersection filter test */
#if defined(__INTERSECTION_FILTER__)
        const int geomID = tri.geomID();
        Geometry* geometry = ((Scene*)geom)->get(geomID);
        if (unlikely(geometry->hasOcclusionFilter4()))
        {
          /* calculate hit information */
          const ssef rcpAbsDen = rcp(absDen);
          const ssef u = U*rcpAbsDen;
          const ssef v = V*rcpAbsDen;
          const ssef t = T*rcpAbsDen;
          const int primID = tri.primID();
          valid = runOcclusionFilter4(valid,geometry,ray,u,v,t,Ng,geomID,primID);
        }
#endif
        
        /* update occlusion */
        valid0 &= !valid;
        if (none(valid0)) break;
      }
      return !valid0;
    }
  };
}
