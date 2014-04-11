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

#include "heuristic_binning2.h"

namespace embree
{
  template<int logBlockSize>
  const size_t HeuristicBinning2<logBlockSize>::maxBins;

  template<int logBlockSize>
  void HeuristicBinning2<logBlockSize>::add(atomic_set<PrimRefBlock>& prims)
  {
    atomic_set<PrimRefBlock>::iterator i=prims;
    while (PrimRefBlock* block = i.next()) {
      info(block->base(),block->size());
    }

    new (&mapping) Mapping(pinfo);
    for (size_t i=0; i<mapping.size(); i++) {
      counts[i] = 0;
      geomBounds[i][0] = geomBounds[i][1] = geomBounds[i][2] = empty;
    }

    atomic_set<PrimRefBlock>::iterator j=prims;
    while (PrimRefBlock* block = j.next())
      bin(block->base(),block->size());
  }

  template<int logBlockSize>
  void HeuristicBinning2<logBlockSize>::info(const PrimRef* prims, size_t num)
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

  template<int logBlockSize>
  void HeuristicBinning2<logBlockSize>::bin(const PrimRef* prims, size_t num)
  {
    if (num == 0) return;
    
    size_t i; for (i=0; i<num-1; i+=2)
    {
      /*! map even and odd primitive to bin */
      const BBox3fa prim0 = prims[i+0].bounds(); const Vec3ia bin0 = mapping.bin(prim0); const Vec3fa center0 = Vec3fa(center2(prim0));
      const BBox3fa prim1 = prims[i+1].bounds(); const Vec3ia bin1 = mapping.bin(prim1); const Vec3fa center1 = Vec3fa(center2(prim1));
      
      /*! increase bounds for bins for even primitive */
      const int b00 = bin0.x; counts[b00][0]++; geomBounds[b00][0].extend(prim0); centBounds[b00][0].extend(center0);
      const int b01 = bin0.y; counts[b01][1]++; geomBounds[b01][1].extend(prim0); centBounds[b01][1].extend(center0);
      const int b02 = bin0.z; counts[b02][2]++; geomBounds[b02][2].extend(prim0); centBounds[b02][2].extend(center0);
      
      /*! increase bounds of bins for odd primitive */
      const int b10 = bin1.x; counts[b10][0]++; geomBounds[b10][0].extend(prim1); centBounds[b10][0].extend(center1);
      const int b11 = bin1.y; counts[b11][1]++; geomBounds[b11][1].extend(prim1); centBounds[b11][1].extend(center1);
      const int b12 = bin1.z; counts[b12][2]++; geomBounds[b12][2].extend(prim1); centBounds[b12][2].extend(center1);
    }
    
    /*! for uneven number of primitives */
    if (i < num)
    {
      /*! map primitive to bin */
      const BBox3fa prim0 = prims[i].bounds(); const Vec3ia bin0 = mapping.bin(prim0); const Vec3fa center0 = Vec3fa(center2(prim0));
      
      /*! increase bounds of bins */
      const int b00 = bin0.x; counts[b00][0]++; geomBounds[b00][0].extend(prim0); centBounds[b00][0].extend(center0);
      const int b01 = bin0.y; counts[b01][1]++; geomBounds[b01][1].extend(prim0); centBounds[b01][1].extend(center0);
      const int b02 = bin0.z; counts[b02][2]++; geomBounds[b02][2].extend(prim0); centBounds[b02][2].extend(center0);
    }
  }
  
  template<int logBlockSize>
  void HeuristicBinning2<logBlockSize>::best(Split& split)
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
    if (split.pos == 0) return;
    
    /* calculate geometry info from binning data */
    size_t pos = split.pos, dim = split.dim;
    size_t numLeft = 0, numRight = 0;
    BBox3fa lcentBounds = empty, rcentBounds = empty;
    BBox3fa lgeomBounds = empty, rgeomBounds = empty;
    for (size_t i=0; i<pos; i++) {
      numLeft += counts[i][dim];
      lcentBounds.extend(centBounds[i][dim]);
      lgeomBounds.extend(geomBounds[i][dim]);
    }
    for (size_t i=pos; i<mapping.size(); i++) {
      numRight += counts[i][dim];
      rcentBounds.extend(centBounds[i][dim]);
      rgeomBounds.extend(geomBounds[i][dim]);
    }
    assert(numLeft + numRight == pinfo.size());
    new (&split.linfo) PrimInfo(numLeft ,lgeomBounds,lcentBounds);
    new (&split.rinfo) PrimInfo(numRight,rgeomBounds,rcentBounds);
  }
  
  template<int logBlockSize>
  void HeuristicBinning2<logBlockSize>::reduce(const HeuristicBinning2 binners[], size_t num, HeuristicBinning2& binner_o)
  {
    binner_o = binners[0];
    for (size_t tid=1; tid<num; tid++) {
      const HeuristicBinning2& binner = binners[tid];
      for (size_t bin=0; bin<binner.mapping.size(); bin++) 
      {
        for (size_t dim=0; dim<3; dim++) {
          binner_o.counts    [bin][dim] += binner.counts[bin][dim];
          binner_o.geomBounds[bin][dim].extend(binner.geomBounds[bin][dim]);
          binner_o.centBounds[bin][dim].extend(binner.centBounds[bin][dim]);
        }
      }
    }
  }

  template<int logBlockSize>
  void HeuristicBinning2<logBlockSize>::Split::split(size_t thread, PrimRefAlloc* alloc, 
                                                     atomic_set<PrimRefBlock>& prims, 
                                                     atomic_set<PrimRefBlock>& lprims, Split& lsplit,
                                                     atomic_set<PrimRefBlock>& rprims, Split& rsplit)
  {
    HeuristicBinning2 lheuristic(this->linfo);
    HeuristicBinning2 rheuristic(this->rinfo);
    atomic_set<PrimRefBlock>::item* lblock = lprims.insert(alloc->malloc(thread));
    atomic_set<PrimRefBlock>::item* rblock = rprims.insert(alloc->malloc(thread));
    
    while (atomic_set<PrimRefBlock>::item* block = prims.take()) 
    {
      for (size_t i=0; i<block->size(); i++) 
      {
        const PrimRef& prim = block->at(i); 
                
        if (left(prim)) 
        {
          if (likely(lblock->insert(prim))) continue; 
          lheuristic.bin(lblock->base(),lblock->size());
          lblock = lprims.insert(alloc->malloc(thread));
          lblock->insert(prim);
        } 
        else 
        {
          if (likely(rblock->insert(prim))) continue;
          rheuristic.bin(rblock->base(),rblock->size());
          rblock = rprims.insert(alloc->malloc(thread));
          rblock->insert(prim);
        }
      }
      alloc->free(thread,block);
    }
    lheuristic.bin(lblock->base(),lblock->size()); linfo = this->linfo; lheuristic.best(lsplit); 
    rheuristic.bin(rblock->base(),rblock->size()); rinfo = this->rinfo; rheuristic.best(rsplit);
  }

  /*! explicit template instantiations */
  template class HeuristicBinning2<0>;
  template class HeuristicBinning2<2>;
  template class HeuristicBinning2<3>;
}
