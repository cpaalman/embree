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

#include "bvh4_statistics.h"

namespace embree
{
  BVH4Statistics::BVH4Statistics (BVH4* bvh) : bvh(bvh)
  {
    numNodes = depth = 0;
    for (size_t i=0; i<4; i++) numLeaves[i] = numPrimBlocks[i] = numPrims[i] = 0;
    bvhSAH = leafSAH = 0.0f;
    statistics(bvh->root,safeArea(bvh->bounds),depth);
    bvhSAH /= area(bvh->bounds);
    leafSAH /= area(bvh->bounds);
    assert(depth <= BVH4::maxDepth);
  }

  size_t BVH4Statistics::bytesUsed()
  {
    size_t bytesNodes = numNodes*sizeof(Node);
    size_t bytesPrims = 0;
    for (size_t i=0; i<4; i++)
      bytesPrims += numPrimBlocks[i]*(bvh->primTy[i] ? bvh->primTy[i]->bytes : 0);
    size_t numVertices = bvh->numVertices;
    size_t bytesVertices = numVertices*sizeof(Vec3fa); 
    return bytesNodes+bytesPrims+bytesVertices;
  }

  std::string BVH4Statistics::str()  
  {
    std::ostringstream stream;
    size_t bytesNodes = numNodes*sizeof(Node);
    size_t bytesPrim[4]; memset(bytesPrim,0,sizeof(bytesPrim));
    size_t bytesPrims = 0;
    size_t numLeavesTotal = 0;
    for (size_t i=0; i<4; i++) {
      bytesPrim[i] = 0;
      numLeavesTotal += numLeaves[i];
      if (bvh->primTy[i] == NULL) continue;
      bytesPrim[i] = numPrimBlocks[i]*bvh->primTy[i]->bytes;
      bytesPrims += bytesPrim[i];
    }
    size_t numVertices = bvh->numVertices;
    size_t bytesVertices = numVertices*sizeof(Vec3fa); 
    size_t bytesTotal = bytesNodes+bytesPrims+bytesVertices;
    size_t bytesTotalAllocated = bvh->bytesAllocated();
    if (g_benchmark) std::cout << "BENCHMARK_TRIANGLE_ACCEL " << bvhSAH << " " << bytesTotal << std::endl;
    stream.setf(std::ios::fixed, std::ios::floatfield);
    stream << "  primitives = " << bvh->numPrimitives << ", vertices = " << bvh->numVertices << std::endl;
    stream.precision(4);
    stream << "  sah = " << bvhSAH << ", leafSAH = " << leafSAH;
    stream.setf(std::ios::fixed, std::ios::floatfield);
    stream.precision(1);
    stream << ", depth = " << depth << std::endl;
    stream << "  used = " << bytesTotal/1E6 << " MB, allocated = " << bytesTotalAllocated/1E6 << " MB, perPrimitive = " << double(bytesTotal)/double(bvh->numPrimitives) << " B" << std::endl;
    stream.precision(1);
    stream << "  nodes = "  << numNodes << " "
           << "(" << bytesNodes/1E6  << " MB) "
           << "(" << 100.0*double(bytesNodes)/double(bytesTotal) << "% of total) "
           << "(" << 100.0*(numNodes-1+numLeavesTotal)/(BVH4::N*numNodes) << "% used)" 
           << std::endl;
    for (size_t i=0; i<4; i++) {
      if (!bvh->primTy[i]) continue;
      stream << "  " << bvh->primTy[i]->name << " leaves = " << numLeaves[i] << " "
             << "(" << bytesPrim[i]/1E6  << " MB) "
             << "(" << 100.0*double(bytesPrim[i])/double(bytesTotal) << "% of total) "
             << "(" << 100.0*double(numPrims[i])/double(bvh->primTy[i]->blockSize*numPrimBlocks[i]) << "% used)" 
             << std::endl;
    }
    stream << "  vertices = " << numVertices << " "
           << "(" << bytesVertices/1E6 << " MB) " 
           << "(" << 100.0*double(bytesVertices)/double(bytesTotal) << "% of total) "
           << "(" << 100.0*12.0f/float(sizeof(Vec3fa)) << "% used)" 
           << std::endl;
    return stream.str();
  }

  void BVH4Statistics::statistics(NodeRef node, const float A, size_t& depth)
  {
    if (node.isNode())
    {
      numNodes++;
      BVH4::Node* n = node.getNode();
      for (size_t i=0; i<BVH4::N; i++) {
        if (n->child(i) == BVH4::emptyNode) {
          for (; i<BVH4::N; i++) {
            if (n->child(i) != BVH4::emptyNode)
              throw std::runtime_error("invalid node");
          }
          break;
        }
      }
    }    

    if (node.isUANode())
    {
      depth = 0;
      size_t cdepth = 0;
      BVH4::UANode* n = node.getUANode();
      bvhSAH += A*BVH4::travCost;
      for (size_t i=0; i<BVH4::N; i++) {
        statistics(n->child(i),safeArea(n->bounds(i)),cdepth); 
        depth=max(depth,cdepth);
      }
      depth++;
    }
    else if (node.isUUNode())
    {
      depth = 0;
      size_t cdepth = 0;
      BVH4::UUNode* n = node.getUUNode();
      bvhSAH += A*BVH4::travCost;
      for (size_t i=0; i<BVH4::N; i++) {
        statistics(n->child(i),safeArea(n->extend(i)),cdepth); 
        depth=max(depth,cdepth);
      }
      depth++;
    }
    else if (node.isLeaf())
    {
      depth = 0;
      size_t num,ty; const char* tri = node.getLeaf(num,ty);
      if (!num) return;
      if (ty >= 4) throw std::runtime_error("invalid leaf type");
      
      numLeaves[ty]++;
      numPrimBlocks[ty] += num;
      for (size_t i=0; i<num; i++) {
        numPrims[ty] += bvh->primTy[ty]->size(tri+i*bvh->primTy[ty]->bytes);
      }
      float sah = A * bvh->primTy[ty]->intCost * num;
      bvhSAH += sah;
      leafSAH += sah;
    }
    else
      throw std::runtime_error("unsupported node type");
  }
}
