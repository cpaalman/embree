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

#include "../common/tutorial/tutorial_device.h"

/* scene data */
RTCScene g_scene = NULL;
Vec3f* colors = NULL;

/* render function to use */
renderPixelFunc renderPixel;

/* shadow ray structure that includes total transparency along the ray */
struct RTCRay2
{
  RTCRay ray;
  float transparency;
};

/* 3D procedural transparency */
inline float transparencyFunction(RTCRay2& ray)
{
  Vec3f h = add(ray.ray.org,mul(ray.ray.dir,ray.ray.tfar));
  float v = abs(sin(4.0f*h.x)*cos(4.0f*h.y)*sin(4.0f*h.z));
  float T = clamp((v-0.1f)*3.0f,0.0f,1.0f);
  return T;
}

/* intersection filter function */
void intersectionFilter(void* ptr, RTCRay2& ray)
{
  float T = transparencyFunction(ray);
  if (T >= 1.0f) ray.ray.geomID = -1;
  else ray.transparency = T;
}

/* occlusion filter function */
void occlusionFilter(void* ptr, RTCRay2& ray)
{
  float T = transparencyFunction(ray);
  T *= ray.transparency;
  ray.transparency = T;
  if (T != 0.0f) {
    ray.ray.geomID = -1;
  }
}

/* adds a cube to the scene */
unsigned int addCube (RTCScene scene_i)
{
  /* create a triangulated cube with 12 triangles and 8 vertices */
  unsigned int mesh = rtcNewTriangleMesh (scene_i, RTC_GEOMETRY_STATIC, 12, 8);

  /* set vertices */
  Vertex* vertices = (Vertex*) rtcMapBuffer(scene_i,mesh,RTC_VERTEX_BUFFER); 
  vertices[0].x = -1; vertices[0].y = -1; vertices[0].z = -1; 
  vertices[1].x = -1; vertices[1].y = -1; vertices[1].z = +1; 
  vertices[2].x = -1; vertices[2].y = +1; vertices[2].z = -1; 
  vertices[3].x = -1; vertices[3].y = +1; vertices[3].z = +1; 
  vertices[4].x = +1; vertices[4].y = -1; vertices[4].z = -1; 
  vertices[5].x = +1; vertices[5].y = -1; vertices[5].z = +1; 
  vertices[6].x = +1; vertices[6].y = +1; vertices[6].z = -1; 
  vertices[7].x = +1; vertices[7].y = +1; vertices[7].z = +1; 
  rtcUnmapBuffer(scene_i,mesh,RTC_VERTEX_BUFFER); 

  /* create triangle color array */
  colors = new Vec3f[12];

  /* set triangles and colors */
  int tri = 0;
  Triangle* triangles = (Triangle*) rtcMapBuffer(scene_i,mesh,RTC_INDEX_BUFFER);
  
  // left side
  colors[tri] = Vec3f(1,0,0); triangles[tri].v0 = 0; triangles[tri].v1 = 2; triangles[tri].v2 = 1; tri++;
  colors[tri] = Vec3f(1,0,0); triangles[tri].v0 = 1; triangles[tri].v1 = 2; triangles[tri].v2 = 3; tri++;

  // right side
  colors[tri] = Vec3f(0,1,0); triangles[tri].v0 = 4; triangles[tri].v1 = 5; triangles[tri].v2 = 6; tri++;
  colors[tri] = Vec3f(0,1,0); triangles[tri].v0 = 5; triangles[tri].v1 = 7; triangles[tri].v2 = 6; tri++;

  // bottom side
  colors[tri] = Vec3f(0.5f);  triangles[tri].v0 = 0; triangles[tri].v1 = 1; triangles[tri].v2 = 4; tri++;
  colors[tri] = Vec3f(0.5f);  triangles[tri].v0 = 1; triangles[tri].v1 = 5; triangles[tri].v2 = 4; tri++;

  // top side
  colors[tri] = Vec3f(1.0f);  triangles[tri].v0 = 2; triangles[tri].v1 = 6; triangles[tri].v2 = 3; tri++;
  colors[tri] = Vec3f(1.0f);  triangles[tri].v0 = 3; triangles[tri].v1 = 6; triangles[tri].v2 = 7; tri++;

  // front side
  colors[tri] = Vec3f(0,0,1); triangles[tri].v0 = 0; triangles[tri].v1 = 4; triangles[tri].v2 = 2; tri++;
  colors[tri] = Vec3f(0,0,1); triangles[tri].v0 = 2; triangles[tri].v1 = 4; triangles[tri].v2 = 6; tri++;

  // back side
  colors[tri] = Vec3f(1,1,0); triangles[tri].v0 = 1; triangles[tri].v1 = 3; triangles[tri].v2 = 5; tri++;
  colors[tri] = Vec3f(1,1,0); triangles[tri].v0 = 3; triangles[tri].v1 = 7; triangles[tri].v2 = 5; tri++;

  rtcUnmapBuffer(scene_i,mesh,RTC_INDEX_BUFFER);

  /* set intersection filter for the cube */
  rtcSetIntersectionFilterFunction(scene_i,mesh,(RTCFilterFunc)&intersectionFilter);
  rtcSetOcclusionFilterFunction   (scene_i,mesh,(RTCFilterFunc)&occlusionFilter);

  return mesh;
}

/* adds a ground plane to the scene */
unsigned int addGroundPlane (RTCScene scene_i)
{
  /* create a triangulated plane with 2 triangles and 4 vertices */
  unsigned int mesh = rtcNewTriangleMesh (scene_i, RTC_GEOMETRY_STATIC, 2, 4);

  /* set vertices */
  Vertex* vertices = (Vertex*) rtcMapBuffer(scene_i,mesh,RTC_VERTEX_BUFFER); 
  vertices[0].x = -10; vertices[0].y = -2; vertices[0].z = -10; 
  vertices[1].x = -10; vertices[1].y = -2; vertices[1].z = +10; 
  vertices[2].x = +10; vertices[2].y = -2; vertices[2].z = -10; 
  vertices[3].x = +10; vertices[3].y = -2; vertices[3].z = +10;
  rtcUnmapBuffer(scene_i,mesh,RTC_VERTEX_BUFFER); 

  /* set triangles */
  Triangle* triangles = (Triangle*) rtcMapBuffer(scene_i,mesh,RTC_INDEX_BUFFER);
  triangles[0].v0 = 0; triangles[0].v1 = 2; triangles[0].v2 = 1;
  triangles[1].v0 = 1; triangles[1].v1 = 2; triangles[1].v2 = 3;
  rtcUnmapBuffer(scene_i,mesh,RTC_INDEX_BUFFER);

  return mesh;
}

/* called by the C++ code for initialization */
extern "C" void device_init (int8* cfg)
{
  /* initialize ray tracing core */
  rtcInit(cfg);

  /* create scene */
  g_scene = rtcNewScene(RTC_SCENE_STATIC,RTC_INTERSECT1);

  /* add cube */
  addCube(g_scene);

  /* add ground plane */
  addGroundPlane(g_scene);

  /* commit changes to scene */
  rtcCommit (g_scene);

  /* set start render mode */
  renderPixel = renderPixelStandard;
}

/* task that renders a single screen tile */
Vec3fa renderPixelStandard(int x, int y, const Vec3fa& vx, const Vec3fa& vy, const Vec3fa& vz, const Vec3fa& p)
{
  float weight = 1.0f;
  Vec3f color = Vec3f(0.0f);

  /* initialize ray */
  RTCRay2 primary;
  primary.ray.org = p;
  primary.ray.dir = normalize(add(mul(x,vx), mul(y,vy), vz));
  primary.ray.tnear = 0.0f;
  primary.ray.tfar = inf;
  primary.ray.geomID = -1;
  primary.ray.primID = -1;
  primary.ray.mask = -1;
  primary.ray.time = 0;
  primary.transparency = 0.0f;

  while (true)
  {
    /* intersect ray with scene */
    rtcIntersect(g_scene,primary.ray);
    
    /* shade pixels */
    if (primary.ray.geomID == -1) 
      break;

    float opacity = 1.0f-primary.transparency;
    Vec3f diffuse = colors[primary.ray.primID];
    Vec3f La = mul(diffuse,0.5f);
    color = add(color,mul(weight*opacity,La));
    Vec3f lightDir = normalize(Vec3f(-1,-1,-1));
      
    /* initialize shadow ray */
    RTCRay2 shadow;
    shadow.ray.org = add(primary.ray.org,mul(primary.ray.tfar,primary.ray.dir));
    shadow.ray.dir = neg(lightDir);
    shadow.ray.tnear = 0.001f;
    shadow.ray.tfar = inf;
    shadow.ray.geomID = RTC_INVALID_GEOMETRY_ID;
    shadow.ray.primID = RTC_INVALID_GEOMETRY_ID;
    shadow.ray.mask = -1;
    shadow.ray.time = 0;
    shadow.transparency = 1.0f;
    
    /* trace shadow ray */
    rtcOccluded(g_scene,shadow.ray);
    
    /* add light contribution */
    if (shadow.ray.geomID) {
      Vec3f Ll = mul(diffuse,shadow.transparency*clamp(-dot(lightDir,normalize(primary.ray.Ng)),0.0f,1.0f));
      color = add(color,mul(weight*opacity,Ll));
    }

    /* shoot transmission ray */
    weight *= primary.transparency;
    primary.ray.org = p;
    primary.ray.dir = normalize(add(mul(x,vx), mul(y,vy), vz));
    primary.ray.tnear = 1.001f*primary.ray.tfar;
    primary.ray.tfar = inf;
    primary.ray.geomID = -1;
    primary.ray.primID = -1;
    primary.ray.mask = -1;
    primary.ray.time = 0;
    primary.transparency = 0.0f;

/*    primary.ray.tfar = inf;
    primary.ray.geomID = -1;
    primary.ray.primID = -1;
    primary.transparency = 0.0f;*/
  }
  return color;
}
  
/* task that renders a single screen tile */
void renderTile(int taskIndex, int* pixels,
                     const int width,
                     const int height, 
                     const float time,
                     const Vec3f& vx, 
                     const Vec3f& vy, 
                     const Vec3f& vz, 
                     const Vec3f& p,
                     const int numTilesX, 
                     const int numTilesY)
{
  const int tileY = taskIndex / numTilesX;
  const int tileX = taskIndex - tileY * numTilesX;
  const int x0 = tileX * TILE_SIZE_X;
  const int x1 = min(x0+TILE_SIZE_X,width);
  const int y0 = tileY * TILE_SIZE_Y;
  const int y1 = min(y0+TILE_SIZE_Y,height);

  for (int y = y0; y<y1; y++) for (int x = x0; x<x1; x++)
  {
    /* calculate pixel color */
    Vec3f color = renderPixel(x,y,vx,vy,vz,p);

    /* write color to framebuffer */
    unsigned int r = (unsigned int) (255.0f * clamp(color.x,0.0f,1.0f));
    unsigned int g = (unsigned int) (255.0f * clamp(color.y,0.0f,1.0f));
    unsigned int b = (unsigned int) (255.0f * clamp(color.z,0.0f,1.0f));
    pixels[y*width+x] = (b << 16) + (g << 8) + r;
  }
}

/* called by the C++ code to render */
extern "C" void device_render (int* pixels,
                    const int width,
                    const int height,
                    const float time,
                    const Vec3f& vx, 
                    const Vec3f& vy, 
                    const Vec3f& vz, 
                    const Vec3f& p)
{
  const int numTilesX = (width +TILE_SIZE_X-1)/TILE_SIZE_X;
  const int numTilesY = (height+TILE_SIZE_Y-1)/TILE_SIZE_Y;
  launch_renderTile(numTilesX*numTilesY,pixels,width,height,time,vx,vy,vz,p,numTilesX,numTilesY); 
  rtcDebug();
}

/* called by the C++ code for cleanup */
extern "C" void device_cleanup ()
{
  rtcDeleteScene (g_scene);
  delete[] colors;
  rtcExit();
}
