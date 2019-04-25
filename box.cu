/* 
 * Copyright (c) 2018, NVIDIA CORPORATION. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *  * Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *  * Neither the name of NVIDIA CORPORATION nor the names of its
 *    contributors may be used to endorse or promote products derived
 *    from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ``AS IS'' AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
 * OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <optix.h>
#include <optixu/optixu_math_namespace.h>
#include <optixu/optixu_matrix_namespace.h>
#include <optixu/optixu_aabb_namespace.h>

using namespace optix;

// Box attributes
rtDeclareVariable(float3, boxmin, , );
rtDeclareVariable(float3, boxmax, , );

// Current ray state
rtDeclareVariable(optix::Ray, ray, rtCurrentRay, );

// Atribute to be used in closest_hit_radiance0
/**
 * . Attributes MUST be modified between rtPotentialIntersection and rtReportIntersection calls.
**/
rtDeclareVariable(float3, shadingnormal, attribute shadingnormal, );
rtDeclareVariable(float2, tvalue, attribute tvalue, );

static __device__ float3 boxnormal (float t, float3 t0, float3 t1)
{
  float3 neg = make_float3(t == t0.x ? 1 : 0,
                           t == t0.y ? 1 : 0,
                           t == t0.z ? 1 : 0);

  float3 pos = make_float3(t == t1.x ? 1 : 0,
                           t == t1.y ? 1 : 0,
                           t == t1.z ? 1 : 0);
  
  return pos - neg;
}

/**
 * . Attribute Variables may only be set in between these two function calls. This
 * ensures that the Variables always represent the values of the closest hit yet found.
 * If the Variables are read inside the Closest Hit Program it is guaranteed, that they
 * represent the values of the closest intersection.
 *
 * . Typically Attribute Variables are used to communicate intersection specific information 
 * to the Closest Hit or Any Hit Programs such as surface normal vectors or texture coordinates.
**/
RT_PROGRAM void box_intersect (int primIdx)
{
  float3 t0   = (boxmin - ray.origin) / ray.direction;
  float3 t1   = (boxmax - ray.origin) / ray.direction;
  float3 near = fminf(t0, t1);
  float3 far  = fmaxf(t0, t1);
  float tmin  = fmaxf(near);
  float tmax  = fminf(far);


  if(tmin <= tmax)
  {
    bool check_second = true;
    // returns true if the t-value lies inside the allowed range of the ray
    if (rtPotentialIntersection(tmin))
    {
       tvalue = make_float2(tmin, tmax);
       shadingnormal = boxnormal(tmin, t0, t1);
       // is called with the material index 
       // what material should be used on this part of the geometry
       if (rtReportIntersection(0))
         check_second = false;
    } 
    if(check_second)
    {
      if (rtPotentialIntersection(tmax))
      {
        tvalue = make_float2(tmin, tmax);
        shadingnormal = boxnormal(tmax, t0, t1);
        rtReportIntersection(0);
      }
    }
  }
}

/**
 * . Geometry objects also need a Bounding Box Program assigned. It must return an
 * axis aligned bounding box that fully encloses the primitive at the primitive index
 * given as an argument. It is used by the Acceleration Structures while traversing
 * and for building the acceleration tree. While a fast implementation is desirable, the
 * accuracy of the bounding box is also important to build good quality acceleration
 * trees. Accurate means that the bounding box should be as small as possible while
 * still fully enclosing the primitive.
**/
RT_PROGRAM void box_bounds (int, float result[6])
{
  optix::Aabb* aabb = (optix::Aabb*)result;
  aabb->set(boxmin, boxmax);
}
