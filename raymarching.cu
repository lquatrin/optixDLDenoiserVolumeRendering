#include "tutorial.h"
#include "random.h"

// Internally providade variables from OptiX
rtDeclareVariable(uint2,               launch_index, rtLaunchIndex, );
rtDeclareVariable(PerRayData_radiance, prd_radiance, rtPayload    , );
rtDeclareVariable(optix::Ray,          ray,          rtCurrentRay , );

// Atribute set on box and parallelogram
rtDeclareVariable(float3, shadingnormal, attribute shadingnormal, );
rtDeclareVariable(float2,        tvalue, attribute tvalue       , );

rtDeclareVariable(float,    scene_epsilon, , );
rtDeclareVariable(rtObject, top_object   , , );

// Pinhole camera implementation
rtDeclareVariable(float3,        eye, , );
rtDeclareVariable(float3,        U, , );
rtDeclareVariable(float3,        V, , );
rtDeclareVariable(float3,        W, , );
rtDeclareVariable(float3,        bad_color, , );
rtDeclareVariable(Matrix3x3, normal_matrix, , );
rtDeclareVariable(float3, VolSizes, , );

rtDeclareVariable(unsigned int, frame_number, , );

rtBuffer<float4, 2> out_render_buffer;
rtBuffer<float4, 2> out_albedo_buffer;
rtBuffer<float4, 2> out_normal_buffer;

rtTextureSampler<float, 3> TexVolume;
rtTextureSampler<float4, 1> TexTransferFunction;
rtTextureSampler<float4, 3> TexGradientVolume;

RT_PROGRAM void RayGenerationProgram ()
{
  size_t2 screen = out_render_buffer.size();

  float2 d = make_float2(launch_index) / make_float2(screen) * 2.f - 1.f;
  float3 ray_origin = eye;
  float3 ray_direction = normalize(d.x*U + d.y*V + W);

  optix::Ray ray(ray_origin, ray_direction, RADIANCE_RAY_TYPE, scene_epsilon );

  unsigned int seed = tea<16>(screen.x*launch_index.y + launch_index.x, frame_number);

  PerRayData_radiance prd;
  prd.importance = 1.f;
  prd.depth = 0;
  prd.alpha_ray = 0.0f;
  prd.normal_output = make_float3(0.0f, 0.0f, 0.0f);
  prd.normal_total_weight = 0.0f;
  prd.albedo_output = make_float3(0.0f, 0.0f, 0.0f);
  prd.id_noise = frame_number % 1024;
  prd.seed = seed;

  // trace a ray giving the object hierarchy group that will be traced
  rtTrace(top_object, ray, prd);

  // Only update the current state of the output buffer
  if (frame_number > 1)
  {
    //*
    float a = 1.0f / (float)frame_number;

    // update output normalized color
    float3 normalized_color = prd.result * prd.alpha_ray + make_float3(0.0f, 0.0f, 0.0f) * (1.0f - prd.alpha_ray);
    float4 curr_color = make_float4(normalized_color.x, normalized_color.y, normalized_color.z, prd.alpha_ray);
    float4 old_color_f4 = out_render_buffer[launch_index];
    out_render_buffer[launch_index] = lerp(old_color_f4, curr_color, a);
    
    // update output normal
    float4 curr_normal = make_float4(prd.normal_output / prd.normal_total_weight, prd.alpha_ray);
    float4 old_normal_f4 = out_normal_buffer[launch_index];
    out_normal_buffer[launch_index] = lerp(old_normal_f4, curr_normal, a);

    // update output albedo
    float4 curr_albedo = make_float4(prd.albedo_output.x, prd.albedo_output.y, prd.albedo_output.z, prd.alpha_ray);
    float4 old_albedo_f4 = out_albedo_buffer[launch_index];
    out_albedo_buffer[launch_index] = lerp(old_albedo_f4, curr_albedo, a);
    // */
  }
  // First frame!
  else
  {
    // output normalized color
    float3 normalized_color = prd.result * prd.alpha_ray + make_float3(0.0f, 0.0f, 0.0f) * (1.0f - prd.alpha_ray);
    out_render_buffer[launch_index] = make_float4(normalized_color.x, normalized_color.y, normalized_color.z, prd.alpha_ray);

    // The current normal is the amount of normals averaged by their weights
    float3 normals_res = prd.normal_output / prd.normal_total_weight;
    out_normal_buffer[launch_index] = make_float4(normals_res.x, normals_res.y, normals_res.z, prd.alpha_ray);

    // Color without the gradient shading
    out_albedo_buffer[launch_index] = make_float4(prd.albedo_output.x, prd.albedo_output.y, prd.albedo_output.z, prd.alpha_ray);
  }
}

/**
 * If the ray misses and finds no intersection. we may implement the background.
**/
// Returns solid color for miss rays
rtDeclareVariable(float3, bg_color, , );
RT_PROGRAM void miss ()
{
  prd_radiance.result = bg_color;
}
 
const int NUMBER_OF_CASTED_RAYS = 1;
const float DO_INITIAL_STEP = 2.0f;
const float DO_MAX_DISTANCE = 80.0f;
float EvaluateDirectionalOcclusion (float3 Wpos)
{
  float3 VolBBMin = -VolSizes * 0.5f;
  float3 VolBBMax =  VolSizes * 0.5f;

  float S_Vt = 0.0;
  float S_Wt = 0.0;

  int iray = 0;
  while (iray < NUMBER_OF_CASTED_RAYS)
  {
    // get noise values for direction [0...1, 
    float3 sn = make_float3(rnd(prd_radiance.seed), (rnd(prd_radiance.seed) * 2.0f) - 1.0f, (rnd(prd_radiance.seed) * 2.0f) - 1.0f);

    float3 fake_light_direction = normalize(-W * sn.z + U * sn.y + V * sn.x);
    
    // Let's evaluate shadow
    float Vt = 1.0f;
    float s = DO_INITIAL_STEP;
    while (s < DO_MAX_DISTANCE)
    {
      float d = min(0.5f, DO_MAX_DISTANCE - s);
    
      // Check if we're out of the bounding box
      float3 spos = Wpos + fake_light_direction * s;

      if (!(spos.x > VolBBMin.x && spos.x < VolBBMax.x &&
            spos.y > VolBBMin.y && spos.y < VolBBMax.y &&
            spos.z > VolBBMin.z && spos.z < VolBBMax.z))
        break;
          
      float3 txpos = (spos + (VolSizes * 0.5f)) / VolSizes;
      float st = tex1D(TexTransferFunction, tex3D(TexVolume, txpos.x, txpos.y, txpos.z)).w;

      Vt = Vt * exp(-st * d);
      if (Vt < 0.05) break;

      s = s + d;
    }

    float wtc = dot(normalize(-W), fake_light_direction);

    S_Vt = S_Vt + Vt * wtc;
    S_Wt = S_Wt + wtc;
    
    iray = iray + 1;
    prd_radiance.id_noise = prd_radiance.id_noise + 1;
  }
  
  return S_Vt / S_Wt;
}

const bool apply_phong = true;
float3 ShadeWithPhong (float3 Li, float ka, float kd, float ks, float3 pos, float expt, float3 Wpos, float ta)
{
  float3 WorldLightingPos = eye;
  float4 tex_gt_n = tex3D(TexGradientVolume, pos.x, pos.y, pos.z);
  
  if (apply_phong)
  {
    float3 gradient_normal = make_float3(tex_gt_n.x, tex_gt_n.y, tex_gt_n.z);
    
    if (gradient_normal.x != 0 && gradient_normal.y != 0 && gradient_normal.z != 0)
    {
      gradient_normal = normalize(gradient_normal);
    
      float3 light_direction = normalize(WorldLightingPos - Wpos);
      float3 eye_direction = normalize(eye - Wpos);
      float3 halfway_vector = normalize(eye_direction + light_direction);
      
      float dot_diff = max(0.0f, dot(gradient_normal, light_direction));
      float dot_spec = max(0.0f, dot(halfway_vector, gradient_normal));

      float3 Is = make_float3(1, 1, 1);
      float Nshininess = 20.0f;
      
      float clr_norm = (1.0 / (ka + kd + ks));
      Li = clr_norm * (Li * ka * ta + Li * kd * dot_diff  + ks * Is * pow(dot_spec, Nshininess));
    }
  }
  else
  {
    Li = Li * ta;
  }

  /////////////////////////////////////////////////////////////////////////////
  // Accumulate gradient
  float3 normal_g = make_float3(tex_gt_n.x, tex_gt_n.y, tex_gt_n.z);
  float3 normal_eyespace = (length(normal_g) > 0.f) ? normalize(normal_matrix * normal_g) : make_float3(0., 0., 1.);
  prd_radiance.normal_output += normal_eyespace * (1.0 - expt);
  prd_radiance.normal_total_weight += (1.0 - expt);

  return Li;
}

// Returns shading normal as the surface shading result
/**
 * . The intersection details such as texture coordinates should be 
 * communicated to the Closest Hit Program via Attribute Variables 
 * that were calculated in the intersection program.
 *
 * . rtTrace may be recursively called to implement shadows or reflection.
 *
 * . result may be stored in the payload of the ray.
 *
 * . Apply ray marching here!
**/
const float Ka = 0.60f;
const float Kd = 0.15f;
const float Ks = 0.25f;
RT_PROGRAM void closest_hit_radiance0 ()
{
  float3 r0 = ray.origin;
  float3 rd = ray.direction;

  float step_pos = rnd(prd_radiance.seed);// ((frame_number * 10) % 100) * 0.01;
  float tmin_pos_epsilon = 0.0001f;
  float3 pt = (r0 + rd * (tvalue.x + tmin_pos_epsilon));

  float integration_step = 0.5f;

  float3 L = make_float3(0.0f);
  float3 Lalbedo = make_float3(0.0f);
  float T = 1.0f;

  float s = tvalue.x;
  while (s < tvalue.y)
  {
    float d = min(integration_step, tvalue.y - s);

    float3 Wpos = r0 + rd * (s + d * step_pos);
    float3 pos = (Wpos + (VolSizes * 0.5f)) / VolSizes;

    float4 Ltex = tex1D(TexTransferFunction, tex3D(TexVolume, pos.x, pos.y, pos.z));
    float3 Li = make_float3(Ltex.x, Ltex.y, Ltex.z);
    float expt = exp(-Ltex.w * d);
    
    if (Ltex.w > 0.0f)
    {
      float S_Vt = EvaluateDirectionalOcclusion(Wpos);
      Li = ShadeWithPhong(Li, Ka, Kd, Ks, pos, expt, Wpos, S_Vt);
    }

    /////////////////////////////////////////////////////////////////////////////
    // Accumulate albedo
    Lalbedo = Lalbedo + T * (1.0 - expt) * make_float3(Ltex.x, Ltex.y, Ltex.z);

    L = L + T * (1.0 - expt) * Li;
    T = T * expt;
    if ((1.0 - T) > 0.95) break;

    s = s + d;
  }

  prd_radiance.result = L;
  prd_radiance.alpha_ray = 1.0 - T;

  prd_radiance.albedo_output = Lalbedo;
}

// any hit program
/**
 * . as opposed to the Closest Hit Program, it might be
 * called multiple times for a single ray cast.
 *
 * . The intersections for which the program is executed may not be ordered 
 * along the ray, but eventually all intersections can be enumerated by 
 * calling rtIgnoreIntersection on each of them.
 *
 * . The Any Hit Program can be used if the application requires to perform 
 * actions at each surface intersection. Rays can also be terminated with 
 * rtTerminateRay, then * the trace ends and no further executions of the Any 
 * Hit Program may be invoked. This can be used if only the knowledge whether
 * the ray hits anything or not is needed, for example in shadow rays.
**/
// Set pixel to solid color upon failur
RT_PROGRAM void exception ()
{
  out_render_buffer[launch_index] = make_float4(bad_color.x, bad_color.y, bad_color.z, 1.0f);
}
