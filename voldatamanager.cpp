#include "voldatamanager.h"
#include <sutil.h>
#include "random.h"

VolumeDataManager::VolumeDataManager ()
{
  vdata = nullptr;
  vtrft = nullptr;
}

VolumeDataManager::~VolumeDataManager ()
{
  DestroyVolume();
  DestroyTransferFunction();
}

float VolumeDataManager::GetNormalizedScalarValue (int x, int y, int z)
{
  return float(GetVolume()->SampleVolume(x, y, z)) / GetVolume()->GetMaxDensity();
}

void VolumeDataManager::ReadVolume (std::string filename)
{
  vr::VolumeRenderingReader vreader;
  DestroyVolume();
  vdata = vreader.ReadVolume(filename);
}

void VolumeDataManager::ReadTransferFunction (std::string filename)
{
  DestroyTransferFunction();
  vtrft = vr::ReadTransferFunction(filename);
}

optix::int3 VolumeDataManager::GetVolumeSizes ()
{
  return optix::make_int3(vdata->GetWidth(), vdata->GetHeight(), vdata->GetDepth());
}

optix::float3 VolumeDataManager::GetScaledBoundingBox ()
{
  return optix::make_float3((float)vdata->GetWidth()  * vdata->GetScaleX(),
                            (float)vdata->GetHeight() * vdata->GetScaleY(),
                            (float)vdata->GetDepth()  * vdata->GetScaleZ()
  );
}

vr::Volume* VolumeDataManager::GetVolume ()
{
  return vdata;
}

vr::TransferFunction* VolumeDataManager::GetTransferFunction ()
{
  return vtrft;
}

void VolumeDataManager::SetVolumeTextureSamplerToOptix (RTcontext context)
{
  optix::int3 vol_sizes = GetVolumeSizes();
  optix::float3 vol_scaledsizes = GetScaledBoundingBox();

  rtBufferCreate(context, RT_BUFFER_INPUT, &tex_buffer);
  rtBufferSetFormat(tex_buffer, RT_FORMAT_FLOAT);
  rtBufferSetSize3D(tex_buffer, vol_sizes.x, vol_sizes.y, vol_sizes.z);
  rtBufferMap(tex_buffer, &tex_data_ptr);

  float* tex_data = (float*)tex_data_ptr;
  int width = vol_sizes.x;
  int height = vol_sizes.y;
  int depth = vol_sizes.z;
  for (int z = 0; z < depth; z++)
    for (int y = 0; y < height; y++)
      for (int x = 0; x < width; x++)
        tex_data[x + (y * width) + (z * width * height)] = GetNormalizedScalarValue(x, y, z);

  rtBufferUnmap(tex_buffer);

  rtTextureSamplerCreate(context, &tex_sampler);
  rtTextureSamplerSetWrapMode(tex_sampler, 0, RT_WRAP_CLAMP_TO_EDGE);
  rtTextureSamplerSetWrapMode(tex_sampler, 1, RT_WRAP_CLAMP_TO_EDGE);
  rtTextureSamplerSetWrapMode(tex_sampler, 2, RT_WRAP_CLAMP_TO_EDGE);
  rtTextureSamplerSetFilteringModes(tex_sampler, RT_FILTER_LINEAR, RT_FILTER_LINEAR, RT_FILTER_NONE);
  rtTextureSamplerSetIndexingMode(tex_sampler, RT_TEXTURE_INDEX_NORMALIZED_COORDINATES);
  rtTextureSamplerSetReadMode(tex_sampler, RT_TEXTURE_READ_NORMALIZED_FLOAT);
  rtTextureSamplerSetMaxAnisotropy(tex_sampler, 1.0f);
  rtTextureSamplerSetMipLevelCount(tex_sampler, 1);
  rtTextureSamplerSetArraySize(tex_sampler, 1);
  rtTextureSamplerSetBuffer(tex_sampler, 0, 0, tex_buffer);
  rtContextDeclareVariable(context, "TexVolume", &tex_sampler_var);
  rtVariableSetObject(tex_sampler_var, tex_sampler);
}

void VolumeDataManager::SetTransferFunctionTextureSamplerToOptix (RTcontext context)
{
  rtBufferCreate(context, RT_BUFFER_INPUT, &tex_tf_buffer);
  rtBufferSetFormat(tex_tf_buffer, RT_FORMAT_FLOAT4);
  rtBufferSetSize1D(tex_tf_buffer, (int)GetVolume()->GetMaxDensity() + 1);
  rtBufferMap(tex_tf_buffer, &tex_tf_data_ptr);

  
  float* tex_tf_data = (float*)tex_tf_data_ptr;
  for (int i = 0; i <= (int)GetVolume()->GetMaxDensity(); i++)
  {
    glm::vec4 v = GetTransferFunction()->Get(i, GetVolume()->GetMaxDensity());
    tex_tf_data[i * 4 + 0] = v.x;
    tex_tf_data[i * 4 + 1] = v.y;
    tex_tf_data[i * 4 + 2] = v.z;
    tex_tf_data[i * 4 + 3] = GetTransferFunction()->GetExt(i, GetVolume()->GetMaxDensity());
  }

  rtBufferUnmap(tex_tf_buffer);

  rtTextureSamplerCreate(context, &tex_tf_sampler);
  rtTextureSamplerSetWrapMode(tex_tf_sampler, 0, RT_WRAP_CLAMP_TO_EDGE);
  rtTextureSamplerSetFilteringModes(tex_tf_sampler, RT_FILTER_LINEAR, RT_FILTER_LINEAR, RT_FILTER_NONE);
  rtTextureSamplerSetIndexingMode(tex_tf_sampler, RT_TEXTURE_INDEX_NORMALIZED_COORDINATES);
  rtTextureSamplerSetReadMode(tex_tf_sampler, RT_TEXTURE_READ_NORMALIZED_FLOAT);
  rtTextureSamplerSetMaxAnisotropy(tex_tf_sampler, 1.0f);
  rtTextureSamplerSetMipLevelCount(tex_tf_sampler, 1);
  rtTextureSamplerSetArraySize(tex_tf_sampler, 1);
  rtTextureSamplerSetBuffer(tex_tf_sampler, 0, 0, tex_tf_buffer);
  rtContextDeclareVariable(context, "TexTransferFunction", &tex_tf_sampler_var);
  rtVariableSetObject(tex_tf_sampler_var, tex_tf_sampler);
}

void VolumeDataManager::SetGradientTextureSamplerToOptix (RTcontext context)
{
  optix::int3 vol_sizes = GetVolumeSizes();
  optix::float3 vol_scaledsizes = GetScaledBoundingBox();
  int width = vol_sizes.x;
  int height = vol_sizes.y;
  int depth = vol_sizes.z;
  
  rtBufferCreate(context, RT_BUFFER_INPUT, &tex_gd_buffer);
  rtBufferSetFormat(tex_gd_buffer, RT_FORMAT_FLOAT4);
  rtBufferSetSize3D(tex_gd_buffer, vol_sizes.x, vol_sizes.y, vol_sizes.z);
  rtBufferMap(tex_gd_buffer, &tex_gd_data_ptr);
  
  glm::vec3* gvalues = vr::GenerateSobelFeldmanGradientTextureT(GetVolume());
  
  float* tex_gd_data = (float*)tex_gd_data_ptr ;
  for (int x = 0; x < width * height * depth; x++)
  {
    tex_gd_data[x * 4 + 0] = gvalues[x].x;
    tex_gd_data[x * 4 + 1] = gvalues[x].y;
    tex_gd_data[x * 4 + 2] = gvalues[x].z;
    tex_gd_data[x * 4 + 3] = gvalues[x].z;
  }
  delete[] gvalues;
  
  rtBufferUnmap(tex_gd_buffer);
  
  rtTextureSamplerCreate(context, &tex_gd_sampler);
  rtTextureSamplerSetWrapMode(tex_gd_sampler, 0, RT_WRAP_CLAMP_TO_EDGE);
  rtTextureSamplerSetWrapMode(tex_gd_sampler, 1, RT_WRAP_CLAMP_TO_EDGE);
  rtTextureSamplerSetWrapMode(tex_gd_sampler, 2, RT_WRAP_CLAMP_TO_EDGE);
  rtTextureSamplerSetFilteringModes(tex_gd_sampler, RT_FILTER_LINEAR, RT_FILTER_LINEAR, RT_FILTER_NONE);
  rtTextureSamplerSetIndexingMode(tex_gd_sampler, RT_TEXTURE_INDEX_NORMALIZED_COORDINATES);
  rtTextureSamplerSetReadMode(tex_gd_sampler, RT_TEXTURE_READ_NORMALIZED_FLOAT);
  rtTextureSamplerSetMaxAnisotropy(tex_gd_sampler, 1.0f);
  rtTextureSamplerSetMipLevelCount(tex_gd_sampler, 1);
  rtTextureSamplerSetArraySize(tex_gd_sampler, 1);
  rtTextureSamplerSetBuffer(tex_gd_sampler, 0, 0, tex_gd_buffer);
  rtContextDeclareVariable(context, "TexGradientVolume", &tex_gd_sampler_var);
  rtVariableSetObject(tex_gd_sampler_var, tex_gd_sampler);
}

void VolumeDataManager::DestroyVolume ()
{
  if (vdata != nullptr)
  {
    delete vdata;
    vdata = nullptr;
  }
}

void VolumeDataManager::DestroyTransferFunction ()
{
  if (vtrft != nullptr)
  {
    delete vtrft;
    vtrft = nullptr;
  }
}