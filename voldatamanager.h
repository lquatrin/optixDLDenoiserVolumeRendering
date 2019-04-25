#ifndef VOLUME_DATA_MANAGER_H
#define VOLUME_DATA_MANAGER_H

#include <optixu/optixpp_namespace.h>
#include <optixu/optixu_math_stream_namespace.h>

#include "reader.h"
#include "volume.h"
#include "transferfunction.h"

class VolumeDataManager {
public:
  VolumeDataManager ();
  ~VolumeDataManager ();

  float GetNormalizedScalarValue (int x, int y, int z);

  void ReadVolume (std::string filename);
  void ReadTransferFunction (std::string filename);

  optix::int3 GetVolumeSizes ();
  optix::float3 GetScaledBoundingBox ();

  vr::Volume* GetVolume ();
  vr::TransferFunction* GetTransferFunction ();

  void SetVolumeTextureSamplerToOptix (RTcontext context);
  void SetTransferFunctionTextureSamplerToOptix (RTcontext context);
  void SetGradientTextureSamplerToOptix (RTcontext context);


  RTbuffer tex_gd_buffer;
  RTtexturesampler tex_gd_sampler;
  void* tex_gd_data_ptr;
  RTvariable tex_gd_sampler_var;
protected:

private:
  void DestroyVolume ();

  vr::Volume* vdata;
  RTbuffer  tex_buffer;
  RTtexturesampler tex_sampler;
  void* tex_data_ptr;
  RTvariable tex_sampler_var;

  void DestroyTransferFunction ();

  vr::TransferFunction* vtrft;
  RTbuffer  tex_tf_buffer;
  RTtexturesampler tex_tf_sampler;
  void* tex_tf_data_ptr;
  RTvariable tex_tf_sampler_var;
};

#endif