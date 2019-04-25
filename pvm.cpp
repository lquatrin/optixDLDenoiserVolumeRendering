#include "pvm.h"
#include <fstream>
#include <iostream>
#include <cmath>
#include <cassert>

#include "glm/glm.hpp"
#include "glm/exponential.hpp"
#include "glm/ext.hpp"
#include "glm/common.hpp"

#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <cerrno>

#define DDS_MAXSTR (256)

#define DDS_BLOCKSIZE (1<<20)
#define DDS_INTERLEAVE (1<<24)

#define DDS_RL (7)

Pvm::Pvm (const char *file_name)
{
  std::string filename(file_name);

  DDSV3 ddsloader;
  unsigned char* raw = ddsloader.readPVMvolume(filename.c_str(),
                                               &width, &height, &depth,
                                               &components,
                                               &scalex, &scaley, &scalez);
  
  printf("Sizes: [%d %d %d]\n", width, height, depth);
  printf("Scale: [%.4f %.4f %.4f]\n", scalex, scaley, scalez);
  printf("Components: %d\n", components);

  pvm_data = PostProcessData(raw);
  
  free(raw);
}

Pvm::~Pvm ()
{
  if (pvm_data) delete[] pvm_data;
  pvm_data = NULL;
}

float* Pvm::GetData ()
{
  return pvm_data;
}

float* Pvm::GenerateNormalizeData ()
{
  float* custom_data = new float[width*height*depth];

  float max_density_value = pow(2, components * 8) - 1;
  // min is = 0

  // Reescale between 0 ~ 1 using max_density_value
  for (int i = 0; i < width*height*depth; i++)
    custom_data[i] = pvm_data[i] / max_density_value;

  return custom_data;
}

float* Pvm::GenerateReescaledMinMaxData (bool normalized, float* fmin, float* fmax)
{
  float max_density_value = pow(2, components * 8) - 1;
  float min = max_density_value;
  float max = 0;

  for (int i = 0; i < width*height*depth; i++)
  {
    min = glm::min(min, pvm_data[i]);
    max = glm::max(max, pvm_data[i]);
  }
  
  if (fmin) *fmin = min;
  if (fmax) *fmax = max;

  float* custom_data = new float[width*height*depth];
  
  // First we reescale between 0 ~ 1 using min and max value found
  for (int i = 0; i < width*height*depth; i++)
    custom_data[i] = (pvm_data[i] - min) / (max - min);
  
  // if normalized is FALSE, reescale by the max_density_value
  if (!normalized) 
  {
    for (int i = 0; i < width*height*depth; i++)
      custom_data[i] = custom_data[i] * max_density_value;
  }

  return custom_data;
}

void Pvm::GetScale (float* sx, float* sy, float* sz)
{
  *sx = scalex;
  *sy = scaley;
  *sz = scalez;
}

void Pvm::GetDimensions (unsigned int* _width, unsigned int* _height, unsigned int* _depth)
{
  *_width  = width;
  *_height = height;
  *_depth  = depth;
}

int Pvm::GetComponents ()
{
  return components;
}

void Pvm::TransformData (unsigned int dst_bytes_per_voxel)
{
  float old_max_density_value = pow(2, components * 8) - 1;
  float new_max_density_value = pow(2, dst_bytes_per_voxel * 8) - 1;

  // First we reescale between 0~1
  for (int i = 0; i < width*height*depth; i++)
  {
    pvm_data[i] = (pvm_data[i] / old_max_density_value) * new_max_density_value;
  }

  components = dst_bytes_per_voxel;
}

float* Pvm::PostProcessData (unsigned char* data)
{
  float* prc_data = new float[width*height*depth];

  for (int i = 0; i < width*height*depth; i++) 
  {
    if (components == 1)
    {
      prc_data[i] = (float)data[i];
    }
    else
    {
      unsigned short v1 = data[i * 2];
      unsigned short v2 = data[i * 2 + 1];
      
      prc_data[i] = (float)(v1 * 256 + v2); // == v1 << 8 | v2
      assert(v1 << 8 | v2 == v1 * 256 + v2);
    }
  }

  return prc_data;
}

////////////////////////////////////
// DDSV3 class methods
////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////
// PVM/DDS reader (8/16 bits):

// Copyright:
// (c) by Stefan Roettger, licensed under GPL 2+
// The Volume Library: http://www9.informatik.uni-erlangen.de/External/vollib/
// V^3 Package code: https://code.google.com/p/vvv/

// read a compressed PVM volume
unsigned char* DDSV3::readPVMvolume (const char* filename,
                                     unsigned int* width, unsigned int* height, unsigned int* depth,
                                     unsigned int* components,
                                     float* scalex, float* scaley, float* scalez,
                                     unsigned char** description,
                                     unsigned char** courtesy,
                                     unsigned char** parameter,
                                     unsigned char** comment)
{
  unsigned char* data;
  unsigned char* ptr;
  unsigned int bytes, numc;

  int version = 1;

  unsigned char* volume;

  float sx = 1.0f, sy = 1.0f, sz = 1.0f;

  unsigned int len1 = 0, len2 = 0, len3 = 0, len4 = 0;

  if ((data = readDDSfile(filename, &bytes)) == NULL)
    if ((data = readRAWfile(filename, &bytes)) == NULL) return(NULL);

  if (bytes<5) return(NULL);

  if ((data = (unsigned char*)realloc(data, bytes + 1)) == NULL) ERRORMSG();
  data[bytes] = '\0';

  if (strncmp((char*)data, "PVM\n", 4) != 0)
  {
    if (strncmp((char*)data, "PVM2\n", 5) == 0) version = 2;
    else if (strncmp((char*)data, "PVM3\n", 5) == 0) version = 3;
    else return(NULL);

    ptr = &data[5];
    if (sscanf_s((char*)ptr, "%d %d %d\n%g %g %g\n", width, height, depth, &sx, &sy, &sz) != 6) ERRORMSG();
    if (*width < 1 || *height < 1 || *depth < 1 || sx <= 0.0f || sy <= 0.0f || sz <= 0.0f) ERRORMSG();
    ptr = (unsigned char*)strchr((char*)ptr, '\n') + 1;
  }
  else
  {
    ptr = &data[4];
    while (*ptr == '#')
      while (*ptr++ != '\n');

    if (sscanf_s((char*)ptr, "%d %d %d\n", width, height, depth) != 3) ERRORMSG();
    if (*width < 1 || *height < 1 || *depth < 1) ERRORMSG();
  }

  if (scalex != NULL && scaley != NULL && scalez != NULL)
  {
    *scalex = sx;
    *scaley = sy;
    *scalez = sz;
  }

  ptr = (unsigned char*)strchr((char*)ptr, '\n') + 1;
  if (sscanf_s((char *)ptr, "%d\n", &numc) != 1) ERRORMSG();
  if (numc<1) ERRORMSG();

  if (components != NULL) *components = numc;
  else if (numc != 1) ERRORMSG();

  ptr = (unsigned char*)strchr((char*)ptr, '\n') + 1;
  if (version == 3) len1 = strlen((char*)(ptr + (*width)*(*height)*(*depth)*numc)) + 1;
  if (version == 3) len2 = strlen((char*)(ptr + (*width)*(*height)*(*depth)*numc + len1)) + 1;
  if (version == 3) len3 = strlen((char*)(ptr + (*width)*(*height)*(*depth)*numc + len1 + len2)) + 1;
  if (version == 3) len4 = strlen((char*)(ptr + (*width)*(*height)*(*depth)*numc + len1 + len2 + len3)) + 1;
  if ((volume = (unsigned char*)malloc((*width)*(*height)*(*depth)*numc + len1 + len2 + len3 + len4)) == NULL) ERRORMSG();
  if (data + bytes != ptr + (*width)*(*height)*(*depth)*numc + len1 + len2 + len3 + len4) ERRORMSG();

  memcpy(volume, ptr, (*width)*(*height)*(*depth)*numc + len1 + len2 + len3 + len4);
  free(data);


  if (description != NULL)
    if (len1 > 1) *description = volume + (*width)*(*height)*(*depth)*numc;
    else *description = NULL;
  //if (description != NULL) {
  //  if (len1>1) *description = volume + (*width)*(*height)*(*depth)*numc;
  //}
  //else {
  //  //*description=NULL;
  //}

  if (courtesy != NULL)
    if (len2 > 1) *courtesy = volume + (*width)*(*height)*(*depth)*numc + len1;
    else *courtesy = NULL;
  //if (courtesy != NULL) {
  //  if (len2>1) *courtesy = volume + (*width)*(*height)*(*depth)*numc + len1;
  //}
  //else {
  //  //*courtesy=NULL;
  //}

  if (parameter != NULL)
    if (len3 > 1) *parameter = volume + (*width)*(*height)*(*depth)*numc + len1 + len2;
    else *parameter = NULL;
  //if (parameter != NULL) {
  //  if (len3>1) *parameter = volume + (*width)*(*height)*(*depth)*numc + len1 + len2;
  //}
  //else {
  //  //*parameter=NULL;
  //}

  if (comment != NULL)
    if (len4 > 1) *comment = volume + (*width)*(*height)*(*depth)*numc + len1 + len2 + len3;
    else *comment = NULL;
  //if (comment != NULL) {
  //  if (len4>1) *comment = volume + (*width)*(*height)*(*depth)*numc + len1 + len2 + len3;
  //}
  //else {
  //  //*comment=NULL;
  //}

  return(volume);
}

// read a possibly compressed PNM image
unsigned char* DDSV3::readPNMimage (const char* filename,
                                    unsigned int* width,
                                    unsigned int* height,
                                    unsigned int* components)
{
  const int maxstr = 100;

  char str[maxstr];

  unsigned char *data, *ptr1, *ptr2;
  unsigned int bytes;

  int pnmtype, maxval;
  unsigned char* image;

  if ((data = readDDSfile(filename, &bytes)) == NULL)
    if ((data = readRAWfile(filename, &bytes)) == NULL) return(NULL);

  if (bytes<4) return(NULL);

  memcpy(str, data, 3);
  str[3] = '\0';

  if (sscanf_s(str, "P%1d\n", &pnmtype) != 1) return(NULL);

  ptr1 = data + 3;
  while (*ptr1 == '\n' || *ptr1 == '#')
  {
    while (*ptr1 == '\n')
      if (++ptr1 >= data + bytes) ERRORMSG();
    while (*ptr1 == '#')
      if (++ptr1 >= data + bytes) ERRORMSG();
      else
        while (*ptr1 != '\n')
          if (++ptr1 >= data + bytes) ERRORMSG();
  }

  ptr2 = ptr1;
  while (*ptr2 != '\n' && *ptr2 != ' ')
    if (++ptr2 >= data + bytes) ERRORMSG();
  if (++ptr2 >= data + bytes) ERRORMSG();
  while (*ptr2 != '\n' && *ptr2 != ' ')
    if (++ptr2 >= data + bytes) ERRORMSG();
  if (++ptr2 >= data + bytes) ERRORMSG();
  while (*ptr2 != '\n' && *ptr2 != ' ')
    if (++ptr2 >= data + bytes) ERRORMSG();
  if (++ptr2 >= data + bytes) ERRORMSG();

  if (ptr2 - ptr1 >= maxstr) ERRORMSG();
  memcpy(str, ptr1, ptr2 - ptr1);
  str[ptr2 - ptr1] = '\0';

  if (sscanf_s(str, "%d %d\n%d\n", width, height, &maxval) != 3) ERRORMSG();

  if (*width<1 || *height<1) ERRORMSG();

  if (pnmtype == 5 && maxval == 255) *components = 1;
  else if (pnmtype == 5 && (maxval == 32767 || maxval == 65535)) *components = 2;
  else if (pnmtype == 6 && maxval == 255) *components = 3;
  else ERRORMSG();

  if ((image = (unsigned char*)malloc((*width)*(*height)*(*components))) == NULL) ERRORMSG();
  if (data + bytes != ptr2 + (*width)*(*height)*(*components)) ERRORMSG();

  memcpy(image, ptr2, (*width)*(*height)*(*components));
  free(data);

  return(image);
}

// decode a Differential Data Stream
void DDSV3::DDS_decode (unsigned char* chunk, unsigned int size,
                        unsigned char** data, unsigned int* bytes,
                        unsigned int block)
{
  unsigned int skip, strip;

  unsigned char *ptr1, *ptr2;

  unsigned int cnt, cnt1, cnt2;
  int bits, act;

  DDS_initbuffer();

  DDS_clearbits();
  DDS_loadbits(chunk, size);

  skip = DDS_readbits(2) + 1;
  strip = DDS_readbits(16) + 1;

  ptr1 = ptr2 = NULL;
  cnt = act = 0;

  while ((cnt1 = DDS_readbits(DDS_RL)) != 0)
  {
    bits = DDS_decode(DDS_readbits(3));

    for (cnt2 = 0; cnt2<cnt1; cnt2++)
    {
      if (strip == 1 || cnt <= strip) act += DDS_readbits(bits) - (1 << bits) / 2;
      else act += *(ptr2 - strip) - *(ptr2 - strip - 1) + DDS_readbits(bits) - (1 << bits) / 2;

      while (act<0) act += 256;
      while (act>255) act -= 256;

      if ((cnt&(DDS_BLOCKSIZE - 1)) == 0) {
        if (ptr1 == NULL)
        {
          if ((ptr1 = (unsigned char*)malloc(DDS_BLOCKSIZE)) == NULL) MEMERROR();
          ptr2 = ptr1;
        }
        else
        {
          if ((ptr1 = (unsigned char*)realloc(ptr1, cnt + DDS_BLOCKSIZE)) == NULL) MEMERROR();
          ptr2 = &ptr1[cnt];
        }
      }

      *ptr2++ = act;
      cnt++;
    }
  }

  if (ptr1 != NULL)
    if ((ptr1 = (uint8_t *)realloc(ptr1, cnt)) == NULL) MEMERROR();

  DDS_interleave(ptr1, cnt, skip, block);

  *data = ptr1;
  *bytes = cnt;
}

// interleave a byte stream
void DDSV3::DDS_interleave (unsigned char* data,
                            unsigned int bytes,
                            unsigned int skip,
                            unsigned int block)
{
  DDS_deinterleave(data, bytes, skip, block, true);
}

// deinterleave a byte stream
void DDSV3::DDS_deinterleave (unsigned char* data,
                              unsigned int bytes,
                              unsigned int skip,
                              unsigned int block,
                              bool restore)
{
  unsigned int i, j, k;

  unsigned char *data2, *ptr;

  if (skip <= 1) return;

  if (block == 0)
  {
    if ((data2 = (unsigned char*)malloc(bytes)) == NULL) MEMERROR();

    if (!restore)
      for (ptr = data2, i = 0; i<skip; i++)
        for (j = i; j<bytes; j += skip) *ptr++ = data[j];
    else
      for (ptr = data, i = 0; i<skip; i++)
        for (j = i; j<bytes; j += skip) data2[j] = *ptr++;

    memcpy(data, data2, bytes);
  }
  else
  {
    if ((data2 = (unsigned char*)malloc((bytes < skip*block) ? bytes : skip*block)) == NULL) MEMERROR();

    if (!restore)
    {
      for (k = 0; k < bytes / skip / block; k++)
      {
        for (ptr = data2, i = 0; i < skip; i++)
          for (j = i; j<skip*block; j += skip) *ptr++ = data[k*skip*block + j];

        memcpy(data + k*skip*block, data2, skip*block);
      }

      for (ptr = data2, i = 0; i<skip; i++)
        for (j = i; j<bytes - k*skip*block; j += skip) *ptr++ = data[k*skip*block + j];

      memcpy(data + k*skip*block, data2, bytes - k*skip*block);
    }
    else
    {
      for (k = 0; k<bytes / skip / block; k++)
      {
        for (ptr = data + k*skip*block, i = 0; i<skip; i++)
          for (j = i; j<skip*block; j += skip) data2[j] = *ptr++;

        memcpy(data + k*skip*block, data2, skip*block);
      }

      for (ptr = data + k*skip*block, i = 0; i<skip; i++)
        for (j = i; j<bytes - k*skip*block; j += skip) data2[j] = *ptr++;

      memcpy(data + k*skip*block, data2, bytes - k*skip*block);
    }
  }

  free(data2);
}

// read a Differential Data Stream
unsigned char* DDSV3::readDDSfile (const char *filename, unsigned int *bytes)
{
  char DDS_ID[] = "DDS v3d\n";
  char DDS_ID2[] = "DDS v3e\n";

  int version = 1;

  FILE *file;
  errno_t err;

  int cnt;

  unsigned char *chunk, *data;
  unsigned int size;


  if ((err = fopen_s(&file, filename, "rb")) != 0) return(NULL);

  for (cnt = 0; DDS_ID[cnt] != '\0'; cnt++)
  {
    if (fgetc(file) != DDS_ID[cnt])
    {
      fclose(file);
      version = 0;
      break;
    }
  }

  if (version == 0)
  {
    if ((err = fopen_s(&file, filename, "rb")) != 0) return(NULL);

    for (cnt = 0; DDS_ID2[cnt] != '\0'; cnt++)
    {
      if (fgetc(file) != DDS_ID2[cnt])
      {
        fclose(file);
        return(NULL);
      }
    }

    version = 2;
  }

  if ((chunk = readRAWfiled(file, &size)) == NULL) IOERROR();

  fclose(file);

  DDS_decode(chunk, size, &data, bytes, version == 1 ? 0 : DDS_INTERLEAVE);

  free(chunk);

  return(data);
}

// read from a RAW file
unsigned char* DDSV3::readRAWfiled (FILE *file, unsigned int *bytes)
{
  unsigned char* data;
  unsigned int cnt, blkcnt;

  data = NULL;
  cnt = 0;

  do
  {
    if (data == NULL)
    {
      if ((data = (unsigned char*)malloc(DDS_BLOCKSIZE)) == NULL) MEMERROR();
    }
    else
      if ((data = (unsigned char*)realloc(data, cnt + DDS_BLOCKSIZE)) == NULL) MEMERROR();

    blkcnt = fread(&data[cnt], 1, DDS_BLOCKSIZE, file);
    cnt += blkcnt;
  } while (blkcnt == DDS_BLOCKSIZE);

  if (cnt == 0)
  {
    free(data);
    return(NULL);
  }

  if ((data = (unsigned char*)realloc(data, cnt)) == NULL) MEMERROR();

  *bytes = cnt;

  return(data);
}

// read a RAW file
unsigned char* DDSV3::readRAWfile (const char *filename, unsigned int *bytes)
{
  FILE* file;
  errno_t err;

  unsigned char* data;

  if ((err = fopen_s(&file, filename, "rb")) != 0) return(NULL);

  data = readRAWfiled(file, bytes);

  fclose(file);

  return(data);
}