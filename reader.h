#ifndef VOLREND_READER_H
#define VOLREND_READER_H

#include "volume.h"

#include <iostream>
#include <cstring>
#include <cstdlib>
#include <string>

#include "pvm.h"
#include "transferfunction.h"

namespace vr
{
  double* ReadVolvisRaw (std::string volfilename, size_t size, int w, int h, int d);
  Volume* ReadRawFile (std::string filepath);
  Volume* ReadSynFile (std::string filename);

  TransferFunction* ReadTransferFunction (std::string file);
  TransferFunction* ReadTransferFunction_tf1d (std::string file);

  class VolumeRenderingReader
  {
  public:
    VolumeRenderingReader () {}
    
    ~VolumeRenderingReader () {}

    Volume* ReadVolume (std::string filepath)
    {
      Volume* ret = NULL;
      
      int found = filepath.find_last_of('.');
      std::string extension = filepath.substr(found + 1);

      printf("------------Reading Volume Model------------\n");
      if (extension.compare("pvm") == 0) {
        ret = readpvm(filepath);
      } else if (extension.compare("raw") == 0) {
        ret = ReadRawFile(filepath);
      } else if (extension.compare("syn") == 0) {
        ret = ReadSynFile(filepath);
      } else {
        printf("Extension not found\n");
      }
      printf("--------------------------------------------\n");

      return ret;
    }

  protected:
    Volume* readpvm (std::string filename)
    {
      printf("Read pvm:\n");
      
      double* scalar_values = NULL;
      unsigned int width, height, depth, components;
      float scalex, scaley, scalez;
      vr::DataStorageSize data_tp;
      
      //////////////////////////////////////////////
      // PVM Class
      Pvm fpvm(filename.c_str());

      fpvm.GetDimensions(&width, &height, &depth);
      components = fpvm.GetComponents();
      scalar_values = new double[width*height*depth];
      fpvm.GetScale(&scalex, &scaley, &scalez);

      // Get Volume Data
      float* vol_data = 
        //fpvm.GetData()
        fpvm.GenerateReescaledMinMaxData()
      ;

      for (int i = 0; i < width*height*depth; i++)
        scalar_values[i] = (double)vol_data[i];
       
      assert(components > 0);

      //GLubyte - 8 bits
      if (components == 1)
        data_tp = vr::DataStorageSize::_8_BITS;
      //GLushort - 16 bits
      else if (components == 2)
        data_tp = vr::DataStorageSize::_16_BITS;
      
      delete[] vol_data;
      //////////////////////////////////////////////

      Volume* ret = new Volume(width, height, depth, scalar_values, data_tp);
      ret->SetName(filename);
      ret->SetScale(scalex, scaley, scalez);

      delete[] scalar_values;

      return ret;
    }

  private:

  };
}
#endif