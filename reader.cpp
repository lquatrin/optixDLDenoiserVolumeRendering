#include "reader.h"

#include <cstdio>
#include <fstream>
#include <string>

#include "rawloader.h"

#include "transferfunction.h"
#include "transferfunction1d.h"

namespace vr
{
  double* ReadVolvisRaw(std::string volfilename, size_t bytes_per_value, int w, int h, int d)
  {
    double* scalar_values = NULL;
    IRAWLoader rawLoader = IRAWLoader(volfilename, bytes_per_value, w*h*d, bytes_per_value);

    //GLushort - 16 bits
    if (bytes_per_value == sizeof(unsigned short))
    {
      scalar_values = new double[w*h*d];
      unsigned short *b = new unsigned short[w*h*d];
      memcpy(b, rawLoader.GetData(), sizeof(unsigned short)*w*h*d);

      for (int i = 0; i < w*h*d; i++)
      {
        scalar_values[i] = (double)b[i];
      }
      // TODO: NORMALIZE

      delete[] b;
    }
    //GLubyte - 8 bits
    else if (bytes_per_value == sizeof(unsigned char))
    {
      scalar_values = new double[w*h*d];
      unsigned char *b = new unsigned char[w*h*d];
      memcpy(b, rawLoader.GetData(), sizeof(unsigned char)*w*h*d);

      for (int i = 0; i < w*h*d; i++)
        scalar_values[i] = (double)b[i];
      // TODO: NORMALIZE

      delete[] b;
    }

    return scalar_values;
  }

  Volume* ReadRawFile(std::string filepath)
  {
    Volume* ret = NULL;

    printf("Started  -> Read Volume From .raw File\n");
    printf("  - File .raw Path: %s\n", filepath.c_str());

    std::ifstream file(filepath.c_str());
    if (file.is_open())
    {
      int foundinit = filepath.find_last_of('\\');
      std::string filename = filepath.substr(foundinit + 1);
      printf("  - File .raw: %s\n", filename.c_str());

      int foundfp = filename.find_last_of('.');
      filename = filename.substr(0, foundfp);

      int foundsizes = filename.find_last_of('.');
      std::string t_filesizes = filename.substr(foundsizes + 1,
        filename.size() - foundsizes);

      filename = filename.substr(0, filename.find_last_of('.'));

      int foundbytesize = filename.find_last_of('.');
      std::string t_filebytesize = filename.substr(foundbytesize + 1,
        filename.size() - foundbytesize);

      int fw, fh, fd;
      int bytesize;

      // Read the Volume Sizes
      int foundd = t_filesizes.find_last_of('x');
      fd = atoi(t_filesizes.substr(foundd + 1,
        t_filesizes.size() - foundd).c_str());

      t_filesizes = t_filesizes.substr(0, t_filesizes.find_last_of('x'));

      int foundh = t_filesizes.find_last_of('x');
      fh = atoi(t_filesizes.substr(foundh + 1,
        t_filesizes.size() - foundh).c_str());

      t_filesizes = t_filesizes.substr(0, t_filesizes.find_last_of('x'));

      int foundw = t_filesizes.find_last_of('x');
      fw = atoi(t_filesizes.substr(foundw + 1,
        t_filesizes.size() - foundw).c_str());
      // Byte Size
      bytesize = atoi(t_filebytesize.c_str());

      double* scalar_values = ReadVolvisRaw(filepath, (size_t)bytesize, fw, fh, fd);
      ret = new Volume(fw, fh, fd, scalar_values, GetStorageSizeType((size_t)bytesize));
      ret->SetName(filepath);

      printf("  - Volume Name     : %s\n", filepath.c_str());
      printf("  - Volume Size     : [%d, %d, %d]\n", fw, fh, fd);
      printf("  - Volume Byte Size: %d\n", bytesize);

      file.close();
      printf("Finished -> Read Volume From .raw File\n");
    }
    else
      printf("Finished -> Error on opening .raw file\n");

    return ret;
  }
  
  Volume* ReadSynFile (std::string filename)
  {
    Volume* ret = NULL;

    printf("Started  -> Read Volume From .syn File\n");
    printf("  - File .syn Path: %s\n", filename.c_str());

    std::string fileifstream = "";
    fileifstream.append(filename);

    std::ifstream file(fileifstream.c_str());
    if (file.is_open())
    {
      int width, height, depth;
      file >> width >> height >> depth;
      std::cout << "  - Volume Size: [" << width << ", " << height << ", " << depth << "]" << std::endl;

      ret = new Volume(width, height, depth);

      int new_data = 0;
      while (file >> new_data)
      {
        if (new_data == 1)
        {
          int x0, y0, z0, x1, y1, z1, v;
          file >> x0 >> y0 >> z0 >> x1 >> y1 >> z1 >> v;
          for (int x = x0; x < x1; x++)
          {
            for (int y = y0; y < y1; y++)
            {
              for (int z = z0; z < z1; z++)
              {
                ret->SetSampleVolume(x, y, z, (double)v);
              }
            }
          }
        }
        else
        {
          int xt, yt, zt, v;
          file >> xt >> yt >> zt >> v;
          ret->SetSampleVolume(xt, yt, zt, (double)v);
        }
      }

      printf("lqc: Finished -> Read Volume From .syn File\n");
    }
    else
      printf("lqc: Finished -> Error on opening .syn file\n");

    return ret;
  }

TransferFunction* ReadTransferFunction(std::string file)
{
  TransferFunction* ret = NULL;

  int found = file.find_last_of('.');
  std::string extension = file.substr(found + 1);

  printf("---------Reading Transfer Function----------\n");
  printf(" - File: %s\n", file.c_str());

  if (extension.compare("tf1d") == 0)
    ret = ReadTransferFunction_tf1d(file);
  //else if (extension.compare("tf1dnorm") == 0)
  //  ret = ReadTransferFunction_tf1dnorm(file);
  //else if (extension.compare ("tfg1d") == 0)
  //  ret = ReadTransferFunction_tfg1d (file);
  //else if (extension.compare ("tfgersa") == 0)
  //  ret = ReadTransferFunction_tfgersa (file);

  printf("--------------------------------------------\n");
  return ret;
}

TransferFunction* ReadTransferFunction_tf1d(std::string file)
{
  std::ifstream myfile(file);
  if (myfile.is_open())
  {
    std::string interpolation;
    std::getline(myfile, interpolation);
    printf(" - Interpolation: %s\n", interpolation.c_str());

    TransferFunction1D* tf = NULL;
    int init;
    myfile >> init;

    if (init == 2)
    {
      int max_density;
      myfile >> max_density;

      int extuse;
      myfile >> extuse;

      tf = new TransferFunction1D(max_density);
      tf->SetExtinctionCoefficientInput(extuse == 1);
    }
    else if (init == 1)
    {
      int max_density;
      myfile >> max_density;

      tf = new TransferFunction1D(max_density);
    }
    else
    {
      tf = new TransferFunction1D();
    }

    int cpt_rgb_size;
    myfile >> cpt_rgb_size;
    double r, g, b, a;
    int isovalue;
    for (int i = 0; i < cpt_rgb_size; i++)
    {
      myfile >> r >> g >> b >> isovalue;
      tf->AddRGBControlPoint(TransferControlPoint(r, g, b, isovalue));
    }

    int cpt_alpha_size;
    myfile >> cpt_alpha_size;
    for (int i = 0; i < cpt_alpha_size; i++)
    {
      myfile >> a >> isovalue;
      tf->AddAlphaControlPoint(TransferControlPoint(a, isovalue));
    }

    myfile.close();


    if (interpolation.compare("linear") == 0)
      tf->m_interpolation_type = TFInterpolationType::LINEAR;
    else if (interpolation.compare("cubic") == 0)
      tf->m_interpolation_type = TFInterpolationType::CUBIC;

    int foundname = file.find_last_of('\\');
    std::string tfname = file.substr(foundname + 1);
    tf->SetName(file);

    return tf;
  }
  return NULL;
}
}