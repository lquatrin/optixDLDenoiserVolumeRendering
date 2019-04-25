#include "volume.h"

#include <iostream>
#include <string>
#include <cstdlib>
#include <fstream>

namespace vr
{
  Volume::Volume ()
    : m_width(0), m_height(0), m_depth(0)
  {
    m_scalar_values = NULL;
    data_storage_size = DataStorageSize::_8_BITS;

    scalex = 1.0f;
    scaley = 1.0f;
    scalez = 1.0f;
  }

  Volume::Volume (unsigned int width, unsigned int height, unsigned int depth)
    : m_width(width), m_height(height), m_depth(depth)
  {
    m_scalar_values = NULL;
    data_storage_size = DataStorageSize::_8_BITS;
    m_scalar_values = new double[m_width*m_height*m_depth];

    for (int i = 0; i < width * height * depth; i++)
      m_scalar_values[i] = 0.0;

    scalex = 1.0f;
    scaley = 1.0f;
    scalez = 1.0f;

    m_pmax = glm::dvec3(width, height, depth);
  }

  Volume::Volume (unsigned int width, unsigned int height, unsigned int depth, double* scalars, DataStorageSize storage_size)
  : m_width(width), m_height(height), m_depth(depth)
  {
    m_scalar_values = NULL;
    data_storage_size = storage_size;
    m_scalar_values = new double[m_width*m_height*m_depth];

    if (scalars != NULL)
      for (int i = 0; i < width * height * depth; i++)
        m_scalar_values[i] = (double)scalars[i];
    
    scalex = 1.0f;
    scaley = 1.0f;
    scalez = 1.0f;

    m_pmax = glm::dvec3(width, height, depth);
  }

  Volume::~Volume ()
  {
    delete[] m_scalar_values;
  }

  int Volume::GetWidth ()
  {
    return m_width;
  }

  int Volume::GetHeight ()
  {
    return m_height;
  }

  int Volume::GetDepth ()
  {
    return m_depth;
  }

  void Volume::SetScale (float sx, float sy, float sz)
  {
    scalex = sx;
    scaley = sy;
    scalez = sz;
  }

  float Volume::GetScaleX ()
  {
    return scalex;
  }
  
  float Volume::GetScaleY ()
  {
    return scaley;
  }
  
  float Volume::GetScaleZ ()
  {
    return scalez;
  }

  glm::vec3 Volume::GetScale ()
  {
    return glm::vec3(GetScaleX(),
                     GetScaleY(),
                     GetScaleZ());
  }
  
  glm::dvec3 Volume::GetAnchorMin ()
  {
    return m_pmin;
  }

  glm::dvec3 Volume::GetAnchorMax ()
  {
    return m_pmax;
  }

  void Volume::SetAnchors (glm::dvec3 pmin, glm::dvec3 pmax)
  {
    m_pmin = pmin;
    m_pmax = pmax;
  }

  int Volume::SampleVolume (int x, int y, int z)
  {
    if (x < 0 || y < 0 || z < 0)
      return 0.0;

    if (x > m_width  - 1 ||
        y > m_height - 1 || 
        z > m_depth  - 1)
      return 0.0;

    //x = glm::clamp<int>(x, 0, m_width  - 1);
    //y = glm::clamp<int>(y, 0, m_height - 1);
    //z = glm::clamp<int>(z, 0, m_depth  - 1);

    return (int)m_scalar_values[x + (y * m_width) + (z * m_width * m_height)];
  }

  int Volume::SampleVolume (double x, double y, double z)
  {
    if ((x >= m_pmin.x && x <= m_pmax.x) && (y >= m_pmin.y && y <= m_pmax.y) && (z >= m_pmin.z && z <= m_pmax.z))
      return -1;

    int px = (int)x;
    int py = (int)y;
    int pz = (int)z;

    return (int)m_scalar_values[px + (py * m_width) + (pz * m_width * m_height)];
  }

  int Volume::SampleVolume (int id)
  {
    return (int)m_scalar_values[id];
  }

  void Volume::SetSampleVolume (int x, int y, int z, double value)
  {
    m_scalar_values[x + (y * GetWidth()) + (z * GetWidth() * GetHeight())] = value;
  }

  double Volume::InterpolatedValue (double px, double py, double pz, glm::dvec3* pmin, glm::dvec3* pmax)
  {
    glm::dvec3 bbmin = pmin != NULL ? (*pmin) : m_pmin;
    glm::dvec3 bbmax = pmax != NULL ? (*pmax) : m_pmax;

    double x = ((px - bbmin.x) / (bbmax.x - bbmin.x)) * (double)(m_width - 1);
    double y = ((py - bbmin.y) / (bbmax.y - bbmin.y)) * (double)(m_height - 1);
    double z = ((pz - bbmin.z) / (bbmax.z - bbmin.z)) * (double)(m_depth - 1);

    int x0 = (int)x; int x1 = x0 + 1;
    int y0 = (int)y; int y1 = y0 + 1;
    int z0 = (int)z; int z1 = z0 + 1;

    if (x0 == (double)(m_width - 1))
    {
      x1 = (int)x0;
      x0 = x1 - 1;
    }

    if (y0 == (double)(m_height - 1))
    {
      y1 = (int)y0;
      y0 = y1 - 1;
    }

    if (z0 == (double)(m_depth - 1))
    {
      z1 = (int)z0;
      z0 = z1 - 1;
    }
    
    double xd = (x - (double)x0) / (double)(x1 - x0);
    double yd = (y - (double)y0) / (double)(y1 - y0);
    double zd = (z - (double)z0) / (double)(z1 - z0);
        
    // X interpolation
    double c00 = SampleVolume(x0, y0, z0)*(1.0 - xd) + SampleVolume(x1, y0, z0)*xd;
    double c10 = SampleVolume(x0, y1, z0)*(1.0 - xd) + SampleVolume(x1, y1, z0)*xd;
    double c01 = SampleVolume(x0, y0, z1)*(1.0 - xd) + SampleVolume(x1, y0, z1)*xd;
    double c11 = SampleVolume(x0, y1, z1)*(1.0 - xd) + SampleVolume(x1, y1, z1)*xd;

    // Y interpolation
    double c0 = c00*(1.0 - yd) + c10*yd;
    double c1 = c01*(1.0 - yd) + c11*yd;

    // Z interpolation
    double c = c0*(1.0 - zd) + c1*zd;

    return c;
  }
   
  double Volume::InterpolatedValue (glm::dvec3 pos)
  {
    return InterpolatedValue (pos.x, pos.y, pos.z);
  }

  /**
   * |----------
   * |
   * |
   * |
   * |----------
  **/
  double Volume::InterpolatedValueTextureBased (glm::dvec3 pos)
  {
    glm::dvec3 gpos = pos;
    
    // transform between [0,0,0] and [w,h,d]
    gpos = gpos / glm::dvec3(GetScale());

    double x = gpos.x - 0.5;
    int x0 = glm::clamp((int)glm::floor(x), 0, GetWidth() - 1);
    int x1 = glm::clamp((int)glm::ceil(x) , 0, GetWidth() - 1);

    double y = gpos.y - 0.5;
    int y0 = glm::clamp((int)glm::floor(y), 0, GetHeight() - 1);
    int y1 = glm::clamp((int)glm::ceil(y) , 0, GetHeight() - 1);

    double z = gpos.z - 0.5;
    int z0 = glm::clamp((int)glm::floor(z), 0, GetDepth() - 1);
    int z1 = glm::clamp((int)glm::ceil(z) , 0, GetDepth() - 1);
    

    double xd = 0.0;
    if (x1 - x0 != 0)
      xd = (x - (double)x0) / (double)(x1 - x0);

    double yd = 0.0;
    if (y1 - y0 != 0)
      yd = (y - (double)y0) / (double)(y1 - y0);

    double zd = 0.0;
    if (z1 - z0 != 0)
      zd = (z - (double)z0) / (double)(z1 - z0);

    // X interpolation
    double c00 = SampleVolume(x0, y0, z0)*(1.0 - xd) + SampleVolume(x1, y0, z0)*xd;
    double c10 = SampleVolume(x0, y1, z0)*(1.0 - xd) + SampleVolume(x1, y1, z0)*xd;
    double c01 = SampleVolume(x0, y0, z1)*(1.0 - xd) + SampleVolume(x1, y0, z1)*xd;
    double c11 = SampleVolume(x0, y1, z1)*(1.0 - xd) + SampleVolume(x1, y1, z1)*xd;

    // Y interpolation
    double c0 = c00 * (1.0 - yd) + c10 * yd;
    double c1 = c01 * (1.0 - yd) + c11 * yd;

    // Z interpolation
    double c = c0 * (1.0 - zd) + c1 * zd;

    return c;
  }

  float Volume::TrilinearScalarFunction (glm::dvec3 pos, glm::dvec3 rayeye, glm::dvec3 raydirection)
  {
    double fx = pos.x, fy = pos.y, fz = pos.z;

    fx = (fx - m_pmin.x) / abs (m_pmax.x - m_pmin.x) * m_width;
    fy = (fy - m_pmin.y) / abs (m_pmax.y - m_pmin.y) * m_height;
    fz = (fz - m_pmin.z) / abs (m_pmax.z - m_pmin.z) * m_depth;

    int x = (int)fx;
    int y = (int)fy;
    int z = (int)fz;

    if (x == m_width) --x;
    if (y == m_height) --y;
    if (z == m_depth) --z;

    //////////////////////////////////////////////////////////////////////////////
    //f(x, y, z) = c0 + c1*x + c2*y + c3*z + c4*x*y + c5*y*z + c6*x*z + c7*x*y*z//
    //////////////////////////////////////////////////////////////////////////////

    //cell coeficients
    double c[8];
    //scalar values of the cells
    double s[8];
    int m = 8, n = 8;
    int px[8] = { x, x + 1, x, x + 1, x, x + 1, x, x + 1 };
    int py[8] = { y, y, y, y, y + 1, y + 1, y + 1, y + 1 };
    int pz[8] = { z, z, z + 1, z + 1, z, z, z + 1, z + 1 };
    for (int i = 0; i < n; i++)
      s[i] = (double)SampleVolume(px[i], py[i], pz[i]);

    //mallocs
    double **a = (double**)malloc(sizeof(double*)* m);
    for (int i = 0; i < m; i++)
    {
      a[i] = (double*)malloc(sizeof(double)* n);
      a[i][0] = 1;
      a[i][1] = (double)px[i];
      a[i][2] = (double)py[i];
      a[i][3] = (double)pz[i];
      a[i][4] = (double)px[i] * py[i];
      a[i][5] = (double)py[i] * pz[i];
      a[i][6] = (double)px[i] * pz[i];
      a[i][7] = (double)px[i] * py[i] * pz[i];
    }
    double *svdw = (double*)malloc(sizeof(double)* n);
    double **v = (double**)malloc(sizeof(double*)* n);
    for (int i = 0; i < n; i++)
      v[i] = (double*)malloc(sizeof(double)* n);

    //Make SVD on Matrix a
    //lqc::dcksvd (a, m, n, svdw, v);

    //////////////////////////////////////////////////////////////
    //f(t) = w3*(t*t*t) + w2*(t*t) + w1*t + w0, t[tback, tfront]//
    //////////////////////////////////////////////////////////////

    //TODO
    double t = 0.0;

    glm::dvec3 e = rayeye;
    glm::dvec3 d = raydirection;

    float w[4];
    w[0] = c[0] + c[1] * e.x + c[2] * e.y + c[4] * e.x*e.y + c[3] * e.z
      + c[6] * e.x*e.z + c[5] * e.y*e.z + c[7] * e.x*e.y*e.z
      + c[7] * d.x*d.y*e.z;
    w[1] = c[1] * d.x + c[2] * d.y + c[3] * d.z + c[4] * d.y*e.x + c[6] * d.z*e.x
      + c[4] * d.x*e.y + c[5] * d.z*e.y + c[7] * d.z*e.x*e.y + c[6] * d.x*e.z + c[5] * d.y*e.z
      + c[7] * d.y*e.x*e.z + c[7] * d.x*e.y*e.z;
    w[2] = c[4] * d.x*d.y + c[6] * d.x*d.z + c[5] * d.y*d.z + c[7] * d.y*d.z*e.x + c[7] * d.x*d.z*e.y;
    w[3] = c[7] * d.x*d.y*d.z;

    double f_t = w[3] * (t*t*t) + w[2] * (t*t) + w[1] * t + w[0];

    //frees
    for (int i = 0; i < m; i++)
      free (a[i]);
    free (a);

    for (int i = 0; i < n; i++)
      free (v[i]);
    free (v);

    free (svdw);

    // return f_t;
    return 0.0f;
  }

  bool Volume::IsOutOfBoundary (int x, int y, int z)
  {
    return !((x >= 0 && x < m_width)
          && (y >= 0 && y < m_height)
          && (z >= 0 && z < m_depth));
  }

  unsigned long long Volume::CheckSum ()
  {
    unsigned long long csum = 0;
    for (int i = 0; i < m_width * m_height * m_depth; i++)
      csum += (unsigned long long)m_scalar_values[i];
    return csum;
  }

  // https://en.wikipedia.org/wiki/Sobel_operator  
  glm::vec3* GenerateSobelFeldmanGradientTextureT(Volume* vol)
  {
    int width = vol->GetWidth();
    int height = vol->GetHeight();
    int depth = vol->GetDepth();

    int n = 1;
    glm::dvec3* gradients = new glm::dvec3[width * height * depth];
    for (int z = 0; z < depth; z++)
    {
      for (int y = 0; y < height; y++)
      {
        for (int x = 0; x < width; x++)
        {
          glm::dvec3 sg(0.0);
          for (int v1 = -1; v1 <= 1; v1++)
          {
            for (int v2 = -1; v2 <= 1; v2++)
            {
              sg.z += ((double)vol->SampleVolume(x + v1, y + v2, z - 1) / vol->GetMaxDensity()) * (4.0 / pow(2.0, glm::abs(v1) + glm::abs(v2)))
                + ((double)vol->SampleVolume(x + v1, y + v2, z + 1) / vol->GetMaxDensity()) * (-4.0 / pow(2.0, glm::abs(v1) + glm::abs(v2)));

              sg.y += ((double)vol->SampleVolume(x + v1, y - 1, z + v2) / vol->GetMaxDensity()) * (4.0 / pow(2.0, glm::abs(v1) + glm::abs(v2)))
                + ((double)vol->SampleVolume(x + v1, y + 1, z + v2) / vol->GetMaxDensity()) * (-4.0 / pow(2.0, glm::abs(v1) + glm::abs(v2)));

              sg.x += ((double)vol->SampleVolume(x - 1, y + v2, z + v1) / vol->GetMaxDensity()) * (4.0 / pow(2.0, glm::abs(v1) + glm::abs(v2)))
                + ((double)vol->SampleVolume(x + 1, y + v2, z + v1) / vol->GetMaxDensity()) * (-4.0 / pow(2.0, glm::abs(v1) + glm::abs(v2)));
            }
          }

          // not normalized (for tests...)
          gradients[x + (y * width) + (z * width * height)] = sg;
        }
      }
    }

    glm::vec3* gradients_values = new glm::vec3[width * height * depth];
    for (int k = 0; k < depth; k++)
    {
      for (int j = 0; j < height; j++)
      {
        for (int i = 0; i < width; i++)
        {
          gradients_values[i + (j * width) + (k * width * height)] = gradients[i + (j * width) + (k * width * height)];
        }
      }
    }

    delete[] gradients;
    return gradients_values;
  }

}