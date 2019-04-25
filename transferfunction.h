#ifndef VOLREND_TRANSFERFUNCTION_H
#define VOLREND_TRANSFERFUNCTION_H

#include <iostream>
#include <vector>

#include "glm/glm.hpp"

namespace vr
{
  enum TFInterpolationType
  {
    LINEAR,
    CUBIC
  };

  enum TFFormatType
  {
    LQC
  };
  
  class TransferControlPoint
  {
  public:
    TransferControlPoint (double r, double g, double b, int isovalue);
    TransferControlPoint (double alpha, int isovalue);

    glm::vec3 operator -(const TransferControlPoint& v)  
    {  
      glm::vec4 ret;
      ret = this->m_color - v.m_color;

      return glm::vec3 (ret.x, ret.y, ret.z);
    }

    glm::vec3 operator +(const TransferControlPoint& v)
    {  
      glm::vec4 ret;
      ret = this->m_color + v.m_color;

      return glm::vec3 (ret.x, ret.y, ret.z);
    }

    glm::vec3 operator +(const glm::vec3& v)
    {  
      return glm::vec3 (this->m_color.x + v.x, this->m_color.y + v.y, this->m_color.z + v.z);
    }

    glm::vec4 m_color;
    int m_isoValue;
  };
    
  class Cubic
  {
  public:
    Cubic ();
    Cubic (glm::vec4 in1, glm::vec4 in2, glm::vec4 in3, glm::vec4 in4);
    ~Cubic ();

    glm::vec4 GetPointOnSpline (double s);
    static Cubic* CalculateCubicSpline (int n, std::vector<TransferControlPoint> v);

    glm::vec4 a, b, c, d;
  };
  Cubic* CalculateCubicSpline (int n, std::vector<TransferControlPoint> v);

  class TransferFunction
  {
  public:
    TransferFunction() {}
    ~TransferFunction() {}

    virtual const char* GetNameClass() = 0;

    virtual glm::vec4 Get (double value, double max_input_value = -1.0) = 0;

    virtual float GetOpc (double value, double max_input_value = -1.0) { return -1.0; }
    virtual float GetExt (double value, double max_input_value = -1.0) { return -1.0; }

    std::string GetName () { return m_name; }
    void SetName (std::string name) { m_name = name; }
    
    //////////////////////////////////////////////////////////////////
    // Interface from:
    // . A Coherent Projection Approach for Direct Volume Rendering
    // . Jane Wilhelms and Allen Van Gelder
    // . 1991
    // Opacity goes from 0 to 1 when extinction goes from 0 to ~15 using floats
    double ExtinctionToMaterialOpacity (float extinction)
    {
      return (1.0 - glm::exp(-extinction));
    }
    
    double MaterialOpacityToExtinction (float opacity)
    {
      return (glm::log(1.0 / (1.0 - opacity)));
    }
    //////////////////////////////////////////////////////////////////


  protected:
    std::string m_name;
  private:

  };
}

#endif