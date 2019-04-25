/**
 * Transfer function based on tutorial from Graphics Runner:
 * . http://graphicsrunner.blogspot.com/2009/01/volume-rendering-101.html 
 * . http://graphicsrunner.blogspot.com/2009/01/volume-rendering-102-transfer-functions.html
 *
 * Format:
 * |interpolation_type   // can be "linear" or "cubic"
 * |type                 // 0, 1 or 2
 * Next Line based on the type in the previous line
 * |max_density[for type = 1, 2] extinction_input[for type = 2]
 * |number_of_rgb_values
 * |r g b isovalue
 * |r g b isovalue
 * |...
 * |number_of_a_or_t_values
 * |a/t isovalue
 * |a/t isovalue
 * |...
**/

#ifndef VOLREND_TRANSFERFUNCTION1D_H
#define VOLREND_TRANSFERFUNCTION1D_H

#include "transferfunction.h"

#include <vector>
#include <iostream>

#include "glm/glm.hpp"

namespace vr
{
  class TransferFunction1D : public TransferFunction
  {
  public:
    TransferFunction1D (int max_value = 255);
    ~TransferFunction1D ();

    virtual const char* GetNameClass ();
    virtual glm::vec4 Get (double value, double max_data_value = -1.0);
    
    virtual float GetOpc (double value, double max_input_value = -1.0);
    virtual float GetExt (double value, double max_input_value = -1.0);

    void SetExtinctionCoefficientInput (bool s);

    void AddRGBControlPoint (TransferControlPoint rgb);
    void AddAlphaControlPoint (TransferControlPoint alpha);
    void ClearControlPoints ();

    glm::vec3 GetRGBPointOnSpline (float s);
    glm::vec3 GetAlphaPointOnSpline (float s);

    //If we don't have a file with the values of the TF, we need to compute the TF
    void Build (TFInterpolationType type);

    void PrintControlPoints ();
    void PrintTransferFunction ();

    bool Save ();
    bool Load ();

    TFInterpolationType m_interpolation_type;
    bool m_built;
  private:
    void BuildLinear ();

    std::vector<TransferControlPoint> m_cpt_rgb;
    std::vector<TransferControlPoint> m_cpt_alpha;
    glm::dvec4* m_transferfunction;
    glm::vec3* m_gradients;

    int max_density;
    bool extinction_coef_type;
  };

}

#endif