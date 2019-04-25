#include "transferfunction.h"

#include <fstream>
#include <cstdlib>

namespace vr
{
  TransferControlPoint::TransferControlPoint (double r, double g, double b, int isovalue)
  {
    m_color.x = r;
    m_color.y = g;
    m_color.z = b;
    m_color.w = 1.0f;
    m_isoValue = isovalue;
  }

  TransferControlPoint::TransferControlPoint (double alpha, int isovalue)
  {
    m_color.x = 0.0f;
    m_color.y = 0.0f;
    m_color.z = 0.0f;
    m_color.w = alpha;
    m_isoValue = isovalue;
  }

  Cubic::Cubic ()
  {
  }

  Cubic::Cubic (glm::vec4 ain, glm::vec4 bin, glm::vec4 cin, glm::vec4 din)
    : a(ain), b(bin), c(cin), d(din)
  {
  }

  Cubic::~Cubic ()
  {
  }

  glm::vec4 Cubic::GetPointOnSpline (double s)
  {
    return (((d * float(s)) + c) * float(s) + b) * float(s) + a;
  }

  Cubic* Cubic::CalculateCubicSpline (int n, std::vector<TransferControlPoint> v)
  {
    glm::vec4* gamma = new glm::vec4[n + 1];
    glm::vec4* delta = new glm::vec4[n + 1];
    glm::vec4* D = new glm::vec4[n + 1];
    int i;
    /* We need to solve the equation
    * taken from: http://mathworld.wolfram.com/CubicSpline.html
    [2 1       ] [D[0]]   [3(v[1] - v[0])  ]
    |1 4 1     | |D[1]|   |3(v[2] - v[0])  |
    |  1 4 1   | | .  | = |      .         |
    |    ..... | | .  |   |      .         |
    |     1 4 1| | .  |   |3(v[n] - v[n-2])|
    [       1 2] [D[n]]   [3(v[n] - v[n-1])]

    by converting the matrix to upper triangular.
    The D[i] are the derivatives at the control points.
    */

    //this builds the coefficients of the left matrix
    gamma[0] = glm::vec4(0);
    gamma[0].x = 1.0 / 2.0;
    gamma[0].y = 1.0 / 2.0;
    gamma[0].z = 1.0 / 2.0;
    gamma[0].w = 1.0 / 2.0;
    for (i = 1; i < n; i++)
    {
      gamma[i] = glm::vec4(1) / ((float(4.0) * glm::vec4(1)) - gamma[i - 1]);
    }
    gamma[n] = glm::vec4(1) / ((float(2.0) * glm::vec4(1)) - gamma[n - 1]);

    delta[0] = float(3.0) * (v[1].m_color - v[0].m_color) * gamma[0];
    for (i = 1; i < n; i++)
    {
      delta[i] = (float(3.0) * (v[i + 1].m_color - v[i - 1].m_color) - delta[i - 1]) * gamma[i];
    }
    delta[n] = (float(3.0) * (v[n].m_color - v[n - 1].m_color) - delta[n - 1]) * gamma[n];

    D[n] = delta[n];
    for (i = n - 1; i >= 0; i--)
    {
      D[i] = delta[i] - gamma[i] * D[i + 1];
    }

    // now compute the coefficients of the cubics 
    Cubic* C = new Cubic[n];
    for (i = 0; i < n; i++)
    {
      C[i] = Cubic(v[i].m_color, D[i], 3.0f * (v[i + 1].m_color- v[i].m_color) - 2.0f * D[i] - D[i + 1],
        2.0f * (v[i].m_color - v[i + 1].m_color) + D[i] + D[i + 1]);
    }
    return C;
  }

  Cubic* CalculateCubicSpline (int n, std::vector<TransferControlPoint> v)
  {
    glm::vec4* gamma = new glm::vec4[n + 1];
    glm::vec4* delta = new glm::vec4[n + 1];
    glm::vec4* D = new glm::vec4[n + 1];
    int i;
   // We need to solve the equation
   // taken from: http://mathworld.wolfram.com/CubicSpline.html
   // [2 1       ] [D[0]]   [3(v[1] - v[0])  ]
   // |1 4 1     | |D[1]|   |3(v[2] - v[0])  |
   // |  1 4 1   | | .  | = |      .         |
   // |    ..... | | .  |   |      .         |
   // |     1 4 1| | .  |   |3(v[n] - v[n-2])|
   // [       1 2] [D[n]]   [3(v[n] - v[n-1])]

    //by converting the matrix to upper triangular.
    //The D[i] are the derivatives at the control points.
    

    //this builds the coefficients of the left matrix
    gamma[0] = glm::vec4(0);
    gamma[0].x = 1.0 / 2.0;
    gamma[0].y = 1.0 / 2.0;
    gamma[0].z = 1.0 / 2.0;
    gamma[0].w = 1.0 / 2.0;
    for (i = 1; i < n; i++)
    {
      gamma[i] = glm::vec4(1) / ((4.0f * glm::vec4(1)) - gamma[i - 1]);
    }
    gamma[n] = glm::vec4(1) / ((2.0f * glm::vec4(1)) - gamma[n - 1]);

    delta[0] = 3.0f * (v[1].m_color - v[0].m_color) * gamma[0];
    for (i = 1; i < n; i++)
    {
      delta[i] = (3.0f * (v[i + 1].m_color - v[i - 1].m_color) - delta[i - 1]) * gamma[i];
    }
    delta[n] = (3.0f * (v[n].m_color - v[n - 1].m_color) - delta[n - 1]) * gamma[n];

    D[n] = delta[n];
    for (i = n - 1; i >= 0; i--)
    {
      D[i] = delta[i] - gamma[i] * D[i + 1];
    }

    // now compute the coefficients of the cubics 
    Cubic* C = new Cubic[n];
    for (i = 0; i < n; i++)
    {
      C[i] = Cubic(v[i].m_color, D[i], 3.0f * (v[i + 1].m_color- v[i].m_color) - 2.0f * D[i] - D[i + 1],
        2.0f * (v[i].m_color - v[i + 1].m_color) + D[i] + D[i + 1]);
    }
    return C;
  }
}