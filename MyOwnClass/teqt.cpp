#include <iostream>
#include <cstdlib>
#include <math.h>

class Point2D
{
public:
  float X{ 0 };
  float Y{ 0 };

  Point2D() = default;

  Point2D(float x, float y)
    : X{x}, Y{y}
  {
  }
};

class Equation
{
private:
  double X;
  double Y;

  bool cutOffState = false;

protected:

public:
  double m_a;
  double m_b;
  double m_r;
  double m_x;
  double m_y;
  double m_slope;
  
  double D{ 0 };
  double a{ 0 };
  double b{ 0 };
  double c{ 0 };

  Equation();
  Equation(float x, float y);
  Equation(Point2D &point);

  ~Equation(){ };

  float getX();
  float getY();

  void setCircle(Point2D &point, float r);
  float slope(Point2D &point);

  bool isCutOff();
  float getCutDistance();

  void pritResult();
};

Equation::Equation()
  : X{0}, Y{0}, m_a{0}, m_b{0}, m_r{0}
{
}

Equation::Equation(float x, float y)
  : X{x}, Y{y}, m_a{0}, m_b{0}, m_r{0}
{
}

Equation::Equation(Point2D &point)
  : X{point.X}, Y{point.Y}, m_a{0}, m_b{0}, m_r{0}
{
}

float Equation::getX()
{
  return X;
}

float Equation::getY()
{
  return Y;
}

void Equation::setCircle(Point2D &point, float r)
{
  m_a = point.X;
  m_b = point.Y;
  m_r = r;
}

float Equation::slope(Point2D &point)
{
  m_x = point.X;
  m_y = point.Y;

  m_slope = (Y - point.Y) / (X - point.X);
  
  return m_slope;
}

bool Equation::isCutOff()
{
  if (m_r == 0.0)
  {
    return false;
  }

  a = 1 + pow(m_slope, 2);
  b = 2*m_slope*m_y - 2*m_slope*m_b - 2*m_a - 2*(pow(m_slope, 2))*m_x;
  c = pow(m_a, 2) + pow(m_b, 2) - pow(m_r, 2) + pow(m_y, 2) + (pow(m_slope, 2))*pow(m_x, 2) - 2*m_slope*m_y*m_x - 2*m_b*m_y;
  
  D = pow(b, 2) - 4*a*c;

  if (D > 0)
  {
    cutOffState = true;

    return true;
  }

  return false;
}

float Equation::getCutDistance()
{
  float dist = -1.0;
  if (cutOffState == true)
  {
    dist = ( fabs(-m_slope*m_a + m_b + (m_slope*m_x - m_y)) / sqrt((pow(m_slope, 2) + pow(1, 2))) );
  
    return dist;
  }

  return dist;
}

void Equation::pritResult()
{
  std::cout << "a, b, and c value " << a << " " << b << " " << c << "\n";
  std::cout << "determinant is " << D << "\n";
}

int main()
{
  float r = 2.23606798;

  Point2D initial_p(0, 5);
  Point2D sec_p(5, -5);
  Point2D center(0, 0);

  Equation eqt(initial_p);

  eqt.setCircle(center, r);

  float m = eqt.slope(sec_p);
  float state = eqt.isCutOff();

  std::cout << "slope is " << m << "\n";
  std::cout << "cut off state " << state << "\n";
  std::cout << "cut off distance " << eqt.getCutDistance() << "\n";
  eqt.pritResult();
}