#include <iostream>
#include <cstdlib>
#include <math.h>
#include <vector>

namespace alg
{
  bool valueInside(float value, float min, float max)
  {
    return (value > min && value < max);
  }
};

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
  double X{ 0 };
  double Y{ 0 };

  float m_x1{ 0 };
  float m_y1{ 0 };
  float m_x2{ 0 };
  float m_y2{ 0 };

  double m_a{ 0 };
  double m_b{ 0 };
  double m_r{ 0 };
  double m_x{ 0 };
  double m_y{ 0 };

  bool cutOffState = false;
  bool determinanCircletState = false;
  bool distanceCutState = false;

protected:

public:
  struct equationComp
  {
    double D{ 0 };
    double a{ 0 };
    double b{ 0 };
    double c{ 0 };
    double dist{ -1.0 };
    double GCD{ 1.0 };
    double slope{ 0 };

    double addt1{ 0 };
    double addt2{ 0 };
  };

  equationComp eqt1;
  equationComp eqt2;

  Equation() = default;
  Equation(float x, float y);
  Equation(Point2D &point);

  ~Equation(){ };

  float getX();
  float getY();

  void setCircle(Point2D &point, float r);
  void setLinier(Point2D &point);
  float getSlope(float x, float y);

  bool isCutOff();
  float getCutDistance();
  float getCircleDeterminant();

  void printResult();
  void printEquation();

  float getTangentPoint();
  float getDeterminant(float a, float b, float c);

  float findGCD(std::vector<double> &var);
  float countGCD(float a, float b);
};

Equation::Equation(float x, float y)
  : X{x}, Y{y}
{
}

Equation::Equation(Point2D &point)
  : X{point.X}, Y{point.Y}
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

void Equation::setLinier(Point2D &point)
{
  m_x = point.X;
  m_y = point.Y;
}

float Equation::getSlope(float x, float y)
{  
  return (Y - y) / (X - x);
}

bool Equation::isCutOff()
{
  if (determinanCircletState == false)
    getCircleDeterminant();
  if (distanceCutState == false)
    getCutDistance();

  if (eqt1.D > 0 && !alg::valueInside(eqt1.dist, m_r - 0.0001, m_r + 0.0001))
  {
    cutOffState = true;

    return true;
  }

  return false;
}

float Equation::getCircleDeterminant()
{
  if (m_r == 0.0)
  {
    return eqt1.D;
  }

  eqt1.slope = getSlope(m_x, m_y);

  eqt1.a = 1 + pow(eqt1.slope, 2);
  eqt1.b = 2*eqt1.slope*m_y - 2*eqt1.slope*m_b - 2*m_a - 2*(pow(eqt1.slope, 2))*m_x;
  eqt1.c = pow(m_a, 2) + pow(m_b, 2) - pow(m_r, 2) + pow(m_y, 2) + (pow(eqt1.slope, 2))*pow(m_x, 2) 
    - 2*eqt1.slope*m_y*m_x + 2*eqt1.slope*m_b*m_x - 2*m_b*m_y;
  
  eqt1.D = getDeterminant(eqt1.a, eqt1.b, eqt1.c);

  determinanCircletState = true;
  return eqt1.D;
}

float Equation::getDeterminant(float a, float b, float c)
{
  return pow(b, 2) - 4*a*c;
}

float Equation::getCutDistance()
{
  if (eqt1.D > 0)
  {
    eqt1.dist = ( fabs(-eqt1.slope*m_a + m_b + (eqt1.slope*m_x - m_y)) / sqrt((pow(eqt1.slope, 2) + pow(1, 2))) );
  
    return eqt1.dist;
  }
  else if (determinanCircletState == false)
  {
    getCircleDeterminant();

    getCutDistance();
  }

  distanceCutState == true;
  return eqt1.dist;
}

float Equation::getTangentPoint()
{
  double n = m_b - Y;

  eqt2.a = pow(n, 2) + pow((X - m_a), 2);
  eqt2.b = 2*m_a*( pow(X, 2) - 2*m_a*X + pow(m_a, 2) - pow(n, 2) ) + 2*m_b*( X*(m_b - Y) 
    - m_a*(m_b + Y) + n*(m_a - X) ) - 2*pow(m_r, 2)*( X - m_a );
  eqt2.c = pow(m_a , 2)*( pow(n, 2) + pow((X - m_a), 2) ) + pow(m_b, 2)*( pow((m_b - Y), 2) ) 
    + pow(m_r, 4) + 2*m_a*m_b*( X*(m_b - Y) + m_a*(Y - m_b) ) + 2*pow(m_r, 2)*( m_a*(m_a - X) 
    + m_b*(Y - m_b) ) + 2*n*m_a*m_b*(m_a - X) + 2*n*pow(m_b, 2)*( Y - m_b + n ) + 2*n*m_b*pow(m_r, 2) 
    - pow(m_r*n, 2);

  std::vector<double> eqt2GCD;
  eqt2GCD.resize(3);

  eqt2.a /= pow(n, 2);
  eqt2.b /= pow(n, 2);
  eqt2.c /= pow(n, 2);

  eqt2GCD[0] = fabs(eqt2.a);
  eqt2GCD[1] = fabs(eqt2.b);
  eqt2GCD[2] = fabs(eqt2.c);

  eqt2.GCD = findGCD(eqt2GCD);

  std::cout << "flag " << eqt2.GCD << "\n";

  eqt2.a /= eqt2.GCD;
  eqt2.b /= eqt2.GCD;
  eqt2.c /= eqt2.GCD;

  eqt2.D = getDeterminant(eqt2.a, eqt2.b, eqt2.c);

  m_x1 = ( -(eqt2.b) + sqrt(eqt2.D) ) / 2*eqt2.a;
  m_x2 = ( -(eqt2.b) - sqrt(eqt2.D) ) / 2*eqt2.a;

  m_y1 = ( m_x1*(X - m_a) + m_a*(m_a - X) + m_b*(m_b - Y) - pow(m_r, 2) ) / n;
  m_y2 = ( m_x2*(X - m_a) + m_a*(m_a - X) + m_b*(m_b - Y) - pow(m_r, 2) ) / n;

  // eqt2.addt1 = slope(m_x1, m_y1);
  // eqt2.addt2 = slope(m_x2, m_y2);
}

float Equation::findGCD(std::vector<double> &var)
{
  std::size_t length = var.size();
  float result = var[0];

  for (std::size_t i = 1; i < length; i++)
    result = countGCD(var[i], result);
  
  return result;
}

float Equation::countGCD(float a, float b)
{
  if (a < b) 
    return countGCD(b, a);

  if (fabs(b) < 0.001) 
    return a; 

  else
    return (countGCD(b, a - std::floor(a / b) * b)); 
}

void Equation::printResult()
{
  std::cout << "a, b, and c value " << eqt1.a << " " << eqt1.b << " " << eqt1.c << "\n";
  std::cout << "determinant is " << eqt1.D << "\n";
}

void Equation::printEquation()
{
  std::cout << "circle equation is ";
  printf("(x - (%.2f))^2 + (y - (%.2f))^2 = %.2f\n", m_a, m_b, pow(m_r, 2));
  
  std:: cout << "linier equation (Robot) is ";
  printf("y = %.2fx + (%.2f)\n", eqt1.slope, -(eqt1.slope*m_x) + m_y);

  std:: cout << "tangent equation is ";
  printf("(%.2f)x^2 + (%.2f)x + (%.2f) = 0\n", eqt2.a, eqt2.b, eqt2.c);
  
  std:: cout << "first tangent point is ";
  printf("x = %.2f, y = %.2f\n", m_x1, m_y1);

  std:: cout << "second tangent point is ";
  printf("x = %.2f, y = %.2f\n", m_x2, m_y2);
}

int main()
{
  float r = 20;

  Point2D initial_p(150, 300);
  Point2D sec_p(222.2, 612.45);
  Point2D center(200, 600);

  Equation eqt(initial_p);

  eqt.setCircle(center, r);
  eqt.setLinier(sec_p);

  bool state = eqt.isCutOff();
  eqt.getTangentPoint();

  std::cout << "slope is " << "\n";
  std::cout << "cut off state " << state << "\n";
  std::cout << "cut off distance " << eqt.getCutDistance() << "\n";
  eqt.printResult();
  eqt.printEquation();
}