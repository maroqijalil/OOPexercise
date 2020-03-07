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

  float piValue()
  {
    return 3.14159265359;
  }

  float normalizeRad(float val)
  {
    while (val < -piValue())
      val += 2 * piValue();

    while (val > piValue())
      val -= 2 * piValue();

    return val;
  }

  float direction(float x, float y)
  {
    return normalizeRad(atan2(y, x));
  }

  float rad2Deg()
  {
    return 180.0 / piValue();
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
  std::vector<Point2D> nPoint;
  std::vector<float> nAngel;

  double X{ 0.0 };
  double Y{ 0.0 };

  double m_a{ -1.0 };
  double m_b{ -1.0 };
  double m_r{ 0.0 };
  double m_x{ 0.0 };
  double m_y{ 0.0 };
  double m_th{ 0.0 };

  bool cutOffState = false;
  bool determinanCircletState = false;
  bool distanceCutState = false;
  bool initialPoint = false;

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
  };

  struct couplePoint
  {
    Point2D _1_;
    Point2D _2_;
  };

  equationComp eqt1;
  equationComp eqt2;

  couplePoint intersect;
  couplePoint tangent;

  Equation() = default;
  Equation(float x, float y);
  Equation(const Point2D& point);

  ~Equation(){ };

  float getX();
  float getY();

  void setCircle(const Point2D& point, float r);
  void setLinier(const Point2D& point);
  float getTargetDistance();
  float getSlope(float x, float y);

  bool isCutOff();
  bool getCutState();
  float getCutDistance();
  float getCircleDeterminant();
  void getTransitPoint(Point2D& point);
  void getEdgePoint(Point2D& point);

  void printResult();
  void printEquation();

  float getTangentPoint();
  float getIntersectPoint();
  float getDeterminant(float a, float b, float c);
  float getAngel(float x1, float y1, float x2, float y2);

  float findGCD(std::vector<double> &var);
  float countGCD(float a, float b);
};

Equation::Equation(float x, float y)
  : X{x}, Y{y}
{
}

Equation::Equation(const Point2D& point)
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

float Equation::getAngel(float x1, float y1, float x2, float y2)
{
  float multiplier = 1.0;;

  if (y1 < y2)
    multiplier *= -1.0;

  if ( !(m_a == -1 || m_b == -1) )
  {
    float dx = x1 - x2;
    float dy = y1 - y2;
    float th = 0.0;

    float dist = sqrt( pow(dx, 2) + pow(dy, 2) );

    return (acos( fabs((2*pow(m_r, 2) - pow(dist, 2)) / (2*pow(m_r, 2))) ) * alg::rad2Deg()) * multiplier;
  }

  return 0.0;
}

void Equation::setCircle(const Point2D& point, float r)
{
  m_a = point.X;
  m_b = point.Y;
  m_r = r;
}

void Equation::setLinier(const Point2D& point)
{
  initialPoint = false;

  m_x = point.X;
  m_y = point.Y;

  m_th = getAngel(m_x, m_y, m_a + 20.0, m_b);

  if (m_th > 180.0)
    m_th -= 360.0;
}

float Equation::getTargetDistance()
{
  float dx = X - m_x;
  float dy = Y - m_y;

  return (sqrt(pow(dx, 2) + pow(dy, 2)));
}

float Equation::getSlope(float x, float y)
{  
  return (Y - y) / (X - x);
}

bool Equation::isCutOff()
{
  getCutDistance();
  
  if (!distanceCutState)
    return false;

  getIntersectPoint();
  float targetDistance = getTargetDistance();

  float dx1 = X - intersect._1_.X;
  float dy1 = Y - intersect._1_.Y;
  float dx2 = X - intersect._2_.X;
  float dy2 = Y - intersect._2_.Y; 

  float dist1 = sqrt(pow(dx1, 2) + pow(dy1, 2));
  float dist2 = sqrt(pow(dx2, 2) + pow(dy2, 2));

  std::cout << "distance 1 and 2 is " << dist1 << " " << dist2 << "\n";

  if (alg::valueInside(dist1, targetDistance - 5.0, targetDistance + 5.0) && alg::valueInside(dist2, targetDistance - 5.0, targetDistance + 5.0))
  {
    return false;
  }
  else if (alg::valueInside(dist1, targetDistance - 5.0, targetDistance + 5.0))
  {
    if (dist1 > dist2)
      return true;
    else
      return false;
  }
  else if (alg::valueInside(dist2, targetDistance - 5.0, targetDistance + 5.0))
  {
    if (dist1 < dist2)
      return true;
    else
      return false;
  }
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

  std::vector<double> eqt1GCD;
  eqt1GCD.resize(3);

  eqt1GCD[0] = fabs(eqt1.a);
  eqt1GCD[1] = fabs(eqt1.b);
  eqt1GCD[2] = fabs(eqt1.c);

  eqt1.GCD = findGCD(eqt1GCD);
  std::cout << "flag" << std::endl;

  eqt1.a /= eqt1.GCD;
  eqt1.b /= eqt1.GCD;
  eqt1.c /= eqt1.GCD;

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
  eqt1.dist = ( fabs(-eqt1.slope*m_a + m_b + (eqt1.slope*m_x - m_y)) / sqrt((pow(eqt1.slope, 2) + pow(1, 2))) );

  if (eqt1.dist < m_r)    
    distanceCutState = true;

  return eqt1.dist;
}

float Equation::getTangentPoint()
{
  float dx = X - m_a;
  float dy = Y - m_b;
  float dxr = -dy;
  float dyr = dx;
  float d = sqrt(pow(dx, 2) + pow(dy, 2));
  float rho = m_r / d;
  float ad = pow(rho, 2);
  float bd = rho*sqrt(1 - pow(rho, 2));

  tangent._1_.X = m_a + ad*dx + bd*dxr;
  tangent._1_.Y = m_b + ad*dy + bd*dyr;
  tangent._2_.X = m_a + ad*dx - bd*dxr;
  tangent._2_.Y = m_b + ad*dy - bd*dyr;
}

float Equation::getIntersectPoint()
{
  float dx = m_x - X;
  float dy = m_y - Y;
  float A = pow(dx, 2) + pow(dy, 2);
  float B = 2*( dx*(X - m_a) + dy*(Y - m_b) );
  float C = pow((X - m_a), 2) + pow((Y - m_b), 2) - pow(m_r, 2);

  float det = getDeterminant(A, B, C);

  if (det > 0)
  {
    float t = (-B + sqrt(det)) / (2*A);
    intersect._1_.X = X + t*dx;
    intersect._1_.Y = Y + t*dy;

    t = (-B - sqrt(det)) / (2*A);
    intersect._2_.X = X + t*dx;
    intersect._2_.Y = Y + t*dy;
  }
}

void Equation::getTransitPoint(Point2D& point)
{
  getTangentPoint();

  float dx1 = m_x - tangent._1_.X;
  float dy1 = m_y - tangent._1_.Y;
  float dx2 = m_x - tangent._2_.X;
  float dy2 = m_y - tangent._2_.Y; 

  float dist1 = sqrt(pow(dx1, 2) + pow(dy1, 2));
  float dist2 = sqrt(pow(dx2, 2) + pow(dy2, 2));

  std::cout << "angel 1 and 2 is " << dist1 << " " << dist2 << "\n";

  if (dist1 > dist2)
    point = tangent._2_;
  else
    point = tangent._1_;
  
}

void Equation::getEdgePoint(Point2D& point)
{
  if (point.X == m_x && point.Y == m_y)
  {
    point = point;

    return;
  }
  
  unsigned int size;
  
  float nSize = 0.0;
  float multiplier = 1.0;
  float tg_th = 0.0;
  float th = 0.0;
  float mod_th = 0.0;

  if (!initialPoint)
  {
    tg_th = getAngel(point.X, point.Y, m_a + 20.0, m_b);

    if (m_th < tg_th)
      multiplier *= -1.0;

    float dx = m_x - point.X;
    float dy = m_y - point.Y;
    float dist = sqrt( pow(dx, 2) + pow(dy, 2) );
  
    initialPoint = true;

    th = acos( fabs((2*pow(m_r, 2) - pow(dist, 2)) / (2*pow(m_r, 2))) ) * alg::rad2Deg();
    std::cout << th << " n " << tg_th << "\n";

    size = (int)th / 15;
    mod_th = th - (float)(size*15);

    if (mod_th != 0)
      size++;
    
    nPoint.resize(size);
    nAngel.resize(size);

    nSize = nAngel.size();

    if (mod_th != 0)
      nAngel[--nSize] = tg_th + (multiplier * th);
      
    for (std::size_t i = 0; i < nSize; i++)
    {
      nAngel[i] = tg_th + (multiplier * ((15.0 * (i+1))));

      nPoint[i].X = m_a + 20.0*cos( (nAngel[i]/alg::rad2Deg()) );
      nPoint[i].Y = m_b + 20.0*sin( (nAngel[i]/alg::rad2Deg()) );
    }
  }
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
  
  std::cout << "linier equation (Robot) is ";
  printf("y = %.2fx + (%.2f)\n", eqt1.slope, -(eqt1.slope*m_x) + m_y);

  std::cout << "tangent equation is ";
  printf("(%.2f)x^2 + (%.2f)x + (%.2f) = 0\n", eqt2.a, eqt2.b, eqt2.c);
  
  std::cout << "tangent point is ";
  printf("x = %.2f, y = %.2f\nx = %.2f, y = %.2f\n", tangent._1_.X, tangent._1_.Y, tangent._2_.X, tangent._2_.Y);

  std::cout << "intersect point is ";
  printf("x = %.2f, y = %.2f\nx = %.2f, y = %.2f\n", intersect._1_.X, intersect._1_.Y, intersect._2_.X, intersect._2_.Y);

  std::cout << "distance is " << getTargetDistance() << " target angel is " << m_th <<  "\n";
}

int main()
{
  float r = 20;

  Point2D initial_p(150, 300);
  Point2D sec_p(214.3867, 613.8931);
  Point2D center(200, 600);
  Point2D point;

  Equation eqt(initial_p);

  eqt.setCircle(center, r);
  eqt.setLinier(sec_p);

  bool state = eqt.isCutOff();
  eqt.getTransitPoint(point);

  std::cout << "cut off state " << state << "\n";
  std::cout << "cut off distance " << eqt.getCutDistance() << "\n";
  eqt.printResult();
  eqt.printEquation();
  std::cout << "tangent point is " << point.X << " " << point.Y << "\n";
  eqt.getEdgePoint(point);
  std::cout << "next point is " << point.X << " " << point.Y << "\n";
}