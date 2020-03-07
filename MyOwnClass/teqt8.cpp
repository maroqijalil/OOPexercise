#include <iostream>
#include <cstdlib>
#include <math.h>
#include <vector>
#include <array>

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
  double X{ 0.0 };
  double Y{ 0.0 };

  double m_a{ -1.0 };
  double m_b{ -1.0 };
  double m_r{ 0.0 };
  double m_x{ 0.0 };
  double m_y{ 0.0 };
  double m_th{ 0.0 };

  float err_r{ 0.5 };
  float err_d{ 1.0 };
  float th_step{ 15.0 };

protected:

public:
  struct couplePoint
  {
    Point2D _1_;
    Point2D _2_;
  };

  couplePoint intersect;
  couplePoint tangent;

  Equation() = default;
  Equation(float x, float y);
  Equation(const Point2D& point);

  ~Equation(){ };

  void setPosition(float x, float y);
  void setBall(const Point2D& point, float r);
  void setTarget(const Point2D& point);
  void getTarget(Point2D& point);
  float getDistance(float x1, float y1, float x2, float y2);

  bool isCutOff();
  float getCutDistance();
  void getTransitPoint(Point2D& point);
  void getEdgePoint(Point2D& point);
  void getNextPoint(Point2D& point);

  void countTangentPoint();
  void countIntersectPoint();
  float getDeterminant(float a, float b, float c);
  float getAngel(float x1, float y1, float x2, float y2);
};

Equation::Equation(float x, float y)
  : X{x}, Y{y}
{
}

Equation::Equation(const Point2D& point)
  : X{point.X}, Y{point.Y}
{
}

void Equation::setPosition(float x, float y)
{
  X = x;
  Y = y;
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
    float dist = sqrt( pow(dx, 2) + pow(dy, 2) );

    return (acos( (2*pow(m_r, 2) - pow(dist, 2)) / (2*pow(m_r, 2)) ) * alg::rad2Deg()) * multiplier;
  }

  return 0.0;
}

float Equation::getDistance(float x1, float y1, float x2, float y2)
{
  float dx = x1 - x2;
  float dy = y1 - y2;

  return (sqrt(pow(dx, 2) + pow(dy, 2)));
}

void Equation::setBall(const Point2D& point, float r)
{
  m_a = point.X;
  m_b = point.Y;
  m_r = r;
}

void Equation::setTarget(const Point2D& point)
{
  float dist = getDistance(point.X, point.Y, m_a, m_b);
  float dx = point.X - m_a;
  float dy = point.Y - m_b;
  float th = atan2(dy, dx) * alg::rad2Deg();

  std::cout << th << " " << dist << "\n";
  
  m_x = m_a + ( m_r * cos((th + 180.0) / alg::rad2Deg()) );
  m_y = m_b + ( m_r * sin((th + 180.0) / alg::rad2Deg()) );

  std::cout << m_x << " " << m_y << "\n";

  m_th = getAngel(m_x, m_y, m_a + m_r, m_b);

  if (m_th > 180.0)
    m_th -= 360.0;
}

void Equation::getTarget(Point2D& point)
{
  point.X = m_x;
  point.Y = m_y;
}

float Equation::getDeterminant(float a, float b, float c)
{
  return pow(b, 2) - 4*a*c;
}

float Equation::getCutDistance()
{
  float dx = m_x - X;
  float dy = m_y - Y;
  float A = dy / dx;
  float C = m_y - (A*m_x);

  return ( fabs(A*m_a + (-1.0*m_b) + C) / sqrt((pow(A, 2) + pow(-1.0, 2))) );
}

bool Equation::isCutOff()
{
  if (getDistance(X, Y, m_a, m_b) < m_r)
    return true;
  
  if (getCutDistance() >= m_r)
    return false;

  countIntersectPoint();
  float targetDistance = getDistance(X, Y, m_x, m_y);

  float dx1 = X - intersect._1_.X;
  float dy1 = Y - intersect._1_.Y;
  float dx2 = X - intersect._2_.X;
  float dy2 = Y - intersect._2_.Y; 

  float dist1 = sqrt(pow(dx1, 2) + pow(dy1, 2));
  float dist2 = sqrt(pow(dx2, 2) + pow(dy2, 2));

  if (alg::valueInside(dist1, targetDistance - err_d, targetDistance + err_d) && alg::valueInside(dist2, targetDistance - err_d, targetDistance + err_d))
  {
    return false;
  }
  else if (alg::valueInside(dist1, targetDistance - err_d, targetDistance + err_d))
  {
    if (dist1 > dist2)
      return true;
    else
      return false;
  }
  else if (alg::valueInside(dist2, targetDistance - err_d, targetDistance + err_d))
  {
    if (dist1 < dist2)
      return true;
    else
      return false;
  }
}

void Equation::countTangentPoint()
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

void Equation::countIntersectPoint()
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
  countTangentPoint();

  if (!alg::valueInside(getDistance(X, Y, m_a, m_b), m_r - err_r, m_r + err_r))
  {
    float dx1 = m_x - tangent._1_.X;
    float dy1 = m_y - tangent._1_.Y;
    float dx2 = m_x - tangent._2_.X;
    float dy2 = m_y - tangent._2_.Y;

    // std::cout << tangent._1_.X << ".. " << tangent._1_.Y << std::endl;
    // std::cout << tangent._2_.X << " .." << tangent._2_.Y << std::endl;

    float dist1 = sqrt(pow(dx1, 2) + pow(dy1, 2));
    float dist2 = sqrt(pow(dx2, 2) + pow(dy2, 2));

    if (dist1 > dist2)
      point = tangent._2_;
    else
      point = tangent._1_;
  }
  else 
    point = point;
  
  // std::cout << point.X << " " << point.Y << " ...." << std::endl;
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
  float position_th = 0.0;
  float th = 0.0;
  float target_th = 0.0;

  position_th = getAngel(point.X, point.Y, m_a + m_r, m_b);

  if (m_th < position_th)
    multiplier *= -1.0;

  th = fabs( getAngel(m_x, m_y, point.X, point.Y) );
  th -= th_step;

  if (th >= 0)
  {
    target_th = position_th + (multiplier * th_step);

    point.X = m_a + m_r*cos( (target_th/alg::rad2Deg()) );
    point.Y = m_b + m_r*sin( (target_th/alg::rad2Deg()) );
  }
  else
  {
    target_th = position_th + (multiplier * th);

    point.X = m_a + m_r*cos( (target_th/alg::rad2Deg()) );
    point.Y = m_b + m_r*sin( (target_th/alg::rad2Deg()) );
  }
  // std::cout << point.X << " " << point.Y << std::endl;
}

void Equation::getNextPoint(Point2D& point)
{
  float dist = getDistance(X, Y, m_a, m_b);

  // std::cout << "................... " << dist << std::endl;
  
  if (dist <= m_r || alg::valueInside(dist, m_r - err_r, m_r + err_r))
  {
    if (alg::valueInside(dist, m_r - err_r, m_r + err_r))
    {
      float dx = X - m_a;
      float dy = Y - m_b;
      float change_th = atan2(dy, dx) * alg::rad2Deg();
      
      point.X = m_a + ( m_r * cos(change_th / alg::rad2Deg()) );
      point.Y = m_b + ( m_r * sin(change_th / alg::rad2Deg()) );
      // std::cout << point.X << " " << point.Y << std::endl;

      getEdgePoint(point);
    }
    else
      point = point;    
  }
  else
    getTransitPoint(point);  
}

int main ()
{
  int p = 0;
  int x, y;
  float r = 40.0;
  float th  = -20.0 / alg::rad2Deg();
  float cXB = 200.0;
  float cYB = 300.0;
  
  Point2D center(cXB, cYB);
  Point2D target(140, 310);
  Point2D point;

  Equation eqt;

  for (int i = 0; i <= 40; i++)
  {
    if (p == 0)
    {
      eqt.setPosition(150, 300);
      x = 210;
      y = 280 + i;
    }
    else if (p == 1)
    {
      eqt.setPosition(point.X, point.Y);
      x = 220 - i;
      y = 320;
    }
    else if (p == 2)
    {
      // eqt.setPosition(point.X, point.Y);
      x = 180;
      y = 320 - i;
    }
    else if (p == 3)
    {
      eqt.getTarget(point);
      x = 180 + i;
      y = 280;
    }

    eqt.setBall(center, r);
    eqt.setTarget(target);

    if (eqt.isCutOff())
      eqt.getNextPoint(point);
    else
      eqt.getTarget(point);
    
    std::cout << point.X << " " << point.Y << std::endl;

    if (x == 210 && y == 319)
    {
      std::cout << "1..\n";
      i = 0;
      p = 1;
    }
    else if (x == 180 && y == 320)
    {
      std::cout << "2..\n";
      i = 0;
      p = 2;
    }
    else if (x == 180 && y == 280)
    {
      std::cout << "3..\n";
      i = 0;
      p = 3;
    }
  }
}