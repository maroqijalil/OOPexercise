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
  std::vector<Point2D> nPoint;
  std::vector<float> nAngel;

  unsigned int pQeue = 0u;

  double X{ 0.0 };
  double Y{ 0.0 };

  double m_a{ -1.0 };
  double m_b{ -1.0 };
  double m_r{ 0.0 };
  double m_x{ 0.0 };
  double m_y{ 0.0 };
  double m_th{ 0.0 };

protected:

public:
  struct equationComp
  {
    double D{ 0.0 };
    double a{ 0.0 };
    double b{ 0.0 };
    double c{ 0.0 };
    double dist{ -1.0 };
    double slope{ 0.0 };
  };

  struct couplePoint
  {
    Point2D _1_;
    Point2D _2_;
  };

  equationComp eqt;

  couplePoint intersect;
  couplePoint tangent;

  Equation() = default;
  Equation(float x, float y);
  Equation(const Point2D& point);

  ~Equation(){ };

  void setBall(const Point2D& point, float r);
  void setTarget(const Point2D& point);
  float getDistance(float x1, float y1, float x2, float y2);

  bool isCutOff();
  float getCutDistance();
  void getTransitPoint(Point2D& point);
  void getEdgePoint(Point2D point);
  void getNextPoint(Point2D& point);

  void printResult();
  void printEquation();

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
  // float th = acos( (pow(dist, 2) + pow(dy, 2) - pow(dx, 2)) / (2*dist*dy) ) * alg::rad2Deg();
  float th = atan2(dy, dx) * alg::rad2Deg();

  std::cout << th << " " << dist << "\n";
  
  m_x = m_a + ( m_r * cos((th + 180.0) / alg::rad2Deg()) );
  m_y = m_b + ( m_r * sin((th + 180.0) / alg::rad2Deg()) );

  // std::cout << m_x << " " << m_y << "\n";

  m_th = getAngel(m_x, m_y, m_a + m_r, m_b);

  if (m_th > 180.0)
    m_th -= 360.0;
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

  if (!alg::valueInside(getDistance(point.X, point.Y, m_a, m_b), m_r - 1.0, m_r + 1.0))
  {
    float dx1 = m_x - tangent._1_.X;
    float dy1 = m_y - tangent._1_.Y;
    float dx2 = m_x - tangent._2_.X;
    float dy2 = m_y - tangent._2_.Y; 

    float dist1 = sqrt(pow(dx1, 2) + pow(dy1, 2));
    float dist2 = sqrt(pow(dx2, 2) + pow(dy2, 2));

    if (dist1 > dist2)
      point = tangent._2_;
    else
      point = tangent._1_;
  }
  else 
    point = point;
}

void Equation::getEdgePoint(Point2D point)
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

  tg_th = getAngel(point.X, point.Y, m_a + m_r, m_b);

  if (m_th < tg_th)
    multiplier *= -1.0;

  th = fabs( getAngel(m_x, m_y, point.X, point.Y) );

  size = (int)th / 15;
  mod_th = th - (float)(size*15);

  if (mod_th != 0)
    size++;
  
  nPoint.resize(size);
  nAngel.resize(size);

  nSize = nAngel.size();

  if (mod_th != 0)
  {
    nAngel[--nSize] = tg_th + (multiplier * th);

    nPoint[nSize].X = m_a + m_r*cos( (nAngel[nSize]/alg::rad2Deg()) );
    nPoint[nSize].Y = m_b + m_r*sin( (nAngel[nSize]/alg::rad2Deg()) );
  }
    
  for (std::size_t i = 0; i < nSize; i++)
  {
    nAngel[i] = tg_th + (multiplier * ((15.0 * (i+1))));

    nPoint[i].X = m_a + m_r*cos( (nAngel[i]/alg::rad2Deg()) );
    nPoint[i].Y = m_b + m_r*sin( (nAngel[i]/alg::rad2Deg()) );
  }
}

void Equation::getNextPoint(Point2D& point)
{
  point.X = nPoint[pQeue].X;
  point.Y = nPoint[pQeue].Y;
  
  pQeue++;
}

void Equation::printEquation()
{
  std::cout << "circle equation is ";
  printf("(x - (%.2f))^2 + (y - (%.2f))^2 = %.2f\n", m_a, m_b, pow(m_r, 2));
  
  std::cout << "linier equation (Robot) is ";
  printf("y = %.2fx + (%.2f)\n", eqt.slope, -(eqt.slope*m_x) + m_y);
  
  std::cout << "tangent point is ";
  printf("x = %.2f, y = %.2f\nx = %.2f, y = %.2f\n", tangent._1_.X, tangent._1_.Y, tangent._2_.X, tangent._2_.Y);

  std::cout << "intersect point is ";
  printf("x = %.2f, y = %.2f\nx = %.2f, y = %.2f\n", intersect._1_.X, intersect._1_.Y, intersect._2_.X, intersect._2_.Y);

  std::cout << "distance is " << " target angel is " << m_th <<  "\n";
}

// int main()
// {
//   std::array<Point2D, 5> arrP;

//   arrP[0].X = 150.0;
//   arrP[0].Y = 600.0;
//   arrP[1].X = 219.469;
//   arrP[1].Y = 304.578;
//   arrP[2].X = 219.99;
//   arrP[2].Y = 299.383;
//   arrP[3].X = 219.15;
//   arrP[3].Y = 294.23;
//   arrP[4].X = 218.794;
//   arrP[4].Y = 293.16;
   
//   float r = 20.0;
//   float th  = -20.0 / alg::rad2Deg();
//   float cXB = 200.0;
//   float cYB = 300.0;

//   for (int i = 0; i < 5; i++)
//   {+ 180.0)
//     Point2D initial_p = arrP[i];
//     Point2D sec_p((cXB + (cos(th) * r)), (cYB + (sin(th) * r)));
//     Point2D center(cXB, cYB);
//     Point2D point = initial_p;

//     Equation eqt(initial_p);

//     eqt.setBall(center, r);
//     eqt.setTarget(sec_p);

//     bool state = eqt.isCutOff();
//     eqt.getTransitPoint(point);

//     std::cout << "cut off state " << state << "\n";
//     std::cout << "cut off distance " << eqt.getCutDistance() << "\n";
//     std::cout << "tangent point is " << point.X << " " << point.Y << "\n";
//     eqt.getEdgePoint(point);
//     eqt.getNextPoint(point);
//     std::cout << "next point is " << point.X << " " << point.Y << "\n";
//     eqt.getNextPoint(point);
//     std::cout << "next point is " << point.X << " " << point.Y << "\n";
//   }
// }

int main ()
{
  int p = 0;
  float k;
  float r = 20.0;
  float th  = -20.0 / alg::rad2Deg();
  float cXB = 200.0;
  float cYB = 300.0;
  
  Point2D initial_p(100, 300);
  Point2D sec_p((cXB + (cos(th) * r)), (cYB + (sin(th) * r)));
  Point2D center(cXB, cYB);
  Point2D *point;

  Equation eqt(initial_p);

  for (int i = 0; i < 900; i++)
  {
    if (p == 0)
      point = new Point2D((k = 220 + i), -200);
    else if (p == 1)
      point = new Point2D(k, (-200 + i));
    eqt.setBall(center, r);
    eqt.setTarget(*point);
    if (i == 898 && p == 0){ i = 0; p = 1;}
  }
}