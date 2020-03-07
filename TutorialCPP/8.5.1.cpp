#include <iostream>
#include <string>
#include <cassert>

using namespace std;

class Ball
{
  string m_color{ "black" };
  double m_radius{ 10.0 };

public:
  Ball() = default;

  Ball(string color, double radius)
  {
    m_color = color;
    m_radius = radius;
  }

  Ball(double radius)
  {
    m_radius = radius;
  }

  Ball(string color)
  {
    m_color = color;
  }

  void print()
  {
    cout << "color: " << m_color << ", radius: " << m_radius << "\n";
  }
};

int main ()
{
  Ball def;
	def.print();
 
	Ball blue("blue");
	blue.print();
	
	Ball twenty(20.0);
	twenty.print();
	
	Ball blueTwenty("blue", 20.0);
	blueTwenty.print();
 
	return 0;
}