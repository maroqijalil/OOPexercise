#include <iostream>
#include <string>
#include <cassert>
#include <cstdint>

using namespace std;

class RGBA
{
private:
  uint_fast8_t m_red{ 0 };
  uint_fast8_t m_green{ 0 };
  uint_fast8_t m_blue{ 0 };
  uint_fast8_t m_alpha{ 255 };

public:
  RGBA() = default;

  RGBA(uint_fast8_t red, uint_fast8_t green, uint_fast8_t blue)
  {
    m_red = red;
    m_green = green;
    m_blue = blue;
  }

  void print()
  {
    cout << "r=" << (int)m_red << " g=" << (int)m_green << " b=" << (int)m_blue << " a=" << (int)m_alpha << "\n";
  }
};

int main ()
{
  RGBA teal(0, 127, 127);
	teal.print();

	return 0;
}