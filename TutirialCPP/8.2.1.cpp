#include <iostream>

using namespace std;

class IntPair
{
public:
  int num_1;
  int num_2;

  void set(int num1, int num2)
  {
    num_1 = num1;
    num_2 = num2;
  }

  void print()
  {
    cout << "Pair(" << num_1 << ", " << num_2 << ")" << endl;
  }
};

int main ()
{
  IntPair p1;
  p1.set(1, 1);

  IntPair p2{ 2, 2 };

  p1.print();
  p2.print();
}