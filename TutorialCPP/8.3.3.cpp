#include <iostream>
#include <cassert>

using namespace std;

class Stack
{
  int arr[10];
  int index;

public:
  void reset()
  {
    index = 0;
  }

  bool push(int num)
  {
    if (index == 10)
      return false;
    arr[index] = num;
    return index++;
  }

  int pop()
  {
    assert(index != 0 && "The array is empty");
    index--;
  }

  void print()
  {
    cout << "( ";
    for (int i=0; i<index; i++)
      cout << arr[i] << " ";
    cout << ")" << endl;
  }
};

int main ()
{
  Stack stack;
	stack.reset();
 
	stack.print();
 
	stack.push(5);
	stack.push(3);
	stack.push(8);
	stack.print();
 
	stack.pop();
	stack.print();
 
	stack.pop();
	stack.pop();
 
	stack.print();
  stack.pop();
 
	stack.print();
 
	return 0;
}