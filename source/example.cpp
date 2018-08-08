#include <iostream>
#include "FFInt.hpp"

int main(int argc, char **argv) {
   firefly::FFInt a(3, 7);
   firefly::FFInt b(4, 7);
   a /= b;
   std::cout << a.n << std::endl;
   return 0;
}
