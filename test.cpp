
#define _NTOOLS_USE_INTEGRATION
#include "numericaltools.h"

#include <iostream>
#include <vector>

int main() {
  std::vector<double> mesh(101);
  for(int i = 0; i <= 100; ++i) mesh[i] = 0.01*i;
  std::cout << "0.5 ?= " << NTools::BooleIntegrator::integrate(0.01, mesh) << std::endl;
}
