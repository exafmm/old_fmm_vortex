#include "../misc/constants.h"

double ierfc(double u) {
  double f;
  f = 1/sqrt(pi)*exp(-u*u)-u*erfc(u);
  return f;
}
