// FLAGLET package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

/*! 
 * \file flaglet_about.c
 * Print information about the B3LET package, including version
 * and build numbers. 
 *
 * Usage: flaglet_about
 *
 */

#include <stdio.h>

int main(int argc, char *argv[]) {

  printf("%s\n", "==========================================================");
  printf("%s\n", "  FLAGLET package");
  printf("%s\n", "  Fast 3D Wavelets on the Solid Sphere");
  printf("%s\n", "  By Boris Leistedt & Jason McEwen");

  printf("%s\n", "  See LICENSE.txt for license details.");

  printf("%s%s\n", "  Version: ", FLAGLET_VERSION);
  printf("%s%s\n", "  Build: ", FLAGLET_BUILD);
  printf("%s\n", "==========================================================");

  return 0;

}
