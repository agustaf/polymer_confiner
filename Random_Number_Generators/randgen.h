
#ifndef _RANDGEN_H_INCLUDED_
#define _RANDGEN_H_INCLUDED_
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <time.h>
#include <limits>
namespace randgenbase{
	void init_randgen(unsigned int* init);
	void autoinit_randgen();
	double randgen();
	void advance_randgen(const int count);
	static double gaussrandgen();
}
#endif

