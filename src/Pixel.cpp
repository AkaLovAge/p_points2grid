#include "pct/Pixel.hpp"
#include <float.h>
#include <math.h>
#include "pct/util.hpp"
Pixel::Pixel() {
	sum = 0.0;
	min = FLT_MAX;
	count = 0;
	filled = 0;
}

Pixel::~Pixel() {
	sum = 0.0;
	min = FLT_MAX;
	count = 0;
	filled = 0;
}

int Pixel::is_filled() const 
{
	if (filled == 1)
		return 1;
	else
		return 0;
}

int Pixel::set(float other) 
{
	sum = sum + other;
	count++;
	if (other < min && !compareFloat(other, -9999.0, 1)) {
		min = other;
	}
	if (!filled)
		filled = 1;
	return count;
}

int Pixel::set(const struct Pixel* other) 
{
	if (!other->is_filled())
		return 0;
	if (other->min < min && !compareFloat(other->min,  -9999.0, 1)) {
		min  = other->min;
	}
	sum = sum + other->avg();
	count++;
	if  (!is_filled())
		filled = 1;
	return 1;	
}


float Pixel::avg() const 
{
	if (is_filled())
		return sum / count;
	else 
		return 0.0;
}

	
