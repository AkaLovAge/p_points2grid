#ifndef PIXEL_HPP
#define PIXEL_HPP

typedef struct Pixel {
	int count;
	float sum;
	float min;
	char filled;
	Pixel();
	~Pixel();
	int set(float other);
	int is_filled() const;
	float avg() const;
	int set(const struct Pixel* other);
} Pixel;

#endif
