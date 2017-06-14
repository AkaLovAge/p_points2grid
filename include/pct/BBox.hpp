#ifndef BBOX_HPP
#define BBOX_HPP

#include <iostream>
#include <vector>
#include "pct/Point.hpp"
#include "pct/Grid.hpp"

using namespace std;

class BBox
{
	private:
		struct Point min;
		struct Point max;
	public:
		BBox(): min(0.0, 0.0, 0.0), max(0.0, 0.0, 0.0){};
		BBox(struct Point * _min, struct Point * _max): min(*_min), max(*_max){};
		BBox(const BBox& box): min(box.min), max(box.max){};
		~BBox();
		void update(struct Point * _min, struct Point * _max);
		void updateMin(double _x, double _y, double _z);
		void updateMax(double _x, double _y, double _z);
		struct Point* centroid();
		double area();
		vector<BBox*> subdivide();
		bool validate();
		int getCorners(struct Point* ll, struct Point* lr, struct Point* ur, 
				struct Point* ul);
		int getGridCells(Grid* grid, int cellIdxs[]); // get array of gridcells that intersect corners
		int intersect(const BBox* other);
		void print() const;
};



#endif
