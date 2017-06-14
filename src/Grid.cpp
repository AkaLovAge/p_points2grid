#include <math.h>
#include <time.h>
#include <stdio.h>
#include <string.h>
#include <stdexcept>
#include <mpi.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <float.h>
#include <algorithm>
#include "pct/Point.hpp"
#include "pct/Grid.hpp"
#include "pct/Pixel.hpp"
#include "pct/DTypes.hpp"

#include "gdal.h"
#include "cpl_conv.h" // for CPLMalloc()
#include "cpl_string.h"
#include "ogr_spatialref.h"
using namespace std;

Grid::Grid() {
	cols = 0;
	rows = 0;
	datatype = DT_Unknown;
	data = NULL;
	origin.update(0.0,0.0,0.0);
}
/* Note Currently origin is assuming lower-left corner for index function,
  This is not standard for raster data nor is it recommended. Need to change for
  interoperability
*/
Grid::Grid(const struct Point *_origin, int _cols, int _rows, DType _datatype, double _resX, double _resY) 
{
	origin = *_origin;
	cols = _cols;
	rows = _rows;
	datatype = _datatype;
	data = NULL;
	resX = _resX;
	resY = _resY;
}
/* Copy constructor, note: will not copy the data buffer*/
Grid::Grid(const Grid &other) {
	origin = other.origin;
	cols = other.cols;
	rows = other.rows;
	datatype = other.datatype;
	data = NULL; // NOTE: Doesn't copy data buffer
	resX = other.resX;
	resY = other.resY;
}



Grid::~Grid() {
	cols = 0;
	rows = 0;
	datatype = DT_Unknown;
	dealloc();
	resX = 0.0;
	resY = 0.0;
	origin.update(0.0,0.0,0.0);
	//delete(origin);
}
/** Allocate the memory necessary for the given grid **/
int Grid::alloc() {
	int d_size = sizeof(Pixel);
	if (!d_size > 0) {
		fprintf(stderr, "GridError: Invalid datatype prevented allocation\n");
		return -1;
	}
	if (data != NULL) {
		fprintf(stderr, "GridError: Memory already allocated for grid data\n");
		return -1;
	}
	fprintf(stdout, "Allocating %ix%i %zu size\n", cols, rows, cols*rows*sizeof(Pixel));
	//data = (Pixel*)malloc(sizeof(Pixel) * cols * rows);
	//int i = 0;
	int count = cols * rows;
	//for (i = 0; i < count; i++) {
	//	data[i] = new Pixel();
	//}
	int i = 0;
	data = new Pixel[count];
	for (i = 0; i < count; i++) {
		data[i].count = 0;
		data[i].sum = 0.0;
		data[i].filled = 0;
	}
	printf("Grid allocated successfully\n");
	//data = malloc(sizeof(cols * rows * d_size));
	//return data;
	return 1;
}

void Grid::dealloc() {
	if (data != NULL) {
		delete[] data;
		data = NULL;
	}
}

// Create a geotiff from the raster
int Grid::write(char* outPath, int epsg, int type) {
	// Type 1= Avg, 0 = count
	GDALDatasetH hDataset;
	GDALDataType d_type;
	GDALAllRegister();
	const char *pszFormat = "GTiff";
	GDALDriverH hDriver = GDALGetDriverByName( pszFormat );
	char **papszMetadata;
	char **papszOptions = NULL;
	OGRSpatialReferenceH hSRS = OSRNewSpatialReference( NULL );
	char *pszSRS_WKT = NULL;
	GDALRasterBandH hBand;
	OGRErr err;
	if ( hDriver == NULL )
		exit (1 );

	papszMetadata = GDALGetMetadata( hDriver, NULL );
	if ( CSLFetchBoolean( papszMetadata, GDAL_DCAP_CREATE, FALSE ) )
		printf( "Driver %s supports Create() method.\n", pszFormat );
	if ( CSLFetchBoolean( papszMetadata, GDAL_DCAP_CREATECOPY, FALSE ) )
		printf( "Driver %s supports CreateCopy() method.\n", pszFormat );
	if (type == 1) 
		d_type = GDT_Float32;
	else 
		d_type = GDT_Int32;
	
	hDataset = GDALCreate( hDriver, outPath, cols, rows, 1, d_type, papszOptions);

	double adfGeoTransform[6] = { origin.x, resX, 0, origin.y, 0, resY * -1};
	// Geotransform for the raster
	printf("Setting GeoTranform\n");
	GDALSetGeoTransform(hDataset, adfGeoTransform);
	printf("Import EPSG definition\n");
	err = OSRImportFromEPSG(hSRS, epsg);
	if (err != 0) {
		fprintf(stderr, "Error importing EPSG:%i , ErrCode: %i\n", epsg, err);
	}
	// Create WKT definition of projection
	printf("Translating SRS to WKT\n");
	OSRExportToWkt(hSRS, &pszSRS_WKT);
	OSRDestroySpatialReference(hSRS);
	printf("Setting output projection to %s\n", pszSRS_WKT);
	GDALSetProjection(hDataset, pszSRS_WKT);
	CPLFree( pszSRS_WKT);

	hBand = GDALGetRasterBand(hDataset, 1);
	// Set no data value for output raster
	if (type == 1)
		GDALSetRasterNoDataValue(hBand, -9999.0);
	else
		GDALSetRasterNoDataValue(hBand, 0);
	if (type == 1)
	{
		printf("Allocating avg array\n");
		float* outArr = (float*)malloc(sizeof(float)*cols * rows);
		getSumArr(outArr);
		GDALRasterIO(hBand, GF_Write, 0, 0, cols, rows, outArr, cols, rows, GDT_Float32, 0,0);
		free(outArr);
	}	
	else
	{
		printf("Allocating count array\n");
		int* outArr = (int*)malloc(sizeof(int) * cols * rows);
		getCountArr(outArr);
		int i = 0;
		for (i = 0; i < cellCount(); i++)
			printf("[%i:%i]\n", i, outArr[i]);

		GDALRasterIO(hBand, GF_Write, 0, 0, cols, rows, outArr, cols, rows, d_type, 0,0);
		free(outArr);
	}
	printf("Writing output raster to %s\n", &outPath[0]);
	/*Close the dataset*/
	printf("Closing Raster\n");

	GDALClose( hDataset);

	//free(outArr);
	return 1;
}

int Grid::cellCount() {
	return rows* cols;
}


size_t Grid::getSize() {
	size_t gridSize = (size_t)cols * rows * sizeof(Pixel);
	return gridSize;
}


double Grid::getMaxX() {
	double max = origin.x + (resX * cols);
	return max;
}

double Grid::getMinY() {
	double max = origin.y - (resY * rows);
	return max;
}

int Grid::getCol(double coord) {
	if (coord > getMaxX() || coord < origin.x) {
		return -1;
	}
	int idx = floor((coord - origin.x) / resX);
	
	return idx;
}

int Grid::getRow(double coord) {
	if (coord < getMinY() || coord > origin.y) {
		//printf("%f is out of bounds [%f, %f]\n", coord, getMinY(), origin->y);
		return -1;
	}
	int idx = floor((origin.y - coord) / resY);
	return idx;
}
/** Test if a Point intersects with the grid **/
int Grid::within(const struct Point *pt) {
	int flag = 0;
	if (pt->x >= origin.x && pt->x <= getMaxX()) {
		if (pt->y >= getMinY() && pt->y <= origin.y) {
			flag = 1;
		} else {
			cout << "Point out of Row range" << endl;
		} 
	} else {
		cout << "Point out of Col range" << endl;
	}
	return flag;
}
/* Assumes that the Point is in raster coordinates */
int Grid::set(const struct Point* pt) {
	int i;
	//int count = cols * rows;
	int col, row = 0;
	col = getCol(pt->x);
	row = getRow(pt->y);
	//i = (pt->y * rows) + pt->x;
	//if (i < 0 || i > count) {
	//	return 0;
	//}
	if (col < 0 || col >= cols) {
		return 0;
	}
	if (row < 0 || row >= rows) {
		return 0;
	}
	//printf("Point coord: (%f,%f,%f)\n", pt->x,pt->y,pt->z);
	//printf("Pixel index is: [%i,%i]\n", col, row);
	i = (row * cols) + col;
	// TODO: Check why this fails
	//printf("Looking for pixel at %i/%i\n", i, cols*rows);
	Pixel* cell = &data[i];
	//printf("Cell status: %i\n", cell->filled);
	if (cell->filled != 1) {
		cell->sum = (float)pt->z;
		cell->count = 1;
		cell->filled = 1;
		return 1;
	} else {
		float tmpSum = cell->sum + (float)pt->z;
		int tmpCount = cell->count +1;
		//printf("New sum = %f, count= %i\n", tmpSum, tmpCount);
		cell->sum = tmpSum;
		cell->count = tmpCount;
		return 1;
	}

}

void Grid::getSumArr(float* out) {
//	out = (float*)malloc(sizeof(float) * cols * rows);
	int count = cellCount();
	int i = 0;
	for (i = 0; i < count; i++) {
		if (data[i].filled == 1) {
			out[i] = data[i].avg();
		} else {
			out[i] = -9999.0;
		}
	}
	//return &out[0];
}


void Grid::getCountArr(int* out) {
//	out = (int*)malloc(sizeof(int) * cols * rows);
	int count = cellCount();
	int i = 0;
	for( i = 0; i < count; i++) 
	{
	//	if (data[i].count) {
			//printf("%i -> %i\n", i, data[i].count);
			out[i] = data[i].count;
	//	}
	//	else 
	//	{
	//		out[i] = 0;
	//	}
	}
	//return &out[0];
}

Pixel* Grid::get(int col, int row) {
	//Pixel* pix = NULL;
	if (data == NULL) {
		printf("Error: Array not initialized");
		exit(1);
	} else {
		int vectorIdx = row * rows + col;
		return &data[vectorIdx];
	}
}
// Return the cell id for the given point coordinates
int Grid::getCellIdx(double x, double y) {
	int col = getCol(x);
	int row = getRow(y);
	if (col < 0 || row < 0)
		return -1;
	//printf("Point(%f, %f) ->Grid(%i,%i)\n", x, y, col, row); 
	int idx = (row * cols) + col;
	//printf("Col: %i, Row: %i, Idx: %i\n", col, row, idx);
	return idx;
}

int Grid::getCell(const struct Point *pt, int* idx) {
	int i, j = 0;
	if (!within(pt))  {
	//	cout << "Point outside grid bounds\n";
		return 0;
	}
	j = getCol(pt->x);
	i = getRow(pt->y);
	idx[0] = j;
	idx[1] = i;
	return 1;
}
	
