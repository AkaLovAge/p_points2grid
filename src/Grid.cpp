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
#include <float.h>
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
		data[i].min = FLT_MAX;
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
	if (type > 0 && type < 3)
	{
		printf("Allocating min array\n");
		float* outArr = (float*)malloc(sizeof(float)*cols * rows);
		if (type == 1) 
			getMinArr(outArr);
		else // type is 2 : average
			getSumArr(outArr);
		GDALRasterIO(hBand, GF_Write, 0, 0, cols, rows, outArr, cols, rows, GDT_Float32, 0,0);
		free(outArr);
	} else 
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
	int flag = cell->set(pt->z);

	/*if (cell->filled != 1) {
		cell->sum = (float)pt->z;
		cell->count = 1;
		cell->filled = 1;
		return 1;
	} else {
	//	float tmpSum = cell->sum + (float)pt->z;
	//	int tmpCount = cell->count +1;
		//printf("New sum = %f, count= %i\n", tmpSum, tmpCount);
		cell->sum = cell->sum + (float)pt->z;
		cell->count = cell->count + 1;;
		return 1;
	}*/
	return flag;

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

void Grid::getMinArr(float* out) {
	float min = FLT_MAX;
	int i = 0;
	int count = cellCount();
	for (i = 0; i < count; ++i) {
		if (data[i].filled == 1) {
			out[i] = data[i].min;
		} else {
			out[i] = -9999.0;
		}
	}
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
		return NULL;
	} else {
		if (col < 0 || col > cols || row < 0 || row > rows)
			return NULL;
		int vectorIdx = row * rows + col;
		if (vectorIdx < 0 || vectorIdx > cols * rows)
			return NULL;
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
 
int Grid::fillNull(int radius) {
	Pixel* idx = NULL;
	Pixel* other = NULL;
	int i,j,k,l;
	//double sum = 0.0;
	//int cnt = 0;
	for(i = 0; i < rows; i++)
	{
		for (j = 0; j < cols; j++)
		{
			idx = get(j, i);
			if (!idx) 
				continue;
			if (!idx->is_filled()) 
			{
				//printf("Cell %i,%i null\n", j,i);
				for (k = -radius; k < radius; k++) 
				{
					for (l = -radius; l < radius; l++) 
					{
						//printf("Checking %i,%i\n", l+j, k+i);
						other = get(l+j,k+i);
						if (!other)
							continue;
						if (other->is_filled()) 
						{
							
							idx->set(other);
							//cnt++;
							//printf("Pix(%i,%i) set to %f\n", j, i, idx->avg());
						}
					}
				}
			}
		}
	}
	return 0;
}

/** Array structure { top_ghost, bottom_ghost, left_ghost, right_ghost} */
Pixel* Grid::createGhostCells(int radius) 
{
	//int new_cols = cols + 2 * radius;
	//int new_rows = rows + 2 * radius;
	//int i,j;
	int col_arr_len = cols * radius;
	int row_arr_len = rows * radius;
	Pixel* new_data = (Pixel *)malloc(sizeof(Pixel) * 2 *(col_arr_len + row_arr_len));
	//Pixel* new_data = (Pixel *)malloc(sizeof(Pixel) * new_cols * new_rows);
	
	//for (i = 0 + radius; i < rows; i++)
	//{
	//	for (j = 0 + radius; j < cols; j++) 
	//	{
	//		new_data[i] = get(j, i);
	//	}
	//}
	return new_data;
}

int Grid::fillGhostCells(Pixel* ghost, int radius)
{
	int idx = 0;
	int idx2 = 0;
	int i = 0;
	int col_arr_len = cols * radius;
	int row_arr_len = rows * radius;
	// fill top rows
	memcpy(ghost, data, sizeof(Pixel) * cols * radius);
	// fill bottim rows
	idx = cols * radius;
	idx2 = (cols * rows) - (cols * radius + 1); // -1 to align with index, may not be needed
	memcpy(&ghost[idx], &data[idx2], sizeof(Pixel) * cols * radius);
	// fill left cols
	idx = 0;
	idx2 = cols * 2 * radius;
	for (i = 0; i < rows; i++) {
		idx = col_arr_len * 2 + (i * radius);
		idx2 = i * cols;
		memcpy(&ghost[idx], &data[idx2], sizeof(Pixel) * radius);
	}
	// fill right cols
	for (i = 0; i < rows; i++)
	{
		idx = (col_arr_len * 2 + row_arr_len)  + i * cols - radius;
		idx2 = i * cols - radius;
		memcpy(&ghost[idx], &data[idx2], sizeof(Pixel) * radius);
	}
	return 0;
}

// Send row data to top neighbor
int Grid::putGhostTop(Pixel* ghost, int radius, int target, MPI_Win window)
{
	int idx = cols * radius * sizeof(Pixel);
	int cnt = sizeof(Pixel) * cols * radius;
	int flag = 0;
	flag = MPI_Put(ghost, cnt, MPI_BYTE, target, idx, cnt, MPI_BYTE, window);
	return flag;
}

// Send row data to bottom neighbor
int Grid::putGhostBottom(Pixel* ghost, int radius, int target, MPI_Win window)
{
	int idx = cols * radius;
	int cnt = sizeof(Pixel) * cols * radius;
	int flag = 0;
	flag = MPI_Put(&ghost[idx], cnt, MPI_BYTE, target, 0, cnt, MPI_BYTE, window);
	return flag;
}

int Grid::putGhostLeft(Pixel* ghost, int radius, int target, MPI_Win window)
{
	int idx = cols * radius * 2;
	int idx2 = idx + rows* radius;
	int cnt = sizeof(Pixel) * rows * radius;
	int flag = 0;
	flag = MPI_Put(&ghost[idx], cnt, MPI_BYTE, target, idx2, cnt, MPI_BYTE, window);
	return flag;
}

int Grid::putGhostRight(Pixel* ghost, int radius, int target, MPI_Win window)
{
	int idx = cols * radius * 2 + radius * rows;
	int idx2 = cols * radius * 2;
	int cnt = sizeof(Pixel)* rows * radius;
	int flag = 0;
	flag = MPI_Put(&ghost[idx], cnt, MPI_BYTE, target, idx2, cnt, MPI_BYTE, window);
	return flag;
}
