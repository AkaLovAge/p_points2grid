#include <stdio.h>
#include <stdlib.h>
#include <iomanip>
#include <math.h>
#include "pct/Point.hpp"
#include "pct/Grid.hpp"
#include "pct/DTypes.hpp"
#include "pct/BlockScheme.hpp"

BlockScheme::BlockScheme() {
	cols = 0;
	rows = 0;
	b_cols = 0;
	b_rows = 0;
	datatype = DT_Unknown;
	proc_count = 0;
	block_grid = NULL;
};

BlockScheme::~BlockScheme() {
	cols = 0;
	rows = 0;
	b_cols = 0;
	b_rows = 0;
	datatype = DT_Unknown;
	proc_count = 0;
	if (block_grid != NULL) {
		delete(block_grid);
	}
	block_grid = NULL;
};

void BlockScheme::getBlockPosition(int rank, int* position) {
	position[1] = rank / cols;
	position[0] = rank % cols;
}
// Get the ID of a block based on its col/row location
int BlockScheme::getBlockId(int x, int y) {
	int idx = y * cols + x;
	return idx;
}

// Return the blockID based on the pixel coordinates
int BlockScheme::getBlock(long x, long y) {
	int i = 0;
	int b_x = floor((float)x / (float)b_cols);
	int b_y = floor((float)y / (float)b_rows);
	return getBlockId(b_x, b_y);
}
//TODO: Function to return a grid object from a block index
Grid BlockScheme::getBlockGrid(int blockId, int res) {
	int pos[2] = {0, 0};
	// Make sure the x and y indexes are assigned properly by this function
	getBlockPosition(blockId, &pos[0]);
	int tmp_col,tmp_row;
	Grid grid;	
	struct Point* _origin = new Point();
	double tmp_x,tmp_y,tmp_z;
	double off_x,off_y;
	printf("Block %i: Position [%i,%i]\n", blockId, pos[0],pos[1]);
	off_x = (pos[0] * block_grid->resX);
	off_y = (pos[1] * block_grid->resY);
	printf("Block %i: Offsets [%f, %f]\n", blockId,off_x, off_y);
	tmp_x = block_grid->origin.x + off_x;

	tmp_y = block_grid->origin.y - off_y;
	tmp_z = 0.0;
	_origin->update(tmp_x,tmp_y,tmp_z);
	printf("New grid will have origin (%f,%f,%f)\n", _origin->x,_origin->y,_origin->z);
	grid.origin = *_origin;
	grid.cols = b_cols;
	grid.rows = b_rows;
	grid.datatype = datatype;
	grid.resX = res;
	grid.resY = res;
	//Grid* grid = new Grid(_origin, b_cols, b_rows, datatype, res, res);
	printf("Grid has %i,%i dims, block_dims (%i, %i)\n", grid.cols, grid.rows, b_cols, b_rows);
	printf("Block Grid %i extent (%f,%f,%f,%f)\n", blockId, grid.origin.x,grid.getMinY(),grid.getMaxX(),grid.origin.y);
	printf("Grid created\n");
	delete(_origin);
	return grid;
}
// Get the process rank based on the blockId
int BlockScheme::getBlockRank(int blockId) {
	return blockId % proc_count;
}

int BlockScheme::getBlockCount(int rank) {
	int gridSize = cols * rows;
	int count = (int)gridSize / proc_count;
	int rem = fmod(gridSize, proc_count);
	//if (rem > 0 && rank <= rem) {
	// Try to fix a bug where last process thinks it has 1 more block than it really does
	if (rem > 0 && rank < rem) {
		count++;
	}
	return count;
}

int* BlockScheme::getBlocks(int rank) {
	int gridSize = cols * rows;
	int counter = rank;
	int i = 0;
	int n_blocks = getBlockCount(rank);
	int* blockIds = (int*)malloc(sizeof(int) * n_blocks);
	for (i = 0; i < n_blocks; i++) {	
		blockIds[i] = counter;
		counter = counter + proc_count;
	}
	return blockIds;
};

/** 
 * Alternate constructor: based on grid dimensions and block
 * pixel limit.
 * @param grid: The grid to subdivide
 * @param block_limit: number of cells per block
 */
BlockScheme::BlockScheme(Grid* grid, int block_limit, DType _datatype, int _proc_count) {
	int gridSize = grid->cols * grid->rows;
	//Start with block that is basically square from square root
	b_cols = (int)sqrt(block_limit);
	b_rows = block_limit / b_cols;
	
	int c_rem = fmod(grid->cols, b_cols);
	int r_rem = fmod(grid->rows, b_rows);
	cols = ceil((float)grid->cols / (float)b_cols);
	rows = ceil((float)grid->rows / (float)b_rows);
	b_cols = grid->cols / cols;
	b_rows = grid->rows / rows;
	printf("Leftovers: rows %i, cols %i\n", c_rem, r_rem);
	printf("Block matrix: dims(%ix%i) cols %i, rows %i\n", b_cols, b_rows, cols, rows);
	proc_count = _proc_count;
	datatype = _datatype;
	double res_x = (grid->getMaxX() - grid->origin.x) / cols;
	double res_y = (grid->origin.y - grid->getMinY()) / rows;
	block_grid = new Grid(&grid->origin, cols, rows, _datatype, res_x, res_y);
	printf("BlockScheme has (%i, %i) cols, res(%f,%f)\n", block_grid->cols, block_grid->rows, block_grid->resX, block_grid->resY);

};


