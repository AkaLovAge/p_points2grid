#include <stdio.h>
#include <stdlib.h>
#include <iomanip>
#include <unistd.h>
#include "points2grid/lasfile.hpp"
#include "pct/Grid.hpp"
#include "pct/BBox.hpp"
#include "pct/DTypes.hpp"
#include "pct/Point.hpp"
#include "pct/BlockScheme.hpp"
#include "pct/FileCollection.hpp"
#include "pct/Pixel.hpp"
#include "pct/util.hpp"
#include <float.h>
#include <math.h>
#include <time.h>
#include <mpi.h>
#include <boost/program_options.hpp>
#include <vector>
#include <string.h>


namespace po = boost::program_options;

using namespace std;

int main(int argc, char* argv[]) {
	/****************************************************/
	/*   Parameter declaration                          */
	/****************************************************/
	// Index params
	int i = 0, j = 0,k = 0;
	
	// Directory path params
	char input_path[2056];
	char file_list[1024];
	char scratch_path[2056];
	char output_path[2056];
	char tmp_path[2056];

	// MPI Params
	int world_size, world_rank, mpi_err;
	MPI_Comm world_comm = MPI_COMM_WORLD;
	//MPI_Info info = MPI_INFO_NULL;
	MPI_Status status;
	MPI_Request request;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(world_comm, &world_size);
	MPI_Comm_rank(world_comm, &world_rank);
	MPI_Errhandler_set(world_comm, MPI_ERRORS_RETURN);
	double starttime, endtime;
	// Set memory limit for grid allocations
	size_t buffer_lim = 0;
	size_t max_buffer = 0;
	// MSG Buffers
	int recv_buf[6] = {0,0,0,0,0,0};
	// File specific parameters
	int file_off, file_blk, file_end;
	file_off = file_blk = file_end = 0;
	char ext[5] = ".las";

	// Metadata
	double g_mins[3] = {DBL_MAX,DBL_MAX,DBL_MAX}; // Global min coord
	double g_maxs[3] = {-DBL_MAX,-DBL_MAX,-DBL_MAX}; // Global max coord
	double l_mins[3] = {DBL_MAX,DBL_MAX,DBL_MAX}; // Local min coord
	double l_maxs[3] = {-DBL_MAX,-DBL_MAX,-DBL_MAX}; // Local max coord

	// File specific params
	FileCollection *g_files = NULL; // Global file list
	int g_n_files = 0; // Global file count
	int l_n_files = 0; // Local file count
	BBox* l_file_bbox = NULL; // Local array of file bboxes
	// Grid specific params
	struct Point* origin = new Point();
	DType datatype = DT_Float32;
	struct Grid* g_grid = NULL; // Global grid
	int g_cols = 0; // Global grid column count
	int g_rows = 0; // Global grid row count
	char** block_file_buffers = NULL;
	//vector<vector<string> >block_files;
	//Specify resolution > should be parameterized in prod
	double res = 0.0f;
	// Block scheme params
	struct BlockScheme* blk_scheme = NULL;
	int *l_blks = NULL; // Local block array
	//int *g_blks = NULL; // Global block array
	int l_n_blks = 0; // Local block count
	int g_n_blks = 0; // Global block count
	vector<vector<int> > blocks;
	int *block_hash = NULL;
	int *g_blk_cnts = NULL;
	int *blk_n_files = NULL; // Array holding block specific file counts
	// Las specific params
	las_file las;
	long long g_n_pts = 0; // Global point count
	long long l_n_pts = 0; // Local point count

	po::options_description general("General Options"), desc;

	general.add_options()
		("help", "produce a help message")
		("input_file_path,i",po::value<std::string>(),"required. Name of input directory of LAS files")
		("file_list,f",po::value<std::string>(),"optional. Path to file listing inputs")
		("scratch_path,s", po::value<std::string>(),"required. Name of scratch directory")
		("output_path,o",po::value<std::string>(),"required. Name of output directory")
		("resolution,r", po::value<float>(),"required. Output Resolution")
		("block_size,b",po::value<size_t>(),"required. Size of output blocks");
	desc.add(general);

	po::variables_map vm;

	try {
		po::store(po::parse_command_line(argc, argv, desc), vm);

		if (vm.count("help")) {
			printf("--------------------------------------------------\n");
			printf("           p_points2grid v0.1.1\n");
			printf("--------------------------------------------------\n");
			printf("Usage BlockLas -i inputPath -s scratch_path");
			exit(0);
		}
		po::notify(vm);

		if (vm.count("input_file_path")) {
			strncpy(input_path, vm["input_file_path"].as<std::string>().c_str(), sizeof(input_path));
		} else {
			throw std::logic_error("Input File Path required");
		}
		if (vm.count("file_list")) {
			strncpy(file_list, vm["file_list"].as<std::string>().c_str(), sizeof(file_list));
		}
		if (vm.count("scratch_path")) {
			strncpy(scratch_path, vm["scratch_path"].as<std::string>().c_str(), sizeof(scratch_path));
		} else {
			throw std::logic_error("Scratch Path Required");
		}
		if (vm.count("output_path")) {
			strncpy(output_path, vm["output_path"].as<std::string>().c_str(), sizeof(output_path));
		}	else {
			throw std::logic_error("Output Path Required");
		}
		if(vm.count("resolution")) {
			res = vm["resolution"].as<float>();

			if (res <= 0) {
				throw std::logic_error("Resolution must be greater than 0");
			}
		}
		if(vm.count("block_size")) {
			buffer_lim = vm["block_size"].as<size_t>();

			if (buffer_lim <= 0) {
				throw std::logic_error("Buffer Limit must be greater than 0");
			}
		}
	} catch (std::exception& e) {
		fprintf(stderr, "Error: %s \n", e.what());
		MPI_Finalize();
		exit(1);

	}


	/**************************************************/
	/*       Begin first file scan                    */
	/**************************************************/
	starttime = MPI_Wtime();
	printf("[%i] Beginning file collection from %s\n", world_rank, 
			&input_path[0]);
	
	// Set maximum blck size based on buffer limit
	max_buffer = buffer_lim / sizeof(Pixel);
	g_files = new FileCollection(input_path, &ext[0]);
	// Check the number of available LAS files
	g_n_files = g_files->countFiles();
	// Create tif output dir
	//printf("Creating tmp dir under %s\n", &scratch_path[0]);
	sprintf(&tmp_path[0], "%s/blocks", scratch_path);
	//printf("[%i] Using tmp dir: %s\n", world_rank,tmp_path);
	struct stat st;
	if (world_rank == 0) 
	{
		if (stat(tmp_path, &st) == -1) 
			mkdir(tmp_path, 0700);
		l_n_files = g_files->getMetadata(0, g_n_files);
	} else 
	{
		// Set file Block size
		file_blk = ceil((float)g_n_files /(float)(world_size-1));
		file_off = file_blk * (world_rank-1);
	
	
	
		if (file_off + file_blk > g_n_files) 
		{
			file_end = g_n_files;
			file_blk = g_n_files - file_off;
		} else
		{
			file_end = file_off + file_blk;
		}
		
		// Read subset of file paths from dir
		l_n_files = g_files->getMetadata(file_off, file_end);
		l_file_bbox = new BBox[l_n_files];
		//printf("[%i] Read %i from %i-%i  files\n", world_rank, l_n_files, 
		//		file_off, file_end);
		for (i = file_off; i < file_end; i++) 
		{
		//	printf("[%i] File %i/%i %s\n", world_rank, i, file_end, g_files->fileList[i]);
			las.open(g_files->fileList[i]);
			l_n_pts = l_n_pts + (long long) las.points_count();
			compareMin(l_mins, las.minimums());
			compareMax(l_maxs, las.maximums());
			double *tmp_min = las.minimums();
			double *tmp_max = las.maximums();
			j = i - file_off;
			l_file_bbox[j].updateMin(tmp_min[0], tmp_min[1], tmp_min[2]);
			l_file_bbox[j].updateMax(tmp_max[0], tmp_max[1], tmp_max[2]);
			las.close();
		}
	}

	
	endtime = MPI_Wtime();
	printf("[%i] Metadata gathered in %f seconds\n", world_rank, 
			endtime - starttime);
	MPI_Barrier(world_comm);
	/*****************************************************/
	/*         Gather global min/max point count         */
	/*                 COMMUNICATIONS                    */
	/*****************************************************/
	MPI_Allreduce(&l_mins[0], &g_mins[0], 3, MPI_DOUBLE, MPI_MIN, world_comm);
	MPI_Allreduce(&l_maxs[0], &g_maxs[0], 3, MPI_DOUBLE, MPI_MAX, world_comm);
	MPI_Allreduce(&l_n_pts, &g_n_pts, 1, MPI_INT, MPI_SUM, world_comm);
	// Gather the bounding box values for each file
	//printf("[%i] Gather bbox values\n", world_rank);
	double io_time = MPI_Wtime();
	printf("[%i] Communication finished in %f seconds\n", world_rank, 
			io_time - endtime);

	// Create origin for global grid
	//printf("Max(%f, %f, %f), Min(%f, %f, %f)\n", g_maxs[0],g_maxs[1],
	//		g_maxs[2], g_mins[0], g_mins[1], g_mins[2]);
	origin->update(g_mins[0], g_maxs[1], g_mins[2]);
	g_cols = (int) ceil((g_maxs[0] - g_mins[0]) / res);
	g_rows = (int) ceil((g_maxs[1] - g_mins[1]) / res);
	// Create global grid
	g_grid = new Grid(origin, g_cols, g_rows, datatype, res, res);
	//printf("GLOBAL GRID: YRANGE: [ %f, %f]\n", g_grid->origin->x, g_grid->getMinY());

	/********************************************************/
	/*        Create global block scheme                    */
	/********************************************************/
	blk_scheme = new BlockScheme(g_grid, (long)max_buffer, datatype, 
			world_size-1);
	g_n_blks = blk_scheme->cols * blk_scheme->rows; // Global block count
	block_hash = (int*)malloc(sizeof(int) * g_n_blks);
	if (world_rank == 0) 
	{
		for (i = 0; i < g_n_blks; i++) 
		{
			vector<int> block;
			blocks.push_back(block);
		}
	} 
	else 
	{
		l_n_blks = blk_scheme->getBlockCount(world_rank-1); // Local block count
		
		block_file_buffers = (char**)malloc(sizeof(char*)*l_n_blks);
		// Create list of file descriptros to send to main process
		//for (i = 0; i < l_n_files; i++) {
		//	vector<int> blockMap;
		//	blocks.push_back(blockMap);
		//}
		l_blks = blk_scheme->getBlocks(world_rank -1); // Local block id array
		blk_n_files = (int*)malloc(sizeof(int) * l_n_blks);
		// Set the block hash table for quick lookup
		for (i = 0; i < l_n_blks; i++) {
			block_hash[l_blks[i]] = i;
		}
		//printf("[%i] Block Total: %i,Local: %i \n", world_rank, g_n_blks, 
		//		l_n_blks);
	}
	double block_time = MPI_Wtime();
	printf("[%i] Block scheme generation finished in %f seconds\n", world_rank,
			block_time - io_time);
	
	/***********************************************************/
	/*                 Get file list for block                 */
	/***********************************************************/
	if (world_rank != 0) 
	{
		int blk_ids[4] = {0,0,0,0};
		//int tmp_col, tmp_row, tmp_idx = 0;
		int idx_count = 0;
		int buf_msg[6] = {0,0,0,0,0,0}; // {FileIndex, Unique Count, Block Ids
		Grid* tmp_grid = blk_scheme->block_grid;
		for (i = 0; i < file_blk; i++) 
		{
			idx_count = l_file_bbox[i].getGridCells(tmp_grid, blk_ids);
			//printf("[%i] File %i, %s intersects with %i blocks\n", world_rank, 
			//		file_off + i, g_files->fileList[file_off+i],idx_count);
			//printf("[%i]File %i unique %i corners [%i,%i,%i,%i]\n", world_rank, file_off+i,idx_count,blk_ids[0],blk_ids[1], blk_ids[2],blk_ids[3]); 
			buf_msg[0] = file_off +i;
			buf_msg[1] = idx_count;
			buf_msg[2] = blk_ids[0];
			buf_msg[3] = blk_ids[1];
			buf_msg[4] = blk_ids[2];
			buf_msg[5] = blk_ids[3];
			MPI_Isend(&buf_msg[0], 6, MPI_INT, 0, 1, world_comm, &request);
			//printf("[%i] Sending block %i indexes for file %i\n", world_rank, 
			//		buf_msg[1],buf_msg[0]);
			MPI_Wait(&request, &status);
		}
		printf("[%i] Finished sending file metdata\n", world_rank);
	} 
	else 
	{
		int file_flag = 0;
		//int done_flag = 0;
		int file_cnt = 0;
		int idx_cntr =0;
		while (file_cnt < g_n_files) 
		{
			if (!file_flag) 
			{
				MPI_Iprobe(MPI_ANY_SOURCE, 1, world_comm, &file_flag, &status);
			} 
			else 
			{
				//fprintf(stdout, "Receiving block indexes from %i\n", status.MPI_SOURCE);
				mpi_err = MPI_Recv(recv_buf, 6, MPI_INT, status.MPI_SOURCE, 1, world_comm, &status);
				//printf("Adding %i block indexes for file %i\n", recv_buf[1], recv_buf[0]);
				for (j = 0; j < recv_buf[1]; j++) 
				{
					
					idx_cntr = recv_buf[2 +j];
					//printf("[%i] Adding counter to block %i\n", world_rank, idx_cntr);
					blocks[idx_cntr].push_back(recv_buf[0]);
				}
				file_cnt++;
				file_flag =0;
			}
		}
		//printf("[%i] Finished Reading block metadata\n", world_rank);

	}
	double file_comm_time = MPI_Wtime();
	printf("[%i] Block Comm finished in %f seconds\n", world_rank, 
			file_comm_time - block_time);
	MPI_Barrier(world_comm);
	/***********************************************************/
	/*                 Read blocks                             */
	/***********************************************************/
	// Root process deliver block files to processes
	g_blk_cnts = (int*)malloc(sizeof(int) * g_n_blks);
	int max_block = 0;
	if (world_rank == 0) 
	{
		//int block_buff[2] = {0,0}; // {Index of block, number of files}
		//int dest = 0;
		
		for (i = 0; i < g_n_blks; i++)
		{
			//block_buff[0] = i;
			g_blk_cnts[i] = blocks[i].size();
			if (max_block < g_blk_cnts[i]) 
				max_block = g_blk_cnts[i];
			//block_buff[1] = blocks[i].size();
			//dest = blk_scheme->getBlockRank(i) + 1;
			//MPI_Isend(block_buff, 2, MPI_INT, dest, 1, world_comm, &request);
			//printf("[%i] Sent file %i/%i to %i\n", world_rank,
			//		i, g_n_blks, dest);
			//printf("[%i]Sent file count %i for block %i to %i\n", world_rank,
			//		block_buff[1], block_buff[0], dest);
			//MPI_Wait(&request, &status);
		}
	}
	MPI_Barrier(world_comm);
	MPI_Bcast(g_blk_cnts, g_n_blks, MPI_INT, 0, world_comm);
	MPI_Bcast(&max_block, 1, MPI_INT, 0, world_comm);
	/*else
	{
		int block_buff[2] = {0,0};
		printf("[%i]Waiting for %i block counts\n", world_rank, l_n_blks);
		for (i = 0; i < l_n_blks; i++)
		{
			MPI_Recv(block_buff,2, MPI_INT, 0, 1, world_comm, &status);
			printf("[%i] Received %i/%i files\n", world_rank, i, l_n_blks);
			printf("[%i] Block %i has %i files to read\n", world_rank,
				block_buff[0], block_buff[1]);
			j = block_hash[block_buff[0]];
			blk_n_files[j] = block_buff[1];
		}
	}*/
	MPI_Barrier(world_comm);
	printf("[%i] Finished broadcast of block file counts\n", world_rank);
	char* pathBuffer = (char*)malloc(sizeof(char)*1024*max_block);
	if (world_rank == 0) 
	{
		//Should we send 2d char or int array to distribute across processors?
		for (i = 0; i < g_n_blks; i++) {
			//printf("Sending paths for block %i\n", i);
			memset(pathBuffer,0,max_block*1024);
			for (j = 0; j < g_blk_cnts[i]; j++) {
				int filecntr = blocks[i][j];
			//	printf("Reading filename for fileindex %i\n", filecntr);
				//printf("Copying value from %s\n", g_files->fileList[filecntr]);
				int counter= j*1024;
				strncpy(&pathBuffer[counter],g_files->fileList[filecntr], 1024);
			}
			int dest = blk_scheme->getBlockRank(i) +1;
			//printf("Sending pathbuffer to %i\n", dest);
			MPI_Send(pathBuffer, max_block*1024, MPI_CHAR, dest, 1, world_comm);

		}
	} 
	else
	{
		g_files->clear(); // Release memory held by file paths
		printf("[%i] Receiving data for %i blocks\n", world_rank, l_n_blks);
		for (i = 0; i < l_n_blks; i++) 
		{
			block_file_buffers[i] = (char*)malloc(sizeof(char) * 1024*
					g_blk_cnts[l_blks[i]]);
			memset(block_file_buffers[i], 0, 1024*g_blk_cnts[l_blks[i]]);
			MPI_Recv(pathBuffer, max_block*1024, MPI_CHAR, 0, 1, world_comm, &status);
			//int counter = g_blk_cnts[l_blks[i]];
			printf("[%i]Copying paths for block%i  %i files\n", world_rank, l_blks[i], g_blk_cnts[l_blks[i]]);
			for (j = 0; j < g_blk_cnts[l_blks[i]]; j++) 
			{
				int counter = j* 1024;
				strncpy(&block_file_buffers[i][counter], &pathBuffer[counter], 1024);
				printf("[%i]Found buffer %i/%i %s\n",world_rank, j, g_blk_cnts[l_blks[i]], &block_file_buffers[i][counter]);
			}
		}
	}
	free(pathBuffer);
			//	MPI_Send(
	
	double block_count_time = MPI_Wtime();
	printf("[%i] Block Count finished in %f seconds\n", world_rank,
			block_count_time - file_comm_time);
	MPI_Barrier(world_comm);
	/***********************************************************/
	/*                Write                                    */
	/***********************************************************/
	// Global process writes out density raster
	if (world_rank == 0) 
	{
		char dens_tif[1024];
		sprintf(dens_tif, "%s/density.tif", scratch_path);
		blk_scheme->block_grid->alloc();
		for (i = 0; i < g_n_blks; i++)
		{

			// FOr testing only
			//blk_scheme->block_grid->data[i].count = i;
			//blk_scheme->block_grid->data[i].count = blk_scheme->getBlockRank(i);
			blk_scheme->block_grid->data[i].count = g_blk_cnts[i];
			//printf("[%i] Block %i -> %i files\n", world_rank, i, g_blk_cnts[i]);
		}
		blk_scheme->block_grid->write(dens_tif, 3443, 0);
		blk_scheme->block_grid->dealloc();
		sprintf(dens_tif, "%s/proc.tif", scratch_path);
		blk_scheme->block_grid->alloc();
		for ( i =0; i < g_n_blks; i++)
		{
			blk_scheme->block_grid->data[i].count = blk_scheme->getBlockRank(i);
		}
		blk_scheme->block_grid->write(dens_tif, 3443,0);
		blk_scheme->block_grid->dealloc();
		sprintf(dens_tif, "%s/block.vrt", output_path);
		//blk_scheme->buildVRT(dens_tif);

	}
	else 
	{
		int block_idx = 0;
		Grid blk_grid;
		struct Point* p = new Point();
		int n_pts = 0;
		char outPath[1024];
		char path_buf[1024];
		Pixel* pixels = NULL;
		printf("[%i] Allocating grids, size is (%i,%i) %zu\n", world_rank, 
				blk_scheme->b_cols, blk_scheme->b_rows, 
				sizeof(Pixel)*blk_scheme->b_cols*blk_scheme->b_rows);

		for (i = 0; i < l_n_blks; i++) {
			block_idx = l_blks[i];
			//printf("[%i] Generating grid for block %i/%i\n", world_rank, i, 
			//		l_n_blks);
			blk_grid = blk_scheme->getBlockGrid(block_idx, res);
			printf("[%i]Grid successfully created for block %i, %i files to read\n", world_rank, block_idx, g_blk_cnts[block_idx]);
			//printf("Block %i grid dims: (%i)\n", block_idx, blk_grid->cols
			printf("Block grid cols: %i\n", blk_grid.cols);
			
			blk_grid.alloc();
			//printf("[%i]Grid allocated, reading files\n", world_rank);
			// TODO: Debug the block_file_buffer implementation
			int counter = 0;
			int pnt_cntr = 0;
			int pnt_flag = 0;
			for (j = 0; j < g_blk_cnts[block_idx]; j++) {
				counter = j * 1024;
				//printf("[%i]Reading file %i/%i\n",world_rank, j, g_blk_cnts[block_idx]);
				memset(path_buf, 0, 1024);
				strncpy(path_buf, &block_file_buffers[i][counter], 1024);
				printf("[%i] Reading file %i/%i %s\n", world_rank,j, g_blk_cnts[block_idx],
						path_buf);
				// Check if the path_buf is empty
				if (strlen(path_buf) == 0)
				{
					fprintf(stderr, "[%i] Found empty string for blokc %i %i/%i\n", 
							world_rank, block_idx, j, g_blk_cnts[block_idx]);
					continue; 
				} else {
					las.open(path_buf);
					
					n_pts = (int)las.points_count();
					printf("[%i] %i Points to read\n", world_rank, n_pts);
					pnt_cntr = 0;
					pnt_flag = 0;
					for (k=0; k < n_pts; k++) {
						if (las.getReturnNumber(k) == las.getReturnCount(k)) {
							p->update(las.getX(k),las.getY(k),las.getZ(k));
							pnt_flag = blk_grid.set(p);
							if (pnt_flag)
								pnt_cntr++;
							//else {
							//	fprintf(stderr,"[%i]Check pt: (%f,%f), grid:(%f,%f,%f,%f)\n", world_rank,
							//			p->x, p->y, blk_grid.origin.x,blk_grid.getMaxX(),blk_grid.getMinY(),blk_grid.origin.y);
							//	fprintf(stderr,"[%i]Point should be in block %i not %i\n",world_rank, blk_scheme->block_grid->getCellIdx(p->x,p->y), block_idx);
							//}
						}
					}
					las.close();
					printf("Set %i pixel values in block %i\n", pnt_cntr, block_idx);
					//printf("[%i] Finished reading points for file %s\n", world_rank, path_buf);
				}
				printf("[%i]Continuing to file%i/%i for block %i\n", world_rank, 
						j+1, g_blk_cnts[block_idx], block_idx);
			}
			//blk_grid.fillNull(7);
			printf("[%i]Finished reading, beginning write\n",world_rank);
			sprintf(outPath, "%s/%i.tif", output_path, block_idx);
			/*if (i == 1) {
				for (k = 0; k < blk_grid->cols; k++) {
					printf("Block: %i, Cell %i value:%f\n", i,  k, blk_grid->data[k].avg());
				}
			}*/
			blk_grid.write(outPath, 3443, 1);
			printf("[%i]Wrote output grid to %s\n",world_rank, outPath);
			blk_grid.dealloc();
			printf("[%i] Deleting grid\n", world_rank);
			//delete(blk_grid);
		}
		delete(p);
		//if (blk_grid != NULL) {
		//}

	}
	double write_time = MPI_Wtime();
	printf("[%i] Write finished in %f seconds\n", world_rank, 
			write_time - block_count_time);
	MPI_Barrier(world_comm);
	/***********************************************************/
	/*        Clean up                                         */
	/***********************************************************/
	g_files->clear();
	printf("[%i]Cleaning up\n", world_rank);
	if (l_n_blks > 0){
		for (i = 0; i < l_n_blks; i++){
			if (block_file_buffers[i] != NULL)
				free(block_file_buffers[i]);
		}
	}
	if (block_file_buffers != NULL)
		free(block_file_buffers);
	//free(pathBuffer);
	free(g_blk_cnts);
	free(block_hash);
	if (blk_n_files != NULL) 
		free(blk_n_files);
	delete(origin);
	delete(g_files);
	delete[] l_file_bbox;
	delete(g_grid);
	delete(blk_scheme);
	endtime = MPI_Wtime();
	printf("[%i] Process completed in %f seconds\n", world_rank,
			endtime - starttime);
	MPI_Finalize();
	return 0;
}
