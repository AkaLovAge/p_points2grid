/*
*
COPYRIGHT AND LICENSE

Copyright (c) 2011 The Regents of the University of California.
All rights reserved.

Redistribution and use in source and binary forms, with or
without modification, are permitted provided that the following
conditions are met:

1. Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above
copyright notice, this list of conditions and the following
disclaimer in the documentation and/or other materials provided
with the distribution.

3. All advertising materials mentioning features or use of this
software must display the following acknowledgement: This product
includes software developed by the San Diego Supercomputer Center.

4. Neither the names of the Centers nor the names of the contributors
may be used to endorse or promote products derived from this
software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE REGENTS AND CONTRIBUTORS ``AS IS''
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE REGENTS
OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
*
*
* Based on the notes by Prof. Ramon Arrowsmith(ramon.arrowsmith@asu.edu)
* Authors: Han S Kim (hskim@cs.ucsd.edu), Sriram Krishnan (sriram@sdsc.edu)
*
*/

#include <points2grid/config.h>
#include <points2grid/Interpolation.hpp>
#include <points2grid/Global.hpp>

#include <string.h>
#include <math.h>
#include <float.h>
#include <stddef.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>

#include <points2grid/lasfile.hpp>

#include <boost/scoped_ptr.hpp>


#include <fstream>  // std::ifstream
#include <iostream> // std::cerr
#include <string.h>

/////////////////////////////////////////////////////////////
// Public Methods
/////////////////////////////////////////////////////////////

Interpolation::Interpolation(double x_dist, double y_dist, double radius,
                             int _window_size, int _interpolation_mode = INTERP_AUTO,
                             int _rank = 0, int _process_count = 1) : GRID_DIST_X (x_dist), GRID_DIST_Y(y_dist)
{
    rank = _rank;
    process_count = _process_count;
    data_count = 0;
    radius_sqr = radius * radius;
    window_size = _window_size;
    interpolation_mode = _interpolation_mode;

    min_x = DBL_MAX;
    min_y = DBL_MAX;

    max_x = -DBL_MAX;
    max_y = -DBL_MAX;
    //printf ("mode %i\n", interpolation_mode);
}

Interpolation::~Interpolation()
{
    delete interp;
}

int Interpolation::init(char *inputName, int inputFormat)
{


    //unsigned int i;

    //struct tms tbuf;
    clock_t t0, t1;


    //////////////////////////////////////////////////////////////////////
    // MIN/MAX SEARCHING
    // TODO: the size of data, min, max of each coordinate
    //       are required to implement streaming processing....
    //
    // This code can be eliminated if database can provide these values
    //////////////////////////////////////////////////////////////////////

    //t0 = times(&tbuf);
    t0 = clock();

    if(inputName == NULL)
    {
        cerr << "Wrong Input File Name" << endl;
        return -1;
    }

    //printf("inputName: '%s'\n", inputName);

    if (inputFormat == INPUT_ASCII) {
        FILE *fp;
        char line[1024];
        double data_x, data_y;
        //double data_z;

        if((fp = fopen(inputName, "r")) == NULL)
        {
            cerr << "file open error" << endl;
            return -1;
        }

        // throw the first line away - it contains the header
        fgets(line, sizeof(line), fp);

        // read the data points to find min and max values
        while(fgets(line, sizeof(line), fp) != NULL)
        {
            data_x = atof(strtok(line, ",\n"));
            if(min_x > data_x) min_x = data_x;
            if(max_x < data_x) max_x = data_x;

            data_y = atof(strtok(NULL, ",\n"));
            if(min_y > data_y) min_y = data_y;
            if(max_y < data_y) max_y = data_y;

            data_count++;

        }

        fclose(fp);
    } else { // las input

        las_file las;
        las.open(inputName);

        min_x = las.minimums()[0];
        min_y = las.minimums()[1];
        max_x = las.maximums()[0];
        max_y = las.maximums()[1];

        data_count = las.points_count();
        
        las.close();

    }

    t1 = clock();
    //printf("Min/Max searching time: %10.2f\n", (double)(t1 - t0)/CLOCKS_PER_SEC);


    //t0 = times(&tbuf);

    //////////////////////////////////////////////////////////////////////
    // Intialization Step excluding min/max searching
    //////////////////////////////////////////////////////////////////////

    //cerr << "min_x: " << min_x << ", max_x: " << max_x << ", min_y: " << min_y << ", max_y: " << max_y << endl;

    GRID_SIZE_X = (int)(ceil((max_x - min_x)/GRID_DIST_X)) + 1;
    GRID_SIZE_Y = (int)(ceil((max_y - min_y)/GRID_DIST_Y)) + 1;

    cerr << "GRID_SIZE_X " << GRID_SIZE_X << endl;
    cerr << "GRID_SIZE_Y " << GRID_SIZE_Y << endl;

    if (interpolation_mode == INTERP_AUTO) {
        // if the size is too big to fit in memory,
        // then construct out-of-core structure
        if(GRID_SIZE_X * GRID_SIZE_Y > MEM_LIMIT) {
            interpolation_mode= INTERP_OUTCORE;
        } else {
            interpolation_mode = INTERP_INCORE;
        }
    }

    if (interpolation_mode == INTERP_OUTCORE) {
        cerr << "Using out of core interp code" << endl;;

        interp = new OutCoreInterp(GRID_DIST_X, GRID_DIST_Y, GRID_SIZE_X, GRID_SIZE_Y, radius_sqr, min_x, max_x, min_y, max_y, window_size);
        if(interp == NULL)
        {
            cerr << "OutCoreInterp construction error" << endl;
            return -1;
        }

        cerr << "Interpolation uses out-of-core algorithm" << endl;

    } else if (interpolation_mode == INTERP_MPI){
        cerr << "Using mpi interp code" << endl;

        interp = new MpiInterp(GRID_DIST_X, GRID_DIST_Y, GRID_SIZE_X, GRID_SIZE_Y, radius_sqr, min_x, max_x, min_y, max_y, window_size, rank, process_count);

        cerr << "Interpolation uses mpi algorithm" << endl;
    } else {
        cerr << "Using incore interp code" << endl;

        interp = new InCoreInterp(GRID_DIST_X, GRID_DIST_Y, GRID_SIZE_X, GRID_SIZE_Y, radius_sqr, min_x, max_x, min_y, max_y, window_size);

        cerr << "Interpolation uses in-core algorithm" << endl;
    }
    if(interp->init() < 0)
    {
        cerr << "inter->init() error" << endl;
        return -1;
    }

    //cerr << "Interpolation::init() done successfully" << endl;

    return 0;
}

int Interpolation::interpolation(char *inputName,
                                 char *outputName,
                                 int inputFormat,
                                 int outputFormat,
                                 unsigned int outputType)
{
    int rc;
    //unsigned int i;
    double data_x, data_y;
    double data_z;

    //struct tms tbuf;
    //clock_t t0, t1;

    printf("Interpolation Starts, rank %i\n", rank);

    //t0 = times(&tbuf);

    //cerr << "data_count: " << data_count << endl;

    /*
      if((rc = interp->init()) < 0)
      {
      cerr << "inter->init() error" << endl;
      return -1;
      }
    */

    if (inputFormat == INPUT_ASCII) {
        FILE *fp;
        char line[1024];

        if((fp = fopen(inputName, "r")) == NULL)
        {
            printf("file open error\n");
            return -1;
        }

        // throw the first line away - it contains the header
        fgets(line, sizeof(line), fp);

        // read every point and generate DEM
        while(fgets(line, sizeof(line), fp) != NULL)
        {
            data_x = atof(strtok(line, ",\n"));
            data_y = atof(strtok(NULL, ",\n"));
            data_z = atof(strtok(NULL, ",\n"));

            data_x -= min_x;
            data_y -= min_y;

            //if((rc = interp->update(arrX[i], arrY[i], arrZ[i])) < 0)
            if((rc = interp->update(data_x, data_y, data_z)) < 0)
            {
                cerr << "interp->update() error while processing " << endl;
                return -1;
            }
        }

        fclose(fp);
    } 

    else
    { // input format is LAS

        
        if (interp->getIsReader())
        {
            las_file las;
            las.open (inputName);

            size_t count = las.points_count ();
            size_t index (0);
            while (index < count)
            {
                data_x = las.getX (index);
                data_y = las.getY (index);
                data_z = las.getZ (index);

                data_x -= min_x;
                data_y -= min_y;
                //
                //cerr << "calling update rank "<< rank << endl;
                if ((rc = interp->update (data_x, data_y, data_z)) < 0)
                {
                    cerr << "interp->update() error while processing " << endl;
                    return -1;
                }
                index++;
            }
            interp->comm_done = 1;
            if (process_count > 1)
            {
                int i;
                for (i = 1; i < process_count; i++)
                {
                    if(interp->getWriters()[i])
                    {
                        interp->updateGridPointSend (i, 1, 1, 1, 1);
                    }
                }
            }
        }
        else if (interp->getIsWriter())// rank is a writer
        {
            interp->updateGridPointRecv();
        }
        

    }

    if((rc = interp->finish(outputName, outputFormat, outputType)) < 0)
    {
        cerr << "interp->finish() error" << endl;
        return -1;
    }

    cerr << "Interpolation::interpolation() done successfully, rank " << rank << endl;

    return 0;
}

void Interpolation::setRadius(double r)
{
    radius_sqr = r * r;
}

unsigned int Interpolation::getDataCount()
{
    return data_count;
}

unsigned int Interpolation::getGridSizeX()
{
    return GRID_SIZE_X;
}

unsigned int Interpolation::getGridSizeY()
{
    return GRID_SIZE_Y;
}

