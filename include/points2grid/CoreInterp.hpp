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

#pragma once

#include <points2grid/export.hpp>

class P2G_DLL CoreInterp
{
public:
    CoreInterp() {};
    virtual ~CoreInterp() {};

    virtual int init() = 0;
    virtual int update(double data_x, double data_y, double data_z) = 0;
    virtual int finish(char *outputName, int outputFormat, unsigned int outputType) = 0;
    // mpi adds - start
    virtual void updateGridPointSend(int target_rank, int x, int y, double data_z, double distance){};
    virtual void updateGridPointRecv(){};

    //int comm_done;

    virtual int getIsReader(){return 1;};
    virtual int getIsWriter(){return 1;};
    virtual int* getReaders(){return 0;};
    virtual int* getWriters(){return 0;};
    virtual int getReaderCount(){return 1;};
    virtual int getWriterCount(){return 1;};
    virtual int*
    getReadDone ()
    {
        return 0;
    }
    int
    get_epsg_code () const
    {
        return epsg_code;
    }

    void
    set_epsg_code (int epsgCode)
    {
        epsg_code = epsgCode;
    }

    int
    get_bigtiff () const
    {
        return bigtiff;
    }

    void
    set_bigtiff (int bigtiff)
    {
        this->bigtiff = bigtiff;
    }

    // mpi adds - finish
protected:
    double GRID_DIST_X;
    double GRID_DIST_Y;

    int GRID_SIZE_X; // total size of a grid
    int GRID_SIZE_Y; //

    // for outputting
    double min_x;
    double max_x;
    double min_y;
    double max_y;

    // for DEM filling
    int window_size;

    // for file output
    int bigtiff;
    int epsg_code;


};

