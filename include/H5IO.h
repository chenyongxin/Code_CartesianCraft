// H5IO.h
// Write/Read field data to/from a HDF5 file collectively.
// Code_CartesianCraft (c) Yongxin Chen

#ifndef H5IO_H
#define H5IO_H

#ifdef USE_HDF

#include <mpi.h>
#include <vector>
#include <string>
#include <cstdio>
#include <fstream>
#include <iostream>
#include "hdf5.h"
#include "Grid.h"
#include "MyMPI.h"
#include "Field.h"
#include "VectorField.h"

using std::endl;
using std::vector;
using std::string;

// Write double precision field
void H5write_field(hid_t file_id, 
                   const char *dataset_name,
                   const MyMPI &mympi, 
                   const Grid &grid, 
                   const Fieldd &f);

// Read double precision field
void H5read_field(hid_t file_id,
                  const char *dataset_name,
                  const MyMPI &mympi,
                  const Grid &grid,
                  Fieldd &f);

// Write single precision fluid fields
void H5fluid_fields(const VectorFieldd   &u                 ,
                    const Fieldd         &p                 ,
                    const vector<Fieldd> &other             ,
                    const vector<string> &fieldnames        ,
                    const Grid           &grid              ,
                    const MyMPI          &mympi             ,
                    const char           *prefix   = "./"   ,
                    const char           *filename = "Step" );

#endif // USE_HDF

#endif // H5IO_H
