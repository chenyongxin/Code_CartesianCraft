// MyVTK.h
// Output VTK visualisation with point data
// Code_CartesianCraft (c) Yongxin Chen

#ifndef MYVTK_H
#define MYVTK_H

#include <cstdio>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>

#include "Vec.h"
#include "Grid.h"
#include "Field.h"
#include "parameter.h"
#include "Controller.h"
#include "VectorField.h"

using std::vector;
using std::endl;
using std::ios;
using std::ofstream;
using std::string;

class MyVTK{

  public:

    MyVTK(){}
    ~MyVTK(){}

    // Output a single field with a path prefix
    // Path prefix ends with /.
    void single_field(const Fieldd       &                   ,
                      const Grid         &                   ,
                      const MyMPI        &                   ,
                      const char         *fieldname = "Data" ,
                      const char         *prefix    = "./"   ,
                      const char         *filename  = "fluid") const;

    // Output fluid fields with default velocity and pressure, and optional
    // other fields.
    // Path prefix ends with /.
    // Filename is an integer that indicates a contiguous time step.
    // All the filenames should keep the same digits.
    void fluid_fields(const VectorFieldd   &u                 ,
                      const Fieldd         &p                 ,
                      const vector<Fieldd> &other             ,
                      const vector<string> &fieldnames        ,
                      const Grid           &grid              ,
                      const MyMPI          &mympi             ,
                      const char           *prefix   = "./"   ,
                      const char           *filename = "Step" ) const;

    // Output a slice.
    // Path prefix ends with /.
    // Filename is an integer that indicates a contiguous time step.
    // All the filenames should keep the same digits.
    void slice_fields(const VectorFieldd &u                 ,
                      const Fieldd       &p                 ,
                      const Fieldd       *other             ,
                      char               **fieldnames       ,
                      const int          nfields            ,
                      const Grid         &grid              ,
                      const MyMPI        &mympi             ,
                      const Controller   &controller        ,
                      const char         *prefix   = "./"   ,
                      const char         *filename = "Step" ) const;
};

#endif // MYVTK_H
