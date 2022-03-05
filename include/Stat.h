// Stat.h
// Class holds time average statistics operations
// Code_CartesianCraft (c) Yongxin Chen

#ifndef STAT_H
#define STAT_H

#include "Grid.h"
#include "H5IO.h"
#include "Field.h"
#include "MyMPI.h"
#include "MyVTK.h"
#include "parameter.h"
#include "Controller.h"
#include "VectorField.h"

#include <vector>
#include <string>
#include <fstream>

using std::string;
using std::vector;
using std::ofstream;
using std::ifstream;
using std::ios;

class Stat{

  private:
    VectorFieldd u;           // velocity field
    Fieldd       p;           // pressure
    Fieldd       uiui[3];     // Reynolds normal stress
    Fieldd       uiuj[3];     // Reynolds shear stress
    int          samples=0;   // number of samples
    MyVTK        myvtk;       // Visualisation object

  public:
    Stat(){}
    ~Stat(){}

    void init(const Grid &);                                  // initialise an object
    void do_statistics(const VectorFieldd &, const Fieldd &); // do statistics
    void write(const Grid &, const MyMPI &) const;            // write time average data
    void save(const MyMPI &) const;                           // save for restart
    void read(const MyMPI &);                                 // read for restart
    void save_hdf(const MyMPI &, const Grid &) const;         // save with hdf5 for restart
    void read_hdf(const MyMPI &, const Grid &);               // read with hdf5 for restart
};

void mean1(Fieldd &mf, const Fieldd &f, int n);
void mean2(Fieldd &mf, const Fieldd &f1, const Fieldd &f2, int n);

#endif // STAT_H
