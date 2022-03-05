// MyVTK.cc
// Output VTK visualisation with point data
// Code_CartesianCraft (c) Yongxin Chen

#include "MyVTK.h"

// 0-based number to char converter (from 'a')
char num2char(int n){return (char)('a'+n);}

// Only output single field
void MyVTK::single_field(const Fieldd &f,
                         const Grid   &grid,
                         const MyMPI  &mympi,
                         const char   *fieldname, 
                         const char   *prefix   ,  
                         const char   *filename ) const {

  // Make file name and open file
  ofstream vtr;
  char name[200];
  sprintf(name, "%s%s.x%dx%dx%d.vtr", prefix, filename, mympi.coords(0), mympi.coords(1), mympi.coords(2)); 
  vtr.open(name, ios::binary);

  // Some variables
  int offset = 0;
  char buffer[300];

  // Write head
  vtr << "<VTKFile type=\"RectilinearGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << endl;
  sprintf(buffer, "<RectilinearGrid WholeExtent=\"%d %d %d %d %d %d\">", 
      0, grid.gn(0), 0, grid.gn(1), 0, grid.gn(2));
  vtr << buffer << endl;
  sprintf(buffer, "<Piece Extent=\"%d %d %d %d %d %d\">",		 
                  grid.piece_start(0), grid.piece_end(0),                       
                  grid.piece_start(1), grid.piece_end(1),
                  grid.piece_start(2), grid.piece_end(2));
  vtr << buffer << endl;

  // Write coordinate
  char coords[3] = {'X', 'Y', 'Z'};
  vtr << "<Coordinates>" << endl;
  for(int d=0; d<3; d++){
    vtr << "<DataArray type=\"Float32\" Name=\"" << coords[d]
        << "\" format=\"appended\" offset=\""    << offset << "\" NumberOfComponents=\"1\"/>" << endl;
    offset += sizeof(int)+(grid.n(d)+1)*sizeof(float);
  }
  vtr << "</Coordinates>" << endl;

  // Write point data
  vtr << "<PointData>" << endl;
  vtr << "<DataArray type=\"Float32\" Name=\"" << fieldname <<"\" format=\"appended\" offset=\""
      << offset << "\" NumberOfComponents=\"1\"/>" << endl;
  vtr << "</PointData>" << endl;
  vtr << "</Piece>" << endl;
  vtr << "</RectilinearGrid>" << endl;
  vtr << "<AppendedData encoding=\"raw\">" << endl;
  vtr << "_";

  // Write data
  // Write coordinate
  float out;    // ouptut buffer
  int bsize;    // size in bytes
  for(int d=0; d<3; d++){
    bsize = (grid.n(d)+1)*sizeof(float);
    vtr.write((char*)&bsize, sizeof(int));
    for(int i=grid.start(d); i<grid.end(d)+1; i++){
      out = (float)grid.x(d, i);
      vtr.write((char*)&out, sizeof(float));
    }
  } 

  // Write field
  int size = (f.n()+1).product();      // point data
  bsize = size*sizeof(float);
  double *data = new double[size];
  f.pack(f.start(), f.end()+1, data);  // point data
  vtr.write((char*)&bsize, sizeof(int));
  for(int i=0; i<size; i++){
    out = (float)data[i];
    vtr.write((char*)&out, sizeof(float));
  }

  // Write ending
  vtr << "</AppendedData>" << endl << "</VTKFile>";
  vtr.close();

  // Write .pvtr file
  if(mympi.rank() == 0) {
    sprintf(name, "%s%s.collection.pvtr", prefix, filename);
    ofstream pvtr;
    pvtr.open(name);
    pvtr << "<VTKFile type=\"PRectilinearGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << endl;
    sprintf(buffer, "<PRectilinearGrid WholeExtent=\"%d %d %d %d %d %d\" GhostLevel=\"0\">", \
        0, grid.gn(0), 0, grid.gn(1), 0, grid.gn(2));
    pvtr << buffer << endl;
    pvtr << "<PCoordinates>"  << endl;
    pvtr << "<DataArray type=\"Float32\" Name=\"X\"/>" << endl;
    pvtr << "<DataArray type=\"Float32\" Name=\"Y\"/>" << endl;
    pvtr << "<DataArray type=\"Float32\" Name=\"Z\"/>" << endl;
    pvtr << "</PCoordinates>" << endl;
    pvtr << "<PPointData>"    << endl;
    pvtr << "<DataArray type=\"Float32\" Name=\"" << fieldname << "\" NumberOfComponents=\"1\"/>" << endl;
    pvtr << "</PPointData>"   << endl;

    // Append sources
    Vec3i n, start, end;    // piece dimension, start, end

    // Loop over all the partitions
    for(int k=0; k<mympi.blocks(2); k++){
      for(int j=0; j<mympi.blocks(1); j++){
        for(int i=0; i<mympi.blocks(0); i++){

          // MPI topology coords
          Vec3i ijk(i,j,k); 

          // Determine view in each partition, especially the last one.
          for(int d=0; d<3; d++){
            n(d)      =  grid.gn(d) / mympi.blocks(d);
            start(d)  =  n(d) *  ijk(d);
            end(d)    =  n(d) * (ijk(d)+1);
            if(ijk(d) == mympi.blocks(d)-1){
              n(d)    =  grid.gn(d) - (n(d)*(mympi.blocks(d)-1));
              end(d)  =  start(d) + n(d); 
            }
          }

          // Write file name. The .pvtr file is written in the same path of .vtr files. 
          sprintf(name, "%s.x%dx%dx%d.vtr", filename, ijk(0), ijk(1), ijk(2)); 

          // Write source
          sprintf(buffer, "<Piece Extent=\"%d %d %d %d %d %d\" Source=\"%s\"/>", \
              start(0), end(0), start(1), end(1), start(2), end(2), name);
          pvtr << buffer << endl;
        } 
      }
    } // end each partition

    pvtr << "</PRectilinearGrid>" << endl;
    pvtr << "</VTKFile>";
    pvtr.close();
  }
}

void MyVTK::fluid_fields(const VectorFieldd   &u         ,
                      	 const Fieldd         &p         ,
                         const vector<Fieldd> &other     ,
                         const vector<string> &fieldnames,
                         const Grid           &grid      ,
                         const MyMPI          &mympi     ,
                         const char           *prefix    ,
                         const char           *filename  ) const {
  // Retrieve number of fields
  int nfields = other.size();

  // Make file name and open file
  ofstream vtr;
  char name[200];
  sprintf(name, "%s%s.x%dx%dx%d.vtr", prefix, filename, mympi.coords(0), mympi.coords(1), mympi.coords(2)); 
  vtr.open(name, ios::binary);

  // Some variables
  int offset = 0;
  char buffer[300];

  // Write head
  vtr << "<VTKFile type=\"RectilinearGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << endl;
  sprintf(buffer, "<RectilinearGrid WholeExtent=\"%d %d %d %d %d %d\">", 
      0, grid.gn(0), 0, grid.gn(1), 0, grid.gn(2));
  vtr << buffer << endl;
  sprintf(buffer, "<Piece Extent=\"%d %d %d %d %d %d\">",		 
      grid.piece_start(0), grid.piece_end(0),                       
      grid.piece_start(1), grid.piece_end(1),
      grid.piece_start(2), grid.piece_end(2));
  vtr << buffer << endl;

  // Write coordinate
  char coords[3] = {'X', 'Y', 'Z'};
  vtr << "<Coordinates>" << endl;
  for(int d=0; d<3; d++){
    vtr << "<DataArray type=\"Float32\" Name=\"" << coords[d] << "\" format=\"appended\" offset=\""
        << offset << "\" NumberOfComponents=\"1\"/>" << endl;
    offset += sizeof(int)+(grid.n(d)+1)*sizeof(float);
  }
  vtr << "</Coordinates>" << endl;

  // Write point data
  vtr << "<PointData>" << endl;
  vtr << "<DataArray type=\"Float32\" Name=\"Velocity\" format=\"appended\" offset=\""
      << offset << "\" NumberOfComponents=\"3\"/>" << endl;
  offset += sizeof(int)+3*(u.n()+1).product()*sizeof(float);   // point data
  vtr << "<DataArray type=\"Float32\" Name=\"Pressure\" format=\"appended\" offset=\""
      << offset << "\" NumberOfComponents=\"1\"/>" << endl;
  offset += sizeof(int)+(p.n()+1).product()*sizeof(float);     // point data

  // write other fields if exist
  for(int i=0; i<nfields; i++){
    vtr << "<DataArray type=\"Float32\" Name=\"" << fieldnames[i] << "\" format=\"appended\" offset=\""
        << offset << "\" NumberOfComponents=\"1\"/>" << endl;
    offset += sizeof(int)+(other[i].n()+1).product()*sizeof(float); // point data
  }  
  vtr << "</PointData>" << endl;
  vtr << "</Piece>" << endl;
  vtr << "</RectilinearGrid>" << endl;
  vtr << "<AppendedData encoding=\"raw\">" << endl;
  vtr << "_";

  // Write data
  // Write coordinate
  float out;    // ouptut buffer
  int bsize;    // size in bytes
  for(int d=0; d<3; d++){
    bsize = (grid.n(d)+1)*sizeof(float);
    vtr.write((char*)&bsize, sizeof(int));
    for(int i=grid.start(d); i<grid.end(d)+1; i++){
      out = (float)grid.x(d, i);
      vtr.write((char*)&out, sizeof(float));
    }
  }

  // Write field
  // Velocity
  {
    int size = (u[0].n()+1).product();  // point data
    bsize = 3*size*sizeof(float);
    vtr.write((char*)&bsize, sizeof(int));
    // Point data
    for(int k=p.start(2); k<p.end(2)+1; k++){
      for(int j=p.start(1); j<p.end(1)+1; j++){
        for(int i=p.start(0); i<p.end(0)+1; i++){
          for(int d=0; d<3; d++){
            out = (float)u[d](i,j,k);
            vtr.write((char*)&out, sizeof(float));
          }
        }
      }
    }
  }

  // Pressure
  {
    int size = (p.n()+1).product();   // point data
    bsize = size*sizeof(float);
    double *data = new double[size];
    p.pack(p.start(), p.end()+1, data);  // point data
    vtr.write((char*)&bsize, sizeof(int));
    for(int i=0; i<size; i++){
      out = (float)data[i];
      vtr.write((char*)&out, sizeof(float));
    }
    delete [] data;
  }
  // Other fields
  for(int i=0; i<nfields; i++){
    int size = (other[i].n()+1).product();  // point data
    bsize = size*sizeof(float);
    double *data = new double[size];
    other[i].pack(other[i].start(), other[i].end()+1, data); // point data
    vtr.write((char*)&bsize, sizeof(int));
    for(int i=0; i<size; i++){
      out = (float)data[i];
      vtr.write((char*)&out, sizeof(float));
    }
    delete [] data;
  }

  // Write ending
  vtr << "</AppendedData>" << endl << "</VTKFile>";
  vtr.close();

  // Write .pvtr file
  if(mympi.rank() == 0) {
    sprintf(name, "%s%s.collection.pvtr", prefix, filename);
    ofstream pvtr;
    pvtr.open(name);
    pvtr << "<VTKFile type=\"PRectilinearGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << endl;
    sprintf(buffer, "<PRectilinearGrid WholeExtent=\"%d %d %d %d %d %d\" GhostLevel=\"0\">", \
        0, grid.gn(0), 0, grid.gn(1), 0, grid.gn(2));
    pvtr << buffer << endl;
    pvtr << "<PCoordinates>" << endl;
    pvtr << "<DataArray type=\"Float32\" Name=\"X\"/>" << endl;
    pvtr << "<DataArray type=\"Float32\" Name=\"Y\"/>" << endl;
    pvtr << "<DataArray type=\"Float32\" Name=\"Z\"/>" << endl;
    pvtr << "</PCoordinates>" << endl;
    pvtr << "<PPointData>"    << endl;
    pvtr << "<DataArray type=\"Float32\" Name=\"Velocity\" NumberOfComponents=\"3\"/>" << endl;
    pvtr << "<DataArray type=\"Float32\" Name=\"Pressure\" NumberOfComponents=\"1\"/>" << endl;
    for(int i=0; i<nfields; i++){
      pvtr << "<DataArray type=\"Float32\" Name=\"" << fieldnames[i] << "\" NumberOfComponents=\"1\"/>" << endl;
    }
    pvtr << "</PPointData>" << endl;

    // Append sources
    Vec3i n, start, end;    // piece dimension, start, end

    // Loop over all the partitions
    for(int k=0; k<mympi.blocks(2); k++){
      for(int j=0; j<mympi.blocks(1); j++){
        for(int i=0; i<mympi.blocks(0); i++){

          // MPI topology coords
          Vec3i ijk(i,j,k); 

          // Determine view in each partition, especially the last one.
          for(int d=0; d<3; d++){
            n(d)      =  grid.gn(d) / mympi.blocks(d);
            start(d)  =  n(d) *  ijk(d);
            end(d)    =  n(d) * (ijk(d)+1);
            if(ijk(d) == mympi.blocks(d)-1){
              n(d)    =  grid.gn(d) - (n(d)*(mympi.blocks(d)-1));
              end(d)  =  start(d) + n(d); 
            }
          }

          // Write file name. The .pvtr file is written in the same path of .vtr files. 
          sprintf(name, "%s.x%dx%dx%d.vtr", filename, ijk(0), ijk(1), ijk(2)); 

          // Write source
          sprintf(buffer, "<Piece Extent=\"%d %d %d %d %d %d\" Source=\"%s\"/>", \
              start(0), end(0), start(1), end(1), start(2), end(2), name);
          pvtr << buffer << endl;
        } 
      }
    } // end each partition

    pvtr << "</PRectilinearGrid>" << endl;
    pvtr << "</VTKFile>";
    pvtr.close();
  }
}

void MyVTK::slice_fields(const VectorFieldd &u         ,
                         const Fieldd       &p         ,
                         const Fieldd       *other     ,
                         char              **fieldnames,
                         const int           nfields   ,
                         const Grid         &grid      ,
                         const MyMPI        &mympi     ,
                         const Controller   &controller,
                         const char         *prefix    ,
                         const char         *filename  ) const {

  // Loop over slice normal directions
  char coords[3] = {'X', 'Y', 'Z'};
  for(int norm=0; norm<3; norm++){

    // MPI coordinates
    int mpicoords[2], m = 0;
    for(int d=0; d<3; d++){
      if(d == norm) continue;
      mpicoords[m] = mympi.coords(d);
      m++;
    }

    // Loop over slice
    for(int slice=0; slice<MAX_SLICES; slice++){

      // Check if slice is inside the partition
      if(controller.slice_inside(slice, norm)){

        // Make file name and open file
        ofstream vtr;
        char name[200];
        sprintf(name, "%s%s.%s.%c.%c.x%dx%d.vtr", prefix, filename, "slice", coords[norm], num2char(slice), mpicoords[0], mpicoords[1]);
        vtr.open(name, ios::binary);

        // Some variables
        int offset = 0;
        char buffer[300];

        // Write head
        Vec3i whole_end   = grid.gn()         ; whole_end(norm)   = 0;
        Vec3i piece_start = grid.piece_start(); piece_start(norm) = 0;
        Vec3i piece_end   = grid.piece_end()  ; piece_end(norm)   = 0;
        vtr << "<VTKFile type=\"RectilinearGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << endl;
        sprintf(buffer, "<RectilinearGrid WholeExtent=\"%d %d %d %d %d %d\">", 
            0, whole_end(0), 0, whole_end(1), 0, whole_end(2));
        vtr << buffer << endl;
        sprintf(buffer, "<Piece Extent=\"%d %d %d %d %d %d\">",		 
            piece_start(0), piece_end(0),                       
            piece_start(1), piece_end(1),
            piece_start(2), piece_end(2));
        vtr << buffer << endl;

        // Write coordinate
        vtr << "<Coordinates>" << endl;
        for(int d=0; d<3; d++){
          vtr << "<DataArray type=\"Float32\" Name=\"" << coords[d] <<
            "\" format=\"appended\" offset=\"" << offset << "\" NumberOfComponents=\"1\"/>" << endl;
          offset += sizeof(int)+(piece_end(d)-piece_start(d)+1)*sizeof(float);
        }
        vtr << "</Coordinates>" << endl;

        // Write data
        vtr << "<PointData>" << endl;             // point data
        vtr << "<DataArray type=\"Float32\" Name=\"Velocity\" format=\"appended\" offset=\""
            << offset << "\" NumberOfComponents=\"3\"/>" << endl;
        offset += sizeof(int)+3*(u[0].n()+1).product()/(grid.n(norm)+1)*sizeof(float);  // point data	
        vtr << "<DataArray type=\"Float32\" Name=\"Pressure\" format=\"appended\" offset=\""
            << offset << "\" NumberOfComponents=\"1\"/>" << endl;
        offset += sizeof(int)+(p.n()+1).product()/(grid.n(norm)+1)*sizeof(float);       // point data

        // write other fields if exist
        for(int i=0; i<nfields; i++){
          vtr << "<DataArray type=\"Float32\" Name=\"" << fieldnames[i] << "\" format=\"appended\" offset=\""
              << offset << "\" NumberOfComponents=\"1\"/>" << endl;
          offset += sizeof(int)+(other[i].n()+1).product()/(grid.n(norm)+1)*sizeof(float);  // point data
        }  
        vtr << "</PointData>" << endl;
        vtr << "</Piece>" << endl;
        vtr << "</RectilinearGrid>" << endl;
        vtr << "<AppendedData encoding=\"raw\">" << endl;
        vtr << "_";

        // Write data
        // Write coordinate
        float out;    // ouptut buffer
        int bsize;    // size in bytes
        for(int d=0; d<3; d++){
          bsize = (piece_end(d)-piece_start(d)+1)*sizeof(float);
          vtr.write((char*)&bsize, sizeof(int));
          if(d == norm){
            out = (float)grid.x(d, controller.slice_index(slice, norm));
            vtr.write((char*)&out, sizeof(float));
          }
          else{
            for(int i=grid.start(d); i<grid.end(d)+1; i++){
              out = (float)grid.x(d, i);
              vtr.write((char*)&out, sizeof(float));
            } // end for
          } // end if slice normal direction
        } // end for 3 grid lines

        // Write field
        // Velocity
        {
          int   size  = (u[0].n()+1).product()/(grid.n(norm)+1);  // point data
          Vec3i start = u.start(); start(norm) = controller.slice_index(slice, norm);
          Vec3i end   = u.end()+1  ; end(norm) = start(norm)+1;   // point data
          bsize = 3*size*sizeof(float);
          double **data = new double *[3];
          for(int d=0; d<3; d++) { data[d] = new double [size]; }
          vtr.write((char*)&bsize, sizeof(int));
          // collect 3 components of data
          for(int d=0; d<3; d++){u[d].pack(start, end, data[d]);}
          // output data
          for(int i=0; i<size; i++){
            for(int d=0; d<3; d++){
              out = (float)data[d][i];
              vtr.write((char*)&out, sizeof(float));
            }
          } 
          for(int d=0; d<3; d++) { delete [] data[d]; }
          delete [] data;
        }
        // Pressure
        {
          int   size  = (p.n()+1).product()/(grid.n(norm)+1);     // point data
          Vec3i start = p.start(); start(norm) = controller.slice_index(slice, norm);
          Vec3i end   = p.end()+1; end(norm)   = start(norm)+1;   // point data
          bsize = size*sizeof(float);
          double *data = new double[size];
          p.pack(start, end, data);
          vtr.write((char*)&bsize, sizeof(int));
          for(int i=0; i<size; i++){
            out = (float)data[i];
            vtr.write((char*)&out, sizeof(float));
          }
          delete [] data;
        }
        // Other fields
        for(int i=0; i<nfields; i++){
          int   size  = (other[i].n()+1).product()/(grid.n(norm)+1);      // point data
          Vec3i start = other[i].start(); start(norm) = controller.slice_index(slice, norm);
          Vec3i end   = other[i].end()+1  ; end(norm)   = start(norm)+1;  // point data
          bsize = size*sizeof(float);
          double *data = new double[size];
          other[i].pack(start, end, data);
          vtr.write((char*)&bsize, sizeof(int));
          for(int i=0; i<size; i++){
            out = (float)data[i];
            vtr.write((char*)&out, sizeof(float));
          }
          delete [] data;
        }

        // Write ending
        vtr << "</AppendedData>" << endl << "</VTKFile>";
        vtr.close();

      } // end slice inside partition
    } // end loop over slices
  } // end loop over slice normal directions



  // Write .pvtr file
  for(int norm=0; norm<3; norm++){

    // Loop over slice
    for(int slice=0; slice<MAX_SLICES; slice++){

      // Collective call to check if the slice will be ouptut
      bool with_slice = controller.slice_inside(slice, norm);
      mympi.lor(with_slice);

      if((mympi.rank()==0) && with_slice){

        // Some variables
        int offset = 0;
        char name[200];
        char buffer[300];

        // Write head
        Vec3i whole_end = grid.gn(); whole_end(norm) = 0;

        sprintf(name, "%s%s.%s.%c.%c.collection.pvtr", prefix, filename, "slice", coords[norm], num2char(slice));
        ofstream pvtr;
        pvtr.open(name);
        pvtr << "<VTKFile type=\"PRectilinearGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << endl;
        sprintf(buffer, "<PRectilinearGrid WholeExtent=\"%d %d %d %d %d %d\" GhostLevel=\"0\">", 
            0, whole_end(0), 0, whole_end(1), 0, whole_end(2));
        pvtr << buffer << endl;
        pvtr << "<PCoordinates>" << endl;
        pvtr << "<DataArray type=\"Float32\" Name=\"X\"/>" << endl;
        pvtr << "<DataArray type=\"Float32\" Name=\"Y\"/>" << endl;
        pvtr << "<DataArray type=\"Float32\" Name=\"Z\"/>" << endl;
        pvtr << "</PCoordinates>" << endl;
        pvtr << "<PPointData>"    << endl;
        pvtr << "<DataArray type=\"Float32\" Name=\"Velocity\" NumberOfComponents=\"3\"/>" << endl;
        pvtr << "<DataArray type=\"Float32\" Name=\"Pressure\" NumberOfComponents=\"1\"/>" << endl;
        for(int i=0; i<nfields; i++){
          pvtr << "<DataArray type=\"Float32\" Name=\"" << fieldnames[i] << "\" NumberOfComponents=\"1\"/>" << endl;
        }
        pvtr << "</PPointData>" << endl;

        // Append sources
        Vec3i n, start, end;                                             // piece dimension, start, end
        int ijke[3]={mympi.blocks(0), mympi.blocks(1), mympi.blocks(2)}; // MPI coords
        int ie, je, ke;
        ijke[norm] = 1;                                                  // flatten MPI coords
        ie = ijke[0]; je = ijke[1]; ke = ijke[2];                        // get back  

        // Loop over all the partitions
        for(int k=0; k<ke; k++)
          for(int j=0; j<je; j++)
            for(int i=0; i<ie; i++){

              // MPI topology coords
              Vec3i ijk(i,j,k); 

              // Determine view in each partition, especially the last one.
              for(int d=0; d<3; d++){
                n(d)      =  grid.gn(d) / mympi.blocks(d);
                start(d)  =  n(d) *  ijk(d);
                end(d)    =  n(d) * (ijk(d)+1);
                if(ijk(d) == mympi.blocks(d)-1){
                  n(d)    =  grid.gn(d) - (n(d)*(mympi.blocks(d)-1));
                  end(d)  =  start(d) + n(d); 
                }
              }
              start(norm) = 0;     // flatten in slice normal direction 
              end  (norm) = 0;     // flatten in slice normal direction 

              // MPI coordinates
              int mpicoords[2], m = 0;
              for(int d=0; d<3; d++){
                if(d == norm) { continue; }
                mpicoords[m] = ijk(d);
                m++;
              }

              // Write file name. The .pvtr file is written in the same path of .vtr files. 
              sprintf(name, "%s.%s.%c.%c.x%dx%d.vtr", filename, "slice", coords[norm], num2char(slice), mpicoords[0], mpicoords[1]); 

              // Write source
              sprintf(buffer, "<Piece Extent=\"%d %d %d %d %d %d\" Source=\"%s\"/>", \
                  start(0), end(0), start(1), end(1), start(2), end(2), name);
              pvtr << buffer << endl;
            } // end each partition

        pvtr << "</PRectilinearGrid>" << endl << "</VTKFile>";
        pvtr.close();

      } // end if rank==0 output
    } // end loop over slices
  } // end looop over slice normal directions
}
