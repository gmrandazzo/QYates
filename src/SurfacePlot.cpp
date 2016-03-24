#include "SurfacePlot.h"

Plot::Plot()
{
  setTitle("A Simple SurfacePlot Demonstration");
  
  Polynomial polynomial(this);

  polynomial.setMesh(41,31);
  polynomial.setDomain(-10 ,10,-10, 10);
  polynomial.setMinZ(0);

  polynomial.create();

  setRotation(30,0,15);
  setScale(1,1,1);
  setShift(0.15,0,0);
  setZoom(0.9);

  for (unsigned i=0; i!=coordinates()->axes.size(); ++i)
  {
    coordinates()->axes[i].setMajors(7);
    coordinates()->axes[i].setMinors(4);
  }


  coordinates()->axes[X1].setLabelString("x-axis");
  coordinates()->axes[Y1].setLabelString("y-axis");
  coordinates()->axes[Z1].setLabelString("z-axis");


  setCoordinateStyle(BOX);

  updateData();
  updateGL();
}

int main(int argc, char **argv)
{
    QApplication a(argc, argv);
    Plot plot;
    //a.setMainWidget(&plot);
    plot.resize(800,600);
    plot.show();
    return a.exec();
}

/*
bool Plot::operator()(Plot3D* plot, QString const& fname, QString const& format)
  {     
    FILE* file;
    unsigned int xmesh, ymesh;
    double minx, maxx, miny, maxy;
        
    if ( !collectInfo(file, fname, xmesh, ymesh, minx, maxx, miny, maxy) )
      return false;
        
    /* allocate some space for the mesh */
    double** data = allocateData(xmesh, ymesh);

    for (unsigned int j = 0; j < ymesh; j++) 
    {
      for (unsigned int i = 0; i < xmesh; i++) 
      {
        if (fscanf(file, "%lf", &data[i][j]) != 1) 
        {
          fprintf(stderr, "NativeReader::read: error in data file \"%s\"\n", (const char*)fname.local8Bit());
          return false;
        }

        if (data[i][j] > maxz_)
          data[i][j] = maxz_;
        else if (data[i][j] < minz_)
          data[i][j] = minz_;
      }
    }

    /* close the file */
    fclose(file);

    ((SurfacePlot*)plot)->loadFromData(data, xmesh, ymesh, minx, maxx, miny, maxy);
    deleteData(data,xmesh);

    return true;
  }*/
  