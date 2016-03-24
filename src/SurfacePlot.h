#ifndef SURFACEPLOT_H
#define SURFACEPLOT_H

#include <cmath>
// #include <QApplication>
#include <qwtplot3d/qwt3d_surfaceplot.h>
#include <qwtplot3d//qwt3d_function.h>


using namespace Qwt3D;

class Polynomial : public Function
{
public:

  Polynomial(SurfacePlot* pw)
  :Function(pw)
  {
  }
  
  

  double operator()(double x, double z )
  {
    return ;
  }
};


class Plot : public SurfacePlot
{
public:
    Plot();
};


#endif