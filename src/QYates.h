#ifndef QYates_H
#define QYates_H

#include <QWidget>
#include <QString>
#include <QList>
#include "ui_QYates.h"

#include <scientific.h>

#ifndef WIN32
  typedef unsigned int uint;
#endif

class QYates : public QWidget
{
Q_OBJECT
public:
    QYates();
    virtual ~QYates();
    
private slots:
  void Open();
  void Close();
  void Calculate();
  void cleanWall();
  void HideShowVariable();

private:
  void ReadFile();
    
  double IndVal2PredVal(uint col, double val);
  void GenSignMatrix();

  void GetMaxVal(QList<double> b, 
                double x1_from, double x1_to, uint nx1, 
                double x2_from, double x2_to, uint nx2,
                double x3_from, double x3_to, uint nx3,
                double *y_m, double *x1_m, double *x2_m, double *x3_m);
  
  void PolyGetPlotX1Const(QList< double > b, double x2_from, double x2_to, double x3_from, double x3_to, double k, uint step );
  void PolyGetPlotX2Const(QList< double > b, double x1_from, double x1_to, double x3_from, double x3_to, double k, uint step );
  void PolyGetPlotX3Const(QList< double > b, double x1_from, double x1_to, double x2_from, double x2_to, double k, uint step );
  
//   void FindMaxOnRespSurf(QList<double> b, double x1_min, double x1_max, double x2_min, double x2_max, double x3_min, double x3_max, uint step);
  QList<double> RunSurfExploration();
  void ResponseSurfExploration();
  
  double PlusSummer(uint col_);
  double MinusSummer(uint col_);
  void Yates();
  
  QString file_;
  
  matrix *x, *signmx, *levels;
  dvector *y;
  
  Ui::QYates ui;
};

#endif // QYates_H
