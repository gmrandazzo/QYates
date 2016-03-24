
#include "QYates.h"
#include <QFile>
#include <QFileDialog>
#include <QList>
#include <QVector>
#include <QMessageBox>
#include <QPushButton>
#include <QTextStream>
#include <QString>
#include <QStringList>

#include <QtConcurrentRun>

#include <cmath>

// #include <QDebug>

#include <scientific.h>


void QYates::HideShowVariable()
{
  if(ui.comboBox->currentIndex() == 0 )
    ui.dipvar_log10->show();
  else if(ui.comboBox->currentIndex() == 1 )
    ui.dipvar_log10->hide();
}

void QYates::cleanWall()
{
  ui.log_results->clear();
}

void QYates::Calculate()
{
  ReadFile();
  if(x->row == 0 || x->col == 0 || y->size != x->row)
    return;
  
  if(ui.comboBox->currentIndex() == 0 )
    ResponseSurfExploration();
  else if(ui.comboBox->currentIndex() == 1 )
    Yates();
}

void QYates::Close()
{
  qApp->exit();
}

void QYates::Open()
{
  #ifdef WIN32
  file_ = QFileDialog::getOpenFileName( this, tr("Open File"), QDir::currentPath(), tr("Text Files (*.txt);;All files (*.*)"));
  #else
  file_ = QFileDialog::getOpenFileName( this, tr("Open File"), QDir::currentPath(), tr("Text Files (*.txt);;All files (*.*)"), 0, QFileDialog::DontUseNativeDialog );
  #endif

  ui.filecontrol->setText(file_);

}

void QYates::ReadFile()
{
  QFile file(file_);
  if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
    return;
  QTextStream in(&file);
  QList<QStringList> item;
  while (!in.atEnd()) {
    item.append(in.readLine().split(QRegExp("\\s+"), QString::SkipEmptyParts));
  }
  
  // alloc the memory of data struct
  ResizeMatrix(&x, item.size(), item.last().size()-1);// -1 because the last column is the dipendent variable column
  DVectorResize(&y, item.size());
  
  for(uint i = 0; i < (uint)item.size(); i++) {
    for (uint j = 0; j < (uint)item[i].size()-1; j++ ){
      setMatrixValue(x, i, j, item[i][j].toDouble());
    }
    setDVectorValue(y, i, item[i].last().toDouble()); // the last column is the column of dipendent variable
  }
  
  GenSignMatrix();
}

double QYates::PlusSummer(uint col)
{
  double sum = 0.f;
  int n = 0;
  for (uint i = 0; i < x->row; i++){
    if (getMatrixValue(signmx, i, col) > 0){
      sum += getDVectorValue(y, i) * getMatrixValue(signmx, i, col);
      n++;
    }
  }
  if(n != 0){
    sum /= n;
    return sum;
  }
  else
    return 0;
}

double QYates::MinusSummer(uint col)
{
  double sum = 0.f;
  int n = 0;
  for (uint i = 0; i < x->row; i++){
    if(getMatrixValue(signmx, i, col) < 0){
      sum += getDVectorValue(y, i) * getMatrixValue(signmx, i, col);
      n++;
    }
  }
  if(n != 0){
    sum /= n;
    return sum;
  }
  else
    return 0;
}

// Convert the indipendent value from a number to a Real Predictor Value.
double QYates::IndVal2PredVal(uint col, double val )
{
  // val * sep + centervalue
  return (val * getMatrixValue(levels, col, 3) + getMatrixValue(levels, col, 2));
}

double dist(double d1, double d2){
  return ceil(sqrt((d1-d2)*(d1-d2)));
}

void QYates::GenSignMatrix()
{
  if(ui.comboBox->currentIndex() == 0){
    ResizeMatrix(&signmx, x->row, x->col*x->col+1); // +1 because we need to explorate the surface and we need the composite cetral design and the first is to constant sign
    
    ResizeMatrix(&levels, x->col, 4); // 4 becaus e 0 = min; 1 = max ; 2 = center 3 = semirange
    // fill the zero column with 1
    for(uint i = 0; i < signmx->row; i++){
      setMatrixValue(signmx, i, 0, 1.0);
    }
    
    for(uint j = 0; j < x->col; j++){
      double max, min, near_center, center, semirange;
      
      // research the center value
      MatrixColumnMinMax(x, j, &min, &max);

//       qDebug() << "min: " << min << " max " << max;
      
      center = (min + ((max - min) / 2));
      
      near_center = 0.f;
      uint c = 0;
      for(uint i = 0; i < x->row; i++){
        if(getMatrixValue(x, i, j) > center 
          && getMatrixValue(x, i, j) < max){
          near_center += getMatrixValue(x, i, j);
          c++;
        }
        else{
          continue;
        }
      }
      
      near_center /= c;
//       qDebug() << "Near Center "<< near_center;
      semirange =  near_center - center;
      
      
//       qDebug()<< "center: " << center << "semirange " << semirange;
      
      for(uint i = 0; i < signmx->row; i++){
        setMatrixValue(signmx, i, j+1, (getMatrixValue(x, i, j) - center) / semirange);
      }
      
      setMatrixValue(levels, j, 0, min);
      setMatrixValue(levels, j, 1, max);
      setMatrixValue(levels, j, 2, center);
      setMatrixValue(levels, j, 3, semirange);
    }

    //now create the supplementar column that are maded from the product of generator columns
    for(uint j = x->col+1; j < (x->col+x->col)+1; j++){
      for(uint i = 0; i < x->row; i++){ /*x_j^2*/
        setMatrixValue(signmx, i, j, square(getMatrixValue(signmx, i, j-x->col)));
      }
    }
    
    /* x_ij * x_ik with j != k*/
    for(uint j = 0; j < x->col; j++){
      for(uint k = j+1; k < x->col; k++){
        for(uint i = 0; i < x->row; i++){
          setMatrixValue(signmx, i, j+k+x->col+x->col, getMatrixValue(signmx, i, j+1)*getMatrixValue(signmx, i, k+1));
        }
      }
    }
  }
  else{
    ResizeMatrix(&signmx, x->row, x->col); // +1 because we need to explorate the surface and we need the composite cetral design and the first is to constant sign
    
    ResizeMatrix(&levels, x->col, 4); // 4 becaus e 0 = min; 1 = max ; 2 = center 3 = semirange
    
    for(uint j = 0; j < x->col; j++){
      double max, min, center, near_center, semirange;
      // research the center value
      MatrixColumnMinMax(x, j, &min, &max);

      qDebug() << "min: " << min << " max " << max;
      center = (min + ((max - min) / 2));
      
      near_center = 0.f;
      uint c = 0;
      for(uint i = 0; i < x->row; i++){
        if(getMatrixValue(x, i, j) > center 
          && getMatrixValue(x, i, j) < max){
          near_center += getMatrixValue(x, i, j);
          c++;
        }
        else{
          continue;
        }
      }
      
      near_center /= c;
//       qDebug() << "Near Center "<< near_center;
      semirange =  near_center - center;
      
//       qDebug()<< "center: " << center << "semirange " << semirange;
   
      for(uint i = 0; i < signmx->row; i++){
        setMatrixValue(signmx, i, j, (getMatrixValue(x, i, j) - center) / semirange);
      }
      
      setMatrixValue(levels, j, 0, min);
      setMatrixValue(levels, j, 1, max);
      setMatrixValue(levels, j, 2, center);
      setMatrixValue(levels, j, 3, semirange);
    }
  }
  PrintMatrix(signmx);
}

/*
 * In order to find the coefficient value "b" of the polynomial equation, the procedure applied is this:
 * 
 * b = (Z' * Z')^-1 * Z'* y
 * 
 * were :
 * b = vector of coefficients (b0, b1, b2, b...., bn)
 * Z' = is the sign matrix generated for the Centro Simmetric CCD Design
 * y = Dependent value vector. For each experiment we have a dependent value.
 * 
 * The standard errors of the varuous estimates are the square roots of the corresponding diagoknal elementso of (Z' * Z')^-1 σ^2
 * 
 * For info check the P.Box - Empirical Model-Building and response surfaces  Pag. 304-316
 */

QList<double> QYates::RunSurfExploration()
{
  dvector *b;
  
  initDVector(&b);
  
  if(ui.dipvar_log10->isChecked()){
    dvector *y_;
    NewDVector(&y_, y->size);
    for(uint i = 0; i < y->size; i++){
      setDVectorValue(y_, i, log10(getDVectorValue(y, i)+1)); // +1 in order to prevent the 0 value
    }
    
    OrdinaryLeastSquares(signmx, y_, b);
    
    DelDVector(&y_);
  }
  else{
    OrdinaryLeastSquares(signmx, y, b);
  }
  
  QList<double> bcoeff;
  for(uint i = 0; i < b->size; i++){
    bcoeff.append(b->data[i]);
  }
  
  DelDVector(&b);
  
  return bcoeff;
}

void QYates::GetMaxVal(QList<double> b, 
                       double x1_from, double x1_to, uint nx1, 
                       double x2_from, double x2_to, uint nx2,
                       double x3_from, double x3_to, uint nx3,
                       double *y_m, double *x1_m, double *x2_m, double *x3_m)
{
  // y = b0 + b1x1 + b2x2 + b3x3 + b11x1^2 + b22x2^2 + b33x3^2 + b12x1x2 + b13x1x3 + b23x2x3  
  double x1_range = floor(x1_to - x1_from);
  double x2_range = floor(x2_to - x2_from);
  double x3_range = floor(x3_to - x3_from);
  
  double dx1 = x1_range / nx1;
  double dx2 = x2_range / nx2;
  double dx3 = x3_range / nx3;  
  
  double x1 = dx1;
  for(uint i = 0; i < nx1; i++){
    double x2 = dx2;
    for(uint j = 0; j < nx2; j++){
      double x3 = dx3;
      for(uint k = 0; k < nx3; k++){
        
        double y = b[0] + b[1]*x1 + b[2]*x2 + b[3]*x3 + b[4]*x1*x1 + b[5]*x2*x2 + b[6]*x3*x3
          + b[7]*x1*x2 + b[8]*x1*x3 + b[9]*x2*x3;
        
        if(i == 0 && j == 0 && k == 0){
          (*x2_m) = x2;
          (*x1_m) = x1;
          (*x3_m) = x3;
          (*y_m) = y;
        }
        else{
          if(y > (*y_m)){
            (*x2_m) = x2;
            (*x1_m) = x1;
            (*x3_m) = x3;
            (*y_m) = y;
          }
        }
        
        x3 += dx3;
      }
      x2 += dx2;
    }
    x1 += dx1;
  }
}

void QYates::PolyGetPlotX1Const(QList< double > b, double x2_from, double x2_to, double x3_from, double x3_to, double k, uint step)
{
  QFile file("polysurface_x1_const.txt");
  if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
    return;

  QTextStream out(&file);
  
  
  // y = b0 + b1x1 + b2x2 + b3x3 + b11x1^2 + b22x2^2 + b33x3^2 + b12x1x2 + b13x1x3 + b23x2x3  
  double x2_range = floor(x2_to - x2_from);
  double x3_range = floor(x3_to - x3_from);
  
  double dx2 = x2_range / step;
  double dx3 = x3_range / step;
  
  double incx2 = x2_from;
  for(uint i = 0; i < step; i++){
    double incx3 = x3_from;
    for(uint j = 0; j < step; j++){
      double y = b[0] + 0.5*(b[1]*k + b[2]*incx2 + b[3]*incx3);
/*      double y = b[0] + b[1]*k + b[2]*incx2 + b[3]*incx3 + b[4]*k*k + b[5]*incx2*incx2 + b[6]*incx3*incx3
          + b[7]*k*incx2 + b[8]*k*incx3 + b[9]*incx2*incx3;*/
      out << incx2 << "\t" << incx3 << "\t" << y << "\n";
      incx3 += dx3;
    }
    incx2 += dx2;
  }
}

void QYates::PolyGetPlotX2Const(QList< double > b, double x1_from, double x1_to, double x3_from, double x3_to, double k, uint step)
{
  QFile file("polysurface_x2_const.txt");
  if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
    return;

  QTextStream out(&file);
  
  
  // y = b0 + b1x1 + b2x2 + b3x3 + b11x1^2 + b22x2^2 + b33x3^2 + b12x1x2 + b13x1x3 + b23x2x3  
  double x1_range = floor(x1_to - x1_from);
  double x3_range = floor(x3_to - x3_from);
  double dx1 = x1_range / step;
  double dx2 = x3_range / step;
  
  double incx1 = x1_from;
  for(uint i = 0; i < step; i++){
    double incx3 = x3_from;
    for(uint j = 0; j < step; j++){
      double y = b[0] + 0.5*(b[1]*incx1 + b[2]*k + b[3]*incx3);
      /*double y = b[0] + b[1]*incx1 + b[2]*k + b[3]*incx3 + b[4]*incx1*incx1 + b[5]*k*k + b[6]*incx3*incx3
          + b[7]*incx1*k + b[8]*incx1*incx3 + b[9]*k*incx3;
          */
      out << incx1 << "\t" << incx3 << "\t" << y << "\n";
      incx3 += dx2;
    }
    incx1 += dx1;
  }
}

void QYates::PolyGetPlotX3Const(QList< double > b, double x1_from, double x1_to, double x2_from, double x2_to, double k, uint step )
{
  QFile file("polysurface_x3_const.txt");
  if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
    return;

  QTextStream out(&file);
  
  
  // y = b0 + b1x1 + b2x2 + b3x3 + b11x1^2 + b22x2^2 + b33x3^2 + b12x1x2 + b13x1x3 + b23x2x3  
  double x1_range = floor(x1_to - x1_from);
  double x2_range = floor(x2_to - x2_from);
  double dx1 = x1_range / step;
  double dx2 = x2_range / step;
  
  double incx1 = x1_from;
  for(uint i = 0; i < step; i++){
    double incx2 = x2_from;
    for(uint j = 0; j < step; j++){
      double y = b[0] + 0.5*(b[1]*incx1 + b[2]*incx2 + b[3]*k);
//       double y = b[0] + b[1]*incx1 + b[2]*incx2 + b[3]*k + b[4]*incx1*incx1 + b[5]*incx2*incx2 + b[6]*k*k
//           + b[7]*incx1*incx2 + b[8]*incx1*k + b[9]*incx2*k;
      out << incx1 << "\t" << incx2 << "\t" << y << "\n";
      incx2 += dx2;
    }
    incx1 += dx1;
  }

}


// use this with 15 experiments
void QYates::ResponseSurfExploration()
{ 
  ui.calculate->hide();
  ui.progressBar->show();

  QFuture< QList<double> > future = QtConcurrent::run(this, &QYates::RunSurfExploration);
  
  while(!future.isFinished())
    QApplication::processEvents();
 
  /* first 3 value are the real value x1, x2, x3. The other results are 
   * the coefficient value for the polynomial equation
   * y = b0 + b1x1 + b2x2 + b3x3 + b11x1^2 + b22x2^2 + b33x3^2 + b12x1x2 + b13x1x3 + b23x2x3
   */
  QList<double> result = future.result();
  
  
  qDebug() << "Coefficients";
  for(int i = 0; i < result.size(); i++){
    qDebug() << result[i];
  }
  
  if(ui.systematicMaxSearch->isChecked()){
    double x1_from, x1_to, x2_from, x2_to, x3_from, x3_to, y_m, x1_m, x2_m, x3_m;
    y_m = x1_m = x2_m = x3_m = 0.f;
    
    int step = ui.stepBox->value();
    
    x1_from = x1_to = getMatrixValue(signmx, 0, 0);
    for(uint i = 1; i < signmx->row; i++){
      if(getMatrixValue(signmx, i, 0) < x1_from)
        x1_from = getMatrixValue(signmx, i, 0);
      else if(getMatrixValue(signmx, i, 0) > x1_to)
        x1_to = getMatrixValue(signmx, i, 0);
      else
        continue;
    }
    
    x2_from = x2_to = getMatrixValue(signmx, 0, 1);
    for(uint i = 1; i < signmx->row; i++){
      if(getMatrixValue(signmx, i, 1) < x2_from)
        x2_from = getMatrixValue(signmx, i, 1);
      else if(getMatrixValue(signmx, i, 1) > x2_to)
        x2_to = getMatrixValue(signmx, i, 1);
      else
        continue;
    }
    
    x3_from = x3_to = getMatrixValue(signmx, 0, 2);
    for(uint i = 1; i < signmx->row; i++){
      if(getMatrixValue(signmx, i, 2) < x3_from)
        x3_from = getMatrixValue(signmx, i, 0);
      else if(getMatrixValue(signmx, i, 2) > x3_to)
        x3_to = getMatrixValue(signmx, i, 2);
      else
        continue;
    }
    
    GetMaxVal(result,
                x1_from, x1_to, step, 
                x2_from, x2_to, step,
                x3_from, x3_to, step,
                &y_m, &x1_m, &x2_m, &x3_m);
    
    ui.log_results->append("\n >>> Exploration of Maxima and Ridge System with Second-Order Response Surface\n");
    QString str;
    str += "coord. x:  "+QString::number(x1_m)+"\t";
    str += "First Variable: "+QString::number(IndVal2PredVal(0, x1_m));
    ui.log_results->append(str);
    str.clear();
    
    str += "coord. y:  "+QString::number(x2_m)+"\t";
    str += "Second Variable: "+QString::number(IndVal2PredVal(1, x2_m));
    ui.log_results->append(str);
    str.clear();
    
    
    str += "coord. z:  "+QString::number(x3_m)+"\t";
    str += "Third Variable: "+QString::number(IndVal2PredVal(2, x3_m));
    ui.log_results->append(str);
    str.clear();
    
    ui.log_results->append(QString("\nEstimated Area: %1").arg(y_m));
  }
  else{
    matrix *equation;
    dvector *solution;
    NewMatrix(&equation, 3,4);
    
    /* deriving the polynomial equation: 
    * y = b0 +b1*x1 + b2*x2 + b3*x3 + b11*x1^2 + b22*x2^2 + b33*x3^2 + b12*x1*x2 + b13*x1*x3 + b23*x2*x3
    * we get these tre equation:
    * 
    * 2*b11*x1 + b12*x2 + b13*x3 = -b1
    * 
    * b12*x1 + 2*b22*x2 + b23*x3 = -b2
    * 
    * b13*x1 + b23*x2 + 2*b33*x3 = -b3
    * 
    */

    equation->data[0][0] = 2*result[4]; equation->data[0][1] = result[7]; equation->data[0][2] = result[8]; equation->data[0][3] = -1*result[1];
    equation->data[1][0] = result[7]; equation->data[1][1] = 2*result[5]; equation->data[1][2] = result[9]; equation->data[1][3] = -1*result[2];
    equation->data[2][0] = result[8]; equation->data[2][1] = result[9]; equation->data[2][2] = 2*result[6]; equation->data[2][3] = -1*result[3];
    
    /*
    * box example 
    equation->data[0][0] = 9.38; equation->data[0][1] = 7.13; equation->data[0][2] = 3.27; equation->data[0][3] = 1.50;
    equation->data[1][0] = 7.13; equation->data[1][1] = 12.54; equation->data[1][2] = 2.73; equation->data[1][3] = -2.13;
    equation->data[2][0] = 3.27; equation->data[2][1] = 2.73; equation->data[2][2] = 10.42; equation->data[2][3] = 1.81;
    */
    
    initDVector(&solution);
    SolveLSE(equation, &solution);

    
    /*
    result.append(b->data[0]); // b0
    result.append(b->data[1]); // b1
    result.append(b->data[2]); // b2 
    result.append(b->data[2]); // b3
    */
    
    ui.log_results->append("\n >>> Exploration of Maxima and Ridge System with Second-Order Response Surface\n");
    QString str;
    str += "coord. x:  "+QString::number(solution->data[1])+"\t";
    str += "First Variable: "+QString::number(IndVal2PredVal(0, solution->data[1]));
    ui.log_results->append(str);
    str.clear();
    
    str += "coord. y:  "+QString::number(solution->data[2])+"\t";
    str += "Second Variable: "+QString::number(IndVal2PredVal(1, solution->data[2]));
    ui.log_results->append(str);
    str.clear();
    
    
    str += "coord. z:  "+QString::number(solution->data[3])+"\t";
    str += "Third Variable: "+QString::number(IndVal2PredVal(2, solution->data[3]));
    ui.log_results->append(str);
    str.clear();
    
    
     /* y = b0 +b1*x1 + b2*x2 + b3*x3 + b11*x1^2 + b22*x2^2 + b33*x3^2 + b12*x1*x2 + b13*x1*x3 + b23*x2*x3 */
    double y = result[0] + (result[1]*solution->data[0]) + (result[2]*solution->data[1]) + (result[3]*solution->data[2]) + result[4]*solution->data[0]*solution->data[0] + result[5]*solution->data[1]*solution->data[1] + result[6]*solution->data[2]*solution->data[2] + result[7]*solution->data[0]*solution->data[1] + result[8]*solution->data[0]*solution->data[2] + result[9]*solution->data[1]*solution->data[2];
    
    double area;
    if(ui.dipvar_log10->isChecked())
      area = pow(10, y);
    else
      area = y;

    ui.log_results->append(QString("\nEstimated Area: %1").arg(area));
    
    
    DelDVector(&solution);
    DelMatrix(&equation);
  }
    
  for(uint i = 0; i < 3; i++)
    result.removeFirst();
  /*
  for(uint i = 0; i < (uint)result.size(); i++)
    qDebug() << result[i];  
  */
  
  /*
  double x1_max, x1_min, x2_max, x2_min, x3_max, x3_min;
  
  FindMaxMin(0, &x1_max, &x1_min);
  FindMaxMin(1, &x2_max, &x2_min);
  FindMaxMin(2, &x3_max, &x3_min);
  
  PolyGetPlotX1Const(result, x2_min, x2_max, x3_min, x3_max, 200.0, 60);
  PolyGetPlotX2Const(result, x1_min, x1_max, x3_min, x3_max, 1000.0, 60);
  PolyGetPlotX3Const(result, x1_min, x1_max, x2_min, x2_max, 500, 60);
  */
 
  /*
  PolyGetPlotX1Const(result, 0, 7000, 0, 7000, -334.361, 60);
  PolyGetPlotX2Const(result, -400, 400, 0, 7000, 6641.43, 60);
  PolyGetPlotX3Const(result, -400, 400, 0, 7000, 6327.93, 60);
  */
  
  ui.progressBar->hide();
  ui.calculate->show();
}

// Used for the 2^3 design
void QYates::Yates()
{
  ui.calculate->hide();
  ui.progressBar->show();

  ui.log_results->append("\n >>> Indipendent Variable Effect Results\n");
  QList<double> result;
  for (uint i = 0; i < x->col; i++) {
    double v = 0.f;
    v = MinusSummer(i);
    v += PlusSummer(i);
    v /= 2;
    result.append(v); 
  }
  
  double sum = 0.f;
  foreach(double x, result){
    sum += fabs(x);
  }

  int i =0;
  
  foreach(double x, result){
    QString str;
    if(x < 0)
      str = tr(" var%1: %2 %. To get better result decrease the value of this variable.\n").arg(QString::number(i)).arg((fabs(x)*100)/sum);
    else
      str = tr(" var%1: %2 %. To get better result increase the value of this variable.\n").arg(QString::number(i)).arg((fabs(x)*100)/sum);
    
    ui.log_results->append(str);
    i++;
  }
  
  ui.calculate->show();
  ui.progressBar->hide();
}

QYates::QYates()
{
  ui.setupUi(this);
  
  initMatrix(&x);
  initMatrix(&signmx);
  initMatrix(&levels);
  initDVector(&y);
  
  ui.log_results->append("Welcome to QYates");
  
  connect(ui.buttonBox_2->button(QDialogButtonBox::Open), SIGNAL(clicked()), SLOT(Open()));
  connect(ui.buttonBox->button(QDialogButtonBox::Close), SIGNAL(clicked()), SLOT(Close()));
  connect(ui.calculate, SIGNAL(clicked()), SLOT(Calculate()));
  connect(ui.cleanWall, SIGNAL(clicked()), SLOT(cleanWall()));
  connect(ui.comboBox, SIGNAL(currentIndexChanged(int)), SLOT(HideShowVariable()));
  ui.progressBar->setMinimum(0);
  ui.progressBar->setMaximum(0);
  ui.progressBar->hide();
  
  HideShowVariable();
}

QYates::~QYates()
{
  DelMatrix(&x);
  DelMatrix(&signmx);
  DelMatrix(&levels);
  DelDVector(&y);
}
