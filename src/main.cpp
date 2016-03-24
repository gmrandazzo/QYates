#include <QtGui/QApplication>
#include "QYates.h"


int main(int argc, char** argv)
{
    QApplication app(argc, argv);
    QYates foo;
    foo.show();
    return app.exec();
}
