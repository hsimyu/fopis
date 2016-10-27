#include "tdpic.h"

// Position Class
Position::Position(){}
Position::~Position(){}

void Position::set(double _x, double _y, double _z){
    x = _x;
    y = _y;
    z = _z;
}
double Position::getX() const { return x; }
double Position::getY() const { return y; }
double Position::getZ() const { return z; }
