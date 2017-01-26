#include <tdpic.h>

// Position Class
Position::Position(const double _x, const double _y, const double _z){
    this->setXYZ(_x, _y, _z);
}

Position::Position(const int _i, const int _j, const int _k){
    this->setIJK(_i, _j, _k);
}

Position::Position(const Particle& p){
    this->setXYZ(p.getX(), p.getY(), p.getZ());
}

Position::~Position(){}

void Position::setXYZ(const double _x, const double _y, const double _z){
    x = _x;
    y = _y;
    z = _z;

    i = floor(x) + 1;
    j = floor(y) + 1;
    k = floor(z) + 1;

    this->updateDelta();
}

void Position::setIJK(const int _i, const int _j, const int _k){
    i = _i;
    j = _j;
    k = _k;

    // x, y, zはglue cellなしの座標系にする
    x = static_cast<double>(i - 1);
    y = static_cast<double>(j - 1);
    z = static_cast<double>(k - 1);

    this->updateDelta();
}

void Position::updateDelta(void){
    // delta xyz is automatically updated
    dx1 = x - (i - 1);
    dy1 = y - (j - 1);
    dz1 = z - (k - 1);
    dx2 = 1.0 - dx1;
    dy2 = 1.0 - dy1;
    dz2 = 1.0 - dz1;
}
