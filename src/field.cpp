#include <tdpic.h>

// potential
void Field::setPhi(threeD_array* _phi){
    pPhi = _phi;
}
threeD_array* Field::getPhi(){
    return pPhi;
}

// charge density
void Field::setRho(threeD_array* _rho){
    pRho = _rho;
}
threeD_array* Field::getRho(){
    return pRho;
}

// electric fields
void Field::setEx(threeD_array* _ex){
    pEx = _ex;
}
threeD_array* Field::getEx(){
    return pEx;
}

void Field::setEy(threeD_array* _ey){
    pEy = _ey;
}
threeD_array* Field::getEy(){
    return pEy;
}

void Field::setEz(threeD_array* _ez){
    pEz = _ez;
}
threeD_array* Field::getEz(){
    return pEz;
}

// magnetic fields
void Field::setBx(threeD_array* _bx){
    pBx = _bx;
}
threeD_array* Field::getBx(){
    return pBx;
}

void Field::setBy(threeD_array* _by){
    pBy = _by;
}
threeD_array* Field::getBy(){
    return pBy;
}

void Field::setBz(threeD_array* _bz){
    pBz = _bz;
}
threeD_array* Field::getBz(){
    return pBz;
}

// destructor
Field::~Field(){
    delete [] pPhi;
    delete [] pRho;
    delete [] pEx;
    delete [] pEy;
    delete [] pEz;
    delete [] pBx;
    delete [] pBy;
    delete [] pBz;
}