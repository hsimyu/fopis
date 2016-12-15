#include <tdpic.h>

// potential
void Field::setPhi(threeDArray _phi){
    pPhi = _phi;
}
threeDArray Field::getPhi(){
    return pPhi;
}

// charge density
void Field::setRho(threeDArray _rho){
    pRho = _rho;
}
threeDArray Field::getRho(){
    return pRho;
}

// electric fields
void Field::setEx(threeDArray _ex){
    pEx = _ex;
}
threeDArray Field::getEx(){
    return pEx;
}

void Field::setEy(threeDArray _ey){
    pEy = _ey;
}
threeDArray Field::getEy(){
    return pEy;
}

void Field::setEz(threeDArray _ez){
    pEz = _ez;
}
threeDArray Field::getEz(){
    return pEz;
}

// magnetic fields
void Field::setBx(threeDArray _bx){
    pBx = _bx;
}
threeDArray Field::getBx(){
    return pBx;
}

void Field::setBy(threeDArray _by){
    pBy = _by;
}
threeDArray Field::getBy(){
    return pBy;
}

void Field::setBz(threeDArray _bz){
    pBz = _bz;
}
threeDArray Field::getBz(){
    return pBz;
}

// destructor
Field::~Field(){
    Utils::delete3DArray(pPhi);
    Utils::delete3DArray(pRho);
    Utils::delete3DArray(pEx);
    Utils::delete3DArray(pEy);
    Utils::delete3DArray(pEz);
    Utils::delete3DArray(pBx);
    Utils::delete3DArray(pBy);
    Utils::delete3DArray(pBz);
}
