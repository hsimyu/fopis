#include "tdpic.h"

// potential
void FieldPointers::setPhi(threeD_array* _phi){
    pPhi = _phi;
}
threeD_array* FieldPointers::getPhi(){
    return pPhi;
}

// charge density
void FieldPointers::setRho(threeD_array* _rho){
    pRho = _rho;
}
threeD_array* FieldPointers::getRho(){
    return pRho;
}

// electric fields
void FieldPointers::setEx(threeD_array* _ex){
    pEx = _ex;
}
threeD_array* FieldPointers::getEx(){
    return pEx;
}

void FieldPointers::setEy(threeD_array* _ey){
    pEy = _ey;
}
threeD_array* FieldPointers::getEy(){
    return pEy;
}

void FieldPointers::setEz(threeD_array* _ez){
    pEz = _ez;
}
threeD_array* FieldPointers::getEz(){
    return pEz;
}

// magnetic fields
void FieldPointers::setBx(threeD_array* _bx){
    pBx = _bx;
}
threeD_array* FieldPointers::getBx(){
    return pBx;
}

void FieldPointers::setBy(threeD_array* _by){
    pBy = _by;
}
threeD_array* FieldPointers::getBy(){
    return pBy;
}

void FieldPointers::setBz(threeD_array* _bz){
    pBz = _bz;
}
threeD_array* FieldPointers::getBz(){
    return pBz;
}
