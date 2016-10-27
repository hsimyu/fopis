#include "2dpic.h"

void FieldPointers::setPhi(threeD_array* _phi){
    pPhi = _phi;
}

threeD_array* FieldPointers::getPhi(){
    return pPhi;
}

void FieldPointers::setRho(threeD_array* _rho){
    pRho = _rho;
}

threeD_array* FieldPointers::getRho(){
    return pRho;
}
