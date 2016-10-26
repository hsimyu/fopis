#include <iostream>
#include <boost/format.hpp>
#include "2dpic.h"

using std::cout;
using std::endl;

namespace Utils {
    template <typename T>
    void print2DArray(T** data, int nx, int ny){
        for(int i = 0; i < nx; i++){
            for(int j = 0; j < ny; j++){
                cout << data[i][j] << " ";
            }
            cout << endl;
        }
    }
}
