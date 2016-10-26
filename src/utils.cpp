#include <iostream>
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

    void print3DArray(threeD_array* data){
        double* it = (*data).begin();
        for ( int i = 0 ; i < (*data).shape()[0] ; ++i ) {
            cout << "[" << i << "]" << endl;
            for ( int j = 0 ; j < (*data).shape()[1] ; ++j ) {
                for ( int k = 0 ; k < (*data).shape()[2] ; ++k, ++it ) {
                    cout << *it << " ";
                }
            }
            cout << endl;
        }
    }
}
