#include <iostream>
#include "2dpic.h"

using std::cout;
using std::endl;

namespace Utils {
    void print3DArray(threeD_array* data){
        double* it = data->origin();
        for ( int i = 0 ; i < data->shape()[0] ; ++i ) {
            cout << "[" << i << "]" << endl;
            for ( int j = 0 ; j < data->shape()[1] ; ++j ) {
                for ( int k = 0 ; k < data->shape()[2] ; ++k, ++it ) {
                    cout << *it << " ";
                }
            }
            cout << endl;
        }
    }
}
