#ifndef __TDPIC_UTILS_H_INCLUDED__
#define __TDPIC_UTILS_H_INCLUDED__
#include <picojson.h>
#include "global.hpp"

namespace Utils {
    void printTotalMemory(void);
    std::string prettyMemoryString(double);
    std::string readFile(const std::string&);
    picojson::value::object readJSONFile(const std::string&);

    // 軸方向を指定する文字列をindexに変換する
    int getAxisIndex(const AXIS);
    int getLowOrUpIndex(const AXIS_SIDE);

    float* getTrueEdges2(tdArray const&, tdArray const&, tdArray const&);
    float* getTrueEdges(tdArray const&, const int);
    float* getTrueFaces(const tdArray&, const int);
    void convert1Dto3Darray(double*, const int, const int, const int, tdArray&);
    void clearBoundaryValues(tdArray&, const int, const int, const int);
    void initializeTdarray(tdArray&);
    void createDir(std::string);
}

#endif
