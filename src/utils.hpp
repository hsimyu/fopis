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

    class ProgressManager {
    private:
        std::string activity_name;
        float report_width;
        float next_report_target; // %単位
        float maximum_step;

    public:
        // float に cast 可能な型ならなんでもよい
        template<typename T>
        ProgressManager(const T _max) : 
            next_report_target{10.0},
            report_width{10.0},
            maximum_step{static_cast<float>(_max)},
            activity_name("activity") {}

        //! S&&はユニバーサル参照
        template<typename T, typename S>
        ProgressManager(const T _max, S&& name) : 
            next_report_target{10.0},
            report_width{10.0},
            maximum_step{static_cast<float>(_max)},
            activity_name(std::forward<S>(name)) {}

        ~ProgressManager(void) {
            cout << "[" << activity_name << "] processed." << endl;
        }

        template<typename T>
        float progress(const T step) {
            return 100.0 * static_cast<float>(step) / maximum_step;
        }

        template<typename T>
        void update(const T step) {
            const auto p = progress(step);
            if ( p > next_report_target ) {
                cout << "[" << activity_name << "] Progress: " << p << "%. " << endl;
                next_report_target += report_width;
            }
        }
    };
}

#endif
