#ifndef __TDPIC_UTILS_H_INCLUDED__
#define __TDPIC_UTILS_H_INCLUDED__
#include <string>
#include <picojson.h>
#include "global.hpp"

namespace Utils {
    void printTotalMemory(void);
    std::string prettyMemoryString(double);
    std::string readFile(const std::string&);
    picojson::value::object readJSONFile(const std::string&);

    float* getTrueCells(const tdArray&);
    float* getTrueEdges(const tdArray&, const int);
    float* getTrueFaces(const tdArray&, const int);
    void convert1Dto3Darray(double*, const int, const int, const int, tdArray&);
    void clearBoundaryValues(tdArray&, const int, const int, const int);

    void createDir(std::string);

    class Normalizer {
        protected:
            // static class
            Normalizer();
            ~Normalizer();

        public:
            // static member
            static double x_unit;
            static double t_unit;
            static double m_unit;
            static double e_unit;

            static double normalizeLength(const double);
            static double unnormalizeLength(const double);

            static double normalizeVelocity(const double);
            static double unnormalizeVelocity(const double);

            static double normalizeTime(const double);
            static double unnormalizeTime(const double);

            static double normalizeCharge(const double);
            static double unnormalizeCharge(const double);

            static double normalizeMass(const double);
            static double unnormalizeMass(const double);

            static double normalizeEnergy(const double);
            static double unnormalizeEnergy(const double);

            static double normalizeEpsilon(const double);
            static double unnormalizeEpsilon(const double);
    };
}


#endif
