#ifndef __TDPIC_UTILS_H_INCLUDED__
#define __TDPIC_UTILS_H_INCLUDED__
#include <picojson.h>
#include "global.hpp"
#include <chrono>

namespace Utils {
    void printTotalMemory(void);
    std::string prettyMemoryString(double);

    //! ファイル操作関連
    bool isExistingFile(const std::string& file_name);
    bool isExistingDirectory(const std::string& dir_name);
    std::string readFile(const std::string&);
    picojson::value::object readJSONFile(const std::string&);
    void createDir(std::string);
    std::string extractFileName(const std::string& target);

    //! 保存用のデータ変換系
    std::vector<double> convertPicoJSONArrayToVectorDouble(const picojson::array& pico_array);
    std::vector<std::string> convertPicoJSONArrayToVectorString(const picojson::array& pico_array);
    std::vector<std::string> split(const std::string& target, char delim);

    // 軸方向を指定する文字列をindexに変換する
    int getAxisIndex(const AXIS);
    int getLowOrUpIndex(const AXIS_SIDE);

    void initialize3DArray(tdArray&);
    void initializeRhoArray(std::vector<tdArray>&);

    class ProgressManager {
    private:
        //! 時計
        using Time = std::chrono::high_resolution_clock;
        Time::time_point begin;

        float maximum_step;
        std::string activity_name;
        float next_report_target; // %単位
        float report_width;

    public:
        // float に cast 可能な型ならなんでもよい
        template<typename T>
        ProgressManager(const T _max) : 
            begin{Time::now()},
            maximum_step{static_cast<float>(_max)},
            activity_name("activity"),
            next_report_target{10.0},
            report_width{10.0} {}

        //! S&&はユニバーサル参照
        template<typename T, typename S>
        ProgressManager(const T _max, S&& name) : 
            begin{Time::now()},
            maximum_step{static_cast<float>(_max)},
            activity_name(std::forward<S>(name)),
            next_report_target{10.0},
            report_width{10.0} {}

        ~ProgressManager(void) {
            const auto elapsed_time = std::chrono::duration_cast<std::chrono::seconds>(Time::now() - begin);
            cout << "[" << activity_name << "] Processed. (" << elapsed_time.count() << "sec elapsed.)" << endl;
        }

        template<typename T>
        float progress(const T step) {
            return 100.0 * static_cast<float>(step) / maximum_step;
        }

        template<typename T>
        void update(const T step) {
            const auto p = progress(step);
            if ( p > next_report_target ) {
                const auto elapsed_time = std::chrono::duration_cast<std::chrono::seconds>(Time::now() - begin);
                cout << format("[%s] Progress: %4.2f%%. (%ssec elapsed.)") % activity_name % p % elapsed_time.count() << endl;
                next_report_target += report_width;
            }
        }
    };

    class TimeCounter {
    private:
        // singleton
        static TimeCounter* instance;

        using Time = std::chrono::high_resolution_clock;
        Time::time_point begin_time_point;

        std::map<std::string, std::chrono::milliseconds> elapsed_times;
        std::string current_activity;
        bool print_report_on_destruct;

        void setCurrentActivity(const std::string& activity_name) {
            current_activity = activity_name;
        }

        TimeCounter() : elapsed_times{}, current_activity{}, print_report_on_destruct{false} {}
        ~TimeCounter() {
            if (print_report_on_destruct) {
                cout << "====== Computation Time Report ======" << endl;
                std::chrono::milliseconds total_time(0);
                for(const auto& pairs : elapsed_times) {
                    cout << format("  %s: %s msecs.") % pairs.first % pairs.second.count() << endl;
                    total_time += pairs.second;
                }
                cout << format("Total: %s msecs.") % total_time.count() << endl;
            }
        }

    public:
        static TimeCounter* getInstance() {
            if (instance == nullptr) {
                instance = new TimeCounter;
            }
            return instance;
        }

        static void destroy() {
            delete instance;
        }

        void begin(const std::string& activity_name) {
            setCurrentActivity(activity_name);

            begin_time_point = Time::now();
        }

        void end() {
            if (current_activity != "") {
                const auto elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(Time::now() - begin_time_point);
                constexpr char delim = '/';

                auto activities = split(current_activity, delim);
                std::string nested_activity_name = "";

                for (const auto& act : activities) {
                    nested_activity_name += delim + act;
                    elapsed_times[nested_activity_name] += elapsed_time;
                }
                setCurrentActivity("");
            }
        }

        void switchTo(const std::string& activity_name) {
            end();
            begin(activity_name);
        }

        void enableReport() {
            print_report_on_destruct = true;
        }
    };
}

#endif
