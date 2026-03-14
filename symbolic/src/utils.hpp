#pragma once
#include <iostream>
#include <fstream>
#include <map>
#include <ctime> 
#include <vector>
#include <string>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/string.hpp>

namespace utils{

template <typename T>
T _sum_(std::vector<T> a) {
    T ret{};
    for(auto _a : a) ret+=_a;
    return ret;
}

class _timer_{
public:
    std::map<std::string, std::vector<double>> times;
    std::map<std::string, std::clock_t> starts;
    std::map<std::string, std::clock_t> ends;
    bool exit_info = true;
public:
    _timer_(bool e){exit_info = e;}
    _timer_(){starts["__TIMMER_START___"] = std::clock();};
    ~_timer_(){
        auto start = starts["__TIMMER_START___"]; 
        double time = double(clock() - start) / CLOCKS_PER_SEC;
        if(exit_info)
            std::cout << "*************** timer exists for " << time << " secends. **********************" << std::endl;
    };
public:
    void tik(std::string name){
        starts[name] = clock();
        ends[name] = clock();
    }

    void tok(std::string name){
        ends[name] = clock();
        times[name].push_back( double(ends[name] - starts[name]) / CLOCKS_PER_SEC);
    }

    void summary(){
        std::cout << "********************************************************************" << std::endl;
        for(auto [name, time] : times){
            std::cout << "\t[" << name << "]: t_all = " << _sum_(time) << "(s),\t" ;
            std::cout << "\t\t n_all = " << time.size() << ", t_avg = " << _sum_(time) / time.size() << " s/r." << std::endl;
        }
    }
};


template <typename T>
void save_binary(const T& obj, const std::string& filename) {
    std::ofstream ofs(filename, std::ios::binary);
    if(!ofs){throw std::runtime_error("Cannot open file for writing!");}
    boost::archive::binary_oarchive oa(ofs); oa << obj;
    ofs.close();
}


template <typename T>
T load_binary(const std::string& filename) {
    T ret;
    std::ifstream ifs(filename, std::ios::binary);
    if(!ifs){throw std::runtime_error("Cannot open file for reading!");}
    boost::archive::binary_iarchive ia(ifs); ia >> ret; ifs.close();
    return ret;
}

template<typename real>
void save(std::string filename, const std::vector<real>& d){
    std::ofstream os(filename);
    for(auto _d : d) os << _d << " ";
    os.close();
}
}
