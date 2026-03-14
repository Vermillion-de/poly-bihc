#pragma once
#include <iostream>
#include <fstream>
#include <map>
#include <ctime> 
#include <vector>
#include <string>

namespace utils{
template<typename A, typename B>
void save(std::vector<double>& data, std::map<A, B> m){
    data.push_back(m.size());
    for(auto [a, b] : m){
        save(data, a);
        save(data, b);
    }
}

void save(std::vector<double>& data, int d){
    data.push_back(d);
}

void save(std::vector<double>& data, double d){
    data.push_back(d);
}

};