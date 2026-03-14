#include <map>
#include <string>
#include <vector>
#include <functional>

using std::map;
using std::string;
using std::vector;
using std::function;

template<typename T>
vector<T> erase(vector<T> v, T d){
    vector<T> ret;
    for(auto _v : v){
        if (_v != d)
            ret.push_back(_v);
    }
    return ret;
}

template<typename T>
class subsec_manager{
public:
void add_object(string name, T obj){
    data[name] = obj;
    if(std::find(names.begin(), names.end(), name) == names.end())
        names.push_back(name);
}
void apply(function<void(string, T)> f){
    for(auto name : names){
        f(name, data[name]);
    }
}   
void clear(){
    data.clear();
    names.clear();
}
private:
    map<string, T> data;
    vector<string> names;
};

template<typename T>
class sec_manager{
public:
void add_object(string sec, string name, T obj){
    data[sec][name] = obj;
    if(std::find(sec_names.begin(), sec_names.end(), sec) == sec_names.end())
        sec_names.push_back(sec);
    if(std::find(subsec_names[sec].begin(), subsec_names[sec].end(), name) == subsec_names[sec].end())
        subsec_names[sec].push_back(name);
}

void apply(function<void(string, string, T)> f){
    for(auto sec : sec_names){
        for(auto subsec : subsec_names[sec])
            f(sec, subsec, data[sec][subsec]);
    }
}

void apply(string sec, function<void(string, T)> f){
  for(auto subsec : subsec_names[sec])
      f(subsec, data[sec][subsec]);
}

void apply(string sec, string subsec, function<void(T)> f){
    f(data[sec][subsec]);
}

void clear(){
    data.clear();
    sec_names.clear();
    subsec_names.clear();
}

void clear(string sec, string subsec){
    data[sec].erase(subsec);
    subsec_names[sec] = erase(subsec_names[sec], subsec);
    if(subsec_names[sec].size() == 0) {
        data.erase(sec);
        sec_names = erase(sec_names, sec);
        subsec_names.erase(sec);
    }
}

private:
    map<string, map<string, T> > data;
    vector<string> sec_names;
    map<string, vector<string> > subsec_names;
};

vector<string> split(string data, char key){
    vector<string> ret;
    ret.push_back({});
    for(auto s : data){
        if(s != key)
            ret.back().push_back(s);
        else if(ret.back().size() != 0)
            ret.push_back({});
    }
    if(ret.size() == 1) return {"", ret[0]};
    if(ret.size() == 0) return {"", ""};
    return ret;
}

template<typename T>
void unique_append(vector<T>& l, T d)
{
    for(int i = 0; i < l.size(); i++)
        if(l[i] == d) return;
    l.push_back(d);
}

