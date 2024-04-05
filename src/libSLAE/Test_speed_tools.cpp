#include "Test_speed_tools.hpp"

void print_to_file(const test_pair& to_write, std::string filename){
    std::fstream out(filename + ".txt", std::ios_base::app);
    for(auto i = 0u; i < to_write.first.size(); i++){
        out << to_write.first[i] << std::endl;
    }
    for(auto i = 0u; i < to_write.second.size(); i++){
        out << to_write.second[i] << std::endl;
    }
    out.close();
}