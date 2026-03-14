#include "types.hpp"
#include "surface.hpp"

int main(int argc, char** argv){
    std::cout << argv[0] << " : " 
              << argv[1] << " : "  
              << argv[2] << " : "  
              << argv[3] << " : "  << std::endl;
    if(argc < 4){
        std::cout << "Usage (transfer bezier to linear): ./transfer --bzp_to_mesh [INPUT_CAGE] [TARGET_CAGE]" << std::endl;
        return 0;
    } 
    string arg(argv[1]);
    if(arg == "--bzp_to_mesh"){
        string input_file = string(argv[2]); 
        string output_file = string(argv[3]); 
        std::cout << "--> Transfer " << input_file << " to " << output_file << std::endl;
        BezierLoader cage(input_file);
        MeshLoader cage_new = cage.to_mesh();
        cage_new.save(output_file);
    }
}