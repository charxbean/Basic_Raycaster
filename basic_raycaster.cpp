#include <fstream>
#include <iostream>
#include <string>


//read text file from command line

//create PPM output file like last assignment

/*
You need to provide build / run instructions and input files.
 Points will be deducted on future assignments.

 -- Check what that means????
*/

int main(int argc, const char * argv[]){
    
    //check for one input file name
    if(argc < 2){
        std::cerr << "Format: program inputFile.txt" << std::endl;
        return -1;
    }
    else if(argc > 2){
        std::cerr << "Format: program inputFile.txt" << std::endl;
        return -1;
    }

    std::string filename = argv[1];

    //open the input file stream
    std::ifstream inputFile(filename);
    if(!inputFile.is_open()){
        std::cerr << "Error: unable to open file: " << filename << std::endl;
        return -1;
    }

    //open the output file stream, set the output file name as input filename.ppm
    std::string newFilename = filename.erase(filename.length() - 4, 4);
    std::ofstream fout(newFilename + ".ppm"); 
    if(fout.fail()) return -1; //error opening file


    //set up variables to be used
    std::string eye;
    int eye_i;
    int eye_j;
    int eye_k;

    std::string viewdir;
    int viewdir_i;
    int viewdir_j;
    int viewdir_k;

    std::string updir;
    int updir_i;
    int updir_j;
    int updir_k;

    std::string h_fov;
    int hfov;

    std::string imsize;
    int width;
    int height;

    std::string bkgcolor;
    int r;
    int g;
    int b;

    std::string mtlcolor;
    int m_r;
    int m_g;
    int m_b;

    std::string sphere; //rename components
    int one;
    int two;
    int three;
    int four;

    std::string line, val1, val2, val3, val4;
    //check if input file is in correct format

    std::string deliminater = " ";
    while(std::getline(inputFile, line)){
        val1 = line.substr(0, line.find(deliminater));
        val2 = line.substr(0, line.find(deliminater) + 1);
        val3 = line.substr(line.find(deliminater) + 2);

        std::cout << line << std::endl;
        std::cout << val1 << std::endl;
        std::cout << val2 << std::endl;

        if(val1 == "eye"){
            eye_i = std::stoi(val2);
            eye_j = std::stoi(val3);
            eye_k = std::stoi(val4);
        }
        else{
            std::cout << "wtf" << line << std::endl;
            return -1;
        }
    }

    std::cout << val1 << val2 << val3 << val4 << std::endl;
    
    inputFile.close();
    return 0;
}