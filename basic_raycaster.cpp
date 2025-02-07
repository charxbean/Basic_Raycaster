#include <fstream>
#include <iostream>
#include <string>
#include <sstream>


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
    float eye_i;
    float eye_j;
    float eye_k;

    std::string viewdir;
    float viewdir_i;
    float viewdir_j;
    float viewdir_k;

    std::string updir;
    float updir_i;
    float updir_j;
    float updir_k;

    std::string h_fov;
    float hfov;

    std::string imsize;
    int width;
    int height;

    std::string bkgcolor;
    float r;
    float g;
    float b;

    std::string mtlcolor;
    int m_r;
    int m_g;
    int m_b;

    std::string sphere; //rename components
    float one;
    float two;
    float three;
    float four;

    std::string line;
    //check if input file is in correct format

    while(std::getline(inputFile, line)){
        std::istringstream ss(line);
        std::string word;
        float num1;
        float num2;
        float num3;

        if(ss >> word >> num1 >> num2 >> num3){
            eye_i = num1;
            eye_j = num2;
            eye_k = num3;
            fout << num1 << std::endl;
        }
    }


    /*
     std::string delimiter = ",";
    while(std::getline(inputFile, line)){
        val1 = line.substr(0, line.find(delimiter));
        val2 = line.substr(line.find(delimiter) + 2, line.find(delimiter));
        val3 = line.substr(line.find(delimiter) + line.find(delimiter) + 2, line.find(delimiter));
        val4 = line.substr(line.find(delimiter)+ line.find(delimiter) + line.find(delimiter) + 3, line.find(delimiter));

        
        std::cout << line << std::endl;
        std::cout << val1 << std::endl;
        std::cout << val2 << std::endl;
        std::cout << val3 << std::endl;
        std::cout << val4 << std::endl;
        

        if(val1 == "eye"){
            eye_i = std::stod(val2);
            eye_j = std::stod(val3);
            eye_k = std::stod(val4);
        }
        else{
            std::cout << "wtf" << line << std::endl;
            return -1;
        }
    }

    std::cout << val1 << " " << val2  << " " << val3  << " " << val4 << std::endl;
    */
   
    
    inputFile.close();
    return 0;
}