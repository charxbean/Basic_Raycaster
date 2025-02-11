#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <cctype>
#include <algorithm>

//read text file from command line

//create PPM output file like last assignment

/*
You need to provide build / run instructions and input files.
 Points will be deducted on future assignments.

 -- Check what that means????
*/

//vector type for: eye, viewdir, updir
typedef struct{
    float i, j, k;
} VectorType;

//h_fov;
float vfov;

//imsize;
int px_width;
int px_height;

//arbitrary d value
int d = 5;

//color type for bg color and material color
typedef struct{
    float r, g, b;
} ColorType;

//sphere defined by it's center point (x, y, z) and a radius
    //has m as it's material property -> an index into an array of material properties
typedef struct{
    float x, y, z;
    float r;
    int m;
} SphereType;

//equation of a ray = (x, y, z) - t*(dx, dy, dz)
typedef struct{
    VectorType origin;
    VectorType intersection;
} Raytype;

//Returns a vector that is the cross product of vectors v1 and v2
VectorType crossProduct(VectorType v1, VectorType v2){
    float i = (v1.k * v2.i) - (v2.k * v1.i);
    float j = (v1.j * v2.k) - (v2.j * v1.k);
    float k = (v1.i * v2.j) - (v2.i * v1.j);
    VectorType newVector;
    newVector.i = i;
    newVector.j = j;
    newVector.k = k;

    return newVector;
}

//returns a new vector that is the sum of vectors v1 and v2
VectorType vectorAdd(VectorType v1, VectorType v2){
    VectorType newVec;
    newVec.i = v1.i + v2.i;
    newVec.j = v1.j + v2.j;
    newVec.k = v1.k + v2.k;
    
    return newVec;
}

//returns a new vector that is the difference of vectors v1 and v2
VectorType vectorSubtract(VectorType v1, VectorType v2){
    VectorType newVec;
    newVec.i = v1.i - v2.i;
    newVec.j = v1.j - v2.j;
    newVec.k = v1.k - v2.k;
    
    return newVec;
}

//returns a new vector that is vector v1 multiplied by vector v2
VectorType vectorMultiply(VectorType v1, VectorType v2){
    VectorType newVec;
    newVec.i = v1.i * v2.i;
    newVec.j = v1.j * v2.j;
    newVec.k = v1.k * v2.k;
    
    return newVec;
}

//returns a new vector that is v1 multiplied by the scalar
VectorType vectorScalar(VectorType v1, float scalar){
    VectorType newVec;
    newVec.i = v1.i * scalar;
    newVec.j = v1.j * scalar;
    newVec.k = v1.k * scalar;

    return newVec;
}

//returns a new vector that is v1 divided by the scalar
VectorType vectorDivide(VectorType v1, float scalar){
    VectorType newVec;
    newVec.i = v1.i / scalar;
    newVec.j = v1.j / scalar;
    newVec.k = v1.k / scalar;
    
    return newVec;
}

//normalizes the input vector v1
void normalize(VectorType& v1){
    float length = std::sqrt(pow(v1.i, 2) + pow(v1.j, 2) + pow(v1.k, 2));
    v1.i = v1.i/length;
    v1.j = v1.j / length;
    v1.k = v1.k/length;
}

float find_t(float t1, float t2){
    if(t1 > 0 && t2 > 0){
        if(t1 > t2){
            return t1;
        }
        else if(t2 > t1){
            return t2;
        }else{
            return t1;
        }
    }
    else if(t1 < 0 && t2 > 0){
        return t2;
    }else if(t2 < 0 && t1 > 0){
        return t1;
    }else{
        return -1;
    }
}
//prints a vector to the terminal
void printVector(VectorType v1){
    std::cout << "<" << v1.i << ", " << v1.j << ", " << v1.k << ">" << std::endl;
}

int main(int argc, const char * argv[]){
    
    //check for one input file name
    if(argc < 2){
        std::cerr << "Format: executable file.txt" << std::endl;
        return -1;
    }
    else if(argc > 2){
        std::cerr << "Format: executable file.txt" << std::endl;
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

    VectorType eye;
    VectorType viewdir;
    VectorType updir;
    ColorType bkgcolor;
    ColorType mtlcolor;

    //array of material colors
    std::vector<ColorType> materialArray;

    //the array of pixel colors to be output to the PPM file
    std::vector<ColorType> outputColors;
    
    //array of spheres
    std::vector<SphereType> sphereArray;


    //read each line of the input file and put file data into respective variables/structs
    std::string line;
    while (std::getline(inputFile, line)) {

        // Skip empty lines
        if (std::all_of(line.begin(), line.end(), ::isspace)) {
            continue;
        }

        std::istringstream ss(line);
        std::string word;
        float num1, num2, num3, num4;

        ss >> word; //Read the first word of every line

        if (word.empty()) {
            std::cerr << "Error: Empty keyword in line: " << line << std::endl;
            continue;
        }

        if (word == "sphere") {
            if (ss >> num1 >> num2 >> num3 >> num4) {
                SphereType sphere{num1, num2, num3, num4, static_cast<int>(materialArray.size() - 1)};
                sphereArray.push_back(sphere);
                fout << word << " " << num1 << " " << num2 << " " << num3 << " " << num4 << std::endl;
            } else {
                std::cerr << "Error parsing sphere line: " << line << std::endl;
            }
        }
        else if (word == "eye" || word == "viewdir" || word == "updir" || word == "bkgcolor" || word == "mtlcolor") {
            if (ss >> num1 >> num2 >> num3) {
                if (word == "eye")         eye = {num1, num2, num3};
                if (word == "viewdir")     viewdir = {num1, num2, num3};
                if (word == "updir")       updir = {num1, num2, num3};
                if (word == "bkgcolor")    bkgcolor = {num1, num2, num3};
                if (word == "mtlcolor") {
                    mtlcolor = {num1, num2, num3};
                    materialArray.push_back(mtlcolor);
                }
                fout << word << " " << num1 << " " << num2 << " " << num3 << std::endl;
            } else {
                std::cerr << "Error parsing 3-float line: " << line << std::endl;
            }
        }
        else if (word == "imsize") {
            if (ss >> num1 >> num2) {
                px_width = num1;
                px_height = num2;
                fout << word << " " << num1 << " " << num2 << std::endl;
            } else {
                std::cerr << "Error parsing imsize line: " << line << std::endl;
            }
        }
        else if (word == "vfov") {
            if (ss >> num1) {
                if (word == "vfov") vfov = num1;
                fout << word << " " << num1 << std::endl;
            } else {
                std::cerr << "Error parsing vfov line: " << line << std::endl;
            }
        }
        else {
            std::cerr << "Unrecognized keyword: " << word << " in line: " << line << std::endl;
        }
    }

    
    
    //find vectors u and v
    VectorType u = crossProduct(viewdir, updir);
    normalize(u);
    VectorType v = crossProduct(u, viewdir);

    //calculate the width and height in 3d world coordinates
    int height = 2*d*std::tan((1/2) * vfov);
    int width = (px_width / px_height) * height;

    //Find the 4 corners of the viewing window
    VectorType ul;
    ul = vectorAdd(eye, vectorScalar(viewdir, d)); //e + d*n
    ul = vectorSubtract(ul, vectorScalar(u, (width/2))); //(w/2) * u
    ul = vectorAdd(ul, vectorScalar(v, (height/2))); //(h/2) * v
    printVector(ul);

    VectorType ur;
    ur = vectorAdd(eye, vectorScalar(viewdir, d));
    ur = vectorAdd(ur, vectorScalar(u, (width/2)));
    ur = vectorAdd(ur, vectorScalar(v, (height/2)));

    VectorType ll;
    ll = vectorAdd(eye, vectorScalar(viewdir, d));
    ll = vectorSubtract(ll, vectorScalar(u, (width/2)));
    ll = vectorSubtract(ll, vectorScalar(v, (height/2)));

    VectorType lr;
    lr = vectorAdd(eye, vectorScalar(viewdir, d));
    lr = vectorAdd(lr, vectorScalar(u, (width/2)));
    lr = vectorSubtract(lr, vectorScalar(v, (height/2)));

    VectorType h_change = vectorDivide(vectorSubtract(ur, ul), px_width);
    VectorType v_change = vectorDivide(vectorSubtract(ll, ul), px_height);

    //iterate through each pixel (i, j)
    //viewing window location corresponding to pixel (i, j) == ul + i * ^h + j *^v
    for(int i = 0; i < width -1; i++){
        for(int j = 0; j < height -1; j++){
            //the point where each ray should pass through the viewing window correspoindng to each pixel
            VectorType intersect_point = vectorAdd(vectorAdd(ul, vectorScalar(h_change, i)), vectorScalar(v_change, j));
            Raytype ray;
            ray.origin = eye;
            ray.intersection = intersect_point;
            
            //calculate A, B and C for each sphere
            for(int i = 0; i < sphereArray.size(); i++){
                float xc = sphereArray[i].x;
                float yc = sphereArray[i].y;
                float zc = sphereArray[i].z;
                float r = sphereArray[i].r;
                
                float A = pow(ray.intersection.i, 2) + pow(ray.intersection.j, 2) + pow(ray.intersection.k, 2);
                float B = 2*(ray.intersection.i*(ray.origin.i - xc) + ray.intersection.j * (ray.origin.j - yc) + ray.intersection.k*(ray.origin.k - zc));
                float C = pow((ray.origin.i - xc), 2) + pow((ray.origin.j - yc), 2) + pow((ray.origin.k - zc), 2) - pow(r, 2);
                
                float determinant = (pow(B, 2) - 4*A*C);
                if(determinant > 0){
                    float t1 = -1 * B + (std::sqrt(determinant))/2*A;
                    float t2 = -1 * B - (std::sqrt(determinant))/2*A;
                    float t = find_t(t1, t2);
                }
                else{
                    float t = -1;
                }
                //how to check closest t-value between different spheres>
                //have an array of rays, an array of t values (-1 if intersection wasn't found)

                //identifier to knwo which object was intersected: the sphere's index
            }

        }
    }


    //double for loop to iterate through each pixel (i, j)
    //for each pixel, find the equation of the ray that goes through the pixel at (ul + i * ^h + j *^v)
    //for each ray, check if it intersects with a sphere by checking the discriminant
    //if the discriminant is positive then find the intersecting points using the quadratic equation and choose the closest intersection
    //determine which sphere it is by it's position and get it's material color from that sphere's struct
    //put that color into the array corresponding to that pixel

    //at the end, output the color data in the array to the output file
    
    inputFile.close();
    return 0;
}

