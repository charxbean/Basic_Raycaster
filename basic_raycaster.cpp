#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <cctype>
#include <algorithm>

//vector type for: eye, viewdir, updir
typedef struct{
    float i, j, k;
} VectorType;

//vfov
float vfov;

//imsize;
int px_width;
int px_height;

//arbitrary d value
int d = 500;

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
    VectorType origin; //(x, y, z)
    VectorType intersection; //(dx, dy, dz)
} Raytype;

//initialize variables for the input
VectorType eye;
VectorType viewdir;
VectorType updir;
ColorType bkgcolor;
ColorType mtlcolor;

//array of material colors
std::vector<ColorType> materialArray;

//array of spheres
std::vector<SphereType> sphereArray;

//Returns a vector that is the cross product of vectors v1 and v2
VectorType crossProduct(VectorType v1, VectorType v2){
    float i = (v1.j * v2.k) - (v2.j * v1.k);
    float j = (v2.i * v1.k) - (v1.i * v2.k);
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

//normalizes the vector v1
void normalize(VectorType& v1){
    float length = std::sqrt(pow(v1.i, 2) + pow(v1.j, 2) + pow(v1.k, 2));
    v1.i = v1.i/length;
    v1.j = v1.j / length;
    v1.k = v1.k/length;
}

//finds and returns the closest positive t value beteween t1 and t2, returns -1 if both are invalid
float find_t(float t1, float t2){
    if(t1 > 0 && t2 > 0){
        if(t1 < t2){
            return t1;
        }
        else if(t2 < t1){
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

//takes the index of a sphere in the sphere array and returns the material color for that sphere
ColorType shadeRay(int sphere_index){
    int mat_index = sphereArray[sphere_index].m;
    ColorType color = materialArray[mat_index];
    return color;
}

//checks each object in the scene for a ray intersection and determines the closest valid intersection
//either returns the color of the object it intersects or returns the background color if no object is intersected
ColorType traceRay(Raytype ray){

    //variables to track the closest t value and corresponding sphere/object
    float t_closest = -1;
    int sphere_index;
    
    //iterate through each object in the scene to check for ray intersections
    for(int k = 0; k < sphereArray.size(); k++){
        float xc = sphereArray[k].x;
        float yc = sphereArray[k].y;
        float zc = sphereArray[k].z;
        float r = sphereArray[k].r;
        
        //calculate A, B and C for each sphere
        float A = pow(ray.intersection.i, 2) + pow(ray.intersection.j, 2) + pow(ray.intersection.k, 2);
        float B = 2*(ray.intersection.i*(ray.origin.i - xc) + ray.intersection.j * (ray.origin.j - yc) + ray.intersection.k*(ray.origin.k - zc));
        float C = pow((ray.origin.i - xc), 2) + pow((ray.origin.j - yc), 2) + pow((ray.origin.k - zc), 2) - pow(r, 2);
        

        //calculate the determinant
        float determinant = (pow(B, 2) - 4*A*C);
        
        //if the determinant is positive then calculate both t values of this object
        if(determinant > 0){
            float t1 = (-B + (std::sqrt(determinant)))/(2.0*A);
            float t2 = (-B - (std::sqrt(determinant)))/(2.0*A);
            
            float t = find_t(t1, t2); //find the closest positive t value out of the 2 calculated
            t_closest = find_t(t, t_closest); //keep track of the closest positive t of all objects
            if(t_closest == t){
                sphere_index = k; //keep track of the sphere index corresponding to the closest t value
            }
            
        } //if the determinant is 0 then calculate the t value of this object
        else if(determinant == 0){
            float t = (-B + (std::sqrt(determinant)))/(2.0*A);
            t_closest = find_t(t, t_closest); //keep track of the closest positive t of all objects
            if(t_closest == t){
                sphere_index = k; //keep track of the sphere index corresponding to the closest t value
            }
        }
        else{ //if no valid t value is found, t = -1
            float t = -1;
        }
        
    }

    //after iterating through every object, use the closest t value to determine the pixel color
    //if no valid t value was found, set the pixel to the background color
        if(t_closest == -1){
        return bkgcolor;
    }
    else{
        return shadeRay(sphere_index);
    }
    return bkgcolor;
}

//Main function
int main(int argc, const char * argv[]){
    
    //check for exactly one input file name as a command line arg, otherwise error
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
    if(fout.fail()){
        std::cerr << "Error: unable to open output file " << newFilename << std::endl;
        return -1; 
    } 

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
                //set sphere material to the index of the last material in the material array
                SphereType sphere = {num1, num2, num3, num4, static_cast<int>(materialArray.size() - 1)};
                sphereArray.push_back(sphere); //add spheres to the spere array
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
                    materialArray.push_back(mtlcolor); //add new material to the material array
                }
            } else {
                std::cerr << "Error parsing " << word << " line: " << line << std::endl;
            }
        }
        else if (word == "imsize") {
            if (ss >> num1 >> num2) {
                px_width = num1;
                px_height = num2;
            } else {
                std::cerr << "Error parsing imsize line: " << line << std::endl;
            }
        }
        else if (word == "vfov") {
            if (ss >> num1) {
               vfov = num1;
            } else {
                std::cerr << "Error parsing vfov line: " << line << std::endl;
            }
        }
        else {
            std::cerr << "Unrecognized keyword: " << word << " in line: " << line << std::endl;
        }
    }

    //create the header for the ppm file
    fout << "P3\n";
    fout << px_width << " # image width\n" << px_height << " # image height\n";
    fout << "255\n";
    
    //PRELIMINARY MATH

    //find vectors u and v
    VectorType u = crossProduct(viewdir, updir);
    normalize(u);
    VectorType v = crossProduct(u, viewdir);

    //calculate the width and height in 3d world coordinates
    float radian_vfov = (vfov * M_PI)/180;
    int height = (2*d*std::tan(0.5 * radian_vfov));
    int width = (px_width / px_height) * height;

    //Find the 4 corners of the viewing window
    VectorType ul;
    ul = vectorAdd(eye, vectorScalar(viewdir, d)); //e + d*n
    ul = vectorSubtract(ul, vectorScalar(u, (width/2.0))); //(w/2.0) * u
    ul = vectorAdd(ul, vectorScalar(v, (height/2.0))); //(h/2.0) * v

    VectorType ur;
    ur = vectorAdd(eye, vectorScalar(viewdir, d));
    ur = vectorAdd(ur, vectorScalar(u, (width/2.0)));
    ur = vectorAdd(ur, vectorScalar(v, (height/2.0)));
    
    VectorType ll;
    ll = vectorAdd(eye, vectorScalar(viewdir, d));
    ll = vectorSubtract(ll, vectorScalar(u, (width/2.0)));
    ll = vectorSubtract(ll, vectorScalar(v, (height/2.0)));
   
    VectorType lr;
    lr = vectorAdd(eye, vectorScalar(viewdir, d));
    lr = vectorAdd(lr, vectorScalar(u, (width/2.0)));
    lr = vectorSubtract(lr, vectorScalar(v, (height/2.0)));

    //find change in h and v
    VectorType h_change = vectorDivide(vectorSubtract(ur, ul), px_width-1);
    
    VectorType v_change = vectorDivide(vectorSubtract(ll, ul), px_height-1);

    //iterate through each pixel (j, i) where i = 0 to px_height-1, and j = 0 to px_width-1.
    for(int i = 0; i < px_height; i++){
        for(int j = 0; j < px_width; j++){

            //the point where each ray should pass through the viewing window correspoindng to each pixel
            VectorType intersect_point = vectorAdd(vectorAdd(ul, vectorScalar(h_change, j)), vectorScalar(v_change, i));

            //initialize a ray for this pixel
            Raytype ray;
            ray.origin = eye;
            ray.intersection = intersect_point;
            normalize(ray.intersection);

            //use trace ray to find the color for each pixel
            ColorType color = traceRay(ray);
            fout << color.r * 255 << " " << color.g * 255 << " " << color.b * 255 << std::endl;
        }
    }

    inputFile.close();
    return 0;
}