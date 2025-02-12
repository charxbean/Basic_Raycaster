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

//finds and returns the closest positive t value beteween t1 and t2, returns -1 if both are invalid
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

//iterate through each object in the scene and check for a ray intersection
//returns the pixel color at the intersection point, or the background color if no intersection is found
ColorType TraceRay(Raytype ray, std::vector<SphereType> &sphereArray){

    for(int i = 0; i < sphereArray.size(); i++){
        float xc = sphereArray[i].x;
        float yc = sphereArray[i].y;
        float zc = sphereArray[i].z;
        float r = sphereArray[i].r;
        
        float A = pow(ray.intersection.i, 2) + pow(ray.intersection.j, 2) + pow(ray.intersection.k, 2);
        float B = 2*(ray.intersection.i*(ray.origin.i - xc) + ray.intersection.j * (ray.origin.j - yc) + ray.intersection.k*(ray.origin.k - zc));
        float C = pow((ray.origin.i - xc), 2) + pow((ray.origin.j - yc), 2) + pow((ray.origin.k - zc), 2) - pow(r, 2);
        
        //check if the determinant is negative
        float determinant = (pow(B, 2) - 4*A*C);
        if(determinant > 0){
            float t1 = -1 * B + (std::sqrt(determinant))/2*A;
            float t2 = -1 * B - (std::sqrt(determinant))/2*A;
            float t = find_t(t1, t2);
        }
        else{
            float t = -1; //no "valid" t value
        }

        
        //identifier to knwo which object was intersected: the sphere's index
    }
        ColorType color;
        color.r = .1;
        color.g = .2;
        color.b = .3;
        //how to check closest t-value between different spheres>
        //have an array of rays, an array of t values (-1 if intersection wasn't found)
        return color;
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
                //fout << word << " " << num1 << " " << num2 << " " << num3 << " " << num4 << std::endl;
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
                //fout << word << " " << num1 << " " << num2 << " " << num3 << std::endl;
            } else {
                std::cerr << "Error parsing 3-float line: " << line << std::endl;
            }
        }
        else if (word == "imsize") {
            if (ss >> num1 >> num2) {
                px_width = num1;
                px_height = num2;
                //fout << word << " " << num1 << " " << num2 << std::endl;
            } else {
                std::cerr << "Error parsing imsize line: " << line << std::endl;
            }
        }
        else if (word == "vfov") {
            if (ss >> num1) {
                if (word == "vfov") vfov = num1;
                //fout << word << " " << num1 << std::endl;
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
    //printVector(u);
    //printVector(v);

    //calculate the width and height in 3d world coordinates
    float radian_vfov = (vfov * M_PI)/180;
    int height = (2*d*std::tan(0.5 * radian_vfov));
    int width = (px_width / px_height) * height;

    //std::cout << height << std::endl;
    //std::cout << width << std::endl;

    //Find the 4 corners of the viewing window
    VectorType ul;
    ul = vectorAdd(eye, vectorScalar(viewdir, d)); //e + d*n
    ul = vectorSubtract(ul, vectorScalar(u, (width/2.0))); //(w/2.0) * u
    ul = vectorAdd(ul, vectorScalar(v, (height/2.0))); //(h/2.0) * v
    //printVector(ul);

    VectorType ur;
    ur = vectorAdd(eye, vectorScalar(viewdir, d));
    ur = vectorAdd(ur, vectorScalar(u, (width/2.0)));
    ur = vectorAdd(ur, vectorScalar(v, (height/2.0)));
    //printVector(ur);
    VectorType ll;
    ll = vectorAdd(eye, vectorScalar(viewdir, d));
    ll = vectorSubtract(ll, vectorScalar(u, (width/2.0)));
    ll = vectorSubtract(ll, vectorScalar(v, (height/2.0)));
    //printVector(ll);
    VectorType lr;
    lr = vectorAdd(eye, vectorScalar(viewdir, d));
    lr = vectorAdd(lr, vectorScalar(u, (width/2.0)));
    lr = vectorSubtract(lr, vectorScalar(v, (height/2.0)));

    //find change in h and v
    VectorType h_change = vectorDivide(vectorSubtract(ur, ul), px_width);
    //printVector(h_change);
    VectorType v_change = vectorDivide(vectorSubtract(ll, ul), px_height);
    //printVector(v_change);

    //iterate through each pixel (i, j)
    for(int i = 0; i < px_height; i++){
        for(int j = 0; j < px_width; j++){

            //the point where each ray should pass through the viewing window correspoindng to each pixel
            VectorType intersect_point = vectorAdd(vectorAdd(ul, vectorScalar(h_change, j)), vectorScalar(v_change, i));
            //initialize a ray for this pixel
            Raytype ray;
            ray.origin = eye;
            ray.intersection = intersect_point;
            normalize(ray.intersection);

            //variables to track the closest t value and corresponding object
            float t_closest = -1;
            int sphere_index;
            
            //for every pixel, iterate through each object in the scene to check for ray intersections
            for(int k = 0; k < sphereArray.size(); k++){
                float xc = sphereArray[k].x;
                float yc = sphereArray[k].y;
                float zc = sphereArray[k].z;
                float r = sphereArray[k].r;
                
                //calculate A, B and C for each sphere
                float A = pow(ray.intersection.i, 2) + pow(ray.intersection.j, 2) + pow(ray.intersection.k, 2);
                float B = 2*(ray.intersection.i*(ray.origin.i - xc) + ray.intersection.j * (ray.origin.j - yc) + ray.intersection.k*(ray.origin.k - zc));
                float C = pow((ray.origin.i - xc), 2) + pow((ray.origin.j - yc), 2) + pow((ray.origin.k - zc), 2) - pow(r, 2);

                //check for a positive determinant, if positive then calculate the t values of this object
                float determinant = (pow(B, 2) - 4*A*C);
                if(determinant > 0){
                    float t1 = -1 * B + (std::sqrt(determinant))/2.0*A;
                    float t2 = -1 * B - (std::sqrt(determinant))/2.0*A;
                    float t = find_t(t1, t2);
                    t_closest = find_t(t, t_closest); //keep track of the closest positive t of all objects
                    if(find_t(t, t_closest) == t){
                        sphere_index = k; //track which object is currently being intercepted first by the ray
                    }
                    //std::cout << "t: " << t << std::endl;
                    //std::cout << "t_closest: " << t_closest<<  std::endl;
                    //std::cout << "find_t: " << find_t(t, t_closest) << std::endl;
                    
                }
                else{
                    float t = -1;
                }
                
            }

            //std::cout << t_closest << std::endl;
            //if t_closest = -1 no valid intersection was found, use the background color for this pixel
            //if t_closest != -1 then use the material color of the intersected sphere for this pixel
            if(t_closest == -1){
                fout << bkgcolor.r*255 << " " << bkgcolor.g*255 << " " << bkgcolor.b*255 << "\n";
            }
            else{
                int mat_index = sphereArray[sphere_index].m;
                ColorType color = materialArray[mat_index];
                fout << color.r*255 << " " << color.g*255 << " " << color.b*255 << "\n";
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

