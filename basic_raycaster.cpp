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

typedef struct{
    ColorType Od;
    ColorType Os;
    float ka, kd, ks;
    int n;
}MaterialColor;

//sphere defined by it's center point (x, y, z) and a radius
//has m as it's material property -> an index into an array of material properties
typedef struct{
    VectorType center;
    float r;
    int m;
} SphereType;

//equation of a ray = (x, y, z) - t*(dx, dy, dz)
typedef struct{
    VectorType origin; //(x, y, z)
    VectorType intersection; //(dx, dy, dz)
} Raytype;

//type 1 = point light 
//type 0 = directional light
typedef struct{
    VectorType position;
    int type;
    float intensity;
}LightType;

typedef struct{
    ColorType depth_color;
    float a_min, a_max, dist_near, dist_far;
}DepthCue;

typedef struct{
    LightType light;
    float c1, c2, c3;
}AttlightType;

//initialize variables for the input
VectorType eye;
VectorType viewdir;
VectorType updir;
ColorType bkgcolor;
MaterialColor mtlcolor;
LightType light;
DepthCue depth;
AttlightType attlight;

//array of material colors
std::vector<MaterialColor> materialArray;

//array of spheres
std::vector<SphereType> sphereArray;

//array of light sources
std::vector<LightType> lightArray;

//array of L vectors that corresponds to each light in the light array
std::vector<VectorType> L_array;

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

//Computes the dot product of 2 vectors
float dotProduct(VectorType v1, VectorType v2){
    float dot = (v1.i * v2.i) + (v1.j * v2.j) + (v1.k * v2.k);
    return dot;
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

//Multiplies each channel in a color by a float value to adjust the intensity
ColorType colorMultiply(ColorType color, float intensity){
    ColorType newColor;
    newColor.r = color.r * intensity;
    newColor.g = color.g * intensity;
    newColor.b = color.b * intensity;
    return newColor;
}

//add the rgb of 2 colors together
ColorType colorAdd(ColorType c1, ColorType c2){
    ColorType newColor;
    newColor.r = c1.r + c2.r;
    newColor.g = c1.g + c2.g;
    newColor.b = c1.b + c2.b;
    return newColor;
}

//gets the length of a vector
float vectorLength(VectorType v1){
    float length = std::sqrt(pow(v1.i, 2) + pow(v1.j, 2) + pow(v1.k, 2));
    return length;
}

//normalizes the vector v1
void normalize(VectorType& v1){
    float length = vectorLength(v1);
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

//returns the max of 2 float values
float max(float a, float b){
    if(a > b){
        return a;
    }
    else{
        return b;
    }
}

//prints a vector to the terminal
void printVector(VectorType v1){
    std::cout << "<" << v1.i << ", " << v1.j << ", " << v1.k << ">" << std::endl;
}

void printColor(ColorType color){
    std::cout << "r: " << color.r << ", g: " << color.g << ", b: " << color.b << std::endl;
}

//Initialize vectors N, L and H which will be computed for each sphere/object
VectorType vector_N;
VectorType vector_L;
VectorType vector_H;
VectorType vector_V;

//takes the index of a sphere in the sphere array and returns the material color for that sphere
ColorType shadeRay(int sphere_index, Raytype ray, float t){
    int mat_index = sphereArray[sphere_index].m;
    MaterialColor material = materialArray[mat_index];
    SphereType sphere = sphereArray[sphere_index];
    VectorType intersection = vectorAdd(ray.origin, vectorScalar(ray.intersection, t));

    //Calculate illumination 
    //calculate N, L and H
    vector_N = vectorDivide(vectorSubtract(intersection, sphere.center), sphere.r);
    normalize(vector_N);

    //calculate L for each point light source
    for(int k = 0; k < lightArray.size(); k++){
        if(lightArray[k].type == 1){
            //calculate light_position - surface position
            VectorType new_L = vectorSubtract(lightArray[k].position, intersection);
            //calculate vector L
            vector_L = vectorDivide(new_L, vectorLength(new_L));
            normalize(vector_L);
            L_array.insert(L_array.begin() + k, vector_L);
        }
    }
    //Define V
    vector_V = vectorDivide(vectorSubtract(ray.origin, ray.intersection), vectorLength(vectorSubtract(ray.origin, ray.intersection)));
    normalize(vector_V);
    
    //initialize intensity
    ColorType illumination = colorMultiply(material.Od, material.ka); //ka + Od

    //for each light source, compute it's intensity and add to the total intensity value
    for(int k = 0; k < lightArray.size(); k++){
        //define H for each light
        vector_H = vectorDivide(vectorAdd(L_array[k], vector_V), 2); //H = (L + V) /2
        vector_H = vectorDivide(vector_H, vectorLength(vector_H)); // H = H / ||H||
        normalize(vector_H);

        ColorType diffuse = colorMultiply(material.Od, material.kd); //kd * od
        float max_NL = max(0.0, dotProduct(vector_N, L_array[k]));
        diffuse = colorMultiply(diffuse, max_NL);

        ColorType specular = colorMultiply(material.Os, material.ks); //ks * Os
        float max_NH = max(0.0, dotProduct(vector_N, vector_H));
        specular = colorMultiply(specular, pow(max_NH, material.n)); //(ks * Os) * (N dot H)^n
        
        //intensity += IL(diffuse + specular)
        illumination = colorAdd(illumination, colorMultiply(colorAdd(diffuse, specular), lightArray[k].intensity));
        
    }
    //printColor(illumination);
    return illumination;
}


//checks each object in the scene for a ray intersection and determines the closest valid intersection
//either returns the color of the object it intersects or returns the background color if no object is intersected
ColorType traceRay(Raytype ray){

    //variables to track the closest t value and corresponding sphere/object
    float t_closest = -1;
    int sphere_index;
    
    //iterate through each object in the scene to check for ray intersections
    for(int k = 0; k < sphereArray.size(); k++){
        float xc = sphereArray[k].center.i;
        float yc = sphereArray[k].center.j;
        float zc = sphereArray[k].center.k;
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
        return shadeRay(sphere_index, ray, t_closest);
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
        float num1, num2, num3, num4, num5, num6, num7, num8, num9;
        int num10;

        ss >> word; //Read the first word of every line

        if (word.empty()) {
            std::cerr << "Error: Empty keyword in line: " << line << std::endl;
            continue;
        }

        if (word == "mtlcolor"){
            if (ss >> num1 >> num2 >> num3 >> num4 >> num5 >> num6 >> num7 >> num8 >> num9 >> num10){
                ColorType Od = {num1, num2, num3};
                ColorType Os = {num4, num5, num6};
                mtlcolor = {Od, Os, num7, num8, num9, num10};
                materialArray.push_back(mtlcolor); //add new material to the material array
            }
            else{
                std::cerr << "Error parsing mtcolor line: " << line << std::endl;
            }
        }
        else if (word == "attlight"){
            if(ss >> num1 >> num2 >> num3 >> num4 >> num5 >> num6 >> num7 >> num8){
                if(num4 != 1 && num4 != 0){
                    std::cerr << "Error in attlight line, input 4: light type must be 0 (directional) or 1 (point)" << std::endl;
                }
                LightType a_light = {num1, num2, num3, static_cast<int>(num4), num5};
                attlight = {a_light, num6, num7, num8};
            }
            else{
                std::cerr << "Error parsing attlight line: " << line << std::endl;
            }
        }
        else if (word == "depthcueing"){
            if(ss >> num1 >> num2 >> num3 >> num4 >> num5 >> num6 >> num7){
                ColorType depth_color = {num1, num2, num3};
                depth = {depth_color, num4, num5, num6, num7};
            }
            else{
                std::cerr << "Error parsing depthcueing line: " << line << std::endl;
            }
        }
        else if (word == "light"){
            if (ss >> num1 >> num2 >> num3 >> num4 >> num5){
                if(num4 != 1 && num4 != 0){
                    std::cerr << "Error in light line, input 4: light type must be 0 (directional) or 1 (point)" << std::endl;
                }
                VectorType position = {num1, num2, num3}; 
                light = {position, static_cast<int>(num4), num5};
                lightArray.push_back(light);
            }
            else{
                std::cerr << "Error parsing light line: " << line << std::endl;
            }
        }
        else if (word == "sphere") {
            if (ss >> num1 >> num2 >> num3 >> num4) {
                //set sphere material to the index of the last material in the material array
                VectorType center = {num1, num2, num3};
                SphereType sphere = {center, num4, static_cast<int>(materialArray.size() - 1)};
                sphereArray.push_back(sphere); //add spheres to the spere array
            } else {
                std::cerr << "Error parsing sphere line: " << line << std::endl;
            }
        }
        else if (word == "eye" || word == "viewdir" || word == "updir" || word == "bkgcolor") {
            if (ss >> num1 >> num2 >> num3) {
                if (word == "eye")         eye = {num1, num2, num3};
                if (word == "viewdir")     viewdir = {num1, num2, num3};
                if (word == "updir")       updir = {num1, num2, num3};
                if (word == "bkgcolor")    bkgcolor = {num1, num2, num3};
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

    //calculate the vector L for Blinn model of the directional lights
    for(int k = 0; k < lightArray.size(); k++){
        if(lightArray[k].type == 0){
            VectorType L_inv = {-vector_L.i, -vector_L.j, -vector_L.k};
            vector_L = vectorDivide(L_inv, vectorLength(L_inv));
            normalize(vector_L);
            //insert this L vector into the L array at the same index the light is at in light array
            L_array.insert(L_array.begin() + k, vector_L);
        }
    }

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
            fout << std::round(color.r * 255) << " " << std::round(color.g * 255) << " " << std::round(color.b * 255) << std::endl;
        }
    }

    inputFile.close();
    return 0;
}