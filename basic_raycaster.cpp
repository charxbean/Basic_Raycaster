#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <cctype>
#include <algorithm>

//vector type
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

//color type
typedef struct{
    float r, g, b;
} ColorType;

//material type 
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

//depth cue 
typedef struct{
    ColorType depth_color;
    float a_min, a_max, dist_near, dist_far;
}DepthCue;

typedef struct{
    int a, b, c;
}Triangle;

//return type of the trace_ray function
typedef struct{
    int sphere_index; //intersected sphere
    Raytype ray; //ray that intersects with the sphere
    float t; //the t value of the intersection
}Intersection;

//initialize variables for the input
VectorType eye;
VectorType viewdir;
VectorType updir;
ColorType bkgcolor;
MaterialColor mtlcolor;
LightType light;
DepthCue depth;

//struct to check if all neccessary input was given
typedef struct{
    bool eye_ = false;
    bool viewdir_= false;
    bool updir_= false;
    bool bkgcolor_= false;
    bool vfov_= false;
    bool px_width_= false;
    bool px_height_= false;
}ValidInput;

ValidInput valid_input;
//boolean to check if other inputs were given
bool sphere_ = false;
bool material_ = false;
bool depth_cue = false;

//a_dc value for depth cue
float a_dc;

//array of material colors
std::vector<MaterialColor> materialArray;

//array of spheres
std::vector<SphereType> sphereArray;

//array of light sources
std::vector<LightType> lightArray;

//array of L vectors that corresponds to each light in the light array
std::vector<VectorType> L_array;

//array of input vertices
std::vector<VectorType> vertex_array;

//array of triangles
std::vector<Triangle> triangle_array;


//Returns a vector that is the cross product of vectors v1 and v2
VectorType crossProduct(const VectorType& v1, const VectorType& v2) {
    return {
        v1.j * v2.k - v1.k * v2.j,  // x
        v1.k * v2.i - v1.i * v2.k,  // y
        v1.i * v2.j - v1.j * v2.i   // z
    };
}

//Computes the dot product of 2 vectors
float dotProduct(const VectorType v1, const VectorType v2) {
    float dot = (v1.i * v2.i) + (v1.j * v2.j) + (v1.k * v2.k);
    return dot;
}

//returns a new vector that is the sum of vectors v1 and v2
VectorType vectorAdd(const VectorType v1, const VectorType v2){
    VectorType newVec;
    newVec.i = v1.i + v2.i;
    newVec.j = v1.j + v2.j;
    newVec.k = v1.k + v2.k;
    
    return newVec;
}

//returns a new vector that is v1 - v2
VectorType vectorSubtract(const VectorType v1, const VectorType v2){
    VectorType newVec;
    newVec.i = v1.i - v2.i;
    newVec.j = v1.j - v2.j;
    newVec.k = v1.k - v2.k;
    
    return newVec;
}

//returns a new vector that is vector v1 multiplied by vector v2
VectorType vectorMultiply(const VectorType v1, const VectorType v2){
    VectorType newVec;
    newVec.i = v1.i * v2.i;
    newVec.j = v1.j * v2.j;
    newVec.k = v1.k * v2.k;
    
    return newVec;
}

//returns a new vector that is v1 multiplied by the scalar
VectorType vectorScalar(const VectorType v1, const float scalar){
    VectorType newVec;
    newVec.i = v1.i * scalar;
    newVec.j = v1.j * scalar;
    newVec.k = v1.k * scalar;

    return newVec;
}

//returns a new vector that is v1 divided by the scalar
VectorType vectorDivide(const VectorType v1, const float scalar){
    VectorType newVec;
    if(scalar == 0){
        newVec = {0, 0, 0};
        return newVec;
    }
    newVec.i = v1.i / scalar;
    newVec.j = v1.j / scalar;
    newVec.k = v1.k / scalar;
    
    return newVec;
}

//gets the length of a vector
float vectorLength(const VectorType& v) {
    return std::sqrt(v.i * v.i + v.j * v.j + v.k * v.k);
}

//normalizes the vector v1
bool normalize(VectorType& v1) {
    float length = vectorLength(v1);
    
    if (length == 0.0f) {
        // The vector is zero; normalization is not possible.
        return false;
    }

    v1.i /= length;
    v1.j /= length;
    v1.k /= length;
    return true; // Successfully normalized
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

bool inside_triangle(float b, float g, float a){
    if((a >= 0 && a <= 1) && (b >= 0 && b <= 1) && (g >= 0 && g <= 1)){
        if(a + b + g == 1){
            return true;
        }
        else{
            return false;
        }
    }
    return false;
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

//returns true if two colors are equal, false if they are not equal
bool colorEquals(ColorType c1, ColorType c2){
    if(c1.r == c2.r && c1.g == c2.g && c1.b == c2.b){
        return true;
    }
    else{
        return false;
    }
}

//make sure a color's rgb values stay between 0 and 1
void clampColor(ColorType &c){
    if(c.r > 1){
        c.r = 1;
    }
    else if(c.r < 0){
        c.r = 0;
    }

    if(c.g > 1){
        c.g = 1;
    }
    else if(c.g < 0){
        c.g = 0;
    }

    if(c.b > 1){
        c.b = 1;
    }
    else if(c.b < 0){
        c.b = 0;
    }
}

//prints a vector to the terminal
void printVector(VectorType v1){
    std::cout << "<" << v1.i << ", " << v1.j << ", " << v1.k << ">" << std::endl;
}

//prints a color to the terminal
void printColor(ColorType color){
    std::cout << "r: " << color.r << ", g: " << color.g << ", b: " << color.b << std::endl;
}

//Initialize vectors N, L and H which will be computed for each sphere/object
VectorType vector_N;
VectorType vector_L;
VectorType vector_H;
VectorType vector_V;

//initialize traceRay so it can be used inside ShadeRay
Intersection traceRay(Raytype ray, int self_index);

//shadeRay function takes the index of a sphere in the sphere array, the ray that intersected it, and the t value from the intersection
//returns the color for that pixel computed using the Blinn-Phong method
ColorType shadeRay(int sphere_index, Raytype ray, float t){
    
    //if t < 0 then no object intersection was found, return the background color
    if(t < 0){
        return bkgcolor;
    }

    //find the material that corresponds with the given sphere
    int mat_index = sphereArray[sphere_index].m;
    MaterialColor material = materialArray[mat_index];
    SphereType sphere = sphereArray[sphere_index];

    //the intersection point of the sphere from the eye
    VectorType intersection = vectorAdd(ray.origin, vectorScalar(ray.intersection, t));

    //Calculate illumination 
    //calculate vector N
    vector_N = vectorDivide(vectorSubtract(intersection, sphere.center), sphere.r);
    normalize(vector_N);

    //calculate vector L for each point light source
    for(int k = 0; k < lightArray.size(); k++){
        if(lightArray[k].type == 1){
            //(light_position - surface position)
            VectorType new_L = vectorSubtract(lightArray[k].position, intersection);
            //calculate vector L
            vector_L = vectorDivide(new_L, vectorLength(new_L));
            normalize(vector_L);
            //add into vector L array corresponding to the index of this light in the light array
            L_array.insert(L_array.begin() + k, vector_L);
        }
    }

    //Define V
    vector_V = vectorDivide(vectorSubtract(ray.origin, ray.intersection), vectorLength(vectorSubtract(ray.origin, ray.intersection)));
    normalize(vector_V);
    
    //initialize intensity to the ambient light
    ColorType illumination = colorMultiply(material.Od, material.ka); //ka + Od

    //for each light source, compute it's individual intensity and add to the total intensity value
    for(int k = 0; k < lightArray.size(); k++){
        //define vector H for each light source
        vector_H = vectorDivide(vectorAdd(L_array[k], vector_V), 2); //H = (L + V) /2
        vector_H = vectorDivide(vector_H, vectorLength(vector_H)); // H = H / ||H||
        normalize(vector_H);

        //calculate the diffuse color
        ColorType diffuse = colorMultiply(material.Od, material.kd); //kd * od
        float max_NL = max(0.0, dotProduct(vector_N, L_array[k])); //clamp diffuse at 0
        diffuse = colorMultiply(diffuse, max_NL);

        //calculate specular color
        ColorType specular = colorMultiply(material.Os, material.ks); //ks * Os
        float max_NH = max(0.0, dotProduct(vector_N, vector_H)); //clamp specular at 0
        specular = colorMultiply(specular, pow(max_NH, material.n)); //(ks * Os) * (N dot H)^n
        
        //distance from the sphere to the light source
        float light_dist = vectorLength(vectorSubtract(lightArray[k].position, intersection));
        
        //Shadows - for each light, cast a ray to it and determine if it's in shadow(0) or not (1)
        int shadow_flag;
        Raytype shadow_ray;
        shadow_ray.origin = intersection; //the surface intersection
        shadow_ray.intersection = vectorSubtract(lightArray[k].position, shadow_ray.origin); //light position from the sphere surface
        normalize(shadow_ray.intersection);
        
        //use traceRay to send a ray from the surface to the light source
        Intersection trace_ray = traceRay(shadow_ray, sphere_index);
        if(trace_ray.t < 0){
            shadow_flag = 1; //no intersection found, no shadow
        }
        else{
            float blocking_dist = trace_ray.t;
            //if light is directional, the object blocks the light
            if(lightArray[k].type == 0){
                shadow_flag = 0;
            }//if t is greater than the distance from the sphere to a point light source, the object doesn't block the light
            else if(blocking_dist > light_dist){
                shadow_flag = 1;
            }
            else{//point light and the object is blocking the light
                shadow_flag = 0;
            }
        }

        //if there is no shadow, add the illumination from this light source to the total
        if(shadow_flag == 1){
            ColorType illum_L = colorMultiply(colorAdd(diffuse, specular), lightArray[k].intensity);
            illumination = colorAdd(illumination, illum_L);
            clampColor(illumination);
        }

    }

    //if depth cue was input, compute depth cue
    if(depth_cue == true){
        if(t <= depth.dist_near){
            a_dc = depth.a_max;
        }
        else if(t > depth.dist_far){
            a_dc = depth.a_min;
        }
        else{
            a_dc = depth.a_max - (depth.a_max - depth.a_min) * ((depth.dist_far - t)/(depth.dist_far - depth.dist_near));
        }
        //compute illumination with depth cue
        illumination = colorAdd(colorMultiply(illumination, a_dc), colorMultiply(depth.depth_color, (1-a_dc)));
    }

    clampColor(illumination);
    return illumination;
}

ColorType oldShadeRay(int sphere_index, int t){
    //std::cout << sphere_index << " " << t << std::endl;
    if(t < 0){
        return bkgcolor;
    }
    else{
        int mat_index = sphereArray[sphere_index].m;
        MaterialColor material = materialArray[mat_index];
        ColorType color = material.Od;
        return color;
    }
}

//checks each sphere in the scene for a ray intersection and determines the closest valid intersection
//skips the input self index if it's a valid index number
//returns a struct containting the closest t-value, the intersected sphere index and the ray
Intersection traceRay(Raytype ray, int self_index){

    //variables to track the closest t value and corresponding sphere/object
    float t_closest = -1;
    int sphere_index = -1;
    
    //iterate through each object in the scene to check for ray intersections
    for(int k = 0; k < sphereArray.size(); k++){
        //when tracing a shadow ray, if the sphere is itself, don't look for an intersection
        if(self_index == k){    
            continue;
        }
        float xc = sphereArray[k].center.i;
        float yc = sphereArray[k].center.j;
        float zc = sphereArray[k].center.k;
        float r = sphereArray[k].r;
        
        //calculate A, B and C for each sphere
        float A = 1.0;
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
            float t = (-B /(2.0*A));
            t_closest = find_t(t, t_closest); //keep track of the closest positive t of all objects
            if(t_closest == t){
                sphere_index = k; //keep track of the sphere index corresponding to the closest t value
            }
        }
        else{ //if no valid t value is found, t = -1
            float t = -1;
        }
    }

    Intersection traceRay_intersection = {sphere_index, ray, t_closest};
    return traceRay_intersection;
}

//checks each triangle in the sccene for a ray intersection

//goes through each triangle in the triangle array
Intersection trace_polygon(Raytype ray){
    //p0, p1 and p2 are a b and c of the traingle type
    //define e1, e2 and n
    float t_closest = -1;
    int triangle_index = -1;

    //n = p1 - p0 X p2 - p0
    VectorType p0;
    VectorType p1;
    VectorType p2;

    VectorType e1;
    VectorType e2;

    for(int k = 0; k < triangle_array.size(); k++){
        //Define p0, p1 and p2 of the triangle
        p0 = vertex_array[triangle_array[k].a];
        p1 = vertex_array[triangle_array[k].b];
        p2 = vertex_array[triangle_array[k].c];

        e1 = vectorSubtract(p1, p0);
        e2 = vectorSubtract(p2, p0);

        //Define vector n
        VectorType n = crossProduct(e1, e2);

        //Define A, B, C and D
        float A = n.i;
        float B = n.j;
        float C = n.k;
        float D = -1 * ((A*p0.i) + (B*p0.j) + (C*p0.k));

        //find the denominator of the possible ray intersection
        //xd, yx and zd = ray intersection
        float denominator = (A*ray.intersection.i) + (B*ray.intersection.j) + (C * ray.intersection.k);
        if(denominator == 0){
            float t = -1;
        }
        else{
            float t = (-1 * (A * ray.origin.i + B*ray.origin.j + C*ray.origin.k))/denominator;
            t_closest = find_t(t, t_closest);
            if(t_closest == t){
                triangle_index = k;
            }
        }
    }

    //if the ray didn't intersect with any plane, end here
    if(t_closest == -1){
        return {-1, ray, -1};
    }

    //if there was an intersection, check if it intersects inside the triangle
    //MODIFY SO YOU HAVE THE RIGHT P0, E1 and E2!!!!
    p0 = vertex_array[triangle_array[triangle_index].a];
    p1 = vertex_array[triangle_array[triangle_index].b];
    p2 = vertex_array[triangle_array[triangle_index].c];
    e1 = vectorSubtract(p1, p0);
    e2 = vectorSubtract(p2, p0);

    VectorType p = {(ray.origin.i + (t_closest * ray.intersection.i)), (ray.origin.j + (t_closest * ray.intersection.j)), (ray.origin.k + (t_closest * ray.intersection.k))};
    VectorType ep = vectorSubtract(p, p0);

    float d11 = dotProduct(e1, e1);
    float d22 = dotProduct(e2, e2);
    float d12 = dotProduct(e1, e2);
    float d1p = dotProduct(e1, ep);
    float d2p = dotProduct(e2, ep);
    
    float determinant = ((d11 * d22) - (d12*d12));
    if(determinant == 0){
        return {-1, ray, -1};
        std::cerr << "erm" << std::endl;
    }
    else{
        float beta = (d11*d1p - d12 * d2p)/determinant;
        float gamma = (d11 * d2p - d12 * d1p)/determinant;
        float alpha = 1 - (beta + gamma);
        //intersection with plane was found AND ray is inside triangle
        if(inside_triangle(beta, gamma, alpha)){
            Intersection polygon_intersection = {triangle_index, ray, t_closest};
            return polygon_intersection;
        }
        else{
            //no intersection with triangle, t and index are -1
            return {-1, ray, -1};
        }
    }

    return {-1, ray, -1};
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
                material_ = true;
            }
            else{
                std::cerr << "Error parsing mtcolor line: " << line << std::endl;
            }
        }
        else if (word == "depthcueing"){
            if(ss >> num1 >> num2 >> num3 >> num4 >> num5 >> num6 >> num7){
                ColorType depth_color = {num1, num2, num3};
                depth = {depth_color, num4, num5, num6, num7};
                depth_cue = true;
            }
            else{
                std::cerr << "Error parsing depthcueing line: " << line << std::endl;
            }
        }
        else if (word == "light"){
            if (ss >> num1 >> num2 >> num3 >> num4 >> num5){
                //check that the 4th input number is either 0 (directional light) or 1 (point light)
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
                sphere_ = true;
            } else {
                std::cerr << "Error parsing sphere line: " << line << std::endl;
            }
        }
        else if (word == "eye" || word == "viewdir" || word == "updir" || word == "bkgcolor" || word == "v" || word == "f") {
            if (ss >> num1 >> num2 >> num3) {
                if (word == "eye"){
                    eye = {num1, num2, num3};
                    valid_input.eye_ = true;
                }         
                if (word == "viewdir"){
                    viewdir = {num1, num2, num3};
                    valid_input.viewdir_ = true;
                }     
                if (word == "updir"){
                    updir = {num1, num2, num3};
                    valid_input.updir_ = true;
                }       
                if (word == "bkgcolor"){
                    bkgcolor = {num1, num2, num3};
                    valid_input.bkgcolor_ = true;
                }
                if(word == "v"){ 
                    VectorType v = {num1, num2, num3};
                    vertex_array.push_back(v);
                }
                if(word == "f"){ //make sure first triangle index starts at 1, not 0??
                    Triangle f = {static_cast<int>(num1), static_cast<int>(num2), static_cast<int>(num3)};
                    triangle_array.push_back(f);
                }
            } else {
                std::cerr << "Error parsing " << word << " line: " << line << std::endl;
            }
        }
        else if (word == "imsize") {
            if (ss >> num1 >> num2) {
                px_width = num1;
                px_height = num2;
                valid_input.px_height_ = true;
                valid_input.px_width_ = true;
            } else {
                std::cerr << "Error parsing imsize line: " << line << std::endl;
            }
        }
        else if (word == "vfov") {
            if (ss >> num1) {
               vfov = num1;
               valid_input.vfov_ = true;
            } else {
                std::cerr << "Error parsing vfov line: " << line << std::endl;
            }
        }
        else {
            std::cerr << "Unrecognized keyword: " << word << " in line: " << line << std::endl;
        }
    }

    //check that all necessary input was given
    if(!valid_input.eye_ || !valid_input.updir_|| !valid_input.viewdir_|| !valid_input.vfov_){
        std::cerr << "Error: one or more viewer inputs (eye, view direction, updirection or vfov) was not given" << std::endl;
        return -1;
    }
    else if(!valid_input.bkgcolor_){
        std::cerr << "Error: No background color input" << std::endl;
        return -1;
    }
    else if(!valid_input.px_height_ || !valid_input.px_width_){
        std::cerr << "Error: imsize width and height were not properly input" << std::endl;
        return -1;
    }

    printVector(vertex_array[3]);
    std::cout << triangle_array[3].a << std::endl;

    //create the header for the ppm file
    fout << "P3\n";
    fout << px_width << " # image width\n" << px_height << " # image height\n";
    fout << "255\n";

    //if there is no input object, just return the background color
    if(!sphere_){
        for(int i = 0; i < px_height; i++){
            for(int j = 0; j < px_width; j++){
                fout << std::round(bkgcolor.r * 255) << " " << std::round(bkgcolor.g * 255) << " " << std::round(bkgcolor.b * 255) << std::endl;
            }
        }
        inputFile.close();
        return 0;
    }//check that if there is a sphere, there is also a material color
    else if(sphere_ && !material_){
        std::cerr << "a sphere was input without a material color" << std::endl;
        inputFile.close();
        return -1;
    }
    
    //PRELIMINARY MATH

    //find vectors u and v
    VectorType u = crossProduct(viewdir, updir);
    normalize(u);
    VectorType v = crossProduct(u, viewdir);

    //calculate the width and height in 3d world coordinates
    float radian_vfov = (vfov * M_PI)/180;
    int height = (2*d*std::tan(0.5 * radian_vfov));
    float aspect = (float)px_width / (float)px_height;
    int width = aspect * height;

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
            VectorType L_inv = vectorScalar(lightArray[k].position, -1.0);
            if(L_inv.i == 0 && L_inv.j == 0 && L_inv.k == 0){ //if zero vector, dont divide by 0
                L_array.insert(L_array.begin() + k, lightArray[k].position);
            }
            else{
                vector_L = vectorDivide(L_inv, vectorLength(L_inv));
                 //insert this L vector into the L array at the same index the light is at in light array
                L_array.insert(L_array.begin() + k, vector_L);
            }
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
            if(!normalize(ray.intersection)){
                std::cerr << "Unable to normalize" << std::endl;
            }
    
            //use trace_ray and shade ray to determine the color for each pixel
            Intersection trace_ray = traceRay(ray, -1);
            ColorType color = shadeRay(trace_ray.sphere_index, trace_ray.ray, trace_ray.t);  
            fout << std::round(color.r * 255) << " " << std::round(color.g * 255) << " " << std::round(color.b * 255) << std::endl;
            
        }
    }

    inputFile.close();
    return 0;
}