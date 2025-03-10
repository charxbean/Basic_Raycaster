#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <cctype>
#include <algorithm>

#define _USE_MATH_DEFINES

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
//and t as it's texture -> also an index into an array of textures
typedef struct{
    VectorType center;
    float r;
    int m;
    int t;
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
    int n1, n2, n3;
    int t1, t2, t3;
    int material, texture;
}Triangle;

//return type of the trace_ray and trace_polygon functions
typedef struct{
    int shape_index; //index of the intersected object in it's object array
    int shape; //shape type, 0 for sphere, 1 for triangle
    Raytype ray; //ray that intersects with the sphere
    float t; //the t value of the intersection
}Intersection;

//defines a point P's beta, gamma and alpha
typedef struct{
    float beta, gamma, alpha;
}BGA;

//contains the 2D array of a texture files colors and the texture image width/height
typedef struct{
    int success; //0 if the texture file was read successfully, -1 otherwise
    int width, height;
    std::vector<std::vector<ColorType>> array;
}TextureFile;

//struct of texture coordinate float values u and v
typedef struct{
    float u, v;
}TextureCoords;

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
bool triangle_ = false;
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

//array of normal vectors
std::vector<VectorType> normal_array;

//array of texture coordinates (u,v)
std::vector<TextureCoords> texture_coords;

//array of TextureFile types which store texture 2D arrays and width/height
std::vector<TextureFile> texture_array;


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

//Check the input barycentric coordinates to determine if they meet the conditions for being inside a triangle
//if conditions are met, return true, otherwise return false
bool inside_triangle(float b, float g, float a){

    if((a + b + g > .9) && (a + b + g < 1.09)){
        if((a >= 0 && a <= 1) && (b >= 0 && b <= 1) && (g >= 0 && g <= 1)){
            return true;
        }
        else{
            return false;
        }
    }
    else{
        return false;
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

BGA get_bga(int triangle_index, Raytype ray, float t){
    VectorType p0 = vertex_array[(triangle_array[triangle_index].a)-1];
    VectorType p1 = vertex_array[(triangle_array[triangle_index].b)-1];
    VectorType p2 = vertex_array[(triangle_array[triangle_index].c)-1];

    VectorType e1 = vectorSubtract(p1, p0);
    VectorType e2 = vectorSubtract(p2, p0);

    VectorType p = {(ray.origin.i + (t * ray.intersection.i)), (ray.origin.j + (t * ray.intersection.j)), (ray.origin.k + (t * ray.intersection.k))};
    VectorType ep = vectorSubtract(p, p0);

    float d11 = dotProduct(e1, e1);
    float d22 = dotProduct(e2, e2);
    float d12 = dotProduct(e1, e2);
    float d1p = dotProduct(e1, ep);
    float d2p = dotProduct(e2, ep);
    
    float determinant = ((d11 * d22) - (d12*d12));
    if(determinant == 0){
        return {-1, -1, -1};
        std::cerr << "erm" << std::endl;
    }
    else{
        float beta = (d22*d1p - d12*d2p)/determinant;
        float gamma = (d11 * d2p - d12 * d1p)/determinant;
        float alpha = 1 - (beta + gamma);
        return {beta, gamma, alpha};
    }
}

VectorType interpolate_t_normal(int triangle_index, Raytype ray, float t){
    if(triangle_array[triangle_index].n1 == -1){
        return {-1, -1, -1};
        std::cerr << "Trying to interpolate normals, but the triangle has no input normal values" << std::endl;
    }
    VectorType interpolated_n;
    BGA bga = get_bga(triangle_index, ray, t);
    VectorType n0 = normal_array[(triangle_array[triangle_index].n1) -1];
    VectorType n1 = normal_array[(triangle_array[triangle_index].n2) -1];
    VectorType n2 = normal_array[(triangle_array[triangle_index].n3) -1];
        //n = n0*alpha + n1*beta + n2*gamma
    interpolated_n = vectorAdd(vectorAdd(vectorScalar(n0, bga.alpha), vectorScalar(n1, bga.beta)), vectorScalar(n2, bga.gamma));
    normalize(interpolated_n);
    return interpolated_n;
}

//initialize traceRay and tracePolygon so they can be used inside ShadeRay
Intersection traceRay(Raytype ray, int self_index);
Intersection tracePolygon(Raytype ray, int self_index);

//shadeRay function takes the index of a sphere in the sphere array, the ray that intersected it, and the t value from the intersection
//returns the color for that pixel computed using the Blinn-Phong method
ColorType shadeRay(int shape_index, int shape, Raytype ray, float t){
    
    //if t < 0 then no object intersection was found, return the background color
    if(t < 0){
        return bkgcolor;
    }

    //initialize material color
    MaterialColor material;
    ColorType texture_color = {-1, -1, -1};

    //the intersection point of the surface from the eye
    VectorType intersection = vectorAdd(ray.origin, vectorScalar(ray.intersection, t));

    if(shape == 0){ //the shape is a sphere
        //find the material that corresponds with the given sphere
        int mat_index = sphereArray[shape_index].m;
        material = materialArray[mat_index];
        SphereType sphere = sphereArray[shape_index];

        //calculate normal vector N
        vector_N = vectorDivide(vectorSubtract(intersection, sphere.center), sphere.r);
        normalize(vector_N);

        //if the sphere has a texture, calculate u and v to get the texture color
        if(sphere.t != -1){ 
            TextureFile texture = texture_array[sphere.t];
            float phi = std::acos((intersection.k - sphere.center.k)/sphere.r);
            float theta = std::atan2((intersection.j - sphere.center.j), (intersection.i - sphere.center.i));

            float v = phi/M_PI;
            float u;
            if(theta > 0){
                u = theta/(2*M_PI);
            }
            else{
                u = (theta + 2*M_PI)/(2*M_PI);
            }

            int i = std::round(u*(texture.width-1));
            int j = std::round(v*(texture.height-1));

            texture_color = texture.array[i][j];
        }
    }
    else if(shape == 1){ //the shape is a triangle
        Triangle triangle = triangle_array[shape_index];
        material = materialArray[triangle.material];

        if(triangle.n1 == -1){//calculate the normal of the triangle(FLAT SHADING)
            VectorType e1 = vectorSubtract(vertex_array[(triangle.b)-1], vertex_array[(triangle.a)-1]);
            VectorType e2 = vectorSubtract(vertex_array[(triangle.c)-1], vertex_array[(triangle.a)-1]);

            vector_N = crossProduct(e1, e2);
            return material.Od;
        }
        else if(triangle.n1 > 0){ //calculate normal for SMOOTH SHADING
            vector_N = interpolate_t_normal(shape_index, ray, t);
        }
        
        //TEXTURE MAPPING
        //if a texture and texture coords were input, find u,v to lookup Od color
        if(triangle.t1 > 0){
            TextureFile texture = texture_array[triangle.texture];
            BGA bga = get_bga(shape_index, ray, t);
            float u0 = texture_coords[(triangle.t1)-1].u;
            float u1 = texture_coords[(triangle.t2)-1].u;
            float u2 = texture_coords[(triangle.t3)-1].u;

            float v0 = texture_coords[(triangle.t1)-1].v;
            float v1 = texture_coords[(triangle.t2)-1].v;
            float v2 = texture_coords[(triangle.t3)-1].v;

            float u = bga.alpha*u0 + bga.beta*u1 + bga.gamma*u2;
            float v = bga.alpha*v0 + bga.beta*v1 + bga.gamma*v2;

            int i = std::round(u*(texture.width-1));
            int j = std::round(v*(texture.height-1));

            //change the od color to the pixel at (i,j) in the texture image
            texture_color = texture.array[i][j];
        }
    }
    else{
        std::cerr << "Attempting to shade an undefined shape" << std::endl;
        return bkgcolor;
    }

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
    
    //initialize illumination to the ambient light
    //if theres a texture, switch od with texture color
    ColorType illumination;
    if(texture_color.r == -1){ //no texture
        illumination = colorMultiply(material.Od, material.ka); //ka + Od
    }
    else{
        illumination = colorMultiply(texture_color, material.ka); //ka + texture_color
    }
    

    //for each light source, compute it's individual intensity and add to the total intensity value
    for(int k = 0; k < lightArray.size(); k++){
        //define vector H for each light source
        vector_H = vectorDivide(vectorAdd(L_array[k], vector_V), 2); //H = (L + V) /2
        vector_H = vectorDivide(vector_H, vectorLength(vector_H)); // H = H / ||H||
        normalize(vector_H);

        //calculate the diffuse color, if theres a texture, Od = texture color
        ColorType diffuse;
        if(texture_color.r == -1){ //no texture
            diffuse = colorMultiply(material.Od, material.kd); //kd * od
        }
        else{ //texture
            diffuse = colorMultiply(texture_color, material.kd); //kd * texture_color
        }
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
        
        //send a ray from the surface of the object to the light source
        Intersection trace_ray;
        if(shape == 0){ //use trace_ray if the object is a sphere
            trace_ray = traceRay(shadow_ray, shape_index);
        }
        else if(shape == 1){ //use trace_polygon if the object is a triangle
            trace_ray = tracePolygon(shadow_ray, shape_index);
        }
        else{
            std::cerr << "Attempting to trace undefined shape" << std::endl;
            trace_ray = {-1, 0, shadow_ray, -1};
        }

        //use the t value from trace_ray to determine if there is a shadow or not
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

    Intersection traceRay_intersection = {sphere_index, 0, ray, t_closest};
    return traceRay_intersection;
}

//checks each triangle in the sccene for a ray intersection

//goes through each triangle in the triangle array
Intersection tracePolygon(Raytype ray, int self_index){
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
        if(self_index == k){
            continue;
        }
        //std::cout << k << std::endl;
        //Define p0, p1 and p2 of the triangle
        p0 = vertex_array[(triangle_array[k].a)-1];
        //printVector(p0);
        p1 = vertex_array[(triangle_array[k].b)-1];
        p2 = vertex_array[(triangle_array[k].c)-1];

        e1 = vectorSubtract(p1, p0);
        //printVector(e1);
        e2 = vectorSubtract(p2, p0);

        //Define vector n
        VectorType n = crossProduct(e1, e2);
        //printVector(n);

        //Define A, B, C and D
        float A = n.i;
        float B = n.j;
        float C = n.k;
        float D = -1 * ((A*p0.i) + (B*p0.j) + (C*p0.k));
       // printVector(n);
        //printVector(ray.origin);
        
        //find the denominator of the possible ray intersection
        //xd, yx and zd = ray intersection
        float denominator = (A*ray.intersection.i) + (B*ray.intersection.j) + (C * ray.intersection.k);
        //printf("%f\n", denominator);
        if(denominator == 0){
            float t = -1;
        }
        else{
            float t = (-1 * (A * ray.origin.i + B*ray.origin.j + C*ray.origin.k))/denominator;
            
            VectorType p = {(ray.origin.i + (t * ray.intersection.i)), (ray.origin.j + (t * ray.intersection.j)), (ray.origin.k + (t * ray.intersection.k))};
            VectorType ep = vectorSubtract(p, p0);
        
            float d11 = dotProduct(e1, e1);
            float d22 = dotProduct(e2, e2);
            float d12 = dotProduct(e1, e2);
            float d1p = dotProduct(e1, ep);
            float d2p = dotProduct(e2, ep);
            
            float determinant = ((d11 * d22) - (d12*d12));
            if(determinant == 0){
                return {-1, 1, ray, -1};
                std::cerr << "erm" << std::endl;
            }
            else{
                float beta = (d22*d1p - d12*d2p)/determinant;
                float gamma = (d11 * d2p - d12 * d1p)/determinant;
                float alpha = 1 - (beta + gamma);
        
                //intersection with plane was found AND ray is inside triangle
                if(inside_triangle(beta, gamma, alpha)){
                    t_closest = find_t(t, t_closest);
                    if(t_closest == t){
                        triangle_index = k;
                    }
                    //Intersection polygon_intersection = {triangle_index, 1, ray, t};
                    //return polygon_intersection;
                }
            }
        }
    }
    //if the ray didn't intersect with any plane, end here
    if(t_closest == -1){
        return {-1, 1, ray, -1};
    }
    else{
        return {triangle_index, 1, ray, t_closest};
    }
}

//reads the data in the input texture file into a 2D array of colors, returns the 2D color array
std::vector<std::vector<ColorType>>readTextureFile(int width, int height, std::ifstream& file){

    std::vector<std::vector<ColorType>> array(height, std::vector<ColorType>(width));
    for(int i = 0; i < width ; i++){
        for(int j = 0; j < height; j++){
            float num1, num2, num3;
            if (file >> num1 >> num2 >> num3) {
                array[i][j] = {num1 / 255.0f, num2 / 255.0f, num3 / 255.0f}; // Normalize RGB
            } else {
                std::cerr << "Error reading texture pixel at (" << i << ", " << j << ")" << std::endl;
                return {}; //empty array to signal failure
            }
        }
    }
    return array;
}

//opens the texture ppm file and checks for a correct header
//returns a TextureFile type containting a 2D array and a success int: -1 on failure, 0 on success
TextureFile openTextureFile(std::string filename){
    std::ifstream file(filename);
    if(!file.is_open()){
        std::cerr << "Error: unable to open file: " << filename << std::endl;
        std::vector<std::vector<ColorType>> array;
        file.close();
        return {-1, -1, -1, array};
    }
    std::string line;
    std::getline(file, line);
    std::istringstream ss(line);
    std::string format;
    int width, height, max_color;

    if(ss >> format >> width >> height >> max_color){
        if(format == "P3" && max_color == 255){
            //std::cout << width << " " << height << std::endl;
            std::vector<std::vector<ColorType>> array = readTextureFile(width, height, file);
            file.close();
            return {0, width, height, array};
        }
        else{
            std::cerr << "Incorrect header for texture ppm file" << line << std::endl;
            std::vector<std::vector<ColorType>> array;
            file.close();
            return {-1, -1, -1, array};
        }
    }
    else{
        std::cerr << "Incorrect header for texture ppm file" << line << std::endl;
        std::vector<std::vector<ColorType>> array;
        file.close();
        return {-1, -1, -1, array};
    }
    
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
        std::string textureFile;

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
                if(texture_array.size() > 0){
                    SphereType sphere = {center, num4, static_cast<int>(materialArray.size() - 1), static_cast<int>(texture_array.size() - 1)};
                    sphereArray.push_back(sphere); //add spheres to the spere array
                    sphere_ = true;
                }
                else{
                    SphereType sphere = {center, num4, static_cast<int>(materialArray.size() - 1), -1};
                    sphereArray.push_back(sphere); //add spheres to the spere array
                    sphere_ = true;
                }
                
            } else {
                std::cerr << "Error parsing sphere line: " << line << std::endl;
            }
        }
        else if (word == "eye" || word == "viewdir" || word == "updir" || word == "bkgcolor" || word == "v" || word == "vn") {
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
                if(word == "vn"){
                    VectorType n = {num1, num2, num3};
                    normal_array.push_back(n);
                }
            }else {
                std::cerr << "Error parsing " << word << " line: " << line << std::endl;
            }
        }
        else if(word == "f"){
            int triangle_type = 0;
            std::string token;
            std::vector<int> indices;
            Triangle f;
            int end_count = 0;

            while (ss >> token) {
                int count = 0;
                std::stringstream tokenStream(token);
                std::string index;
                
                while (std::getline(tokenStream, index, '/')) {
                    if (!index.empty()) {
                        indices.push_back(std::stoi(index));
                        triangle_type ++;
                    }
                    count ++;
                }
                end_count = count;
            }
            if (triangle_type == 3 && end_count == 1) { // v v v
                if(materialArray.size() < 1){
                    std::cerr << "Material must be specified before an object in the input file" << std::endl;
                }
                //triangle f = {v1, v2, v3, n1, n2, n3, t1, t2, t3, material index, texture index}
                f = {indices[0], indices[1], indices[2], -1, -1, -1, -1, -1, -1, static_cast<int>(materialArray.size() - 1), -1};
                triangle_array.push_back(f);
                triangle_ = true;
            }
            else if(triangle_type == 6 && end_count == 2){ // v/vt v/vt v/vt
                if(materialArray.size() < 1 || texture_array.size() < 1){
                    std::cerr << "Material and texture must be specified before an object in the input file" << std::endl;
                }
                f = {indices[0], indices[2], indices[4], -1, -1, -1, indices[1], indices[3], indices[5], static_cast<int>(materialArray.size() - 1), static_cast<int>(texture_array.size() - 1)};
                triangle_array.push_back(f);
                triangle_ = true;
            }
            else if(triangle_type == 6 && end_count == 3){ // v//vn v//vn v//vn
                if(materialArray.size() < 1){
                    std::cerr << "Material must be specified before an object in the input file" << std::endl;
                }
                f = {indices[0], indices[2], indices[4], indices[1], indices[3], indices[5], -1, -1, -1, static_cast<int>(materialArray.size() - 1), -1};
                triangle_array.push_back(f);
                triangle_ = true;
            }
            else if(triangle_type == 9 && end_count == 3){ // v/vn/vt v/vn/vt v/vn/vt
                if(materialArray.size() < 1 || texture_array.size() < 1){
                    std::cerr << "Material and texture must be specified before an object in the input file" << std::endl;
                }
                f = {indices[0], indices[3], indices[6], indices[1], indices[4], indices[7], indices[2], indices[5], indices[8], static_cast<int>(materialArray.size() - 1), static_cast<int>(texture_array.size() - 1)};
                triangle_array.push_back(f);
                triangle_ = true;
            }
            else {
                std::cerr << "Error parsing " << word << " line: " << line << std::endl;
            }
        }
        else if (word == "imsize" || word == "vt") {
            if (ss >> num1 >> num2) {
                if(word == "imsize"){
                    px_width = num1;
                    px_height = num2;
                    valid_input.px_height_ = true;
                    valid_input.px_width_ = true;
                }
                if(word == "vt"){
                    TextureCoords t = {num1, num2};
                    texture_coords.push_back(t);
                }
            } else {
                std::cerr << "Error parsing line: " << line << std::endl;
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
        else if(word == "texture"){
            if(ss >> textureFile){
                std::cout << textureFile << std::endl;
                TextureFile texture = openTextureFile(textureFile);
                if(texture.success == 0 && !(texture.array.empty())){
                    texture_array.push_back(texture);
                }
            }
            else{
                std::cerr << "Error parsing texture file " << line << std::endl;
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
    if(!valid_input.bkgcolor_){
        std::cerr << "Error: No background color input" << std::endl;
        return -1;
    }
    if(!valid_input.px_height_ || !valid_input.px_width_){
        std::cerr << "Error: imsize width and height were not properly input" << std::endl;
        return -1;
    }

    //create the header for the ppm file
    fout << "P3\n";
    fout << px_width << " # image width\n" << px_height << " # image height\n";
    fout << "255\n";

    //if there is no input object, just return the background color
    if(!sphere_ && !triangle_){
        for(int i = 0; i < px_height; i++){
            for(int j = 0; j < px_width; j++){
                fout << std::round(bkgcolor.r * 255) << " " << std::round(bkgcolor.g * 255) << " " << std::round(bkgcolor.b * 255) << std::endl;
            }
        }
        std::cout << "No Object in scene" << std::endl;
        inputFile.close();
        return 0;
    }//check that if there is an object, there is also a material color
    else if((sphere_ || triangle_) && !material_){
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
            
            //use trace_ray and trace_polygon to look for object intersections
            Intersection trace_ray = traceRay(ray, -1);
            Intersection trace_polygon = tracePolygon(ray, -1);

            //determine whether a sphere or polygon was intersected closest to the eye
            float t_shape = find_t(trace_polygon.t, trace_ray.t);

            ColorType color;
            if(t_shape == trace_ray.t){
                color = shadeRay(trace_ray.shape_index, trace_ray.shape, trace_ray.ray, trace_ray.t);  
            }
            else if(t_shape = trace_polygon.t){
                color = shadeRay(trace_polygon.shape_index, trace_polygon.shape, trace_polygon.ray, trace_polygon.t);  
            }
            fout << std::round(color.r * 255) << " " << std::round(color.g * 255) << " " << std::round(color.b * 255) << std::endl;
        }
    }

    inputFile.close();
    return 0;
}