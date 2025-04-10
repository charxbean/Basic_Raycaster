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
int d = 1;

//color type, r, g, b
typedef struct{
    float r, g, b;
} ColorType;

//Defines a material color
typedef struct{
    ColorType Od;
    ColorType Os;
    float ka, kd, ks;
    int n;
    float alpha, aeda, F0;
}MaterialColor;

//Defines a sphere
typedef struct{
    VectorType center;
    float r;
    int m; //index into material array
    int t; //index into texture array
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
    std::vector<std::vector<ColorType>> array; //2D array of texture colors
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

float bkg_aeda;
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

/**
 * Check the input barycentric coordinates to determine if they meet the conditions for being inside a triangle
 * if conditions are met, return true, otherwise return false
 */
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
VectorType vector_I;

//recursive depth of computeReflection
int rec_depth = 0;
int rec_depth_T = 0;
//reflected color recursively computed in shade_ray
ColorType R_color;
//refracted color recursivley computed in shade_ray
ColorType T_Color;
//total color computed at each pixel in shade_ray
ColorType illumination;

float ni = bkg_aeda;
float nt;

BGA get_bga(int triangle_index, Raytype ray, float t){
    VectorType p0 = vertex_array[(triangle_array[triangle_index].a)-1];
    VectorType p1 = vertex_array[(triangle_array[triangle_index].b)-1];
    VectorType p2 = vertex_array[(triangle_array[triangle_index].c)-1];

    VectorType e1 = vectorSubtract(p1, p0);
    VectorType e2 = vectorSubtract(p2, p0);

    VectorType p = vectorAdd(ray.origin, vectorScalar(ray.intersection, t));
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

//initialize traceRay and tracePolygon so they can be used inside ShadowRay
Intersection traceRay(Raytype ray, int self_index);
Intersection tracePolygon(Raytype ray, int self_index);

/*
* Determine if a given shape is in shadow relative to the given light source
* Returns 1 if the shape is not in shadow, 0 if it is in shadow
*/
int shadowRay(Raytype shadow_ray, int shape_index, int light_index, int shape){
    int shadow_flag;
    //determine shadow ray intersection based on light type
    if(lightArray[light_index].type == 0){ //directional
        shadow_ray.intersection = L_array[light_index];
    }
    else{ //point
        shadow_ray.intersection = vectorSubtract(lightArray[light_index].position, shadow_ray.origin); //light position from the sphere surface
        normalize(shadow_ray.intersection);
    }
    
    //send a ray from the surface of the object to the light source
    Intersection trace_ray;
    Intersection trace_polygon;
    if(shape == 0){ //ignore self index if the object is a sphere
        trace_ray = traceRay(shadow_ray, shape_index);
        trace_polygon = tracePolygon(shadow_ray, -1);
    }
    else if(shape == 1){ //ignore self index if the object is a polygon
        trace_polygon = tracePolygon(shadow_ray, shape_index);
        trace_ray = traceRay(shadow_ray, -1);
    }
    else{
        std::cerr << "Attempting to trace undefined shape" << std::endl;
        trace_ray = {-1, 0, shadow_ray, -1};
    }

    //use the t values from trace_ray and trace_polygon to determine if there is a shadow or not
    if(trace_ray.t < 0 && trace_polygon.t < 0){
        shadow_flag = 1; //no intersection found, no shadow
    }
    else{
        float blocking_dist = find_t(trace_ray.t, trace_polygon.t);
        float light_dist = vectorLength(vectorSubtract(lightArray[light_index].position, shadow_ray.origin));
        //if light is directional, the object blocks the light
        if(lightArray[light_index].type == 0){
            shadow_flag = 0;
        }//if t is greater than the distance from the sphere to a point light source, the object doesn't block the light
        else if(blocking_dist > light_dist){
            shadow_flag = 1;
        }
        else{//point light and the object is blocking the light
            shadow_flag = 0;
        }
    }

    return shadow_flag;
}

//initialize computeReflection and computeRefraction so they can be used inside shadeRay
ColorType computeReflection(Raytype R, VectorType I, VectorType N, MaterialColor material, int rec_depth, int shape_index, int shape);
ColorType computeRefraction(Raytype T, VectorType I, VectorType N, MaterialColor material, int rec_depth, int shape_index, int shape);

/**
 * Takes a struct defining the index of a shape, the shape type, the ray that intersected it and the t value of the intersection
 * Returns the color of the pixel at the intersection point, computed using the Blinn-Phong method
 */
ColorType shadeRay(int shape_index, int shape, Raytype ray, float t){

    //if t < 0 then no object intersection was found, return the background color
    if(t < 0){
        return bkgcolor;
    }

    //initialize material color
    MaterialColor material;
    ColorType texture_color = {-1, -1, -1};
    nt = material.aeda;

    //the intersection point of the surface from the eye
    VectorType intersection = vectorAdd(ray.origin, vectorScalar(ray.intersection, t));

    //Define vector I as the opposite direction of the incoming intersection ray
    vector_I = vectorSubtract(intersection, ray.origin);
    normalize(vector_I);
    vector_I = vectorScalar(ray.intersection, -1);

    //Determine shape normals and texture color (if applicable) based on the shape type
    if(shape == 0){ //the shape is a sphere

        //find the material that corresponds with the given sphere
        int mat_index = sphereArray[shape_index].m;
        material = materialArray[mat_index];
        SphereType sphere = sphereArray[shape_index];

        //calculate normal vector N
        vector_N = vectorDivide(vectorSubtract(intersection, sphere.center), sphere.r);
        normalize(vector_N);
        //if the normal is facing the wrong way, flip it
        //ray coming from inside the sphere
        if(dotProduct(vector_N, vector_I) < 0){
            vector_N = vectorScalar(vector_N, -1); 
        }
        ni = material.aeda;
        nt = bkg_aeda;
        
        //if the sphere has a texture, calculate u and v to get the texture color
        if(sphere.t != -1){ 
            TextureFile *texture = &texture_array[sphere.t];
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

            int i = std::round(u*(texture->width -1));
            int j = std::round(v*(texture->height -1));
            texture_color = texture->array[j][i];
        }
    }
    else if(shape == 1){ //the shape is a triangle
        Triangle& triangle = triangle_array[shape_index];
        material = materialArray[triangle.material];

        //FLAT SHADING: immediately return OD
        if(triangle.n1 == -1){
            VectorType e1 = vectorSubtract(vertex_array[(triangle.b)-1], vertex_array[(triangle.a)-1]);
            VectorType e2 = vectorSubtract(vertex_array[(triangle.c)-1], vertex_array[(triangle.a)-1]);

            vector_N = crossProduct(e1, e2);
            normalize(vector_N);
            if(dotProduct(vector_N, vector_I) < 0){
                vector_N = vectorScalar(vector_N, -1); 
            }
            ni = material.aeda;
            nt = bkg_aeda;
        }
        //SMOOTH SHADING: Calculate normal
        if(triangle.n1 > 0){
            vector_N = interpolate_t_normal(shape_index, ray, t);
            if(dotProduct(vector_N, vector_I) < 0){
                vector_N = vectorScalar(vector_N, -1); 
            }
            ni = material.aeda;
            nt = bkg_aeda;
        }
        
        /**
         * TEXTURE MAPPING
         * if a texture and texture coords were input, find u,v to lookup Od color
         */
        if(triangle.t1 > 0){
            TextureFile* texture = &texture_array[triangle.texture];
            BGA bga = get_bga(shape_index, ray, t);
            if(bga.alpha == -1){
                std::cerr << "Error: bga values are invalid" << std::endl;
            }
            float u0 = texture_coords[(triangle.t1)-1].u;
            float u1 = texture_coords[(triangle.t2)-1].u;
            float u2 = texture_coords[(triangle.t3)-1].u;

            float v0 = texture_coords[(triangle.t1)-1].v;
            float v1 = texture_coords[(triangle.t2)-1].v;
            float v2 = texture_coords[(triangle.t3)-1].v;

            float u = bga.alpha*u0 + bga.beta*u1 + bga.gamma*u2;
            float v = bga.alpha*v0 + bga.beta*v1 + bga.gamma*v2;
            u = std::fmod(u, 1.0f);
            if(u < 0) u += 1.0f;
            
            v = std::fmod(v, 1.0f);
            if(v < 0) v += 1.0f;
            int i = std::round(u*(texture->width-1));
            int j = std::round(v*(texture->height-1));

            //change the od color to the pixel at (i,j) in the texture image
            texture_color = texture->array[j][i];
        }
    }
    else{
        std::cerr << "Attempting to shade an undefined shape" << std::endl;
        return bkgcolor;
    }

    //calculate vector L for each POINT light source
    for(int k = 0; k < lightArray.size(); k++){
        if(lightArray[k].type == 1){
            //(light_position - surface position)
            VectorType new_L = vectorSubtract(lightArray[k].position, intersection);
            //calculate vector L
            vector_L = vectorDivide(new_L, vectorLength(new_L));
            normalize(vector_L);
            //add into vector L array corresponding to the index of this light in the light array
            L_array.at(k) = vector_L;
        }
    }

    //initialize illumination to the ambient light
    //initialize the base diffuse color to Od or texture color
    ColorType diffuse;
    if(texture_color.r == -1){ //no texture
        illumination = colorMultiply(material.Od, material.ka); //ka + Od
        diffuse = colorMultiply(material.Od, material.kd); //kd * od
    }
    else{
        //if theres a texture, switch od with texture color
        illumination = colorMultiply(texture_color, material.ka); //ka + texture_color
        diffuse = colorMultiply(texture_color, material.kd); //kd * texture_color
    }
    //initialize base specular color
    ColorType specular = colorMultiply(material.Os, material.ks); //ks * Os

    //for each light source, compute it's individual illumination and add to the total illumination value
    for(int k = 0; k < lightArray.size(); k++){
        //define vector H for each light source
        vector_H = vectorAdd(vector_L, vector_I);
        normalize(vector_H);

        //calculate diffuse color
        float max_NL = max(0.0, dotProduct(vector_N, L_array[k])); //clamp diffuse at 0
        diffuse = colorMultiply(diffuse, max_NL);

        //calculate specular color
        float max_NH = max(0.0, dotProduct(vector_N, vector_H)); //clamp specular at 0
        specular = colorMultiply(specular, pow(max_NH, material.n)); //(ks * Os) * (N dot H)^n

        //SHADOWS - cast a ray from the object surface to each light source to determine if it's in shadow
        Raytype shadow_ray;
        shadow_ray.origin = intersection; //the surface intersection point

        //returns 1 for no shadow, 0 for in shadow
        int shadow_flag = shadowRay(shadow_ray, shape_index, k, shape);

        //if there is no shadow, add the illumination from this light source to the total
        if(shadow_flag == 1){
            ColorType illum_L = colorMultiply(colorAdd(diffuse, specular), lightArray[k].intensity);
            illumination = colorAdd(illumination, illum_L);
            clampColor(illumination);
        }

    }

    //Add F0 + R_lambda
    //Fxn to compute F0 and R_lambda (recursive)
    //Returns the value of F0 + R_lambda
    if(material.ks > 0){
        Raytype R;
        R.origin = intersection; //current ray/surface point
        float a = dotProduct(vector_N, vector_I);
        VectorType reflection_dir = vectorSubtract(vectorScalar(vector_N, 2 * a), vector_I);
        normalize(reflection_dir);
        R.intersection = reflection_dir;
        //prevent self intersecting messing up the reflection color
        //R.origin = vectorAdd(vectorScalar(vector_N, .001f), intersection);

        rec_depth ++;
        ColorType reflection = computeReflection(R, vector_I,  vector_N, material, rec_depth, shape_index, shape);
        illumination = colorAdd(illumination, reflection);
    }

    if(material.alpha == 1){
        Raytype T;
        T.origin = intersection; //current ray/surface point
        float a = dotProduct(vector_N, vector_I);
        float n = ni/nt;
        float delim = 1 - pow(n, 2) * (1 - pow(a, 2));
        if(delim > 0){
            T.intersection = vectorScalar(vectorScalar(vector_N, -1), std::sqrt(delim));
            T.intersection = vectorAdd(T.intersection, vectorScalar(vectorSubtract(vectorScalar(vector_N, a), vector_I), n));
            normalize(T.intersection);
            rec_depth_T ++;
            ColorType refraction = computeRefraction(T, vector_I, vector_N, material, rec_depth_T, shape_index, shape);
            illumination = colorAdd(illumination, refraction);
        }
        //if delim < 0, total internal reflection occurs, no refraction??
    }

    //if depth cue was input, add depth cue color to the illumination calculation
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

/**
 * Computes the reflective color added to a surface based on the input reflection ray
 * Returns a color type of R_lamda * Fr
 */
ColorType computeReflection(Raytype R, VectorType I, VectorType N, MaterialColor material, int rec_depth, int shape_index, int shape){
    
    //END RECURSION if recursion depth is greater than 10
    if(rec_depth >= 10){
        return R_color;
    }

    //Compute Reflection ray direction: 2(N dot I)N - I
    float a = dotProduct(N, I);
    //std::cout << "a: " << a << std::endl;

    //Compute Fr: F0 + (1 - F0)(1 - (N dot I))^5
    float Fr;
    if(material.F0 > .9){
        Fr = 1;
    }
    else if(material.F0 < 0.01){
        Fr = pow((1-a), 5);
    }
    else{
        Fr = material.F0 + (1 - material.F0) * pow((1 - a), 5);
    }
    if(Fr > 1 || Fr < 0){
        return R_color;
    }
    

    Intersection trace_ray;
    Intersection trace_polygon;
    //find the closest object intersection with the reflection ray
    if(shape == 0){ //sphere
        trace_ray = traceRay(R, shape_index);
        trace_polygon = tracePolygon(R, -1);
    }
    else if(shape == 1){
        trace_ray = traceRay(R, -1);
        trace_polygon = tracePolygon(R, shape_index);
    }

    float t = find_t(trace_ray.t, trace_polygon.t);

    ColorType color = {0, 0, 0}; //initialize color to black

    //Use shade ray to define the pixel color based on the nearest R intersection, if any
    if(t < 0){
        return R_color; //END RECURSION
    }
    else if(t == trace_ray.t){ //Trace_ray.ray and trace_polygon.ray should be R
        color = colorAdd(shadeRay(trace_ray.shape_index, trace_ray.shape, R, trace_ray.t), color);  
    }
    else if(t == trace_polygon.t){
        color = colorAdd(shadeRay(trace_polygon.shape_index, trace_polygon.shape, R, trace_polygon.t), color);  
    }
    else{
        color = colorAdd(shadeRay(trace_ray.shape_index, trace_ray.shape, trace_ray.ray, t), color);  
    }

    color = colorMultiply(color, Fr); //multiply the color by Fr
    clampColor(color);
    R_color = colorAdd(R_color, color); //add the color to the total R_color value 
    clampColor(R_color);

    return R_color;
}

ColorType computeRefraction(Raytype T, VectorType I, VectorType N, MaterialColor material, int rec_depth, int shape_index, int shape){
    
    //END RECURSION if recursion depth is greater than 10
    if(rec_depth >= 10){
        return T_Color;
    }

    float a = dotProduct(N, I);

    //Compute Fr: F0 + (1 - F0)(1 - (N dot I))^5
    float Fr;
    if(material.F0 == 1){
        Fr = 1;
    }
    else if(material.F0 == 0){
        Fr = pow((1-a), 5);
    }
    else{
        Fr = material.F0 + (1 - material.F0) * pow((1 - a), 5);
    }

    if(material.alpha == 0){
        return T_Color = {T_Color.r + (1-Fr), T_Color.g + (1-Fr), T_Color.b + (1-Fr)};
    }
    
    ColorType color = {0, 0, 0}; //initialize color to black

    Intersection trace_ray = traceRay(T, -1);
    Intersection trace_polygon = tracePolygon(T, -1);
    float t = find_t(trace_ray.t, trace_polygon.t);

    if(t < 0){
        return T_Color; //END RECURSION
    }
    else if(t == trace_ray.t){ //Trace_ray.ray and trace_polygon.ray should be T
        color = colorAdd(shadeRay(trace_ray.shape_index, trace_ray.shape, T, trace_ray.t), color);  
    }
    else if(t == trace_polygon.t){
        color = colorAdd(shadeRay(trace_polygon.shape_index, trace_polygon.shape, T, trace_polygon.t), color);  
    }
    else{
        color = colorAdd(shadeRay(trace_ray.shape_index, trace_ray.shape, trace_ray.ray, t), color);  
    }

    float x = (1-Fr)*(1-material.alpha);
    color = colorMultiply(color, x); //multiply the color by Fr
    clampColor(color);
    T_Color = colorAdd(T_Color, color); 
    clampColor(T_Color);
    
    return T_Color;
}

/*
* checks each sphere in the sccene for the closest valid intersection with the given ray
* skips the input self-index if it's a valid index number
* returns a struct containting the closest t-value (or -1), the intersected sphere index and the ray
*/
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


/*
* checks each triangle in the sccene for the closest valid intersection with the given ray
* skips the input self-index if it's a valid index number
* returns a struct containting the closest t-value (Or -1), the intersected triangle index and the ray
*/
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
        
        p0 = vertex_array[(triangle_array[k].a)-1];
        p1 = vertex_array[(triangle_array[k].b)-1];
        p2 = vertex_array[(triangle_array[k].c)-1];

        e1 = vectorSubtract(p1, p0);
        e2 = vectorSubtract(p2, p0);

        //Define normal vector n
        VectorType n = crossProduct(e1, e2);

        //Define A, B, C and D
        float A = n.i;
        float B = n.j;
        float C = n.k;
        float D = -1 * ((A*p0.i) + (B*p0.j) + (C*p0.k));

        //find the denominator of the possible ray intersection
        float denominator = (A*ray.intersection.i) + (B*ray.intersection.j) + (C * ray.intersection.k);

        if(denominator == 0){ //no intersection
            float t = -1;
        }
        else{ //intersection, calculate t
            float t = (-1 * (A * ray.origin.i + B*ray.origin.j + C*ray.origin.k + D))/denominator;
            
            //calculate intersection point p and vector ep
            //VectorType p = {(ray.origin.i + (t * ray.intersection.i)), (ray.origin.j + (t * ray.intersection.j)), (ray.origin.k + (t * ray.intersection.k))};
            VectorType p = vectorAdd(ray.origin, vectorScalar(ray.intersection, t));
            VectorType ep = vectorSubtract(p, p0);
        
            //calculate the determinant and the barycentric coordinates
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
                    //update t-closest to reflect a valid intersection
                    t_closest = find_t(t, t_closest);
                    if(t_closest == t){
                        triangle_index = k;
                    }
                }
            }
        }
    }
    //if the ray didn't intersect with any triangle, return -1
    if(t_closest == -1){
        return {-1, 1, ray, -1};
    }
    else{ //return the closest intersection
        return {triangle_index, 1, ray, t_closest};
    }
}


/**
 * reads the data in the input texture file into a 2D array of colors
 * returns the 2D color array or an empty array on error
 */
std::vector<std::vector<ColorType>>readTextureFile(int width, int height, std::ifstream& file){
    std::vector<std::vector<ColorType>> array(height, std::vector<ColorType>(width));
    for(int i = 0; i < height ; i++){
        for(int j = 0; j < width; j++){
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

/**
 * opens the texture ppm file and checks for a correct header
 * returns a TextureFile type containting a 2D array and a success indicator int: -1 on failure, 0 on success
 */
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
        float num11, num12;
        std::string textureFile;

        ss >> word; //Read the first word of every line

        if (word.empty()) {
            std::cerr << "Error: Empty keyword in line: " << line << std::endl;
            continue;
        }

        //skip # as comments in input file
        if(word == "#"){
            continue;
        }

        if (word == "mtlcolor"){
            if (ss >> num1 >> num2 >> num3 >> num4 >> num5 >> num6 >> num7 >> num8 >> num9 >> num10 >> num11 >> num12){
                ColorType Od = {num1, num2, num3};
                ColorType Os = {num4, num5, num6};

                //clamp alpha value to (0, 1)
                if(num11 > 1){
                    num11 = 1;
                    printf("Alpha value clamped to 1\n");
                }
                else if(num11 < 0){
                    num11 = 0;
                    printf("Alpha value clamped to 0\n");
                }

                //mtcolor = Od, Os, ka, kd, ks, n, alpha, aeda, F0
                mtlcolor = {Od, Os, num7, num8, num9, num10, num11, num12, -1};
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
        else if (word == "sphere" || word == "bkgcolor") {
            if (ss >> num1 >> num2 >> num3 >> num4) {
                if (word == "bkgcolor"){
                    bkgcolor = {num1, num2, num3};
                    bkg_aeda = num4;
                    valid_input.bkgcolor_ = true;
                }
                else if(word == "sphere"){
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
                }
                
            } else {
                std::cerr << "Error parsing line: " << line << std::endl;
            }
        }
        else if (word == "eye" || word == "viewdir" || word == "updir" || word == "v" || word == "vn") {
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
                    std::cout << "texture read a success" << std::endl;
                }
                else{
                    std::cerr << "Error parsing texture file" << line << std::endl;
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
    }
    //check that if there is an object, there is also a material color
    else if((sphere_ || triangle_) && !material_){
        std::cerr << "a sphere was input without a material color" << std::endl;
        inputFile.close();
        return -1;
    }
    
    //PRELIMINARY MATH
    normalize(viewdir);
    normalize(updir);

    //find vectors u and v
    VectorType u = crossProduct(viewdir, updir);
    normalize(u);
    VectorType v = crossProduct(u, viewdir);

    //calculate the width and height in 3d world coordinates
    float radian_vfov = (vfov * M_PI)/180;
    float height = (2*d*std::tan(0.5 * radian_vfov));
    float aspect = (float)px_width / (float)px_height;
    float width = aspect * height;

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

    VectorType h_change = vectorDivide(vectorSubtract(ur, ul), px_width-1);
    VectorType v_change = vectorDivide(vectorSubtract(ll, ul), px_height-1);
    
    //Calculate F0 for each material in the material array
    for(int i = 0; i < materialArray.size(); i++){
        float aeda = materialArray[i].aeda;
        materialArray[i].F0 = (aeda - 1)/(aeda + 1);
        materialArray[i].F0 = pow(materialArray[i].F0, 2);
    }

    //calculate the vector_L for Blinn Phone of each directional light source
    L_array.resize(lightArray.size());
    for(int k = 0; k < lightArray.size(); k++){
        if(lightArray[k].type == 0){
            VectorType L_inv = vectorScalar(lightArray[k].position, -1.0);
            if(L_inv.i == 0 && L_inv.j == 0 && L_inv.k == 0){ //if zero vector, dont divide by 0
                L_array.at(k) = lightArray[k].position;
            }
            else{
                normalize(L_inv);
                vector_L = L_inv;
                //insert this L vector into the L_array at the same index as the light in the light_array
                L_array.at(k) = vector_L;
            }
        }
    }

    //iterate through each pixel (j, i) where i = 0 to px_height-1, and j = 0 to px_width-1.
    for(int i = 0; i < px_height; i++){
        for(int j = 0; j < px_width; j++){
            rec_depth = 0; //initialize the recursion depth to 0
            rec_depth_T = 0; //initialize the recursion depth to 0
            R_color = {0, 0, 0}; //initialize the reflection color to 0
            T_Color = {0, 0, 0}; //initialize the refraction color to 0
            illumination = {0, 0, 0}; //initialize the illumination color to 0

            //the point where each ray should pass through the viewing window correspoindng to each pixel
            VectorType intersect_point = vectorAdd(ul, vectorScalar(h_change, j));
            intersect_point = vectorAdd(intersect_point, vectorScalar(v_change, i));

            VectorType direction = vectorSubtract(intersect_point, eye);
            normalize(direction);

            //initialize a ray for this pixel
            Raytype ray;
            ray.origin = eye;
            ray.intersection = direction;
            
            //use trace_ray and trace_polygon to look for an object intersection
            Intersection trace_ray = traceRay(ray, -1);
            Intersection trace_polygon = tracePolygon(ray, -1);

            //use t values to determine whether a sphere or polygon was intersected closest to the eye
            float t_shape = find_t(trace_polygon.t, trace_ray.t);

            //Use shade ray to define the pixel color based on which object was intersected if any
            ColorType color;

            if(t_shape < 0){ //no intersection
                color = bkgcolor;
            }
            else if(t_shape == trace_ray.t){
                color = shadeRay(trace_ray.shape_index, trace_ray.shape, trace_ray.ray, trace_ray.t);  
            }
            else if(t_shape == trace_polygon.t){
                color = shadeRay(trace_polygon.shape_index, trace_polygon.shape, trace_polygon.ray, trace_polygon.t);  
            }

            //print the color to the output file
            fout << std::round(color.r * 255) << " " << std::round(color.g * 255) << " " << std::round(color.b * 255) << std::endl;
        }
    }

    inputFile.close();
    return 0;
}