//Vec 3 Header
/*
    sphere.h
    Code sourced from Ray Tracing in One Weekend Peter Shirley
    edited by Steve Hollasch and Trevor David Black

    Modified by: Khalil Lyons
*/



#define iostream
#define cmath
#define limits
#define memory
#define cstdlib
#define random
#define vector

struct vec3 sphere
{
    //constructors
    vec3() : e{ 0,0,0 } {}
    vec3(double e0, double e1, double e2) : e{ e0, e1, e2 } {}

    //return functions
    double x() const { return e[0]; }
    double y() const { return e[1]; }
    double z() const { return e[2]; }

    //overload operators for vectors
    vec3 operator-() const { return vec3(-e[0], -e[1], -e[2]); }
    double operator[](int i) const { return e[i]; }
    double& operator[](int i) { return e[i]; }

    vec3& operator+=(const vec3& v) {
        e[0] += v.e[0];
        e[1] += v.e[1];
        e[2] += v.e[2];
        return *this;
    }

    vec3& operator*=(const double t) {
        e[0] *= t;
        e[1] *= t;
        e[2] *= t;
        return *this;
    }

    vec3& operator/=(const double t) {
        return *this *= 1 / t;
    }

    vec3& operator*(const double t) {

        return *this *= t;
    }

    vec3& operator/(const double t) {
        e[0] /= t;
        e[1] /= t;
        e[2] /= t;
        return *this;
    }

    //length functions
    double length() const {
        return sqrt(length_squared());
    }

    double length_squared() const {
        return e[0] * e[0] + e[1] * e[1] + e[2] * e[2];
    }

    double e[3];
};

// Type aliases for vec3
using point3 = vec3;   // 3D point
using color = vec3;    // RGB color





//Ray Header
/*
    ray.h
    Code sourced from Ray Tracing in One Weekend Peter Shirley
    edited by Steve Hollasch and Trevor David Black

    Modified by: Sebastian Meredith
*/



//Ray Class for independent rays
struct ray {

    ray() {}
    ray(const point3& origin, const vec3& direction)
        : orig(origin), dir(direction)
    {}

    point3 origin() const { return orig; }
    vec3 direction() const { return dir; }

    point3 at(double t) const {
        return orig + t * dir;
    }


    point3 orig;
    vec3 dir;
};






//Sphere Header

    /*
    sphere.h
    Code sourced from Ray Tracing in One Weekend Peter Shirley
    edited by Steve Hollasch and Trevor David Black

    Modified by: Sebastian Meredith
*/

//includes

//Class for spheres being able to be hit by rays
struct sphere {

    sphere() {}

    sphere(struct hittable, point3 cen, double r) : center(cen), radius(r) {};

    virtual bool hit(
        const ray& r, double tmin, double tmax, hit_record& rec) const override;


    point3 center;
    double radius;
};
//Determining whether they hit or not
bool sphere::hit(const ray& r, double t_min, double t_max, hit_record& rec) const {
    double oc = r.origin() - center;
    double a = r.direction().length_squared();
    float half_b = dot(oc, r.direction());
    double c = oc.length_squared() - radius * radius;
    double discriminant = half_b * half_b - a * c;

    if (discriminant > 0) {
        float root = sqrt(discriminant);

        double temp = (-half_b - root) / a;
        if (temp < t_max && temp > t_min) {
            rec.t = temp;
            rec.p = r.at(rec.t);
            rec.normal = (rec.p - center) / radius;
            return true;
        }

        temp = (-half_b + root) / a;
        if (temp < t_max && temp > t_min) {
            rec.t = temp;
            rec.p = r.at(rec.t);
            rec.normal = (rec.p - center) / radius;
            return true;
        }
    }

    return false;
}



//RtWeekend Header
/*
    rtweekend.h
    Code sourced from Ray Tracing in One Weekend Peter Shirley
    edited by Steve Hollasch and Trevor David Black

    Modified by: Sebastian Meredith
*/



//Added changes from Chapter 7 


inline double clamp(double x, double min, double max) {
    if (x < min) return min;
    if (x > max) return max;
    return x;
}
}
inline double random_double() {
    // Returns a random real in [0,1).
    return rand() / (RAND_MAX + 1.0);
}
inline double random_double() {
    static uniform_real_distribution<double> distribution(0.0, 1.0);
    static mt19937 generator;
    return distribution(generator);
}

inline double random_double(double min, double max) {
    // Returns a random real in [min,max).
    return min + (max - min) * random_double();
}

// Usings

using shared_ptr;
using make_shared;
using sqrt;

// Constants

const double infinity = numeric_limits<double>::infinity();
const double pi = 3.1415926535897932385;

// Utility Functions

inline double degrees_to_radians(double degrees) {

    return degrees * pi / 180.0;
}


// Common Headers






//Hittable header
/*
    hittable.h
    Code sourced from Ray Tracing in One Weekend Peter Shirley
    edited by Steve Hollasch and Trevor David Black

    Modified by: Sebastian Meredith
*/


//includes

//Creates hit function that will take a ray
struct hit_record {
    point3 p;
    vec3 normal;
    double t;
};

struct hittable {

    virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const = 0;
};




//Hittable list header 
    /*
    hitable_list.h
    Code sourced from Ray Tracing in One Weekend Peter Shirley
    edited by Steve Hollasch and Trevor David Black

    Modified by: Sebastian Meredith
*/


//includes

//Usings
using std::shared_ptr;
using std::make_shared;
//Creating list of objects that the rays can hit
struct hittable_list {

    hittable_list() {}
    hittable_list(struct hittable, shared_ptr<hittable> object) { add(object); }

    void clear() { objects.clear(); }
    void add(shared_ptr<hittable> object) { objects.push_back(object); }

    virtual bool hit(
        const ray& r, double tmin, double tmax, hit_record& rec) const override;


    vector<shared_ptr<hittable>> objects;
};

bool hittable_list::hit(const ray& r, double t_min, double t_max, hit_record& rec) const {
    hit_record temp_rec;
    bool hit_anything = false;
    double closest_so_far = t_max;

    for (const auto& object : objects) {
        if (object->hit(r, t_min, closest_so_far, temp_rec)) {
            hit_anything = true;
            closest_so_far = temp_rec.t;
            rec = temp_rec;
        }
    }

    return hit_anything;
}



//Camera Header 


struct camera {

    camera() {
        float aspect_ratio = 16.0 / 9.0;
        float viewport_height = 2.0;
        float viewport_width = aspect_ratio * viewport_height;
        float focal_length = 1.0;

        origin = point3(0, 0, 0);
        horizontal = vec3(viewport_width, 0.0, 0.0);
        vertical = vec3(0.0, viewport_height, 0.0);
        lower_left_corner = origin - horizontal / 2 - vertical / 2 - vec3(0, 0, focal_length);
    }

    ray get_ray(double u, double v) const {
        return ray(origin, lower_left_corner + u * horizontal + v * vertical - origin);
    }


    point3 origin;
    point3 lower_left_corner;
    vec3 horizontal;
    vec3 vertical;
};


// calcViewport: calculate the viewing plane (viewport) coordinate
//    viewport:       output viewing plane coordinate
//    ndc:            output normalized device coordinate
//    uv:             output screen-space coordinate
//    aspect:         output aspect ratio of screen
//    resolutionInv:  output reciprocal of resolution
//    viewportHeight: input height of viewing plane
//    fragCoord:      input coordinate of current fragment (in pixels)
//    resolution:     input resolution of screen (in pixels)
void calcViewport(out vec2 viewport, out vec2 ndc, out vec2 uv,
    out float aspect, out vec2 resolutionInv,
    in float viewportHeight, in vec2 fragCoord, in vec2 resolution)
{
    // inverse (reciprocal) resolution = 1 / resolution
    resolutionInv = 1.0 / resolution;

    // aspect ratio = screen width / screen height
    aspect = resolution.x * resolutionInv.y;

    // uv = screen-space coordinate = [0, 1) = coord / resolution
    uv = fragCoord * resolutionInv;

    // ndc = normalized device coordinate = [-1, +1) = uv*2 - 1
    ndc = uv * 2.0 - 1.0;

    // viewport: x = [-aspect*h/2, +aspect*h/2), y = [-h/2, +h/2)
    viewport = ndc * (vec2(aspect, 1.0) * (viewportHeight * 0.5));
}


// calcRay: calculate the ray direction and origin for the current pixel
//    rayDirection: output direction of ray from origin
//    rayOrigin:    output origin point of ray
//    viewport:     input viewing plane coordinate (use above function to calculate)
//    focalLength:  input distance to viewing plane
void calcRay(out vec4 rayDirection, out vec4 rayOrigin,
    in vec2 viewport, in float focalLength)
{
    // ray origin relative to viewer is the origin
    // w = 1 because it represents a point; can ignore when using
    rayOrigin = vec4(0.0, 0.0, 0.0, 1.0);

    // ray direction relative to origin is based on viewing plane coordinate
    // w = 0 because it represents a direction; can ignore when using
    rayDirection = vec4(viewport.x, viewport.y, -focalLength, 0.0);
}


// calcColor: calculate the color of a pixel given a ray
//    rayDirection: input ray direction
//    rayOrigin:    input ray origin
vec4 calcColor(in vec4 rayDirection, in vec4 rayOrigin)
{
    //vec3 unit_direction = unit_vector(rayDirection);
    float t = 0.5 * (rayDirection.y + 1.0);
    (1.0 - t)* color(1.0, 1.0, 1.0) + t * color(0.5, 0.7, 1.0);
}
return


// mainImage: process the current pixel (exactly one call per pixel)
//    fragColor: output final color for current pixel
//    fragCoord: input location of current pixel in image (in pixels)
void mainImage(out vec4 fragColor, in vec2 fragCoord)
{

    // Image

    const auto aspect_ratio = 16.0 / 9.0;
    const int image_width = 400;
    const int image_height = static_cast<int>(image_width / aspect_ratio);
    const int samples_per_pixel = 100;

    // World

    hittable_list world;
    world.add(make_shared<sphere>(point3(0, 0, -1), 0.5));
    world.add(make_shared<sphere>(point3(0, -100.5, -1), 100));

    // Camera
    camera cam;

    // Render

    //cout << "P3\n" << image_width << " " << image_height << "\n255\n";

    for (int j = image_height - 1; j >= 0; --j) {
        //std::cerr << "\rScanlines remaining: " << j << ' ' << std::flush;
        for (int i = 0; i < image_width; ++i) {
            color pixel_color(0, 0, 0);
            for (int s = 0; s < samples_per_pixel; ++s) {
                auto u = (i + random_double()) / (image_width - 1);
                auto v = (j + random_double()) / (image_height - 1);
                ray r = cam.get_ray(u, v);
                pixel_color += ray_color(r, world);
            }
            write_color(std::cout, pixel_color, samples_per_pixel);
        }
    }

    //std::cerr << "\nDone.\n";
}


// viewing plane (viewport) info
vec2 viewport, ndc, uv, resolutionInv;
float aspect;
const float viewportHeight = 2.0, focalLength = 1.0;

// ray
vec4 rayDirection, rayOrigin;

// setup
calcViewport(viewport, ndc, uv, aspect, resolutionInv,
    viewportHeight, fragCoord, iResolution.xy);
calcRay(rayDirection, rayOrigin,
    viewport, focalLength);

// color
fragColor = calcColor(rayDirection, rayOrigin);

// TEST COLOR:
//  -> what do the other things calculated above look like?
//fragColor = vec4(viewport, 0.0, 0.0);
//fragColor = vec4(ndc, 0.0, 0.0);
//fragColor = vec4(uv, 0.0, 0.0);
}