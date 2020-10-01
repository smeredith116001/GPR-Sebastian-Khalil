// GLSL STARTER CODE BY DANIEL S. BUCKSTEIN
// asPoint: promote a 3D vector into a 4D vector representing a point (w=1)
//point: input 3D vector

vec4 asPoint(in vec3 point)
{
    return vec4(point, 1.0);
}
// asOffset: promote a 3D vector into a 4D vector representing an offset (w=0)
//    point: input 3D vector
vec4 asOffset(in vec3 offset)
{
    return vec4(offset, 0.0);
}
// calcViewport: calculate the viewing plane (viewport) coordinate
//    viewport:       output viewing plane coordinate
//    ndc:            output normalized device coordinate
//    uv:             output screen-space coordinate
//    aspect:         output aspect ratio of screen
//    resolutionInv:  output reciprocal of resolution
//    viewportHeight: input height of viewing plane
//    fragCoord:      input coordinate of current fragment (in pixels)
//    resolution:     input resolution of screen (in pixels)
void calcViewport(out vec3 viewport, out vec2 ndc, out vec2 uv, out float aspect, out vec2 resolutionInv, in float viewportHeight, in float focalLength, in vec2 fragCoord, in vec2 resolution)
{
    // inverse (reciprocal) resolution = 1 / resolution
    resolutionInv = 1.0 / resolution;
    // aspect ratio = screen width / screen height
    aspect = resolution.x * resolutionInv.y;
    // uv = screen-space coordinate = [0, 1) = coord / resolution
    uv = fragCoord * resolutionInv;
    // ndc = normalized device coordinate = [-1, +1) = uv*2 - 1
    ndc = uv * 2.0 - 1.0;
    // viewport: x = [-aspect*h/2, +aspect*h/2), y = [-h/2, +h/2), z = -f
    viewport = vec3(ndc * vec2(aspect, 1.0) * (viewportHeight * 0.5), -focalLength);
}
// calcRay: calculate the ray direction and origin for the current pixel
//    rayDirection: output direction of ray from origin
//    rayOrigin:    output origin point of ray
//    viewport:     input viewing plane coordinate (use above function to calculate)


//    focalLength:  input distance to viewing plane
void calcRay(out vec4 rayDirection, out vec4 rayOrigin, in vec3 eyePosition, in vec3 viewport)
{
    // ray origin relative to viewer is the origin
    // w = 1 because it represents a point; can ignore when using
    rayOrigin = asPoint(eyePosition);
    // ray direction relative to origin is based on viewing plane coordinate
    // w = 0 because it represents a direction; can ignore when using
    rayDirection = asOffset(viewport - eyePosition);
}
struct sSphere
{
    vec4 center;
    float radius;
};
void initSphere(out sSphere sphere, in vec3 center, in float radius)
{
    sphere.center = asPoint(center);
    sphere.radius = abs(radius);

}
float lenSq(vec2 x)
{
    return dot(x, x);
}
// calcColor: calculate the color of a pixel given a ray
//    rayDirection: input ray direction
//    rayOrigin:    input ray origin
vec4 calcColor(in vec4 rayDirection, in vec4 rayOrigin)
{
    // DUMMY RESULT: OUTPUT RAY DIRECTION AS-IS//  -> what does the ray look like as color?
    //return rayDirection;
    //Scene
    sSphere sphere;
    initSphere(sphere, vec3(0.0, 0.0, -4.0), 0.5);

    //testing proc sphere
    vec3 dp;
    dp.xy = rayDirection.xy - sphere.center.xy;
    float lSq = lenSq(dp.xy);
    float rSq = sphere.radius * sphere.radius;
    if (lSq <= rSq)
        if (length(dp.xy) <= sphere.radius)
        {
            //return vec4(0.0,1.0,1.0,1.0);
            //dx,dy,r
            dp.z = rSq - lSq;
            vec3 position = sphere.center.xyz + vec3(dp.x, dp.y, sqrt(dp.z));
            vec3 normal = //normalize(position - sphere.center.xyz);
                (position - sphere.center.xyz) / sphere.radius;
            return vec4(normal * 0.5 + 0.5, 1.0);
            // this is actually dz_sq
        }
    // BACKGROUND
    const vec3 warm = vec3(0.8, 0.4, 0.2), cool = vec3(0.2, 0.4, 0.8);
    //return vec4(mix(warm, cool, rayDirection.y), 1.0);
    return vec4(0.5);
}
// mainImage: process the current pixel (exactly one call per pixel)
//    fragColor: output final color for current pixel
//    fragCoord: input location of current pixel in image (in pixels)
void mainImage(out vec4 fragColor, in vec2 fragCoord)
{
    // viewing plane (viewport) info
    vec3 viewport;
    vec2 ndc, uv, resolutionInv;
    float aspect;
    const float viewportHeight = 2.0, focalLength = 1.0;
    // ray
    vec4 rayDirection, rayOrigin;
    // setup
    fragColor = vec4(0.0);
    calcViewport(viewport, ndc, uv, aspect, resolutionInv, viewportHeight, focalLength, fragCoord, iResolution.xy);
    calcRay(rayDirection, rayOrigin, vec3(0.0), viewport);
    fragColor += calcColor(rayDirection, rayOrigin);
}
struct pointLight
{
    vec3 center;
    vec4 color;
    float intensity;
};

void initpointLight(out pointLight light, in vec3 center, vec4 color, in float intensity)
{
    light.center = asPoint(center);
    light.color = color;
    light.intenisty = intensity;

}

vec4 LambertanceReflection(out vec4 normal, in vec2 uv)
{
    /*float theta = 6.283185 * uv.x;
        sSphere sphere;
        initSphere(sphere, vec3(0.0, 0.0, -4.0),0.5);
        uv.y = 2.0 * uv.y - 1.0;
        vec3 spherePoint = vec3(sqrt(1.0 - uv.y * uv.y) * vec2(cos(theta), sin(theta)),uv.y);
        vec3 normal = (position - sphere.center.xyz) / sphere.radius;;
        return vec4 (normal * 0.5 + 0.5, 1.0);*/

    pointLight light1;
    float diffIntense;
    float diffCoef;
    float attenIntense;
    vec4 lightPos;

    diffIntense = diffCoef * attenIntense;


    return Id = Kd * Il
}
/*vec3 lambert(vec3 normal, vec2 uv)
{
    float theta = 6.283185 * uv.x;
    float r = sqrt(uv.y);
    vec3 s, t;
vec3 spherePoint = s * cos(theta) * r + t * sin(theta) * r + normal * sqrt(1.0 - uv.y);


    //calcColor(normal, s, t);

    return spherePoint;
}*/

