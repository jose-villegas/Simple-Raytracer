#define PI 3.1415926535897932384626433832795
uniform vec4 iMouse;
uniform vec3 iResolution;

struct Ray {
    vec3 origin;
    vec3 direction;
};

struct Sphere {
    vec3 origin;
    float radius;
};

struct Camera {
    vec3 position;
    vec3 direction;
};

vec2 transform(vec2 p);
float iSphere(Ray r, Sphere s);

void main()
{
    vec4 color = vec4(0.0, 0.0, 0.0, 1.0);
    Camera camera;
    camera.position = vec3(0.0, 0.0, 0.0);
    camera.direction = vec3(-transform(iMouse.xy), -1.0);
    vec3 rayOrigin = camera.position;
    vec3 rayDirection = vec3(transform(gl_FragCoord.xy), 0.0) + camera.direction;
    rayDirection = normalize(rayDirection);
    Ray r;
    r.origin = rayOrigin;
    r.direction = rayDirection;
    Sphere s;
    s.origin = vec3(0.0, 0.0, -2.1);
    s.radius = 0.55;
    Sphere s2;
    s2.origin = vec3(-1.0, 0.0, -3.1);
    s2.radius = 0.35;
    float t = iSphere(r, s);
    color.xyz += t;
    gl_FragColor = color;
}


float iSphere(Ray r, Sphere s)
{
    vec3 oc = r.origin - s.origin;
    float b = 2.0 * dot(oc, r.direction);
    float c = dot(oc, oc) - s.radius * s.radius;
    float h = b * b - 4.0 * c;

    if (h < 0.0) { return -1.0; }

    float t = (-b - sqrt(h)) / 2.0;
    return t;
}

vec2 transform(vec2 p)
{
    p = -1.0 + 2.0 * p / iResolution.xy;
    p.x *= iResolution.x / iResolution.y;
    return p;
}