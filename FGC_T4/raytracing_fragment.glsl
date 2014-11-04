#version 400
uniform vec3 iResolution;
#define PI 3.1415926535897932384626433832795
#define EPSILON 0.000001

struct Material {
    vec3 diffuse;
    vec3 specular;
    float refractiveIndex;
    float reflectiveIndex;
    float alpha;
};

struct Ray {
    vec3 origin;
    vec3 direction;
};

struct Sphere {
    Material mat;
    vec3 center;
    float radius;
};

struct Cube {
    Material mat;
    vec2 facesX;
    vec2 facesY;
    vec2 facesZ;
};

struct Cylinder {
    Material mat;
    vec3 center;
    vec3 axis;
    float radius;
    float min;
    float max;
};

struct Triangle {
    Material mat;
    vec3 point1;
    vec3 point2;
    vec3 point3;
    vec3 normal;
};

struct Camera {
    vec3 position;
    vec3 direction;
};

struct Light {
    vec3 position;
    vec3 color;
};

// Fragment Shader Output
out vec4 fragColor;
// Clip Coordinates Transform
vec2 transform(vec2 p);
// Figures Intersections
float intersection(Ray r, Sphere s, inout vec3 position, inout vec3 normal);
float intersection(Ray r, Cylinder cy, inout vec3 position, inout vec3 normal);
float intersection(Ray r, Triangle tr, inout vec3 position, inout vec3 normal);
// Lighting Model
vec4 phong(Ray r, vec3 n, Material m, Light l, vec3 position);

void main()
{
    vec4 color = vec4(0.1, 0.1, 0.1, 0.0);
    Camera camera;
    camera.position = vec3(0.0, 0.0, 3.0);
    camera.direction = vec3(0.0, 0.0, -1.0);
    vec3 rayOrigin = camera.position;
    vec3 rayDirection = vec3(transform(gl_FragCoord.xy), 0.0) + camera.direction;
    rayDirection = normalize(rayDirection);
    Ray r;
    r.origin = rayOrigin;
    r.direction = rayDirection;
    Sphere s;
    s.center = vec3(0.0, 0.0, 0.0);
    s.radius = 1.0;
    Sphere s2;
    s2.center = vec3(1.5, 0.0, 0.0);
    s2.radius = 0.30;
    Triangle tr;
    tr.point1 = vec3(4.0, 5.0, 6.0);
    tr.point2 = vec3(-1.5, 2.0, -1.0);
    tr.point3 = vec3(1.0, -1.0, -1.0);
    Cylinder cy;
    cy.min = -2.0;
    cy.max = 3.0;
    cy.radius = 1.0;
    cy.axis = vec3(0.0, 0.0, 1.0);
    cy.center = vec3(vec2(-2.0, -2.0), 0.0);
    // Mat -> Azul1
    s.mat.diffuse = vec3(0.1, 0.1, 0.9);
    s.mat.specular = vec3(1.0, 1.0, 1.0);
    s.mat.reflectiveIndex = 0.3;
    s.mat.refractiveIndex = 0.0;
    s.mat.alpha = 0.0;
    // Mat -> Azul
    s2.mat.diffuse = vec3(0.1, 0.1, 0.9);
    s2.mat.specular = vec3(1.0, 1.0, 1.0);
    s2.mat.reflectiveIndex = 0.3;
    s2.mat.refractiveIndex = 0.0;
    s2.mat.alpha = 0.0;
    // Mat -> Rojo1
    tr.mat.diffuse = vec3(1.0, 0.1, 0.1);
    tr.mat.specular = vec3(1.0, 0.0, 0.3);
    tr.mat.reflectiveIndex = 0.3;
    tr.mat.refractiveIndex = 0.0;
    tr.mat.alpha = 0.0;
    // Mat -> Azul1
    cy.mat.diffuse = vec3(0.1, 0.0, 1.0);
    cy.mat.specular = vec3(0.0, .0, 1.0);
    cy.mat.reflectiveIndex = 0.3;
    cy.mat.refractiveIndex = 0.0;
    cy.mat.alpha = 0.0;
    // Light on TOP, only one for TESTING, TODO - MORE LIGHTS
    Light l;
    l.position = vec3(0.0, 5.0, 0.0);
    l.color =  vec3(1.0, 0.0, 0.5);
    vec3 normal, position;
    float t = intersection(r, s, position, normal);
    // TODO - DEPTH TESTING - REFLECTIONS - LIGHT BOUNCES
    //if (t > 0.0) {
    //    color += phong(r, normal, s.mat, l, position);
    //}
    float t2 = intersection(r, s2, position, normal);

    if (t2 > 0.0) {
        color += phong(r, normal, s2.mat, l, position);
    }

    float t4 = intersection(r, cy, position, normal);

    if (t4 > 0.0) {
        color += phong(r, normal, cy.mat, l, position);
    }

    // Unsure
    /* float t3 = intersection(r, tr, position, normal);

     if (t3 > 0.0) {
         color += phong(r, normal, tr.mat, l, position);
     }*/
    fragColor = color;
}


float intersection(Ray r, Sphere s, inout vec3 position, inout vec3 normal)
{
    vec3 oc = r.origin - s.center;
    float b = 2.0 * dot(oc, r.direction);
    float c = dot(oc, oc) - s.radius * s.radius;
    float h = b * b - 4.0 * c;

    if (h < 0.0) { return -1.0; }

    float t = (-b - sqrt(h)) / 2.0;
    position = r.origin + t * r.direction;
    normal = (position - s.center) / s.radius;
    return t;
}

float intersection(Ray r, Triangle tr, inout vec3 position, inout vec3 normal)
{
    // Möller–Trumbore Algorithm
    vec3 e1 = tr.point2 - tr.point1;
    vec3 e2 = tr.point3 - tr.point1;
    vec3 p = (r.direction * e2);
    float det = dot(e1, p);

    if (det > -EPSILON && det < EPSILON) { return -1.0; }

    float invDet = 1.0 / det;
    vec3 tvec = r.origin - tr.point1;
    float u = dot(tvec, p) * invDet;

    if (u < 0.0 || u > 1.0) { return -1.0; }

    vec3 q = tvec * e1;
    float v = dot(r.direction, q) * invDet;

    if (v < 0.0 || u + v > 1.0) { return -1.0; }

    float t  = dot(e2, q) * invDet;

    if (t > EPSILON) {
        position = r.origin + t * r.direction;
        normal = normalize(vec3(e1.y * e2.z - e1.z * e2.y,
                                e1.z * e2.x - e1.x * e2.z,
                                e1.x * e2.y - e1.y * e2.x)); // Placeholder - will pass as a uniform precalculated
        return t;
    }

    return 0.0;
}

float intersection(Ray r, Cylinder cy, inout vec3 position, inout vec3 normal)
{
    vec3 oc = r.origin - cy.center;
    vec3 h = r.direction - dot(r.direction, cy.axis) * cy.axis;
    vec3 q = oc - dot(oc, cy.axis) * cy.axis;
    float A  = dot(h, h);
    float B = 2.0 * dot(h, q);
    float C = dot(q, q) - cy.radius * cy.radius;
    float det = B * B - 4.0 * A * C;
    float tb, tt, t;
    t = tt = tb = -1.0;

    if (det <= 0.0) { return -1.0; }

    t = (-B - sqrt(det)) / 2.0 * A;
    // t = t > 0.0 ? t : (-B + sqrt(det)) / 2.0 * A;
    vec3 p = r.origin + t * r.direction;
    vec3 bCap = cy.axis * cy.min + cy.center;
    vec3 tCap = cy.axis * cy.max + cy.center;
    tb = dot(cy.axis, p - bCap);
    tt = dot(cy.axis, p - tCap);

    if (tb < 0.0 || tt > 0.0) {
        return -1.0;
    }

    vec3 nv = (tCap - bCap) / distance(tCap, bCap);
    float nt = dot(p - tCap, nv);
    vec3 spineP = tCap + nt * nv;
    position = p;
    normal = -normalize(p - spineP);
    return t;
}

vec2 transform(vec2 p)
{
    p = -1.0 + 2.0 * p / iResolution.xy;
    p.x *= iResolution.x / iResolution.y;
    return p;
}

vec4 phong(Ray r, vec3 n, Material m, Light l, vec3 position)
{
    vec3 lightDir = normalize(l.position - position);
    float d = length(lightDir);
    vec3 R = reflect(-lightDir, n);
    vec3 E = normalize(r.origin - position);
    float cosAlpha = clamp(dot(E, R), 0.0, 1.0);
    float cosTheta = clamp(dot(n, lightDir), 0.0, 1.0);
    float att = 1.0 / (1.0 + (0.5 * d) + (0.25 * d * d));
    vec3 ambient = vec3(0.1, 0.1, 0.1);
    vec3 diffuse = m.diffuse * l.color * cosTheta;
    vec3 specular = m.specular * l.color * pow(cosAlpha, 10.0);
    return vec4(ambient + diffuse + specular, 1.0) * att;
}