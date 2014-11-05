#version 400
uniform vec3 iResolution;
#define PI 3.1415926535897932384626433832795
#define MAX_DELTA 1e20
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
bool intersection(Ray r, Sphere s, inout vec3 position, inout vec3 normal);
bool intersection(Ray r, Cylinder cy, inout vec3 position, inout vec3 normal);
bool intersection(Ray r, Triangle tr, inout vec3 position, inout vec3 normal);
// Lighting Model
vec4 phong(Ray r, vec3 n, Material m, Light l, vec3 position);

void main()
{
    vec4 color = vec4(0.1, 0.1, 0.1, 0.0);
    Camera camera;
    camera.position = vec3(0.0, 0.0, 5.0);
    camera.direction = vec3(0.0, 0.0, -1.0);
    Ray r;
    r.origin = camera.position;
    r.direction = normalize(vec3(transform(gl_FragCoord.xy), 0.0) + camera.direction);
    Sphere s;
    s.center = vec3(0.0, 0.0, 0.0);
    s.radius = 1.0;
    Sphere s2;
    s2.center = vec3(1.5, 0.0, 0.0);
    s2.radius = 1.0;
    Triangle tr;
    tr.point1 = vec3(4.0, 5.0, 6.0);
    tr.point2 = vec3(-1.5, 2.0, -1.0);
    tr.point3 = vec3(1.0, -1.0, -1.0);
    Cylinder cy;
    cy.min = -2.0;
    cy.max = 2.0;
    cy.radius = 1.0;
    cy.axis = vec3(0.0, 1.0, 0.0);
    cy.center = vec3(-2.0, .0, 0.0);
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
    // Light on TOP, only one for TESTING, TODO - MORE LIGHTS
    Light l;
    l.position = vec3(0.0, 5.0, 8.0);
    l.color =  vec3(1.0, 0.0, 0.5);
    vec3 normal, position;
    // TODO - DEPTH TESTING - REFLECTIONS - LIGHT BOUNCES

    if (intersection(r, s2, position, normal)) {
        color += phong(r, normal, s2.mat, l, position);
    }

    if (intersection(r, cy, position, normal)) {
        color += phong(r, normal, tr.mat, l, position);
    }

    // Unsure
    /* float t3 = intersection(r, tr, position, normal);

     if (t3 > 0.0) {
         color += phong(r, normal, tr.mat, l, position);
     }*/
    fragColor = color;
}


bool intersection(Ray r, Sphere s, inout vec3 position, inout vec3 normal)
{
    vec3 oc = r.origin - s.center;
    float b = 2.0 * dot(oc, r.direction);
    float c = dot(oc, oc) - s.radius * s.radius;
    float h = b * b - 4.0 * c;

    if (h < 0.0) { return false; }

    float t = (-b - sqrt(h)) / 2.0;
    position = r.origin + t * r.direction;
    normal = (position - s.center) / s.radius;
    return true;
}

bool intersection(Ray r, Triangle tr, inout vec3 position, inout vec3 normal)
{
    // Möller–Trumbore Algorithm
    vec3 e1 = tr.point2 - tr.point1;
    vec3 e2 = tr.point3 - tr.point1;
    vec3 p = (r.direction * e2);
    float det = dot(e1, p);

    if (det > -EPSILON && det < EPSILON) { return false; }

    float invDet = 1.0 / det;
    vec3 tvec = r.origin - tr.point1;
    float u = dot(tvec, p) * invDet;

    if (u < 0.0 || u > 1.0) { return false; }

    vec3 q = tvec * e1;
    float v = dot(r.direction, q) * invDet;

    if (v < 0.0 || u + v > 1.0) { return false; }

    float t  = dot(e2, q) * invDet;

    if (t > EPSILON) {
        position = r.origin + t * r.direction;
        normal = normalize(vec3(e1.y * e2.z - e1.z * e2.y,
                                e1.z * e2.x - e1.x * e2.z,
                                e1.x * e2.y - e1.y * e2.x)); // Placeholder - will pass as a uniform precalculated
        return true;
    }

    return false;
}

vec3 perp(vec3 v)
{
    vec3 b = cross(v, vec3(0, 0, 1));

    if (dot(b, b) < 0.01) {
        b = cross(v, vec3(0, 1, 0));
    }

    return b;
}

bool rayCylinderIntersectionT(Ray r, Cylinder cy, out float t)
{
    vec3 A = cy.axis * cy.min + cy.center;
    vec3 B = cy.axis * cy.max + cy.center;
    float extent = distance(A, B);
    vec3 W = (B - A) / extent;
    vec3 U = perp(W);
    vec3 V = cross(U, W);
    U = normalize(cross(V, W));
    V = normalize(V);
    float rSqr = cy.radius * cy.radius;
    vec3 diff = r.origin - 0.5 * (A + B);
    mat3 basis = mat3(U, V, W);
    vec3 P = diff * basis;
    float dz = dot(W, r.direction);

    if (abs(dz) >= 1.0 - EPSILON) {
        float radialSqrtDist = rSqr - P.x * P.x - P.y * P.y;

        if (radialSqrtDist < 0.0) { return false; }

        t = (dz > 0.0 ? -P.z : P.z) + extent * 0.5;
        return true;
    }

    vec3 D = vec3(dot(U, r.direction), dot(V, r.direction), dz);
    float a0 = P.x * P.x + P.y * P.y - rSqr;
    float a1 = P.x * D.x + P.y * D.y;
    float a2 = D.x * D.x + D.y * D.y;
    float discr = a1 * a1 - a0 * a2;

    if (discr < 0.0) {
        return false;
    }

    if (discr > EPSILON) {
        float root = sqrt(discr);
        float inv = 1.0 / a2;
        t = (-a1 - root) * inv;
        return true;
    }

    t = -a1 / a2;
    return true;
}

bool rayPlaneIntersection(Ray r, vec3 planeNormal, vec3 planePoint, out float t)
{
    float denom = dot(planeNormal, r.direction);
    vec3 po = planePoint - r.origin;
    t = dot(po, planeNormal) / denom;
    return t >= 0.0;
}

bool rayDiskPlaneIntersection(Ray r, vec3 planeNormal, vec3 planePoint, float radius, out float t)
{
    if (rayPlaneIntersection(r, planeNormal, planePoint, t)) {
        vec3 p = r.origin + r.direction * t;
        vec3 v = p - planePoint;
        float d2 = dot(v, v);
        return (sqrt(d2) <= radius);
    }

    return false;
}

bool intersection(Ray r, Cylinder cy, inout vec3 position, inout vec3 normal)
{
    vec3 A = cy.axis * cy.min + cy.center;
    vec3 B = cy.axis * cy.max + cy.center;
    vec3 buttomCapPos, topCapPos, sidePos, sideNormal;
    float tTop, tButtom, tSide, t;
    bool iButtomCap, iTopCap, iSide;
    t = tTop = tButtom = tSide = MAX_DELTA;
    iButtomCap = iTopCap = iSide = false;

    // Test intersection against cylinder caps
    if (rayDiskPlaneIntersection(r, cy.axis, A, cy.radius, t)) {tButtom = t; iButtomCap = true;}

    if (rayDiskPlaneIntersection(r, cy.axis, B, cy.radius, t)) {tTop = t; iTopCap = true;}

    if (rayCylinderIntersectionT(r, cy, t)) {
        vec3 sidePos = r.origin + t * r.direction;

        if (dot(sidePos - B, cy.axis) < 0.0 && dot(sidePos - A, cy.axis) > 0.0) {
            iSide = true;
            tSide = t;
        }
    }

    if (iSide || iTopCap || iButtomCap) {
        t = min(tSide, min(tButtom, tTop));
        position = r.origin + t * r.direction;

        if (t == tButtom || t == tTop) {
            normal = cy.axis;
        } else {
            vec3 v = (B - A) / distance(A, B);
            float tn = dot(position - A, v);
            vec3 spinePoint = A + tn * v;
            normal = normalize(position - spinePoint);
        }

        return true;
    }

    return false;
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
    float att = 1.0 / (1.0 + (0.5 * d) + (0.25 * d * d)); // Light Attenuation
    vec3 ambient = vec3(0.1, 0.1, 0.1);
    vec3 diffuse = m.diffuse * l.color * cosTheta;
    vec3 specular = m.specular * l.color * pow(cosAlpha, 10.0);
    return vec4(ambient + diffuse + specular, m.alpha) * att;
}