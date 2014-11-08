#version 400
uniform vec3 iResolution;
#define PI 3.1415926535897932384626433832795
#define MAX_DELTA 1e20
#define EPSILON 1e-5
out vec4 fragColor;			// Fragment Shader Output

// DO NOT MODIFY THIS SECTION BELOW - SCENE.CPP WRITES HERE
// BEGIN:SCENEWRITER
// -- Structures from Scene.h
struct Material {
    vec3 diffuse;
    vec3 specular;
    float refractiveIndex;
    float reflectiveIndex;
    float alpha;
    int id;
};
struct Ray {
    vec3 origin;
    vec3 direction;
};
struct Sphere {
    Material mat;
    vec3 center;
    float radius;
    int id;
};
struct Cube {
    Material mat;
    vec2 facesX;
    vec2 facesY;
    vec2 facesZ;
    int id;
};
struct Cylinder {
    Material mat;
    vec3 center;
    vec3 axis;
    float radius;
    float min;
    float max;
    int id;
};
struct Triangle {
    Material mat;
    vec3 point1;
    vec3 point2;
    vec3 point3;
    int id;
};
struct Camera {
    vec3 position;
    vec3 direction;
};
struct Light {
    vec3 position;
    vec3 color;
    float intensity;
    int id;
};
struct Intersection {
    vec3 position;
    vec3 normal;
    float dist;
    Material mat;
    int hit_id;
};
// --
const int NUM_BOUNCES = 256;
// -- Input from .yml
// Materials
const int NUM_MATERIALS = 3;
Material materials[NUM_MATERIALS];
// Figures
const int NUM_SPHERES = 4;
const int NUM_TRIANGLES = 1;
const int NUM_CYLINDERS = 2;
Sphere spheres[NUM_SPHERES];
Triangle triangles[NUM_TRIANGLES];
Cylinder cylinders[NUM_CYLINDERS];
// Lights
const int NUM_LIGHTS = 3;
Light lights[NUM_LIGHTS];
// Scene Initializer
void initScene()
{
    // BEGIN:SCENEOBJECTS
    // -- YML MATERIALS
    materials[0] = Material(vec3(0.0, 0.0, 1.0), vec3(0.5, 1.0, 0.5), 0.1, 0.0, 0.7, 1);
    materials[1] = Material(vec3(0.0, 1.0, 0.0), vec3(0.5, 1.0, 0.5), 0.0, 1.33, 0.5, 2);
    materials[2] = Material(vec3(1.0, 0.0, 0.0), vec3(1.0, 1.0, 1.0), 0.2, 0.0, 0.3, 3);
    // --
    // YML SPHERES
    spheres[0] = Sphere(materials[0], vec3(1.0, 0.0, 1.5), 1.0, 1);
    spheres[1] = Sphere(materials[2], vec3(0.0, 0.0, -1.5), 1.0, 2);
    spheres[2] = Sphere(materials[1], vec3(2.0, 3.0, 0.0), 0.75, 3);
    spheres[3] = Sphere(materials[0], vec3(-2.0, -3.0, 0.0), 0.5, 4);
    // YML CYLINDERS
    cylinders[0] = Cylinder(materials[0], vec3(4.5, -1.0, 1.0), vec3(0.0, 1.0, 0.0), 1.0, 2.0, -2.0, 5);
    cylinders[1] = Cylinder(materials[1], vec3(-3.0, 0.0, -1.0), vec3(1.0, 0.0, 0.0), 1.0, 1.0, -1.0, 6);
    // --
    // --
    // YML LIGHTS
    lights[0] = Light(vec3(0.0, 5.0, 8.0), vec3(1.0, 0.0, 0.0), 1.0, 1);
    lights[1] = Light(vec3(-8.0, -10.0, 2.0), vec3(0.0, 1.0, 0.0), 1.7, 2);
    lights[2] = Light(vec3(-5.0, 5.0, 1.0), vec3(0.0, 0.0, 1.0), 1.5, 3);
    // --
    // GARBAGE - NEGATIVE ID DISCARD THESE OBJECTS
    triangles[0].id = -1;
    // --
    // END::SCENEOBJECTS
}
// --
// END:SCENEWRITER
const vec4 color_ambient = vec4(vec3(0.01), 1.0);
const Material ambient = Material(color_ambient.xyz, color_ambient.xyz, 0.0, 0.0, 1.0, 0);

// Some relevant functions
vec2 transform(vec2 p);						// Clip Coordinates Transform
vec3 perp(vec3 v);							// Perpendicular Vector
Material blend(Material m1, Material m2);	// Alpha Blending

// Figures Intersections
bool intersection(Ray sRay, Sphere sph, inout Intersection point);		// Sphere Intersection
bool intersection(Ray sRay, Cylinder cyl, inout Intersection point);	// Finite Cylinder Intersection
bool intersection(Ray sRay, Triangle tri, inout Intersection point);	// Triangle Intersection

// -- Used for finite cilinder intersection
bool rayCylinderIntersectionT(Ray sRay, Cylinder cyl, inout float t);										// Infinite Cylinder Intersection
bool rayPlaneIntersection(Ray sRay, vec3 planeNormal, vec3 planePoint, inout float t);						// Plane Intersection (No Far-Near)
bool rayDiskPlaneIntersection(Ray sRay, vec3 planeNormal, vec3 planePoint, float radius, inout float t);	// Planar Circle Intersection
// --

// Lighting Models
vec4 phong(Ray sRay, Light light, Intersection point);
vec4 cooktorrance(Ray sRay, Light light, Intersection point);

// Ray-Trace
vec4 trace(inout Ray inRay);
void intersectAll(Ray inRay, inout Intersection point);

void main()
{
    vec4 color = vec4(0.1, 0.1, 0.1, 0.0);
    Camera camera;
    camera.position = vec3(0.0, 0.0, 5.0);
    camera.direction = vec3(0.0, 0.0, -1.0);
    Ray sRay;
    sRay.origin = camera.position;
    sRay.direction = normalize(vec3(transform(gl_FragCoord.xy), 0.0) + camera.direction);
    // Init Scene
    initScene();
    // Start Ray-Tracing
    fragColor = trace(sRay);
}

void intersectAll(Ray inRay, inout Intersection point)
{
    for (int i = 0; i < NUM_SPHERES; i++) {
        if (spheres[i].id != -1) { intersection(inRay, spheres[i], point); }
    }

    for (int i = 0; i < NUM_TRIANGLES; i++) {
        if (triangles[i].id != -1) { intersection(inRay, triangles[i], point); }
    }

    for (int i = 0; i < NUM_CYLINDERS; i++) {
        if (cylinders[i].id != -1) { intersection(inRay, cylinders[i], point); }
    }
}

vec4 trace(inout Ray inRay)
{
    vec4 color_final = vec4(0.01);
    Intersection point;
    point.hit_id = 0;
    point.mat = ambient;
    point.normal = vec3(0.0, 0.0, 0.0);
    point.position = vec3(0.0, 0.0, 0.0);
    point.dist = MAX_DELTA;

    for (int i = 0; i < NUM_BOUNCES; i++) {
        intersectAll(inRay, point);

        // Break loop if ray didn't hit anything
        if (point.dist  >= MAX_DELTA) {
            break;
        }

        for (int l = 0; l < NUM_LIGHTS; l++) {
            color_final += cooktorrance(inRay, lights[l], point);
        }

        if (point.mat.reflectiveIndex <= 0.0 && point.mat.refractiveIndex <= 0.0) {
            break;
        }

        if (point.mat.reflectiveIndex > 0.0) {
            inRay.origin = point.position;
            inRay.direction = normalize(reflect(inRay.direction, point.normal));
        }

        if (point.mat.refractiveIndex > 0.0) {
            inRay.origin = point.position;
            inRay.direction = normalize(refract(inRay.direction, point.normal, point.mat.refractiveIndex));
        }

        point.dist = MAX_DELTA;
    }

    return color_final;
}

bool intersection(Ray sRay, Sphere sph, inout Intersection point)
{
    if (point.hit_id == sph.id) { return false; } // Avoid Self Lighting

    vec3 oc = sRay.origin - sph.center;
    float b = 2.0 * dot(oc, sRay.direction);
    float c = dot(oc, oc) - sph.radius * sph.radius;
    float h = b * b - 4.0 * c;

    if (h < 0.0) { return false; }

    float t  = (-b - sqrt(h)) / 2.0;
    point.mat = blend(sph.mat, point.mat);

    if (t > point.dist) { return false; }

    point.dist = t;
    point.position = sRay.origin + point.dist * sRay.direction;
    point.hit_id = sph.id;
    point.normal = (point.position - sph.center) / sph.radius;
    return true;
}

bool intersection(Ray sRay, Triangle tri, inout Intersection point)
{
    if (point.hit_id == tri.id) { return false; }

    // Möller–Trumbore Algorithm
    vec3 e1 = tri.point2 - tri.point1;
    vec3 e2 = tri.point3 - tri.point1;
    vec3 p = cross(sRay.direction, e2);
    float det = dot(e1, p);

    if (det > -EPSILON && det < EPSILON) { return false; }

    float invDet = 1.0 / det;
    vec3 tvec = sRay.origin - tri.point1;
    float u = dot(tvec, p) * invDet;

    if (u < 0.0 || u > 1.0) { return false; }

    vec3 q = cross(tvec, e1);
    float v = dot(sRay.direction, q) * invDet;

    if (v < 0.0 || u + v > 1.0) { return false; }

    float t  = dot(e2, q) * invDet;
    point.mat = blend(tri.mat, point.mat);

    if (t > point.dist) { return false; }

    if (t > EPSILON) {
        point.dist = t;
        point.position = sRay.origin + t * sRay.direction;
        point.hit_id = tri.id;
        point.normal = normalize(cross(e1, e2));
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

bool rayCylinderIntersectionT(Ray sRay, Cylinder cyl, inout float t)
{
    vec3 A = cyl.axis * cyl.min + cyl.center;
    vec3 B = cyl.axis * cyl.max + cyl.center;
    float extent = distance(A, B);
    vec3 W = (B - A) / extent;
    vec3 U = perp(W);
    vec3 V = cross(U, W);
    U = normalize(cross(V, W));
    V = normalize(V);
    float rSqr = cyl.radius * cyl.radius;
    vec3 diff = sRay.origin - 0.5 * (A + B);
    mat3 basis = mat3(U, V, W);
    vec3 P = diff * basis;
    float dz = dot(W, sRay.direction);

    if (abs(dz) >= 1.0 - EPSILON) {
        float radialSqrtDist = rSqr - P.x * P.x - P.y * P.y;

        if (radialSqrtDist < 0.0) { return false; }

        t = (dz > 0.0 ? -P.z : P.z) + extent * 0.5;
        return true;
    }

    vec3 D = vec3(dot(U, sRay.direction), dot(V, sRay.direction), dz);
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

bool rayPlaneIntersection(Ray sRay, vec3 planeNormal, vec3 planePoint, inout float t)
{
    float denom = dot(planeNormal, sRay.direction);
    vec3 po = planePoint - sRay.origin;
    t = dot(po, planeNormal) / denom;
    return t >= 0.0;
}

bool rayDiskPlaneIntersection(Ray sRay, vec3 planeNormal, vec3 planePoint, float radius, inout float t)
{
    if (rayPlaneIntersection(sRay, planeNormal, planePoint, t)) {
        vec3 p = sRay.origin + sRay.direction * t;
        vec3 v = p - planePoint;
        float d2 = dot(v, v);
        return (sqrt(d2) <= radius);
    }

    return false;
}

bool intersection(Ray sRay, Cylinder cyl, inout Intersection point)
{
    if (point.hit_id == cyl.id) { return false; }

    vec3 A = cyl.axis * cyl.min + cyl.center;
    vec3 B = cyl.axis * cyl.max + cyl.center;
    vec3 buttomCapPos, topCapPos, sidePos, sideNormal;
    float tTop, tButtom, tSide, t;
    bool iButtomCap, iTopCap, iSide;
    t = tTop = tButtom = tSide = MAX_DELTA;
    iButtomCap = iTopCap = iSide = false;

    // Test intersection against cylinder caps
    if (rayDiskPlaneIntersection(sRay, cyl.axis, A, cyl.radius, t)) {tButtom = t; iButtomCap = true;}

    if (rayDiskPlaneIntersection(sRay, cyl.axis, B, cyl.radius, t)) {tTop = t; iTopCap = true;}

    if (rayCylinderIntersectionT(sRay, cyl, t)) {
        vec3 sidePos = sRay.origin + t * sRay.direction;

        if (dot(sidePos - B, cyl.axis) > 0.0 && dot(sidePos - A, cyl.axis) < 0.0) {
            iSide = true;
            tSide = t;
        }
    }

    if (iSide || iTopCap || iButtomCap) {
        t = min(tSide, min(tButtom, tTop));
        point.mat = blend(cyl.mat, point.mat);

        if (t > point.dist) { return false; }

        point.position = sRay.origin + t * sRay.direction;
        point.dist = t;
        point.hit_id = cyl.id;

        if (t == tButtom || t == tTop) {
            point.normal = cyl.axis;
        } else {
            vec3 v = (B - A) / distance(A, B);
            float tn = dot(point.position - A, v);
            vec3 spinePoint = A + tn * v;
            point.normal = normalize(point.position - spinePoint);
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

vec4 phong(Ray sRay, Light light, Intersection point)
{
    vec3 lightDirection = normalize(light.position - point.position);
    float d = length(lightDirection);
    vec3 R = reflect(-lightDirection, point.normal);
    vec3 E = normalize(sRay.origin - point.position);
    float cosAlpha = clamp(dot(E, R), 0.0, 1.0);
    float cosTheta = clamp(dot(point.normal, lightDirection), 0.0, 1.0);
    float att = 1.0 / (1.0 + (0.5 * d) + (0.25 * d * d)); // Light Attenuation
    vec3 ambient = vec3(0.01, 0.01, 0.01);
    vec3 diffuse = point.mat.diffuse * cosTheta;
    vec3 specular = point.mat.specular * pow(cosAlpha, 10.0);
    vec3 finalValue = (ambient + diffuse + specular) * light.color * att * light.intensity;
    return vec4(finalValue, point.mat.alpha);
}

vec4 cooktorrance(Ray sRay, Light light, Intersection point)
{
    // Attenuation factor based on the distance between the point on the surface and light position
    vec3 L = light.position - point.position;
    float dist = length(L);
    float A = (1.0 + light.intensity) / dist;
    // Specify surface material values
    vec3 MCOL = point.mat.diffuse + point.mat.specular;
    float R = 0.3;
    float F = 1.0;
    float K = 0.8;
    // Specify some other variables
    L = normalize(L);
    vec3 N = point.normal;
    vec3 V = -sRay.direction;
    vec3 LCOL = light.color * light.intensity;
    float NdotL = max(dot(N, L), 0.0);
    float specular = 0.0;

    if (NdotL > 0.0) {
        // Calculate intermediary values
        vec3 H = normalize(V + L);
        float NdotH = max(dot(N, H), 0.0);
        float NdotV = max(dot(N, V), 0.0);
        float VdotH = max(dot(V, H), 0.0);
        float RR = R * R;
        // Geometric attenuation
        float NH2 = 2.0 * NdotH;
        float g1 = (NH2 * NdotV) / VdotH;
        float g2 = (NH2 * NdotL) / VdotH;
        float geo = min(1.0, min(g1, g2));
        // Roughness (Microfacet/Beckmann distribution function)
        float r1 = 1.0 / (4.0 * RR * pow(NdotH, 4.0));
        float r2 = (NdotH * NdotH - 1.0) / (RR * NdotH * NdotH);
        float rough = r1 * exp(r2);
        // Fresnel (Schlick approximation)
        float fresnel = pow(1.0 - VdotH, 5.0);
        fresnel *= (1.0 - F);
        fresnel += F;
        // Final specular highlight
        specular = (fresnel * geo * rough) / (NdotV * NdotL);
    }

    return vec4(A * LCOL * NdotL * (K * MCOL + specular * (1.0 - K)), point.mat.alpha);
}

Material blend(Material m1, Material m2)
{
    Material final_mat;
    final_mat.diffuse = m2.diffuse * m2.alpha + m1.diffuse * (1.0 - m2.alpha);
    final_mat.specular = m2.specular * m2.alpha + m1.specular * (1.0 - m2.alpha);
    final_mat.refractiveIndex = m1.refractiveIndex;
    final_mat.reflectiveIndex = m1.reflectiveIndex;
    final_mat.alpha = m1.alpha;
    final_mat.id = m1.id;
    return final_mat;
}