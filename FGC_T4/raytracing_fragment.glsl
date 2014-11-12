#version 400
#define NUM_BOUNCES 4
#define PI 3.1415926535897932384626433832795
#define MAX_DELTA 1e20
#define EPSILON 1e-5
#define SUBPIXEL_DIST 0.5
uniform vec3 iResolution;
uniform vec3 cameraPosition = vec3(0.0, 0.0, 10.0);
uniform vec3 cameraDirection = vec3(0.0, 0.0, -1.0);
uniform float enableShadows = 1;
uniform float enableAntiAliasing = 1;

layout(location = 0) out vec4 fragColor;			// Fragment Shader Output

// -- Structures from Scene.h
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
    vec3 color;
    float intensity;
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

// BEGIN:SCENEPARAMS
// END:SCENEPARAMS

const vec4 color_ambient = vec4(vec3(0.1), 1.0);
const Material ambient = Material(color_ambient.xyz, color_ambient.xyz, 0.0, 0.0, color_ambient.w);


// Some relevant functions
float inverseMix(float a, float b, float x);	// mix() inverse
vec2 transform(vec2 p);							// Clip Coordinates Transform
vec3 perp(vec3 v);								// Perpendicular Vector
vec4 blend(in vec4 fg, in vec4 bg);				// Alpha Blending

// Figures Intersections
bool intersectionSphere(Ray sRay, Sphere sph, inout Intersection point);		// Sphere Intersection
bool intersectionCylinder(Ray sRay, Cylinder cyl, inout Intersection point);	// Finite Cylinder Intersection
bool intersectionTriangle(Ray sRay, Triangle tri, inout Intersection point);	// Triangle Intersection

// -- Used for finite cilinder intersection
bool rayCylinderIntersectionT(Ray sRay, Cylinder cyl, inout float t);										// Infinite Cylinder Intersection
bool rayPlaneIntersection(Ray sRay, vec3 planeNormal, vec3 planePoint, inout float t);						// Plane Intersection (No Far-Near)
bool rayDiskPlaneIntersection(Ray sRay, vec3 planeNormal, vec3 planePoint, float radius, inout float t);	// Planar Circle Intersection
// --

// Lighting Models
vec3 phong(Ray sRay, Light light, Intersection point);
vec3 cooktorrance(Ray sRay, Light light, Intersection point);

// Ray-Trace
vec4 trace(in Ray inRay);
void intersectAll(Ray inRay, inout Intersection point);
vec4 sceneInterpolateLerp(Ray inRay);

// Shadow checks if something is in between the ray
bool shadow(Ray inRay, Intersection point, Light light, out Intersection shadowPoint);

// Super Sampling Functions
vec2 calcSubPixel(vec2 p, float x, float y);
void sampleSubPixels(Ray inr, vec2 frag, Camera cam, inout vec4 subPixel[9]);
void superSample(Ray inr, vec2 frag, Camera cam, inout vec4 color);

// Start Ray Trace Per Fragment
void main()
{
    vec4 color = vec4(0.0);
    vec2 fragCoord = gl_FragCoord.xy;
    vec2 uv = transform(fragCoord);
    Camera camera = Camera(cameraPosition, cameraDirection);
    Ray sRay = Ray(camera.position, normalize(vec3(uv, 0.0) + camera.direction), vec3(0.0), 1.0);
    // Init Scene
    initScene();

    if (enableAntiAliasing > 0.0) {
        superSample(sRay, fragCoord, camera, color);
    } else {
        color = trace(sRay);
    }

    fragColor = color;
}

vec4 blend(vec4 fg, vec4 bg)
{
    return fg * fg.a + bg * (1.0 - fg.a);
}

vec2 calcSubPixel(vec2 p, float x, float y)
{
    return transform(vec2(p.x - x, p.y + y));
}

void sampleSubPixels(Ray inr, vec2 frag, Camera cam, inout vec4 subPixel[9])
{
    int count = 0;

    for (int i = -1; i < 2; i++) {
        for (int j = -1; j < 2; j++) {
            inr.direction = normalize(vec3(calcSubPixel(frag, SUBPIXEL_DIST * i, SUBPIXEL_DIST * j), 0.0) + cam.direction);
            subPixel[count++] = trace(inr);
        }
    }
}

void superSample(Ray inr, vec2 frag, Camera cam, inout vec4 color)
{
    vec4 subPixel[9] = vec4[9](
                           vec4(0.0), vec4(0.0), vec4(0.0),
                           vec4(0.0), vec4(0.0), vec4(0.0),
                           vec4(0.0), vec4(0.0), vec4(0.0)
                       );
    sampleSubPixels(inr, frag, cam,  subPixel);

    for (int i = 0; i < 9; i++) { color += subPixel[i] / 9.0; }
}

void intersectAll(Ray inRay, inout Intersection point)
{
    for (int i = 0; i < NUM_SPHERES; i++) {
        if (spheres[i].id != -1) { intersectionSphere(inRay, spheres[i], point); }
    }

    for (int i = 0; i < NUM_CYLINDERS; i++) {
        if (cylinders[i].id != -1) { intersectionCylinder(inRay, cylinders[i], point); }
    }

    for (int i = 0; i < NUM_TRIANGLES; i++) {
        if (triangles[i].id != -1) { intersectionTriangle(inRay, triangles[i], point); }
    }
}

vec4 traceAlpha(in Ray firstRay, in Intersection firstPoint, vec4 fg)
{
    if (firstPoint.mat.alpha >= 1.0) { return vec4(0.0); }

    vec4 alphaColor = vec4(0.0);
    Ray alphaRay = firstRay;
    Intersection alphaPoint = firstPoint, shadowPoint;
    float alphaSdAcc = 1.0;
    vec4 bgColor = vec4(0.0);
    // Change Ray Origin to Closest Intersection
    alphaRay.origin = alphaPoint.position;
    // Intersect with the rest in the same direction
    intersectAll(alphaRay, alphaPoint);

    if (alphaPoint.dist < firstPoint.dist) {
        for (int l = 0; l < NUM_LIGHTS; l++) {
            bgColor.xyz += cooktorrance(alphaRay, lights[l], alphaPoint);

            if (enableShadows > 0.0) {
                if (shadow(alphaRay, alphaPoint, lights[l], shadowPoint)) {
                    alphaSdAcc *= 0.66;
                }
            }
        }

        alphaColor = blend(fg, bgColor * (0.5 + 0.5 * alphaSdAcc));
    } else {
        alphaColor = blend(fg, color_ambient);
    }

    return alphaColor;
}

vec4 trace(in Ray inRay)
{
    Ray firstRay = inRay;
    Intersection currentPoint, shadowPoint, firstPoint;
    currentPoint.hit_id = 0;
    currentPoint.mat = ambient;
    currentPoint.dist = MAX_DELTA;
    float shadowAcc = 1.0;
    bool isShadowed, firstPointSet = false;

    for (int i = 0; i < NUM_BOUNCES; i++) {
        intersectAll(inRay, currentPoint);

        // Break loop if ray didn't hit anything
        if (currentPoint.dist  >= MAX_DELTA) {
            vec4 currentColor = vec4(inRay.color * (0.5 + 0.5 * shadowAcc), firstPoint.mat.alpha);
            return firstPointSet ? blend(currentColor, traceAlpha(firstRay, firstPoint, currentColor)) : currentColor;
        }

        if (!firstPointSet) { firstPointSet = true; firstPoint = currentPoint; }

        for (int l = 0; l < NUM_LIGHTS; l++) {
            if (lights[l].id > 0) { inRay.color += cooktorrance(inRay, lights[l], currentPoint) * inRay.intensity; };

            if (enableShadows > 0.0) {
                if (shadow(inRay, currentPoint, lights[l], shadowPoint)) {
                    shadowAcc *= 0.66;
                }
            }
        }

        if (currentPoint.mat.reflectiveIndex <= 0.0 && currentPoint.mat.refractiveIndex <= 0.0) {
            break;
        }

        if (currentPoint.mat.reflectiveIndex > 0.0) {
            inRay.intensity *= currentPoint.mat.reflectiveIndex;
            inRay.origin = currentPoint.position;
            inRay.direction = normalize(reflect(inRay.direction, currentPoint.normal));
        }

        if (currentPoint.mat.refractiveIndex > 0.0) {
            inRay.origin = currentPoint.position;
            inRay.direction = normalize(refract(inRay.direction, currentPoint.normal, 1.0 / currentPoint.mat.refractiveIndex));
        }

        currentPoint.dist = MAX_DELTA;
        currentPoint.mat = ambient;
    }

    vec4 currentColor = vec4(inRay.color * (0.5 + 0.5 * shadowAcc), firstPoint.mat.alpha);
    vec4 alphaColor = traceAlpha(firstRay, firstPoint, currentColor);
    return blend(currentColor, alphaColor);
}

bool intersectionSphere(Ray sRay, Sphere sph, inout Intersection point)
{
    if (point.hit_id == sph.id) { return false; } // Avoid Self Lighting

    vec3 oc = sRay.origin - sph.center;
    float b = 2.0 * dot(oc, sRay.direction);
    float c = dot(oc, oc) - sph.radius * sph.radius;
    float h = b * b - 4.0 * c;

    if (h < 0.0) { return false; }

    float t  = (-b - sqrt(h)) / 2.0;

    if (t > point.dist) { return false; }

    point.mat = sph.mat;
    point.dist = t;
    point.position = sRay.origin + point.dist * sRay.direction;
    point.hit_id = sph.id;
    point.normal = (point.position - sph.center) / sph.radius;
    return true;
}

bool intersectionTriangle(Ray sRay, Triangle tri, inout Intersection point)
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

    if (t > point.dist) { return false; }

    if (t > EPSILON) {
        point.mat = tri.mat;
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

bool intersectionCylinder(Ray sRay, Cylinder cyl, inout Intersection point)
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
        point.mat = cyl.mat;

        if (t > point.dist) { return false; }

        point.position = sRay.origin + t * sRay.direction;
        point.dist = t;
        point.hit_id = cyl.id;

        if (t == tButtom) {
            point.normal = +cyl.axis;
        } else if (t == tTop) {
            point.normal = -cyl.axis;
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

vec3 phong(Ray sRay, Light light, Intersection point)
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
    return finalValue;
}

vec3 cooktorrance(Ray sRay, Light light, Intersection point)
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

    return A * LCOL * NdotL * (K * MCOL + specular * (1.0 - K));
}

bool shadow(Ray inRay, Intersection point, Light light, out Intersection shadowPoint)
{
    vec3 L = light.position - point.position;
    float distance = length(L);
    L = normalize(L);
    Ray shadowRay;
    Intersection lx = point;
    lx.dist = MAX_DELTA;
    shadowRay.origin = lx.position;
    shadowRay.direction = L;
    intersectAll(shadowRay, lx);

    if (lx.dist < MAX_DELTA && lx.dist < distance) {
        shadowPoint = lx;
        return true;
    }

    return false;
}