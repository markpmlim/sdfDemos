// https://www.shadertoy.com/view/3ssfRr
// Conell Box Soft Shadows
#include <metal_stdlib>

using namespace metal;

constant const int RAY_STEPS = 256;

#define GI 1

struct Intersection
{
    float t;
    float3 color;
    float3 p;
    int object;

    Intersection(float parameter, float3 col, float3 pos, int hitObj) {
        t = parameter;
        color = col;
        p = pos;
        object = hitObj;
    }
};
    
float3 rotateY(float3 p, float a) {
    return float3(cos(a) * p.x + sin(a) * p.z, p.y, -sin(a) * p.x + cos(a) * p.z);
}

float3 rotateX(float3 p, float amt) {
    return float3(p.x, cos(amt) * p.y - sin(amt) * p.z, sin(p.y) + cos(p.z));
}

float3 rotateZ(float3 p, float amt) {
    
    return float3(cos(amt)* p.x - sin(amt)*p.y, cos(amt) * p.y + sin(amt) * p.x, p.z);
}


void raycast(float2 uv, thread float3& dir,
             thread float3& eye, thread float3& ref,
             float2 iResolution) {
    eye = float3(0.0, 5.5, -20.0);
    ref = float3(0.0, 2.5, 0.0);
    
    float len = tan(3.14159 * 0.125) * distance(eye, ref);
    float3 H = normalize(cross(float3(0.0, 1.0, 0.0), ref - eye));
    float3 V = normalize(cross(H, eye - ref));
    V *= len;
    H *= len * iResolution.x / iResolution.y;
    float3 p = ref + uv.x * H + uv.y * V;
    // Normalized ray direction
    dir = normalize(p - eye);
}

float sphere(float3 p, float r, float3 c)
{
    return distance(p, c) - r;
}

float torus(float3 p, float2 t)
{
  float2 q = float2(length(p.xz)-t.x,p.y);
  return length(q) - t.y;
}

float cone(float3 p, float2 c)
{
    // c must be normalized
    float q = length(p.xy);
    return dot(c,float2(q,p.z));
}

// Box with side lengths b
float box(float3 p, float3 b)
{
  return length(max(abs(p) - b, 0.0));
}

float plane( float3 p, float4 n )
{
  // n must be normalized
  return dot(p,n.xyz) + n.w;
}



void sceneMap3D(float3 pos, thread float& t, thread int& obj, float3 lightPos)
{
    t = plane(pos+float3(5.0,2.5,0.0), float4(1.0, 0.0, 0.0, 0.0));
    float t2;
    obj = 0; // 0 is center sphere
    if((t2 = plane(pos+float3(-5.0,2.5,0.0), float4(-1.0, 0.0, 0.0, 0.0))) < t) {
        t = t2;
        obj = 1;
    }
    float t3;
    if((t3 = plane(pos+float3(0.0, -7.5, 0.0), float4(0.0, -1.0, 0.0, 0.0))) < t) {
        t = t3;
        obj = 2;
    }
    float t4;
    if((t4 = plane(pos+float3(0.0,2.5,0.0), float4(0.0,1.0,0.0,0.0))) < t) {
    	t = t4;
        obj = 3;
    }
    float t5;
    if((t5 = plane(pos+float3(0.0,2.5,-5.0), float4(0.0,0.0,-1.0,0.0))) < t) {
        t = t5;
        obj = 4;
    }
    float t6;
    if((t6 = box(rotateY(pos+float3(-2, 1, -0.75),-0.3054), float3(1.5,1.5,1.5))) < t) {
        t = t6;
        obj = 5;
    }
    float t7;
    if((t7 = box(rotateY(pos+float3(2,0,-3), 0.480), float3(1.5,3.0,1.5))) < t) {
        t = t7;
        obj = 6;
    }
}


float shadowMap3D(float3 pos)
{
    float t = box(rotateY(pos+float3(-2, 1, -0.75),-0.3054), float3(1.5,1.5,1.5));
    
    float t7;
    if((t7 = box(rotateY(pos+float3(2,0,-3), 0.480), float3(1.5,3.0,1.5))) < t) {
        t = t7;
    }
    
    return t;
}

float sceneMap3D(float3 pos, float3 lightPos)
{
    thread float t;
    thread int obj;
    sceneMap3D(pos, t, obj, lightPos);
    return t;
}

void march(float3 origin, float3 dir, thread float& t, thread int& hitObj, float3 lightPos)
{
    t = 0.001;
    for(int i = 0; i < RAY_STEPS; ++i)
    {
        float3 pos = origin + t * dir;
    	thread float m;
        sceneMap3D(pos, m, hitObj, lightPos);
        if(m < 0.01)
        {
            return;
        }
        t += m;
    }
    t = -1.0;
    hitObj = -1;
}

float3 computeNormal(float3 pos, float3 lightPos)
{
    float3 epsilon = float3(0.0, 0.001, 0.0);
    return normalize( float3( sceneMap3D(pos + epsilon.yxx, lightPos) - sceneMap3D(pos - epsilon.yxx, lightPos),
                            sceneMap3D(pos + epsilon.xyx, lightPos) - sceneMap3D(pos - epsilon.xyx, lightPos),
                            sceneMap3D(pos + epsilon.xxy, lightPos) - sceneMap3D(pos - epsilon.xxy, lightPos)));
}



float softShadow(float3 lightDir, float3 isect, float max_t, float3 lightPos) {
    float res = 1.0;
    for(float t = 0.5; t < max_t; /**/) {
        float m = shadowMap3D(isect + t * lightDir);
        if(m < 0.0001) {
            return 0.0; 
        }
        res = min(res, 4.0 * m / t);
        t += m;
    }
    return res;
}

float3 computeMaterial(int hitObj, float3 p, float3 n, float3 light, float3 view) {
    float t;
    
    switch(hitObj) {
        case 0:
            // Red wall
            return float3(0.63, 0.065, 0.05) * max(0.0, dot(n, light));
            break;
        case 1:
            // Green wall
            return float3(0.14, 0.45, 0.091) * max(0.0, dot(n, light));
            break;
            
        case 2:
        case 3:
        case 4:
        case 5:
        case 6:
            // floor
            return float3(0.85, 0.81, 0.78) * max(0.0, dot(n, light));
            break;
            
        case -1:
            // Background
            return float3(0.0);
            break;
    }
    return float3(1.0);
}

Intersection sdf3D(float3 dir, float3 eye, float3 lightPos)
{
    float t;
    int hitObj;
    march(eye, dir, t, hitObj, lightPos);

    float3 isect = eye + t * dir;
    float3 nor = computeNormal(isect, lightPos);
    float3 lightDir = normalize(lightPos - isect);
    
    float3 surfaceColor = computeMaterial(hitObj, isect, nor, lightDir, normalize(eye - isect));
    
    float shade = softShadow(lightDir, isect, length(lightPos - isect), lightPos);
    
    return Intersection(t, shade * surfaceColor, isect, hitObj);
}

kernel void
compute(texture2d<float, access::write> output      [[texture(0)]],
        constant float                  &u_time     [[buffer(0)]],
        constant float2                 &u_mouse    [[buffer(1)]],
        uint2                           gid         [[thread_position_in_grid]])
{
    uint width = output.get_width();
    uint height = output.get_height();
    uint column = gid.x;
    uint row = gid.y;
    if ((column >= width) || (row >= height))
    {
        // In case the size of the texture does not match the size of the grid.
        // Return early if the pixel is out of bounds
        return;
    }
    // We don't have to pass the resolution as a parameter.
    float2 u_resolution = float2(width, height);
    float2 fragCoord = float2(gid.x, gid.y);

    // Normalized pixel coordinates (from 0.0 to 1.0)
    float2 uv = fragCoord/u_resolution.xy;
    // [0.0, 1.0] --> [-1.0, 1.0]
    float2 uv2 = 2.0 * uv - float2(1.0);
    uv2.y = -uv2.y;

    float3 lightPos = float3(0.0, 7.45, 0.0);
    
    Intersection aaIsects[4] = {
        Intersection(0, float3(0), float3(0), -1),
        Intersection(0, float3(0), float3(0), -1),
        Intersection(0, float3(0), float3(0), -1),
        Intersection(0, float3(0), float3(0), -1)
    };
    float3 dir, eye, ref;
    int idx = 0;
    for(float i = 0.0; i < 1.0; i += 0.5) {
        for(float j = 0.0; j < 1.0; j += 0.5) {
            raycast(uv2 + float2(i, j) / u_resolution.xy, dir, eye, ref, u_resolution);
            aaIsects[idx++] = sdf3D(dir, eye, lightPos);
        }
    }
    float3 avgColor = float3(0.0);
    for(int i = 0; i < 4; ++i) {
        avgColor += aaIsects[i].color;
    }
    avgColor *= 0.25;
    float4 fragColor = float4(avgColor, 1.0);

    // Output to texture
    output.write(fragColor, gid);
}
