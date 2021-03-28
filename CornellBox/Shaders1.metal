// https://www.shadertoy.com/view/3d2yW3
//  Cornell Box (Raymarching)
#include <metal_stdlib>

using namespace metal;

constant const int RAY_STEPS = 256;
constant const float FOV = 19.5 * 3.14159 / 180.0;

struct Intersection {
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
    
float4x4 scale(float3 s) {
    return float4x4(
        float4(s.x, 0.0, 0.0, 0.0),
        float4(0.0, s.y, 0.0, 0.0),
        float4(0.0, 0.0, s.z, 0.0),
        float4(0.0, 0.0, 0.0, 1.0)
    );
}

float4x4 translate(float3 t){
    return float4x4(
        float4(1.0, 0.0, 0.0, 0.0),
        float4(0.0, 1.0, 0.0, 0.0),
        float4(0.0, 0.0, 1.0, 0.0),
        float4(t.x, t.y, t.z, 1.0)
    );
}

float4x4 rotateX(float theta){
    theta *= 3.14159 / 180.0;
    return float4x4(
        float4(1.,0.,0.,0),
        float4(0.,cos(theta),-sin(theta),0.),
        float4(0.,sin(theta),cos(theta),0.),
        float4(0.,0.,0.,1.));
}

float4x4 rotateY(float theta) {
    theta *= 3.14159 / 180.0;
    return float4x4(
        float4(cos(theta),0.,-sin(theta),0),
        float4(0.,1.,0.,0.),
        float4(sin(theta),0.,cos(theta),0.),
        float4(0.,0.,0.,1.));
}

float4x4 rotateZ(float theta){
    theta *= 3.14159 / 180.0;
    return float4x4(
        float4(cos(theta),-sin(theta),0.,0),
        float4(sin(theta),cos(theta),0.,0.),
        float4(0.,0.,1.,0.),
        float4(0.,0.,0.,1.));
}
    
float sphere(float3 p, float r, float3 c) {
    return distance(p, c) - r;
}


float dot2( float3 v ) {
    return dot(v,v);
}

float udQuad( float3 p, float3 a, float3 b, float3 c, float3 d ) {
  float3 ba = b - a; float3 pa = p - a;
  float3 cb = c - b; float3 pb = p - b;
  float3 dc = d - c; float3 pc = p - c;
  float3 ad = a - d; float3 pd = p - d;
  float3 nor = cross( ba, ad );

  return sqrt(
    (sign(dot(cross(ba,nor),pa)) +
     sign(dot(cross(cb,nor),pb)) +
     sign(dot(cross(dc,nor),pc)) +
     sign(dot(cross(ad,nor),pd))<3.0)
     ?
     min( min( min(
     dot2(ba*clamp(dot(ba,pa)/dot2(ba),0.0,1.0)-pa),
     dot2(cb*clamp(dot(cb,pb)/dot2(cb),0.0,1.0)-pb) ),
     dot2(dc*clamp(dot(dc,pc)/dot2(dc),0.0,1.0)-pc) ),
     dot2(ad*clamp(dot(ad,pd)/dot2(ad),0.0,1.0)-pd) )
     :
     dot(nor,pa)*dot(nor,pa)/dot2(nor) );
}

float4x4 inverse(float4x4 mat) {
    float det = determinant(mat);
    if (abs(det) < 0.005)
        return float4x4(1.0,0.0,0.0,0.0,
                        0.0,1.0,0.0,0.0,
                        0.0,0.0,1.0,0.0,
                        0.0,0.0,0.0,1.0);
    // Calculate 2x2 determinants
    float coef00 = mat[2][2] * mat[3][3] - mat[3][2] * mat[2][3];
    float coef02 = mat[1][2] * mat[3][3] - mat[3][2] * mat[1][3];
    float coef03 = mat[1][2] * mat[2][3] - mat[2][2] * mat[1][3];

    float coef04 = mat[2][1] * mat[3][3] - mat[3][1] * mat[2][3];
    float coef06 = mat[1][1] * mat[3][3] - mat[3][1] * mat[1][3];
    float coef07 = mat[1][1] * mat[2][3] - mat[2][1] * mat[1][3];

    float coef08 = mat[2][1] * mat[3][2] - mat[3][1] * mat[2][2];
    float coef10 = mat[1][1] * mat[3][2] - mat[3][1] * mat[1][2];
    float coef11 = mat[1][1] * mat[2][2] - mat[2][1] * mat[1][2];

    float coef12 = mat[2][0] * mat[3][3] - mat[3][0] * mat[2][3];
    float coef14 = mat[1][0] * mat[3][3] - mat[3][0] * mat[1][3];
    float coef15 = mat[1][0] * mat[2][3] - mat[2][0] * mat[1][3];

    float coef16 = mat[2][0] * mat[3][2] - mat[3][0] * mat[2][2];
    float coef18 = mat[1][0] * mat[3][2] - mat[3][0] * mat[1][2];
    float coef19 = mat[1][0] * mat[2][2] - mat[2][0] * mat[1][2];

    float coef20 = mat[2][0] * mat[3][1] - mat[3][0] * mat[2][1];
    float coef22 = mat[1][0] * mat[3][1] - mat[3][0] * mat[1][1];
    float coef23 = mat[1][0] * mat[2][1] - mat[2][0] * mat[1][1];

    float4 const signA(float4(+1, -1, +1, -1));
    float4 const signB(float4(-1, +1, -1, +1));

    // Calculate the cofactors
    float4 cofac0(coef00, coef00, coef02, coef03);
    float4 cofac1(coef04, coef04, coef06, coef07);
    float4 cofac2(coef08, coef08, coef10, coef11);
    float4 cofac3(coef12, coef12, coef14, coef15);
    float4 cofac4(coef16, coef16, coef18, coef19);
    float4 cofac5(coef20, coef20, coef22, coef23);

    float4 vec0(mat[1][0], mat[0][0], mat[0][0], mat[0][0]);
    float4 vec1(mat[1][1], mat[0][1], mat[0][1], mat[0][1]);
    float4 vec2(mat[1][2], mat[0][2], mat[0][2], mat[0][2]);
    float4 vec3(mat[1][3], mat[0][3], mat[0][3], mat[0][3]);

    float4 Inv0 = signA * (vec1 * cofac0 - vec2 * cofac1 + vec3 * cofac2);
    float4 Inv1 = signB * (vec0 * cofac0 - vec2 * cofac3 + vec3 * cofac4);
    float4 Inv2 = signA * (vec0 * cofac1 - vec1 * cofac3 + vec3 * cofac5);
    float4 Inv3 = signB * (vec0 * cofac2 - vec1 * cofac4 + vec2 * cofac5);

    float4x4 result(Inv0, Inv1, Inv2, Inv3);
    float OneOverDet = 1/det;
    result *= OneOverDet;
    return result;
}

float box(float3 pos, float3 t, float3 r, float3 s) {
    float4x4 worldInverse = inverse(translate(t) * rotateX(r.x) * rotateY(r.y) * rotateZ(r.z) * scale(s));
    float3 p = float3(worldInverse * float4(pos, 1));
    float sFactor = min(min(abs(s.x), abs(s.y)), abs(s.z));
    float dX = 0.0;
    float dY = 0.0;
    float dZ = 0.0;
    if (p.x > 0.5) {
        dX = p.x - 0.5;
    }
    else if (p.x < -0.5) {
        dX = -0.5 - p.x;
    }
    if (p.y > 0.5) {
        dY = p.y - 0.5;
    }
    else if (p.y < -0.5) {
        dY = -0.5 - p.y;
    }
    if (p.z > 0.5) {
        dZ = p.z - 0.5;
    }
    else if (p.z < -0.5) {
        dZ = -0.5 - p.z;
    }
    if (dX == 0.0 && dY == 0.0 && dZ == 0.0) {
        float xmin = min(0.5 - p.x, p.x + 0.5);
        float ymin = min(0.5 - p.y, p.y + 0.5);
        float zmin = min(0.5 - p.z, p.z + 0.5);
        return -min(min(xmin, ymin), zmin) * sFactor;
    }
    else {
        return sqrt(dX * dX + dY * dY + dZ * dZ) * sFactor;
    }
}

void rayCast(float2 uv, thread float3& dir, thread float3& eye, thread float3& ref,
             float2 iResolution) {
    eye = float3(0.0, 5.5, -30.0);
    ref = float3(0.0, 2.5, 0.0);
    
    float len = tan(FOV * 0.5) * distance(eye, ref);
    float3 H = normalize(cross(ref - eye, float3(0.0, 1.0, 0.0)));
    float3 V = normalize(cross(H, ref - eye));
    V *= len;
    H *= len * iResolution.x / iResolution.y;
    float3 p = ref + uv.x * H + uv.y * V;
    dir = normalize(p - eye);
}

float3 computeMaterial(int hitObj, float3 p, float3 d, float3 n) {
    switch(hitObj) {
        case 0: // Left Wall
        return float3(0.63, 0.065, 0.05);
        break;
        case 1: // Back Wall
        return float3(0.85, 0.81, 0.78);
        break;
        case 2: // Right Wall
        return float3(0.14, 0.45, 0.091);
        break;
        case 3: // Floor Wall
        return float3(0.85, 0.81, 0.78);
        break;
        case 4: // Ceiling
        return float3(0.85, 0.81, 0.78);
        break;
        case 5: // Long Cube
        return float3(0.85, 0.81, 0.78);
        break;
        case 6: // Short Cube
        return float3(0.85, 0.81, 0.78);
        break;
        case -1:
        return float3(0.0);
        break;
    }
    return float3(0.0);
}

void sceneMap3D(float3 pos, thread float& t, thread int& obj) {
    t = udQuad(pos, float3(5.0, 7.5, -5.0), float3(5.0, 7.5, 5.0), float3(5.0, -2.5, 5.0), float3(5.0, -2.5, -5.0));
    obj = 0; // Left Wall
    float t2;
    if ((t2 = box(pos, float3(0.0, 2.5, 10.0), float3(0.0, 0.0, 0.0), float3(10.0, 10.0, 10.0))) < t) {
        t = t2;
        obj = 1; // Back Wall
    }
    if ((t2 = udQuad(pos, float3(-5.0, 7.5, -5.0), float3(-5.0, 7.5, 5.0), float3(-5.0, -2.5, 5.0), float3(-5.0, -2.5, -5.0))) < t) {
        t = t2;
        obj = 2; // Right Wall
    }
    if ((t2 = udQuad(pos, float3(5.0, -2.5, -5.0), float3(-5.0, -2.5, -5.0), float3(-5.0, -2.5, 5.0), float3(5.0, -2.5, 5.0))) < t) {
        t = t2;
        obj = 3; // Floor
    }
    if ((t2 = udQuad(pos, float3(5.0, 7.5, -5.0), float3(-5.0, 7.5, -5.0), float3(-5.0, 7.5, 5.0), float3(5.0, 7.5, 5.0))) < t) {
        t = t2;
        obj = 4; // Ceiling
    }
    if ((t2 = box(pos, float3(2.0, 0.0, 3.0), float3(0.0, 27.5, 0.0), float3(3, 6, 3))) < t) {
        t = t2;
        obj = 5;// Long Cube
    }
    if ((t2 = box(pos, float3(-2.0, -1.0, 0.75), float3(0.0, -17.5, 0.0), float3(3, 3, 3))) < t) {
        t = t2;
        obj = 6; // Short Cube
    }
}

float sceneMap3D(float3 pos) {
    float t;
    int o;
    sceneMap3D(pos, t, o);
    return t;
}

void march(float3 origin, float3 dir, thread float& t, thread int& hitObj) {
    t = 0.001;
    for(int i = 0; i < RAY_STEPS; ++i)
    {
        float3 pos = origin + t * dir;
    	float m;
        sceneMap3D(pos, m, hitObj);
        if(m < 0.01)
        {
            return;
        }
        t += m;
    }
    t = -1.0;
    hitObj = -1;
}

float3 computeNormal(float3 pos) {
    float3 epsilon = float3(0.0, 0.001, 0.0);
    return normalize( float3( sceneMap3D(pos + epsilon.yxx) - sceneMap3D(pos - epsilon.yxx),
                            sceneMap3D(pos + epsilon.xyx) - sceneMap3D(pos - epsilon.xyx),
                            sceneMap3D(pos + epsilon.xxy) - sceneMap3D(pos - epsilon.xxy)));
}

Intersection sdf3D(float3 dir, float3 eye) {
    float t;
    int hitObj;
    march(eye, dir, t, hitObj);
    if (t == -1.0) {
        return Intersection(t, float3(0.0), eye + 1000.0 * dir, -1);
    }
    float3 isect = eye + t * dir;
    float3 nor = computeNormal(isect);
    float3 material = computeMaterial(hitObj, isect, dir, nor);
    float3 light = float3(0.0, 7.45, 0.0);
    float lambert = dot(normalize(light - isect), nor);
    float3 col = material * lambert;
   
    return Intersection(t, col, isect, hitObj);
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
    // Normalized pixel coordinates (from 0 to 1)
    float2 uv = fragCoord/u_resolution.xy;
    // [-1, 1]
    uv = 2.0 * uv - float2(1.0);
    
    float3 dir, eye, ref;
    rayCast(uv, dir, eye, ref, u_resolution);
    float3 col = sdf3D(dir, eye).color;
    
    // Output to screen
    float4 fragColor = float4(col,1.0);

    // Output to texture
    output.write(fragColor, gid);
}
