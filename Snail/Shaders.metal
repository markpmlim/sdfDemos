#include <metal_stdlib>

using namespace metal;

//https://www.shadertoy.com/view/ld3Gz2

// Created by inigo quilez - iq/2015
// I share this piece (art and code) here in Shadertoy and through its Public API, only for educational purposes. 
// You cannot use, share or host this piece or modifications of it as part of your own commercial or non-commercial product, website or project.
// You can share a link to it or an unmodified screenshot of it provided you attribute "by Inigo Quilez, @iquilezles and iquilezles.org". 
// If you are a techer, lecturer, educator or similar and these conditions are too restrictive for your needs, please contact me and we'll work it out.

// You can buy a metal print of this shader here:
// https://www.redbubble.com/i/metal-print/Snail-by-InigoQuilez/39845499.0JXQP

//#define ZERO (min(iFrame, 0))

// antialiasing - make AA 2, meaning 4x AA, if you have a fast machine
#define AA 2

constexpr sampler texSampler(mip_filter::linear,
                             mag_filter::linear,
                             min_filter::linear);

// http://www.iquilezles.org/www/articles/distfunctions/distfunctions.htm
float sdSphere(float3 p, float4 s )
{
    // centre=s.xyz, radius=s.w
    return length(p-s.xyz) - s.w;
}

// http://www.iquilezles.org/www/articles/distfunctions/distfunctions.htm
// c=centre, r=radius
float sdEllipsoid( float3 p, float3 c, float3 r )
{
    return (length( (p-c)/r ) - 1.0) * min(min(r.x,r.y),r.z);
}

// http://www.iquilezles.org/www/articles/distfunctions/distfunctions.htm
// c=centre, r=radius
float sdCircle( float2 p, float2 c, float r )
{
    return length(p-c) - r;
}

// http://www.iquilezles.org/www/articles/distfunctions/distfunctions.htm
float sdTorus( float3 p, float2 t )
{
    return length( float2(length(p.xz) - t.x,p.y) ) - t.y;
}

// http://www.iquilezles.org/www/articles/distfunctions/distfunctions.htm
float sdCapsule( float3 p, float3 a, float3 b, float r )
{
	float3 pa = p-a, ba = b-a;
	float h = clamp(dot(pa,ba)/dot(ba,ba),
                    0.0, 1.0 );
	return length( pa - ba*h ) - r;
}

// http://www.iquilezles.org/www/articles/distfunctions/distfunctions.htm
float2 udSegment( float3 p, float3 a, float3 b )
{
	float3 pa = p-a, ba = b-a;
	float h = clamp(dot(pa,ba)/dot(ba,ba),
                    0.0, 1.0 );
	return float2( length( pa - ba*h ), h );
}

// http://research.microsoft.com/en-us/um/people/hoppe/ravg.pdf
// Determinant of 2x2 matrix
float det( float2 a, float2 b )
{
    return a.x*b.y-b.x*a.y;
}

float3 getClosest( float2 b0, float2 b1, float2 b2 )
{
    float a =     det(b0,b2);
    float b = 2.0*det(b1,b0);
    float d = 2.0*det(b2,b1);
    float f = b*d - a*a;
    float2  d21 = b2-b1;
    float2  d10 = b1-b0;
    float2  d20 = b2-b0;
    float2  gf = 2.0*(b*d21+d*d10+a*d20); gf = float2(gf.y,-gf.x);
    float2  pp = -f*gf/dot(gf,gf);
    float2  d0p = b0-pp;
    float ap = det(d0p,d20);
    float bp = 2.0*det(d10,d0p);
    float t = clamp( (ap+bp)/(2.0*a+b+d), 0.0 ,1.0 );
    return float3( mix(mix(b0,b1,t), mix(b1,b2,t),t), t );
}

float4 sdBezier( float3 a, float3 b, float3 c, float3 p )
{
	float3 w = normalize( cross( c-b, a-b ) );  // forward vector
	float3 u = normalize( c-b );                // right vector
	float3 v =          ( cross( w, u ) );      // up vector

	float2 a2 = float2(dot(a-b, u),
                       dot(a-b, v));
	float2 b2 = float2( 0.0 );
	float2 c2 = float2(dot(c-b, u),
                       dot(c-b, v));
	float3 p3 = float3(dot(p-b, u),
                       dot(p-b, v),
                       dot(p-b, w));

	float3 cp = getClosest( a2-p3.xy, b2-p3.xy, c2-p3.xy );

	return float4(sqrt(dot(cp.xy, cp.xy) + p3.z*p3.z),
                  cp.z,
                  length(cp.xy),
                  p3.z);
}

// http://iquilezles.org/www/articles/smin/smin.htm
float smin( float a, float b, float k )
{
	float h = clamp(0.5 + 0.5*(b-a)/k,
                    0.0, 1.0 );
	return mix( b, a, h ) - k*h*(1.0 - h);
}

// http://iquilezles.org/www/articles/smin/smin.htm
float smax( float a, float b, float k )
{
	float h = clamp(0.5 + 0.5*(b-a)/k,
                    0.0, 1.0);
	return mix( a, b, h ) + k*h*(1.0-h);
}

// http://iquilezles.org/www/articles/smin/smin.htm
float3 smax( float3 a, float3 b, float k )
{
	float3 h = clamp(0.5 + 0.5*(b-a)/k,
                     0.0, 1.0);
	return mix( a, b, h ) + k*h*(1.0 - h);
}

//---------------------------------------------------------------------------

float hash1( float n )
{
    return fract(sin(n)*43758.5453123);
}

float3 hash3( float n )
{
    return fract(sin(n + float3(0.0, 13.1, 31.3))*158.5453123);
}

float3 forwardSF( float i, float n)
{
    const float PI  = 3.141592653589793238;
    const float PHI = 1.618033988749894848;
    float phi = 2.0*PI*fract(i/PHI);
    float zi = 1.0 - (2.0*i + 1.0)/n;
    float sinTheta = sqrt(1.0 - zi*zi);
    return float3(cos(phi)*sinTheta,
                  sin(phi)*sinTheta,
                  zi);
}

//---------------------------------------------------------------------------


constant const float pi = 3.1415927;

//---------------------------------------------------------------------------

float mapShell(float3 p, thread float4& matInfo,
               texture2d<float> iChannel2 )
{
    const float sc = 1.0/1.0;
    p -= float3(0.05, 0.12, -0.09);

    p *= sc;

    float3 q = float3x3(-0.6333234236, -0.7332753384, 0.2474039592,
                         0.7738444477, -0.6034162289, 0.1924931824,
                         0.0081370606,  0.3133626215, 0.9495986813) * p;

    const float b = 0.1759;

    float r = length( q.xy );
    float t = atan2( q.y, q.x );

    // https://swiftcoder.wordpress.com/2010/06/21/logarithmic-spiral-distance-field/
    float np = (log(   r)/b-t)/(2.0*pi);
    float nm = (log(0.11)/b-t)/(2.0*pi);
    float n = min(np,nm);

    float ni = floor( n );

    float r1 = exp( b * (t + 2.0*pi*ni));
    float r2 = r1 * 3.019863;

    //-------

    float h1 = q.z + 1.5*r1 - 0.5;
    float d1 = sqrt((r1-r)*(r1-r) + h1*h1) - r1;
    float h2 = q.z + 1.5*r2 - 0.5;
    float d2 = sqrt((r2-r)*(r2-r) + h2*h2) - r2;

    float d, dx, dy;
    if( d1<d2 ) {
        d  = d1;
        dx = r1-r;
        dy = h1;
    }
    else {
        d  = d2;
        dx = r2 - r;
        dy = h2;
    }

    //float di = iChannel2.sample(texSampler, float2(t+r, 0.5)).x;
    float di = iChannel2.sample(texSampler, float2(t+r, 0.5), level(0.0) ).x;
    d += 0.002*di;

    matInfo = float4(dx, dy, r/0.4, t/pi);

    float3 s = q;
    q = q - float3(0.34, -0.1, 0.03);
    q.xy = float2x2( 0.8, 0.6,
                    -0.6, 0.8)*q.xy;
    d = smin( d,  sdTorus( q, float2(0.28,0.05) ), 0.06);
    d = smax( d, -sdEllipsoid(q, float3(0.0,0.0,0.0), float3(0.24,0.36,0.24) ), 0.03 );
    d = smax( d, -sdEllipsoid(s, float3(0.52,-0.0,0.0), float3(0.42,0.23,0.5) ), 0.05 );

    return d/sc;
}

float2 mapSnail(float3 p, thread float4& matInfo,
                float iTime, texture2d<float> iChannel1)
{
    const float3 head = float3(-0.76, 0.6, -0.3);

    float3 q = p - head;

    // body
#if 1
    float4 b1 = sdBezier(float3(-0.13,-0.65,0.0),
                         float3(0.24,0.9+0.1,0.0),
                         head+float3(0.04,0.01,0.0),
                         p);
    float d1 = b1.x;
    d1 -= smoothstep(0.0,0.2,b1.y)*(0.16 - 0.07*smoothstep(0.5,1.0,b1.y));
    b1 = sdBezier(float3(-0.085,0.0,0.0),
                  float3(-0.1,0.9-0.05,0.0),
                  head + float3(0.06,-0.08,0.0),
                  p);
    float d2 = b1.x;
    d2 -= 0.1 - 0.06*b1.y;
    d1 = smin( d1, d2, 0.03 );
    matInfo.xyz = b1.yzw;
#else
    float4 b1 = sdBezier(float3(-0.13,-0.65,0.0),
                         float3(0.24,0.9+0.11,0.0),
                         head+float3(0.05,0.01-0.02,0.0),
                         p);
    float d1 = b1.x;
    d1 -= smoothstep(0.0,0.2,b1.y)*(0.16 - 0.75*0.07*smoothstep(0.5,1.0,b1.y));
    matInfo.xyz = b1.yzw;
    float d2;
#endif
    d2 = sdSphere( q, float4(0.0,-0.06,0.0,0.085) );
    d1 = smin( d1, d2, 0.03 );
    d1 = smin( d1, sdSphere(p,float4(0.05,0.52,0.0,0.13)), 0.07 );

    q.xz = float2x2(0.8,0.6,-0.6,0.8)*q.xz;

    float3 sq = float3( q.xy, abs(q.z) );

    // top antennas
    float3 af = 0.05*sin(0.5*iTime+float3(0.0,1.0,3.0) + float3(2.0,1.0,0.0)*sign(q.z) );
    float4 b2 = sdBezier( float3(0.0), float3(-0.1,0.2,0.2), float3(-0.3,0.2,0.3)+af, sq );
    float d3 = b2.x;
    d3 -= 0.03 - 0.025*b2.y;
    d1 = smin( d1, d3, 0.04 );
    d3 = sdSphere( sq, float4(-0.3,0.2,0.3,0.016) + float4(af,0.0) );
    d1 = smin( d1, d3, 0.01 );    

    // bottom antenas
    float3 bf = 0.02*sin(0.3*iTime+float3(4.0,1.0,2.0) + float3(3.0,0.0,1.0)*sign(q.z) );
    float2 b3 = udSegment( sq, float3(0.06,-0.05,0.0), float3(-0.04,-0.2,0.18)+bf );
    d3 = b3.x;
    d3 -= 0.025 - 0.02*b3.y;
    d1 = smin( d1, d3, 0.06 );
    d3 = sdSphere( sq, float4(-0.04,-0.2,0.18,0.008)+float4(bf,0.0) );
    d1 = smin( d1, d3, 0.02 );

    // bottom
    float3 pp = p-float3(-0.17,0.15,0.0);
    float co = 0.988771078;
    float si = 0.149438132;
    pp.xy = float2x2(co,-si,si,co)*pp.xy;
    d1 = smin( d1, sdEllipsoid( pp, float3(0.0,0.0,0.0), float3(0.084,0.3,0.15) ), 0.05 );
    d1 = smax( d1, -sdEllipsoid( pp, float3(-0.08,-0.0,0.0), float3(0.06,0.55,0.1) ), 0.02 );

    // disp
    float dis = iChannel1.sample(texSampler, 5.0*p.xy, level(0.0) ).x;
    float dx = 0.5 + 0.5*(1.0-smoothstep(0.5,1.0,b1.y));
    d1 -= 0.005*dis*dx*0.5;

    return float2(d1,1.0);
}
    
float mapDrop( float3 p )
{
    p -= float3(-0.26,0.25,-0.02);
    p.x -= 2.5*p.y*p.y;
    return sdCapsule( p, float3(0.0,-0.06,0.0), float3(0.014,0.06,0.0), 0.037 );
}

float mapLeaf( float3 p )
{
    p -= float3(-1.8,0.6,-0.75);

    p = float3x3( 0.671212, 0.366685, -0.644218,
                 -0.479426, 0.877583,  0.000000,
                  0.565354, 0.308854,  0.764842) * p;
 
    p.y += 0.2*exp(-abs(2.0*p.z) );

    float ph = 0.25*50.0*p.x - 0.25*75.0*abs(p.z);// + 1.0*sin(5.0*p.x)*sin(5.0*p.z);
    float rr = sin( ph );
    rr = rr*rr;    
    rr = rr*rr;    
    p.y += 0.005*rr;

    float r = clamp((p.x+2.0)/4.0,
                    0.0, 1.0);
    r = 0.0001 + r*(1.0-r)*(1.0-r)*6.0;

    rr = sin( ph*2.0 );
    rr = rr*rr;    
    rr *= 0.5+0.5*sin( p.x*12.0 );

    float ri = 0.035*rr;

    float d = sdEllipsoid( p, float3(0.0), float3(2.0,0.25*r,r+ri) );

    float d2 = p.y - 0.02;
    d = smax( d, -d2, 0.02 );

    return d;
}

float2 mapOpaque(float3 p, thread float4& matInfo, float iTime,
                 texture2d<float> iChannel1, texture2d<float> iChannel2)
{
    matInfo = float4(0.0);

   	//--------------
    float2 res = mapSnail( p, matInfo, iTime, iChannel1 );

    //---------------
    float4 tmpMatInfo;
    float d4 = mapShell( p, tmpMatInfo, iChannel2);
    if( d4<res.x  ) {
        res = float2(d4, 2.0);
        matInfo = tmpMatInfo;
    }

    //---------------

    float4 b3 = sdBezier(float3(-0.15,-1.5,0.0),
                         float3(-0.1,0.5,0.0),
                         float3(-0.6,1.5,0.0), p);
    d4 = b3.x;
    d4 -= 0.04 - 0.02*b3.y;
    if( d4<res.x  )
    {
        res = float2(d4,3.0);
    }
	
	//----------------------------
    
    float d5 = mapLeaf( p );
    if( d5<res.x )
        res = float2(d5,4.0);

    return res;
}

// http://iquilezles.org/www/articles/normalsSDF/normalsSDF.htm
float3 calcNormalOpaque( float3 pos, float eps, float iTime,
                        texture2d<float> iChannel1, texture2d<float> iChannel2 )
{
    float4 kk;
#if 0
    float2 e = float2(1.0,-1.0)*0.5773*eps;
    return normalize( e.xyy*mapOpaque( pos + e.xyy, kk, iTime, iChannel1, iChannel2 ).x +
					  e.yyx*mapOpaque( pos + e.yyx, kk, iTime, iChannel1, iChannel2 ).x +
					  e.yxy*mapOpaque( pos + e.yxy, kk, iTime, iChannel1, iChannel2 ).x +
					  e.xxx*mapOpaque( pos + e.xxx, kk, iTime, iChannel1, iChannel2 ).x );
#else
    // inspired by tdhooper and klems - a way to prevent the compiler from inlining map() 4 times
    float3 n = float3(0.0);
    for( int i=0; i<4; i++ )
    {
        float3 e = 0.5773*(2.0*float3((((i+3)>>1)&1),((i>>1)&1),(i&1))-1.0);
        n += e*mapOpaque(pos+eps*e,kk, iTime, iChannel1, iChannel2).x;
    }
    return normalize(n);
#endif    
}

//=========================================================================

float mapLeafWaterDrops( float3 p )
{
    p -= float3(-1.8,0.6,-0.75);
    float3 s = p;
    p = float3x3(0.671212, 0.366685, -0.644218,
            -0.479426, 0.877583,  0.000000,
             0.565354, 0.308854,  0.764842)*p;
  
    float3 q = p;
    p.y += 0.2*exp(-abs(2.0*p.z) );
    
    //---------------
    
    float r = clamp((p.x+2.0)/4.0,0.0,1.0);
    r = r*(1.0-r)*(1.0-r)*6.0;
    float d0 = sdEllipsoid( p, float3(0.0), float3(2.0,0.25*r,r) );
    float d1 = sdEllipsoid( q, float3(0.5,0.0,0.2), 1.0*float3(0.15,0.13,0.15) );
    float d2 = sdEllipsoid( q, float3(0.8,-0.07,-0.15), 0.5*float3(0.15,0.13,0.15) );
    float d3 = sdEllipsoid( s, float3(0.76,-0.8,0.6), 0.5*float3(0.15,0.2,0.15) );
    float d4 = sdEllipsoid( q, float3(-0.5,0.09,-0.2), float3(0.04,0.03,0.04) );

    d3 = max( d3, p.y-0.01);
    
    return min( min(d1,d4), min(d2,d3) );
}

float2 mapTransparent( float3 p, thread float4& matInfo )
{
    matInfo = float4(0.0);
    
    float d5 = mapDrop( p );
    float2  res = float2(d5,4.0);

    float d6 = mapLeafWaterDrops( p );
    res.x = min( res.x, d6 );

    return res;
}

// http://iquilezles.org/www/articles/normalsSDF/normalsSDF.htm
float3 calcNormalTransparent( float3 pos, float eps )
{
    float4 kk;
    float2 e = float2(1.0,-1.0)*0.5773*eps;
    return normalize( e.xyy*mapTransparent( pos + e.xyy, kk ).x + 
					  e.yyx*mapTransparent( pos + e.yyx, kk ).x + 
					  e.yxy*mapTransparent( pos + e.yxy, kk ).x + 
					  e.xxx*mapTransparent( pos + e.xxx, kk ).x );
}

//=========================================================================

float calcAO( float3 pos, float3 nor, float iTime,
             texture2d<float> iChannel1, texture2d<float> iChannel2 )
{
    float4 kk;
	float ao = 0.0;
    for( int i=0; i<32; i++ )
    {
        float3 ap = forwardSF( float(i), 32.0 );
        float h = hash1(float(i));
		ap *= sign( dot(ap,nor) ) * h*0.1;
        ao += clamp( mapOpaque( pos + nor*0.01 + ap, kk, iTime, iChannel1, iChannel2 ).x*3.0, 0.0, 1.0 );
    }
	ao /= 32.0;
	
    return clamp( ao*6.0, 0.0, 1.0 );
}

float calcSSS( float3 pos, float3 nor,
              float iTime,
              texture2d<float> iChannel1, texture2d<float> iChannel2 )
{
    float4 kk;
	float occ = 0.0;
    for( int i=0; i<8; i++ )
    {
        float h = 0.002 + 0.11*float(i)/7.0;
        float3 dir = normalize( sin( float(i)*13.0 + float3(0.0,2.1,4.2) ) );
        dir *= sign(dot(dir,nor));
        occ += (h-mapOpaque(pos-h*dir, kk, iTime, iChannel1, iChannel2).x);
    }
    occ = clamp( 1.0 - 11.0*occ/8.0, 0.0, 1.0 );    
    return occ*occ;
}

// http://iquilezles.org/www/articles/rmshadows/rmshadows.htm
float calcSoftShadow( float3 ro, float3 rd, float k,
                     float iTime,
                     texture2d<float> iChannel1, texture2d<float> iChannel2 )
{
    float4 kk;
    float res = 1.0;
    float t = 0.01;
    for( int i=0; i<32; i++ )
    {
        float h = mapOpaque(ro + rd*t, kk,
                            iTime, iChannel1, iChannel2).x;
        res = min(res, smoothstep(0.0, 1.0,
                                  k*h/t));
        t += clamp( h, 0.04, 0.1 );
		if( res<0.01 )
            break;
    }
    return clamp(res,
                 0.0, 1.0);
}

//constant const float3 sunDir = normalize( float3(0.2,0.1,0.02) );
// We cannot normalize vectors in program scope in Metal.
constant const float3 sunDir = float3(0.890871, 0.445435, 0.0890871);

// The textures must be passed as function parameters or as part of a struct.
// Returned values via function parameters must be prefixed with the keyword "thread"
// ro=ray origin, rd=(normalized) ray direction, t=parameter along the ray
// Return the color
float3 shadeOpaque(float3 ro, float3 rd, float t, float m, thread float4& matInfo, float iTime,
                   texture2d<float> iChannel1, texture2d<float> iChannel2, texture2d<float> iChannel3)
{
    float eps = 0.002;

    float3 pos = ro + t*rd;         // point along the ray
    float3 nor = calcNormalOpaque(pos, eps,
                                  iTime, iChannel1, iChannel2);

    float3 mateD = float3(0.0);
    float3 mateS = float3(0.0);
    float2 mateK = float2(0.0);
    float3 mateE = float3(0.0);

    float focc = 1.0;
    float fsha = 1.0;

    if( m<1.5 ) // snail body
    {
        //float dis = texture( iChannel1, 5.0*pos.xy ).x;
        float dis = iChannel1.sample(texSampler, 5.0*pos.xy).x;
        float be = sdEllipsoid( pos, float3(-0.3,-0.5,-0.1), float3(0.2,1.0,0.5) );
        be = 1.0 - smoothstep( -0.01, 0.01, be );
        
        float ff = abs(matInfo.x-0.20);
        
        mateS = 6.0*mix(0.7*float3(2.0,1.2,0.2),    // color0
                        float3(2.5,1.8,0.9),        // color1
                        ff);
        mateS += 2.0*dis;
        mateS *= 1.5;
        mateS *= 1.0 + 0.5*ff*ff;
        mateS *= 1.0-0.5*be;

        mateD = float3(1.0,0.8,0.4);
        mateD *= dis;
        mateD *= 0.015;
        mateD += float3(0.8,0.4,0.3)*0.15*be;

        mateK = float2( 60.0, 0.7 + 2.0*dis );

        float f = clamp( dot( -rd, nor ), 0.0, 1.0 );
        f = 1.0 - pow( f, 8.0 );
        //f = 1.0 - (1.0-f)*(1.0-texture( iChannel2, 0.3*pos.xy ).x);
        f = 1.0 - (1.0-f)*(1.0-iChannel2.sample(texSampler, 0.3*pos.xy ).x);
        mateS *= float3(0.5,0.1,0.0) + f*float3(0.5,0.9,1.0);

        float b = 1.0-smoothstep( 0.25,0.55,abs(pos.y));
        focc = 0.2 + 0.8*smoothstep(0.0, 0.15,
                                    sdSphere(pos, float4(0.05,0.52,0.0, 0.13)));
    }
	else if( m<2.5 ) // shell
    {
        mateK = float2(0.0);

        float tip = 1.0-smoothstep(0.05,0.4, length(pos-float3(0.17,0.2,0.35)) );
        mateD = mix( 0.7*float3(0.2,0.21,0.22), 0.2*float3(0.15,0.1,0.0), tip );

        float2 uv = float2( .5*atan2(matInfo.x,matInfo.y)/3.1416, 1.5*matInfo.w );
        
        //float3 ral = texture( iChannel1, float2(2.0*matInfo.w+matInfo.z*0.5,0.5) ).xxx;
        float3 ral = iChannel1.sample(texSampler, float2(2.0*matInfo.w+matInfo.z*0.5, 0.5) ).xxx;
        mateD *= 0.25 + 0.75*ral;
        
        float pa = smoothstep(-0.2,0.2, 0.3+sin(2.0+40.0*uv.x + 3.0*sin(11.0*uv.x)) );
        float bar = mix(pa,1.0,smoothstep(0.7,1.0,tip));
        bar *= (matInfo.z<0.6) ? 1.0 : smoothstep( 0.17, 0.21, abs(matInfo.w)  );
        mateD *= float3(0.06,0.03,0.0)+float3(0.94,0.97,1.0)*bar;
        
        mateK = float2( 64.0, 0.2 );
        mateS = 1.5*float3(1.0,0.65,0.6) * (1.0-tip);//*0.5;
    }
    else if( m<3.5 ) // plant
    {
        mateD = float3(0.05,0.1,0.0)*0.2;
        mateS = float3(0.1,0.2,0.02)*25.0;
        mateK = float2(5.0,1.0);
        
        float fre = clamp(1.0+dot(nor,rd), 0.0, 1.0 );
        mateD += 0.2*fre*float3(1.0,0.5,0.1);
        
        //float3 te = texture( iChannel2, pos.xy*0.2 ).xyz;
        float3 te = iChannel2.sample(texSampler, pos.xy*0.2 ).xyz;
        mateS *= 0.5 + 1.5*te;
        mateE = 0.5*float3(0.1,0.1,0.03)*(0.2+0.8*te.x);
    }
    else //if( m<4.5 ) // leave
    {
        float3 p = pos - float3(-1.8,0.6,-0.75);
        float3 s = p;
        p = float3x3(0.671212, 0.366685, -0.644218,
                -0.479426, 0.877583,  0.000000,
                 0.565354, 0.308854,  0.764842)*p;

        float3 q = p;
        p.y += 0.2*exp(-abs(2.0*p.z) );

        float v = smoothstep( 0.01, 0.02, abs(p.z));
        
        float rr = sin( 4.0*0.25*50.0*p.x - 4.0*0.25*75.0*abs(p.z) );

        //float3 te = texture( iChannel2, p.xz*0.35 ).xyz;
        float3 te = iChannel2.sample(texSampler, p.xz*0.35 ).xyz;

        float r = clamp((p.x+2.0)/4.0,0.0,1.0);
        r = r*(1.0-r)*(1.0-r)*6.0;
        float ff = length(p.xz/float2(2.0,r));

        mateD = mix( float3(0.07,0.1,0.0), float3(0.05,0.2,0.01)*0.25, v );
        mateD = mix( mateD, float3(0.16,0.2,0.01)*0.25, ff );
        mateD *= 1.0 + 0.25*te;
        mateD *= 0.8;
        
        mateS = float3(0.15,0.2,0.02)*0.8;
        mateS *= 1.0 + 0.2*rr;
        mateS *= 0.8;

        mateK = float2(64.0,0.25);
        
        //---------------------
        
        //nor.xz += v*0.15*(-1.0+2.0*texture( iChannel3, 1.0*p.xz ).xy);
        nor.xz += v*0.15*(-1.0+2.0*iChannel3.sample(texSampler, 1.0*p.xz ).xy);
        nor = normalize( nor );

        float d1 = sdEllipsoid( q, float3( 0.5-0.07, 0.0,  0.20), 1.0*float3(1.4*0.15,0.13,0.15) );
        float d2 = sdEllipsoid( q, float3( 0.8-0.05,-0.07,-0.15), 0.5*float3(1.3*0.15,0.13,0.15) );
        float d4 = sdEllipsoid( q, float3(-0.5-0.07, 0.09,-0.20), 1.0*float3(1.4*0.04,0.03,0.04) );
        float dd = min(d1,min(d2,d4));
        fsha = 0.05 + 0.95*smoothstep(0.0,0.05,dd);
        
        d1 = abs( sdCircle( q.xz, float2( 0.5, 0.20), 1.0*0.15 ));
        d2 = abs( sdCircle( q.xz, float2( 0.8,-0.15), 0.5*0.15 ));
        d4 = abs( sdCircle( q.xz, float2(-0.5,-0.20), 1.0*0.04 ));
        dd = min(d1,min(d2,d4));
        focc *= 0.55 + 0.45*smoothstep(0.0,0.08,dd);
        
        d1 = distance( q.xz, float2( 0.5-0.07, 0.20) );
        d2 = distance( q.xz, float2( 0.8-0.03,-0.15) );
        fsha += (1.0-smoothstep(0.0,0.10,d1))*1.5;
        fsha += (1.0-smoothstep(0.0,0.05,d2))*1.5;    
    }
    
  
    float3 hal = normalize( sunDir-rd );
    float fre = clamp(1.0+dot(nor,rd), 0.0, 1.0 );
    float occ = calcAO( pos, nor, iTime, iChannel1, iChannel2 )*focc;
    float sss = calcSSS( pos, nor, iTime, iChannel1, iChannel2  );
    sss = sss*occ + fre*occ + (0.5+0.5*fre)*pow(abs(matInfo.x-0.2),1.0)*occ;
    
    float dif1 = clamp( dot(nor,sunDir), 0.0, 1.0 );
    float sha = calcSoftShadow( pos, sunDir, 20.0, iTime, iChannel1, iChannel2 );
    dif1 *= sha*fsha;
    float spe1 = clamp( dot(nor,hal), 0.0, 1.0 );

    float bou = clamp( 0.3-0.7*nor.y, 0.0, 1.0 );

    // illumination
    
    float3 col = float3(0.0);
    col += 7.0*float3(1.7,1.2,0.6)*dif1*2.0;            // sun
    col += 4.0*float3(0.2,1.2,1.6)*occ*(0.5+0.5*nor.y); // sky
    col += 1.8*float3(0.1,2.0,0.1)*bou*occ;             // bounce

    col *= mateD;

    col += 0.4*sss*(float3(0.15,0.1,0.05)+float3(0.85,0.9,0.95)*dif1)*(0.05+0.95*occ)*mateS; // sss
    col = pow(col,float3(0.6,0.8,1.0));
    
    col += float3(1.0,1.0,1.0)*0.2*pow( spe1, 1.0+mateK.x )*dif1*(0.04+0.96*pow(fre,4.0))*mateK.x*mateK.y;   // sun lobe1
    col += float3(1.0,1.0,1.0)*0.1*pow( spe1, 1.0+mateK.x/3.0 )*dif1*(0.1+0.9*pow(fre,4.0))*mateK.x*mateK.y; // sun lobe2
	col += 0.1*float3(1.0,max(1.5-0.7*col.y,0.0),2.0)*occ*occ*smoothstep( 0.0, 0.3, reflect( rd, nor ).y )*mateK.x*mateK.y*(0.04+0.96*pow(fre,5.0)); // sky

    col += mateE;

    return col;        
}

float3 shadeTransparent(float3 ro, float3 rd,
                        float t, float m, thread float4& matInfo,
                        float3 col, float depth,
                        texture2d<float> iChannel1)
{
    float3 oriCol = col;
    
    float dz = depth - t;
    float ao = clamp(dz*50.0,0.0,1.0);
    float3  pos = ro + t*rd;
    float3  nor = calcNormalTransparent( pos, 0.002 );
    float fre = clamp( 1.0 + dot( rd, nor ), 0.0, 1.0 );
    float3  hal = normalize( sunDir-rd );
    float3  ref = reflect( -rd, nor );
    float spe1 = clamp( dot(nor,hal), 0.0, 1.0 );
    float spe2 = clamp( dot(ref,sunDir), 0.0, 1.0 );


    float ds = 1.6 - col.y;
    
    col *= mix( float3(0.0,0.0,0.0), float3(0.4,0.6,0.4), ao );

    col += ds*1.5*float3(1.0,0.9,0.8)*pow( spe1, 80.0 );
    col += ds*0.2*float3(0.9,1.0,1.0)*smoothstep(0.4,0.8,fre);
    col += ds*0.9*float3(0.6,0.7,1.0)*smoothstep( -0.5, 0.5, -reflect( rd, nor ).y )*smoothstep(0.2,0.4,fre);
    col += ds*0.5*float3(1.0,0.9,0.8)*pow( spe2, 80.0 );
    col += ds*0.5*float3(1.0,0.9,0.8)*pow( spe2, 16.0 );
    //col += float3(0.8,1.0,0.8)*0.5*smoothstep(0.3,0.6, texture( iChannel1, 0.8*nor.xy ).x)*(0.1+0.9*fre*fre);
    col += float3(0.8,1.0,0.8)*0.5*smoothstep(0.3,0.6, iChannel1.sample(texSampler, 0.8*nor.xy ).x)*(0.1+0.9*fre*fre);
    
    // hide aliasing a bit
    return mix( col, oriCol, smoothstep(0.6,1.0,fre) ); 
}

//--------------------------------------------

float2 intersectOpaque(float3 ro, float3 rd,
                       const float mindist, const float maxdist,
                       thread float4& matInfo, float iTime,
                       texture2d<float> iChannel1, texture2d<float> iChannel2)
{
    float2 res = float2(-1.0);
    
    float t = mindist;
    for( int i=0; i<64; i++ )
    {
        float3 p = ro + t*rd;
        float2 h = mapOpaque( p, matInfo, iTime, iChannel1, iChannel2 );
        res = float2(t,h.y);

        if( h.x<(0.001*t) ||  t>maxdist )
            break;
        
        t += h.x*0.9;
    }
	return res;
}

float2 intersectTransparent(float3 ro, float3 rd,
                            const float mindist, const float maxdist,
                            thread float4& matInfo )
{
    float2 res = float2(-1.0);

    float t = mindist;
    for( int i=0; i<64; i++ )
    {
        float3 p = ro + t*rd;
        float2 h = mapTransparent( p, matInfo );
        res = float2(t,h.y);

        if( h.x<(0.001*t) ||  t>maxdist )
            break;

        t += h.x;
    }
	return res;
}

// ray direction and background texture.
float3 background(float3 d, texture2d<float> iChannel2)
{
    // cheap cubemap
    float3 n = abs(d);
    float2 uv = (n.x>n.y && n.x>n.z) ? d.yz/d.x:
              (n.y>n.x && n.y>n.z) ? d.zx/d.y:
                                     d.xy/d.z;
    
    // fancy blur
    float3  col = float3( 0.0 );
    for( int i=0; i<200; i++ )
    {
        float h = float(i)/200.0;
        float an = 31.0*6.2831*h;
        float2  of = float2( cos(an), sin(an) ) * h;

        //float3 tmp = texture( iChannel2, uv*0.25 + 0.0075*of, 4.0 ).yxz;
        float3 tmp = iChannel2.sample(texSampler, uv*0.25 + 0.0075*of, 4.0 ).yxz;
        col = smax( col, tmp, 0.5 );
    }
    
    return pow(col, float3(3.5,3.0,6.0))*0.2;
}


float3 render(float3 ro, float3 rd, float2 q, float iTime,
              texture2d<float> iChannel1, texture2d<float> iChannel2, texture2d<float> iChannel3)
{
    //-----------------------------

    float3 col = background( rd, iChannel2 );
    
    //-----------------------------

    float mindist = 1.0;
    float maxdist = 4.0;

    float4 matInfo;
    float2 tm = intersectOpaque(ro, rd, mindist, maxdist, matInfo,
                                iTime, iChannel1, iChannel2);
    if( tm.y>-0.5 && tm.x < maxdist )
    {
        col = shadeOpaque(ro, rd, tm.x, tm.y, matInfo, iTime,
                          iChannel1, iChannel2, iChannel3);
        maxdist = tm.x;
    }

    //-----------------------------

    tm = intersectTransparent( ro, rd, mindist, maxdist, matInfo );
    if( tm.y>-0.5 && tm.x < maxdist )
    {
        col = shadeTransparent( ro, rd, tm.x, tm.y, matInfo, col, maxdist, iChannel1 );
    }

    //-----------------------------

    float sun = clamp(dot(rd, sunDir), 0.0,1.0);
    col += 1.0*float3(1.5,0.8,0.7)*pow(sun,4.0);

    //-----------------------------

    col = pow( col, float3(0.45) );

    col = float3(1.05,1.0,1.0)*col*(0.7+0.3*col*max(3.0-2.0*col,0.0)) + float3(0.0,0.0,0.04);

    col *= 0.3 + 0.7*pow(16.0*q.x*q.y*(1.0-q.x)*(1.0-q.y),0.1);

    return clamp( col, 0.0, 1.0 );
}

// ro=ray origin, rt=target
float3x3 setCamera( float3 ro, float3 rt )
{
    // Compute vector from target to origin
    float3 w = normalize(ro - rt);
    float m = sqrt(1.0 - w.y*w.y);
    return float3x3( w.z,   0.0,  -w.x,
                     0.0,   m*m,  -w.z*w.y,
                     w.x*m, w.y*m, w.z*m );
}

typedef struct
{
    float4 renderedCoordinate [[position]]; // clip space
    float2 textureCoordinate;
} TextureMappingVertex;


// Projects provided vertices to corners of offscreen texture.
vertex TextureMappingVertex
vertexShader(unsigned int vertex_id  [[ vertex_id ]])
{
    // Triangle strip in NDC (normalized device coords).
    // The vertices' coord system has (0, 0) at the bottom left.
    float4x4 renderedCoordinates = float4x4(float4(-1.0, -1.0, 0.0, 1.0),   /// (x, y, z, w)
                                            float4( 1.0, -1.0, 0.0, 1.0),
                                            float4(-1.0,  1.0, 0.0, 1.0),
                                            float4( 1.0,  1.0, 0.0, 1.0));
    // The texture coord system has (0, 0) at the upper left
    // The s-axis is +ve right and the t-axis is +ve down
    float4x2 textureCoordinates = float4x2(float2(0.0, 1.0),                /// (s, t)
                                           float2(1.0, 1.0),
                                           float2(0.0, 0.0),
                                           float2(1.0, 0.0));
    TextureMappingVertex outVertex;
    outVertex.renderedCoordinate = renderedCoordinates[vertex_id];
    outVertex.textureCoordinate = textureCoordinates[vertex_id];
    return outVertex;
}

fragment float4
fragmentShader(TextureMappingVertex             mappingVertex   [[stage_in]],
               texture2d<float, access::sample> input1          [[texture(0)]],
               texture2d<float, access::sample> input2          [[texture(1)]],
               texture2d<float, access::sample> input3          [[texture(2)]],
               constant float2                  &u_resolution   [[buffer(0)]],
               constant float                   &u_time         [[buffer(1)]],
               constant float2                  &u_mouse        [[buffer(2)]])
{
    // The pixel's window-relative coordinates are given by TextureMappingVertex's
    // [[position]] attribute. This is identical to the OpenGL gl_FragCoord.
    // cf Table 5.4: Attributes for fragment function input arguments
    float2 fragCoord = mappingVertex.renderedCoordinate.xy;
    float iTime = 0.01*u_time;   // adjust this value slow/speed up the animation

#if AA<2
    // Map [0, w) and [0, h) to [-1.0, 1.0)
    float2  p = (2.0*fragCoord-u_resolution.xy)/u_resolution.y;
    // The origin of Metal's 2D pixel coord system is at the upper-left corner.
    p.y = -p.y;
    float2  q = fragCoord/u_resolution.xy;
    float an = 1.87 - 0.04*(1.0 - cos(0.5*iTime));

    // Init the ray origin
    float3  ro = float3(-0.4,0.2,0.0) + 2.2*float3(cos(an), 0.0, sin(an));
    // Set the target (lookAt)
    float3  ta = float3(-0.6, 0.2, 0.0);
    // Compute the camera matrix
    float3x3  ca = setCamera( ro, ta );
    // Compute the ray direction and normalize it.
    float3  rd = normalize(ca * float3(p, -2.8));

    float3 col = render(ro, rd, q,
                        iTime, input1, input2, input3);
#else
    // Antialiasing using the concept of multi-sampling
    float3 col = float3(0.0);
    for( int m=0; m<AA; m++ ) {
        for( int n=0; n<AA; n++ )
        {
            float2 rr = float2(float(m), float(n))/float(AA);

            // Map [0, w) and [0, h) to [-1.0, 1.0)
            float2 p = (2.0*(fragCoord+rr)-u_resolution.xy)/u_resolution.y;
            // The origin of Metal's 2D pixel coord system is at the upper-left corner.
            p.y = -p.y;
            float an = 1.87 - 0.04*(1.0 - cos(0.5*iTime));
            float2 q = (fragCoord+rr)/u_resolution.xy;

            float3 ro = float3(-0.4,0.2,0.0) + 2.2*float3(cos(an), 0.0, sin(an));
            float3 ta = float3(-0.6,0.2,0.0);
            float3x3 ca = setCamera( ro, ta );
            float3 rd = normalize( ca * float3(p,-2.8) );

            col += render(ro, rd, q,
                          iTime, input1, input2, input3);
        } // inner
    } // outer
    col /= float(AA*AA);
#endif

    float4 fragColor = float4( col, 1.0 );

    // Output the fragment
    return fragColor;
}
