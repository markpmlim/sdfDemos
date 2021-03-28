// https://www.shadertoy.com/view/XsB3Rm

#include <metal_stdlib>

using namespace metal;

// ray marching
constant const int max_iterations = 512;
constant const float stop_threshold = 0.001;
constant const float grad_step = 0.02;
constant const float clip_far = 1000.0;

// math
constant const float PI = 3.14159265359;
constant const float DEG_TO_RAD = PI / 180.0;

// The ray's direction is a 3D vector so the texture must be a cubemap.
constexpr sampler texSampler(mip_filter::linear,
                              mag_filter::linear,
                              min_filter::linear);

float3 mod(float3 x, float3 y)
{
    return x - y * floor(x/y);
}

// iq's distance function
float sdSphere( float3 pos, float r )
{
	return length( pos ) - r;
}

float sdBox( float3 p, float3 b )
{
  float3 d = abs(p) - b;
  return min(max(d.x,max(d.y,d.z)),0.0) + length(max(d,0.0));
}


float sdUnion( float d0, float d1 )
{
    return min( d0, d1 );
}

float sdInter( float d0, float d1 )
{
    return max( d0, d1 );
}

float sdSub( float d0, float d1 )
{
    return max( d0, -d1 );
}

float sdUnion_s( float a, float b, float k )
{
    float h = clamp( 0.5+0.5*(b-a)/k, 0.0, 1.0 );
    return mix( b, a, h ) - k*h*(1.0-h);
}

float sfDisp( float3 p )
{
    return sin(p.x)*sin(p.y)*sin(p.z) ;
}

float3 sdTwist( float3 p, float a )
{
    float c = cos(a*p.y);
    float s = sin(a*p.y);
    float2x2  m = float2x2(c,-s, s,c);
    return float3(m*p.xz, p.y);
}

float3 sdRep( float3 p, float3 c )
{
    return mod(p,c) - 0.5*c;
}

// get distance in the world
float dist_field( float3 p )
{
//  p = sdRep( p, float3( 4.0 ) );
//  p = sdTwist( p, 3.0 );
    
    float d0 = sdBox( p, float3(0.5) );
    float d1 = sdSphere( p, 0.6 );
    
    float d = sdInter( d1, d0 );

    return d;
    //return d + sfDisp( p * 2.5 );
    //return sdUnion_s( d + sfDisp( p * 2.5 * sin( iTime * 1.01 ) ), d1, 0.1 );
}

// get gradient in the world
float3 gradient( float3 pos )
{
	const float3 dx = float3( grad_step, 0.0, 0.0 );
	const float3 dy = float3( 0.0, grad_step, 0.0 );
	const float3 dz = float3( 0.0, 0.0, grad_step );
	return normalize (
		float3(
			dist_field( pos + dx ) - dist_field( pos - dx ),
			dist_field( pos + dy ) - dist_field( pos - dy ),
			dist_field( pos + dz ) - dist_field( pos - dz )			
		)
	);
}

// The origin of Metal's 2D pixel coord system is at left-hand corner.
float3 getTextureColor2(texture2d<float> input, float3 dir) {
    float u = 0.5 + atan2(dir.z, -dir.x) / (2 * PI);
    float v = 0.5 - asin(-dir.y) / PI;
    float4 texColor = input.sample(texSampler, float2(u, v));
    return texColor.rgb;
}

constant float2 invAtan = float2(0.15915, 0.31831);   // 1/2π, 1/π

// Function to use an equirectangular map for environment lookups.
float3 getTextureColor(texture2d<float> input, float3 dir) {
    // tan(θ) = dir.z/-dir.x and sin(φ) = dir.y/1.0
    float2 uv = float2( atan2(dir.z, -dir.x),
                       asin(dir.y));
    // The range of uv.x: [ -π,   π ] --> [-0.5, 0.5]
    // The range of uv.y: [-π/2, π/2] --> [-0.5, 0.5]
    uv *= invAtan;
    uv += 0.5;          // [0, 1] for both u & v coords
    float4 texColor = input.sample(texSampler, uv);
    return texColor.rgb;
}

// Fresnel-Schlick approximation:
// F_Schlick(F_0, n, v) = F_0 + (1 - F_0)(1 - dot(n, v))^5
float3 fresnel( float3 F0, float3 h, float3 l )
{
	return F0 + ( 1.0 - F0 ) * pow( clamp( 1.0 - dot( h, l ), 0.0, 1.0 ), 5.0 );
}

// phong shading
float3 shading(float3 v, float3 n, float3 dir, float3 eye,
               texture2d<float> input )
{
	// ...add lights here...

	float shininess = 16.0;

	float3 final = float3( 0.0 );

    // Both dir and n are normalized vectors.
	float3 ref = reflect( dir, n );

    float3 Ks = float3( 0.5 );  // coefficients
    float3 Kd = float3( 1.0 );

	// light 0
	{
		float3 light_pos   = float3( 20.0, 20.0, 20.0 );
		float3 light_color = float3( 1.0, 0.7, 0.7 );

		float3 vl = normalize( light_pos - v );

		float3 diffuse  = Kd * float3( max( 0.0, dot( vl, n ) ) );
		float3 specular = float3( max( 0.0, dot( vl, ref ) ) );

        float3 F = fresnel( Ks, normalize( vl - dir ), vl );
		specular = pow( specular, float3( shininess ) );

		final += light_color * mix( diffuse, specular, F ); 
	}
	
	// light 1
	{
		float3 light_pos   = float3( -20.0, -20.0, -30.0 );
		float3 light_color = float3( 0.5, 0.7, 1.0 );

		float3 vl = normalize( light_pos - v );

		float3 diffuse  = Kd * float3( max( 0.0, dot( vl, n ) ) );
		float3 specular = float3( max( 0.0, dot( vl, ref ) ) );

        float3 F = fresnel( Ks, normalize( vl - dir ), vl );
		specular = pow( specular, float3( shininess ) );

		final += light_color * mix( diffuse, specular, F );
	}

    //final += texture( u_tex0, ref ).rgb * fresnel( Ks, n, -dir );
    final += getTextureColor(input, ref).rgb * fresnel( Ks, n, -dir );
	return final;
}

// unused
bool ray_vs_aabb(float3 o, float3 dir,
                 float3 bmin, float3 bmax,
                 device float2 &e )
{
    float3 a = ( bmin - o ) / dir;
    float3 b = ( bmax - o ) / dir;

    float3 s = min( a, b );
    float3 t = max( a, b );

    e.x = max( max( s.x, s.y ), max( s.z, e.x ) );
    e.y = max( min( t.x, t.y ), max( t.z, e.y ) );

    return e.x < e.y;
}

// ray marching
// o=ray origin; dir=(normalized) ray direction
// Returned a normalized 3D "normal" vector and depth.
bool ray_marching(float3 o, float3 dir,
                  thread float &depth, thread float3 &n)
{
	float t = 0.0;          // parameter along ray
    float d = 10000.0;      // distance
    float dt = 0.0;         // delta
    for ( int i = 0; i < 128; i++ )
    {
        float3 v = o + dir * t;
        d = dist_field( v );
        if ( d < 0.001 )
        {
            break;
        }
        dt = min( abs(d), 0.1 );
        t += dt;
        if ( t > depth )
        {
            break;
        }
    }

    if ( d >= 0.001 )
    {
        return false;
    }

    t -= dt;
    for ( int i = 0; i < 4; i++ )
    {
        dt *= 0.5;

        float3 v = o + dir * ( t + dt );
        if ( dist_field( v ) >= 0.001 )
        {
            t += dt;
        }
    }

    depth = t;
    // Normalized the normal.
    n = normalize( gradient( o + dir * t ) );
    return true;
}

// Returns normalized ray direction
// fov=vertical field-of-view, size=size of view
float3 ray_dir( float fov, float2 size, float2 pos )
{
	float2 xy = pos - size * 0.5;

    // cot θ = tan (90 - θ) where θ is expressed in radians.
	float cot_half_fov = tan( ( 90.0 - fov * 0.5 ) * DEG_TO_RAD );	
	float z = size.y * 0.5 * cot_half_fov;

    // The z-component of the ray is in the -ve direction of z-axis.
	return normalize( float3( xy, -z ) );
}

// camera rotation : pitch, yaw
float3x3 rotationXY( float2 angle )
{
	float2 c = cos( angle );
	float2 s = sin( angle );

	return float3x3(
		c.y      ,  0.0, -s.y,
		s.y * s.x,  c.x,  c.y * s.x,
		s.y * c.x, -s.x,  c.y * c.x
	);
}

// The parameter "input" is instantiated from a 2:1 equirectangular graphic.
kernel void
compute(texture2d<float, access::write> output      [[texture(0)]],
        texture2d<float>                input       [[texture(1)]],
        constant float                  &u_time     [[buffer(0)]],
        constant float2                 &u_mouse    [[buffer(1)]],
        uint2                           gid         [[thread_position_in_grid]])
{
    uint width = output.get_width();
    uint height = output.get_height();
    uint col = gid.x;
    uint row = gid.y;
    if ((col >= width) || (row >= height))
    {
        // In case the size of the texture does not match the size of the grid.
        // Return early if the pixel is out of bounds
        return;
    }
    // We don't have to pass the resolution as a parameter.
    float2 u_resolution = float2(width, height);
    float2 fragCoord = float2(gid.x, height-gid.y-1);

    // Get default (normalized) ray direction
	float3 dir = ray_dir( 45.0, u_resolution.xy, fragCoord.xy );

    // default ray origin
	float3 eye = float3( 0.0, 0.0, 3.5 );

	// rotate camera
	float3x3 rot = rotationXY( ( u_mouse.xy - u_resolution.xy * 0.5 ).yx * float2( 0.01, -0.01 ) );
	dir = rot * dir;
	eye = rot * eye;

    float4 fragColor = float4(0);
	// The function ray_marching will return "depth" and "n"
    float depth = clip_far;
    float3 n = float3( 0.0 );
	if ( !ray_marching( eye, dir, depth, n ) )
    {
        // On returning from ray marching, n is a normalized vector.
        //fragColor = texture( u_tex0, dir );
        fragColor = float4(getTextureColor(input, dir), 1.0);
	}
    else {
        // The ray hit some object proceed to do shading
        float3 pos = eye + dir * depth;

        // Note: both n and dir are normalized vectors.
        float3 color = shading( pos, n, dir, eye, input );
        fragColor = float4( pow( color, float3(1.0/1.2) ), 1.0 );
    }
    // Output to texture
    output.write(fragColor, gid);
}
