//!HLSL

#ifndef SHOW_REFLECTIONS
#define SHOW_REFLECTIONS 1
#endif

#include "foghelper.hlsl"	
	
struct PerPixelReflect_VO
{
	float4 _position	: POSITION;
	float4 _worldPos	: TEXCOORD0;
#if LOWEND
	float3 _bumpUV1		: TEXCOORD1;	// Use z coordinate for distanceToPlaneAlongNormal, since we can't use .w in LOWEND shaders
	float2 _bumpUV2		: TEXCOORD2;
	float2 _eyeDir2		: TEXCOORD3;
	float2 _halfVector	: TEXCOORD4;	
	float4 _fogValue    : COLOR0;
	float4 _viewAngle   : COLOR1;
#else
	float4 _bumpUV		: TEXCOORD1;
	float4 _eyeVector	: TEXCOORD2;
	float4 _sunDirection : TEXCOORD3;	// w holds viewAngle
	float4 _fogValue    : TEXCOORD4;
#endif
#if SHOW_REFLECTIONS
	float3 _refUVs		: TEXCOORD5;
#else
	float3 _refColor	: TEXCOORD5;
#endif
};

#define PI 3.141592654
void getWaveVariables( out float4 magnitudes, out float4 coss, out float4 sins,
					   in float4 scale, in float4 unscaledMagnitudes, in float4 frequencies, float xPoint, float zPoint, float time,
					   float4 waveDirX, float4 waveDirZ, float4 waveSpeed, float4 offsets )
{

///	height(x, z, t) = sum { M sin[ K ( x cos D + z sin D ) + W t + P]
//	M = amplitude
//	K = 2 pi / wavelength
//	D = direction
//	W = 2 pi * frequency = 2 pi * velocity / wavelength
//	P = phase angle
	
	// wave equation stolen from paper
	float4 kxz = frequencies *(((xPoint * waveDirX) + (zPoint * waveDirZ)) + waveSpeed*time + offsets);
	
	magnitudes = unscaledMagnitudes * scale;

	
	float4 clamped = ( (2*PI)*frac((kxz/(2*PI))+.5) ) - PI;
	sins = sin(clamped);
	coss = cos(clamped);
}


PerPixelReflect_VO ViewspaceOceanVS(
		float4 screenPos : POSITION,
		uniform float4x4 modelViewProjection,
		uniform float4x4 invViewProjection,
		uniform float3 sunDirectionWorld, 
		uniform float3 sunColor,
		uniform float4 waveDirX,
		uniform float4 waveDirZ,
		uniform float3 rampConst,
		uniform float3 eyeWorld,
		uniform float time,
		uniform float4	refPlane,
		uniform float4 bumpScale,
		uniform float4 fogColor,
 		uniform float2 fogParams,	
 		uniform float4		alphaRamp,
		uniform float4x4	worldRef
#ifdef USE_WAVES
		,
		uniform float4  waveFrequency,
		uniform float4  waveMagnitude,
		uniform float4 waveSpeed,
		uniform float4 waveOffset
#endif
#if LOWEND
		, uniform float3	fresnelConsts
#endif
#if !SHOW_REFLECTIONS		
		, uniform float3		horizonColor,
		uniform float3		skyColor
#endif		
		)
		
{

	PerPixelReflect_VO output;
	
	// projection direction screenpos into world space	
	float4 position = screenPos;
	
	float4 worldPos = position;
	worldPos = mul( invViewProjection, worldPos );
	worldPos = worldPos/worldPos.w;
	

	// projection direction			
	float3 worldDir = worldPos.xyz - eyeWorld.xyz;
	worldDir = normalize( worldDir );
	
	if ( (worldDir.y) > -.000001 )
	{
		worldDir.y = -.000001;
	}
			
	// now put back onto the screen with the appropriate depth
	
	// find plane intersection point
	float distanceToPlane = -worldPos.y;
	
	
	float cosAngle = worldDir.y;
	
	distanceToPlane = distanceToPlane/cosAngle;
	


	// now find the actual point
	float4 oceanPoint;
	oceanPoint.xyz = worldPos.xyz + distanceToPlane * worldDir.xyz;
	oceanPoint.w = 1;
	
		
	float3 eyevertdist = oceanPoint.xyz - eyeWorld;
	float3 eyeVec = normalize( eyevertdist );

	// distance from eye to vertex in xz plane
	eyevertdist.y = 0;
	// distance gets clamped to between x (1 value) and y (0 value)
	float distance = length(eyevertdist/eyeVec.y); 
	distance = clamp( distance, rampConst.x, rampConst.y);

	// we are going to scale y so it damps off with distance.
	// this gives us 1 when distance = rampConst.x, and 0 when it = rampConst.y. 
	float	scale;
	scale = (rampConst.y - distance)/rampConst.z;
	
	float3 normal = float3(0,1,0);
	
#ifdef USE_WAVES
	float3 tangent = float3(1,0,0);
	float3 binormal = float3(0,0,1);
	
	float4 magnitudes;
	float4 sins;
	float4 coss;
	
	// add waves
	getWaveVariables( magnitudes,  coss, sins,
					   scale, waveMagnitude, waveFrequency, oceanPoint.x, oceanPoint.z, time,
					   waveDirX, waveDirZ, waveSpeed, waveOffset );
					   	
	
	// take partial derivates
	
	// now the tangent ...
	tangent.y = dot( coss, (magnitudes * waveFrequency *  waveDirX) );
	tangent = normalize( tangent );
	
	// and the binormal
	binormal.y = dot( coss, (magnitudes * waveFrequency * waveDirZ) );
	binormal = normalize( binormal );
	
	// cross product for the normal
	normal = cross( binormal, tangent );
	
	// sum up to get the height
	oceanPoint.y = dot( sins, magnitudes);

   	float3x3 rotation = float3x3(tangent, normal, binormal );
	sunDirectionWorld = mul( rotation, sunDirectionWorld );
#endif

		
#if !LOWEND
	output._sunDirection.xyz = sunDirectionWorld;
#endif	
	
	worldPos = oceanPoint;
	float3 eyeVector = worldPos.xyz - eyeWorld;
	eyeVector = normalize( eyeVector );
	
#ifdef USE_WAVES
	eyeVector = mul( (float3x3)rotation, (float3)eyeVector );	
#endif
	
	// piggy back the alpha
	float viewAngle0 = -(dot(eyeVector, normal.xyz));
	float viewAngle = 1-viewAngle0;
	// moved this to vertex shader, we're out of instruction
	viewAngle = clamp(viewAngle, alphaRamp.x, alphaRamp.y);
	viewAngle = viewAngle/alphaRamp.z; 

#if LOWEND
	float sunDiffuse = saturate(dot( sunDirectionWorld.xyz, normal.xyz ));
	float fresnelViewAngle = saturate( pow(viewAngle0, fresnelConsts.z) );
	output._viewAngle.x = sunDiffuse;
	output._viewAngle.y = normal.y;
	output._viewAngle.z = fresnelViewAngle;
	output._viewAngle.w = viewAngle;
	output._halfVector.xy = normalize( -eyeVector.xyz + sunDirectionWorld.xyz ).xz;	// swizzle here instead of in the pixel shader
	output._eyeDir2.xy = normalize( -eyeVector.xyz - dot(normal, eyeVector) * normal ).xz; 	// Project eye vector to plane and swizzle
#else
	output._sunDirection.w = viewAngle;
	output._eyeVector.xyz = eyeVector;
#endif	
		
	output._worldPos.xyz = worldPos.xyz;
	
	// distance to plane
	float d = dot(oceanPoint.xyz, refPlane.xyz) + refPlane.w;
	
 	// Calculate a reflection vector
	float3 reflectNormal = float3(0,1,0);
	// hack to make the reflections less wonky
  	reflectNormal = reflect(eyeVector.xyz, reflectNormal);
	
	// distance to plane along normal
	float cosPlane = dot(-reflectNormal, refPlane.xyz);
  	float t = d/cosPlane;
#if LOWEND
	// Use z coordinate of bump uvs for distanceToPlaneAlongNormal,
	// since we can't use .w for LOWEND shaders
	output._bumpUV1.z = t;
#else
	output._eyeVector.w = t;
#endif	

#if SHOW_REFLECTIONS
	// calc the ref uvs	
	float4 refPlanePt;
	refPlanePt.xyz = worldPos.xyz - t*reflectNormal.xyz;
	refPlanePt.w = 1;

	float4 refUVs;
	refUVs.x = dot(worldRef[0], refPlanePt );
	refUVs.y = -dot(worldRef[1], refPlanePt );
	refUVs.z = alphaRamp.w;	// piggy back perturbscale on alpha ramp
	refUVs.w = dot(worldRef[3], refPlanePt );
	refUVs.xy /= refUVs.w;
	
	// TO DO: move this into the matrix itself
	refUVs.xy = saturate((refUVs.xy + float2(1,1)) / 2);
	
	output._refUVs.xyz = refUVs.xyz;
#else
	// calc the reflection color
	float fSkyFraction = reflectNormal.y;
	output._refColor = lerp(horizonColor, skyColor, fSkyFraction);
#endif	
	output._worldPos.w = scale;
	
	// put back in bump offset over time
	//  piggy back the speed onto the scale scales are x and z, speed is y and w
	float4 bumpUVs = ((time * waveDirX.x * bumpScale.yyww) + oceanPoint.xzxz * bumpScale.xxzz - floor(eyeWorld.xzxz * bumpScale.xxzz));
#if LOWEND	
	output._bumpUV1.xy = bumpUVs.xy;
	output._bumpUV2.xy = bumpUVs.zw;
#else	
	output._bumpUV.xyzw = bumpUVs;
#endif

			
	worldPos = mul(modelViewProjection, oceanPoint );
	
	output._fogValue = getFogValue(fogColor, fogParams,
								worldPos);


	worldPos = worldPos/worldPos.w;

	output._position = worldPos;

	return output;
	
} 




struct OceanPixelOutput
{
	float4 _color : COLOR;
};

#if LOWEND

// Simplified ocean shader 

float approxDotInverse(float Nx, float Ny, float Hx, float Hy)
{
/*	
	float2 delta = float2(Hx - Nx, Hy - Ny);	// 0 when on center of highlight
	float e = dot(delta,delta);
*/
	float4 N = float4(Nx, Ny, 0, 0);
	float4 H = float4(Hx, Hy, 0, 0);
	float4 D = H - N;	// 0 when on center of highlight
	float e = dot(D,D);
	return e;
}

float approxSpecular(float p,	// specular exponent
					 float Nx, float Ny,	// Nz = sqrt(1 - Nx*Nx - Ny*Ny)
					 float Hx, float Hy)	// Hz = sqrt(1 - Hx*Hx - Hy*Hy)
{
	// Approximate
	//      f = pow( dot(N3,H3), p)
	// with
	//      g = pow( 1 - length2(N2-H2), p/2 )
	//
	// Let a = Nx, b = Ny, c = Hx, d = Hy
	//
	//      f = pow( a*c + b*d + sqrt(1 - (a*a + b*b))*sqrt(1 - (c*c + d*d)), p)
	//
	//      g = pow( 1 - ((a-c)*(a-c) + (b-d)*(b-d)), p/2)
	//        = pow( 1 - (a*a - 2*a*c + c*c + b*b - 2*b*d + d*d), p/2)
	//        = pow( 2*(a*c + b*d) - 1 + (1 - (a*a + b*b)) + (1 - (c*c+ d*d)), p/2 )
	//
	// When length2(a,b) = 0, f = g
	//      f = pow( sqrt(1 - (c*c + d*d)), p) = pow( 1 - (c*c + d*d), p/2)
	//      g = pow( 1 - (c*c + d*d), p / 2)
	
	// When length2(a,b) = 1, qualitatively right when saturation is
	// added so that values are clipped below 0.
	// Highlight is too big? p/2 instead of p
	//      f = pow( a*c + b*d, p )
	//      g = pow( 2*(a*c + b*d) - (c*c + d*d), p/2)
	float e = approxDotInverse(Nx, Ny, Hx, Hy);

	// Power series expansion of pow( 1 - e, q ) around x0 = 1
	//
	// h(x0 - e) = h(x0) - h'(x0) * e / 2!  + h''(x0) * e * e / 3! - h'''(x0) * e * e * e / 4! + ...
	//
	// h(x) = pow( x, q)
	// h'(x) = q * pow( x, q - 1)
	// h''(x) = q * (q - 1) * pow( x, q - 2)
	// h'''(x) = q * (q - 1) * (q - 2) * pow( x, q - 3)
	//
	// For x0 = 1
	// h(1) = 1
	// h'(1) = q
	// h''(1) = q * (q - 1)
	// h'''(1) = q * (q - 1) * (q - 2)
	//
	// h(1-e) = 1 - (q/2)*e + (q*(q-1)/6)*e*e + (q*(q-1)*(q-2)/24)*e*e*e + ...
	// h(1-e) = 1 - e * (q/2)(1 - e * (q-1)/3 * ( 1 - e*(q-2)/4 )) + ...
	
	// q = p/2

	// float result = 1 - e * (p/4) * (1 - e * (p/2-1)/3 * (1 - e * (p/2-2)/4));	// third order
	// float result = 1.f - (p / 4.f ) * e;	// first order
#define FUDGE 8		// fudge factor for ps_1_4 arithmetic troubles
	//float result = 1-FUDGE*p*e;
	float result = 1-FUDGE*e;
	return saturate(result);
}

OceanPixelOutput ViewspaceOcean_LowEnd_PS(
	PerPixelReflect_VO In,
#if SHOW_REFLECTIONS
	uniform sampler2D	reflectionMap,
#endif
	uniform sampler2D	bumpMap,
	uniform sampler2D	bumpMap2,
	uniform float4		refPlane,
	uniform float3		eyeWorld,
	uniform float3		fresnelConsts,
	uniform float3		fresnelColor,
	uniform float3		sunSpecular,
	uniform float4		animationConst,
	uniform float4		oceanAmbient
	)
{
	float2 halfVector	= In._halfVector.xy;
	float viewAngleAlpha = In._viewAngle.w;

	OceanPixelOutput retVal;

	// yz-swizzle already done in halfVector
	float2 normal = tex2D( bumpMap, In._bumpUV1 );
	normal = (normal*2) -  0.99608;
	normal.xy *= animationConst.x;
	float2 normal2 = tex2D( bumpMap2, In._bumpUV2 );
	normal2 = (normal2*2)- 0.99608;
	normal2.xy *= animationConst.y;
	normal.xy += normal2.xy;
	
#if SHOW_REFLECTIONS	
	float3 refUVs = In._refUVs;
	float3 refColor = tex2D(reflectionMap, refUVs.xy);
#else
	float3 refColor = In._refColor;
#endif

	float refl = approxDotInverse( normal.x, normal.y, In._eyeDir2.x, In._eyeDir2.y);
	float3 refColor2 = refColor * refl;

	float fSpecular = approxSpecular( fresnelConsts.y, normal.x, normal.y, halfVector.x, halfVector.y);
	float3 specular = fSpecular * sunSpecular;
	
	retVal._color.xyz = (specular + refColor2);
	retVal._color.w = viewAngleAlpha;

	// TODO: fog
	
	return retVal;
}


#else // !LOWEND

OceanPixelOutput ViewspaceOcean_PerPixelReflect_PS(
	PerPixelReflect_VO In,
#if SHOW_REFLECTIONS
	uniform sampler2D	reflectionMap,
#endif
	uniform sampler2D	bumpMap,
	uniform sampler2D	bumpMap2,
	uniform float4		refPlane,
	uniform float3		eyeWorld,
	uniform float3		fresnelConsts,
	uniform float3		fresnelColor,
	uniform float3		sunSpecular,
	uniform float4		animationConst,
	uniform float4		oceanAmbient
	)
{
	float4 worldPos		= In._worldPos;
	float4 eyeVector	= In._eyeVector;
	float3 sunDirection = In._sunDirection;
	float4 fogValue     = In._fogValue;
#if SHOW_REFLECTIONS	
	float3 refUVs		= In._refUVs;
	float  distanceToPlaneAlongNormal = In._eyeVector.w;
#endif
	float4 bumpUV		= In._bumpUV;
	float viewAngleAlpha = In._sunDirection.w;
	OceanPixelOutput retVal;
	
	// swap y's and z's since we are using y as up 
	float3 normal = tex2D( bumpMap, bumpUV.xy ).xzy;

	normal = (normal*2) -  0.99608;
	normal.xz *= animationConst.x;
	
	// can't normalize because we're out of instructions
	float3 normal2 = tex2D( bumpMap2, bumpUV.zw ).xzy;

	normal2 = (normal2*2)- 0.99608;
	normal2.xz *= animationConst.y;
	
	// this is not quite accurate, but we need to scrimp on instructions
	normal.xz += normal2.xz;
	normal.y = sqrt(1 - (normal.x*normal.x + normal.z*normal.z));
	
	normal = normalize(normal); //Enhanced by Celeste
	
#if SHOW_REFLECTIONS
	// find the perturbed point
	float perturbMagnitude = dot(normal, float3(0,1,0));
	perturbMagnitude = sqrt(1 - (perturbMagnitude*perturbMagnitude));
	float2 perturbAmountScreen = distanceToPlaneAlongNormal * perturbMagnitude;

	// now put into reflection
	refUVs.xy += (perturbAmountScreen * normal.xz)*refUVs.z;
#endif	

	float3 halfAngle = normalize(-eyeVector.xyz+sunDirection.xyz);
	float specular = dot( normal.xyz, halfAngle);
	specular = (fresnelConsts.x * pow( specular, fresnelConsts.y));
	
	float light = saturate(dot( sunDirection.xyz, normal.xyz));

#if SHOW_REFLECTIONS	
	float4 refColor = tex2D(reflectionMap, refUVs.xy);
#else
	float3 refColor = In._refColor;
#endif	
	
	float viewAngle = -(dot(eyeVector, normal.xyz));
	viewAngle = pow( viewAngle, fresnelConsts.z );
	viewAngle = saturate(viewAngle);
	
	fresnelColor = fresnelColor * light * 0.8;	//Enhanced by Celeste - Darken the overall color
	fresnelColor.y = fresnelColor.x;	//Enhanced by Celeste
	fresnelColor.z = fresnelColor.x;	//Enhanced by Celeste
	
	refColor *= 0.8;
	
	fresnelColor.y *= 0.5; //Enhanced by Celeste
	fresnelColor.x = fresnelColor.y; //*= 0.5; //Enhanced by Celeste
	fresnelColor.z *= 0.85; //Enhanced by Celeste
	
	float3 flatColor = lerp(refColor.xyz, fresnelColor, viewAngle);

	flatColor.x *= 0.6;
	flatColor.y *= 0.8;
	
	oceanAmbient *= 0.5; //Enhanced by Celeste
		
	// removing saturate, not enough instructions
	retVal._color.xyz = flatColor.xyz + saturate(specular * sunSpecular) + oceanAmbient * (1 - light);
	retVal._color.x = clamp(retVal._color.x, 0.0, retVal._color.y * 0.9);	 //Enhanced by Celeste

	// fog	
	retVal._color.xyz = lerp(retVal._color.xyz, fogValue.xyz, fogValue.w);
	retVal._color.w = saturate(viewAngleAlpha + 0.2 + fogValue.w); //Enhanced by Celeste

	return retVal;
}

#endif
