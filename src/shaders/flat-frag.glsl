#version 300 es
precision highp float;

uniform vec3 u_Eye, u_Ref, u_Up;
uniform vec2 u_Dimensions;
uniform float u_Time;
uniform float u_Thrust;
uniform float u_Fins;

in vec2 fs_Pos;
out vec4 out_Col;

#define EPSILON 0.000001
#define NORMAL_EPSILON 0.01
#define MAX_DIST 100.0
#define MAX_MARCHING_STEPS 500
#define M_PI 3.1415926535897932384626433832795

#define PLANET_RADIUS 1000.0

mat4 rotationMatrix(vec3 axis, float angle)
{
    axis = normalize(axis);
    float s = -sin(angle);
    float c = cos(angle);
    float oc = 1.0 - c;

    return mat4(oc * axis.x * axis.x + c,           oc * axis.x * axis.y - axis.z * s,  oc * axis.z * axis.x + axis.y * s,  0.0,
                oc * axis.x * axis.y + axis.z * s,  oc * axis.y * axis.y + c,           oc * axis.y * axis.z - axis.x * s,  0.0,
                oc * axis.z * axis.x - axis.y * s,  oc * axis.y * axis.z + axis.x * s,  oc * axis.z * axis.z + c,           0.0,
                0.0,                                0.0,                                0.0,                                1.0);
}

float sdSphere( vec3 p, float s )
{
    return length(p)-s;
}

float sdBox( vec3 p, vec3 b )
{
  vec3 d = abs(p) - b;
  return length(max(d,0.0))
         + min(max(d.x,max(d.y,d.z)),0.0);
}

float sdCylinder(vec3 p, vec3 a, vec3 b, float r)
{
    vec3 pa = p - a;
    vec3 ba = b - a;
    float baba = dot(ba,ba);
    float paba = dot(pa,ba);

    float x = length(pa*baba-ba*paba) - r*baba;
    float y = abs(paba-baba*0.5)-baba*0.5;
    float x2 = x*x;
    float y2 = y*y*baba;
    float d = (max(x,y)<0.0)?-min(x2,y2):(((x>0.0)?x2:0.0)+((y>0.0)?y2:0.0));
    return sign(d)*sqrt(abs(d))/baba;
}

float sdCapsule( vec3 p, vec3 a, vec3 b, float r )
{
	vec3 pa = p-a, ba = b-a;
	float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
	return length( pa - ba*h ) - r;
}

float dot2( in vec2 v ) {
  return dot(v,v);
}

float dot2(in vec3 v) {
  return dot(v,v);
}

float sdCappedCone( in vec3 p, in float h, in float r1, in float r2 )
{
  vec2 q = vec2( length(p.xz), p.y );

  vec2 k1 = vec2(r2,h);
  vec2 k2 = vec2(r2-r1,2.0*h);
  vec2 ca = vec2(q.x-min(q.x,(q.y < 0.0)?r1:r2), abs(q.y)-h);
  vec2 cb = q - k1 + k2*clamp( dot(k1-q,k2)/dot2(k2), 0.0, 1.0 );
  float s = (cb.x < 0.0 && ca.y < 0.0) ? -1.0 : 1.0;
  return s*sqrt( min(dot2(ca),dot2(cb)) );
}

float sdRoundCone(vec3 p, vec3 a, vec3 b, float r1, float r2)
{
    vec3  ba = b - a;
    float l2 = dot(ba,ba);
    float rr = r1 - r2;
    float a2 = l2 - rr*rr;
    float il2 = 1.0/l2;

    vec3 pa = p - a;
    float y = dot(pa,ba);
    float z = y - l2;
    float x2 = dot2( pa*l2 - ba*y );
    float y2 = y*y*l2;
    float z2 = z*z*l2;

    float k = sign(rr)*rr*rr*x2;
    if( sign(z)*a2*z2 > k ) return  sqrt(x2 + z2)        *il2 - r2;
    if( sign(y)*a2*y2 < k ) return  sqrt(x2 + y2)        *il2 - r1;
                            return (sqrt(x2*a2*il2)+y*rr)*il2 - r1;
}

float sdEllipsoid( in vec3 p, in vec3 r )
{
    float k0 = length(p/r);
    float k1 = length(p/(r*r));
    return k0*(k0-1.0)/k1;
}

float udTriangle(vec3 p, vec3 a, vec3 b, vec3 c)
{
    vec3 ba = b - a; vec3 pa = p - a;
    vec3 cb = c - b; vec3 pb = p - b;
    vec3 ac = a - c; vec3 pc = p - c;
    vec3 nor = cross( ba, ac );

    return sqrt(
    (sign(dot(cross(ba,nor),pa)) +
     sign(dot(cross(cb,nor),pb)) +
     sign(dot(cross(ac,nor),pc))<2.0)
     ?
     min( min(
     dot2(ba*clamp(dot(ba,pa)/dot2(ba),0.0,1.0)-pa),
     dot2(cb*clamp(dot(cb,pb)/dot2(cb),0.0,1.0)-pb) ),
     dot2(ac*clamp(dot(ac,pc)/dot2(ac),0.0,1.0)-pc) )
     :
     dot(nor,pa)*dot(nor,pa)/dot2(nor) );
}

float udQuad(vec3 p, vec3 a, vec3 b, vec3 c, vec3 d)
{
    vec3 ba = b - a; vec3 pa = p - a;
    vec3 cb = c - b; vec3 pb = p - b;
    vec3 dc = d - c; vec3 pc = p - c;
    vec3 ad = a - d; vec3 pd = p - d;
    vec3 nor = cross( ba, ad );

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

float sdTriangularPrism(vec3 p, vec3 a1, vec3 b1, vec3 c1, vec3 a2, vec3 b2, vec3 c2) {
  float f1 = udTriangle(p, a1, b1, c1);
  float f2 = udTriangle(p, a2, b2, c2);
  float s1 = udQuad(p, a1, b1, b2, a2);
  float s2 = udQuad(p, b1, c1, c2, b2);
  float s3 = udQuad(p, c1, a1, a2, c2);
  if (abs(f1 + f2 - distance(a1, a2)) < EPSILON) {
    return -min(f1, min(f2, min(s1, min(s2, s3))));
  } else {
    return min(f1, min(f2, min(s1, min(s2, s3))));
  }
}

float random1(vec3 p, vec3 seed) {
  return fract(sin(dot(p + seed, vec3(987.654, 123.456, 531.975))) * 85734.3545);
}

vec2 randvec3(vec2 n, vec2 seed) {
  float x = sin(dot(n + seed, vec2(14.92, 64.42)));
  float y = sin(dot(n + seed, vec2(48.12, 32.42)));
  return fract(334.963f * vec2(x, y));
}

float triangle_wave(float x, float freq, float amplitude) {
  return abs(x * freq - amplitude * floor(x * freq / amplitude) - (0.5 * amplitude));
}

float quinticSmooth(float t) {
  float x = clamp(t, 0.0, 1.0);
  return x * x * x * (x * (x * 6.0  - 15.0) + 10.0);
}

float worleyNoise(vec2 pos) {
  float factor = 8.0;
  vec2 seed = vec2(0.0, 0.0);

  int x = int(floor(pos.x / factor));
  int y = int(floor(pos.y / factor));
  vec2 minWorley = factor * randvec3(vec2(float(x), float(y)), seed) + vec2(float(x) * factor, float(y) * factor);
  float minDist = distance(minWorley, pos);
  for (int i = x - 1; i <= x + 1; i++) {
      for (int j = y - 1; j <= y + 1; j++) {
          vec2 worley = factor * randvec3(vec2(float(i), float(j)), seed) + vec2(float(i) * factor, float(j) * factor);
          if (minDist > distance(pos, worley)) {
              minDist = distance(pos, worley);
              minWorley = worley;
          }
      }
  }
  return clamp(minDist / (factor * 2.0), 0.0, 0.5);
}

float interpRand(float x, float y, float z) {
  vec3 seed = vec3(0.0, 0.0, 0.0);

  float intX = floor(x);
  float fractX = fract(x);
  float intY = floor(y);
  float fractY = fract(y);
  float intZ = floor(z);
  float fractZ = fract(z);

  vec3 c1 = vec3(intX, intY, intZ);
  vec3 c2 = vec3(intX + 1.0, intY, intZ);
  vec3 c3 = vec3(intX, intY, intZ + 1.0);
  vec3 c4 = vec3(intX + 1.0, intY, intZ + 1.0);
  vec3 c5 = vec3(intX, intY + 1.0, intZ);
  vec3 c6 = vec3(intX + 1.0, intY + 1.0, intZ);
  vec3 c7 = vec3(intX, intY + 1.0, intZ + 1.0);
  vec3 c8 = vec3(intX + 1.0, intY + 1.0, intZ + 1.0);

  float v1 = random1(c1, seed);
  float v2 = random1(c2, seed);
  float v3 = random1(c3, seed);
  float v4 = random1(c4, seed);
  float v5 = random1(c5, seed);
  float v6 = random1(c6, seed);
  float v7 = random1(c7, seed);
  float v8 = random1(c8, seed);

  float i1 = mix(v1, v2, quinticSmooth(fractX));
  float i2 = mix(v3, v4, quinticSmooth(fractX));
  float i3 = mix(v5, v6, quinticSmooth(fractX));
  float i4 = mix(v7, v8, quinticSmooth(fractX));
  float j1 = mix(i1, i2, quinticSmooth(fractZ));
  float j2 = mix(i3, i4, quinticSmooth(fractZ));
  return mix(j1, j2, quinticSmooth(fractY));
}

float rocketTexture(float x, float y, float z) {
  float total = 0.0;
  int octaves = 8;
  float persistence = 0.5;
  for (int i = 0; i < octaves; i++) {
    float freq = pow(2.0, float(i));
    float amp = pow(persistence, float(i));
    total += interpRand(x * freq, y * freq, z * freq) * amp;
  }
  return total;
}

float planetTexture(float x, float y, float z) {
  float theta = PLANET_RADIUS * acos(x / PLANET_RADIUS) / 10.0;
  float phi = PLANET_RADIUS * atan(y, z) / 10.0;
  float total = 0.0;
  int octaves = 10;
  float persistence = 0.7;
  for (int i = 0; i < octaves; i++) {
    float freq = pow(2.0, float(i));
    float amp = pow(persistence, float(i));
    total += worleyNoise(vec2(theta * freq, phi * freq)) * amp;
  }
  return total;
}

float smin(float a, float b, float k) {
  float h = clamp(0.5 + 0.5 * (b - a) / k, 0.0, 1.0);
  return mix(b, a, h) - k * h * (1.0 - h);
}

float sdPlanet(vec3 pos) {
  return sdSphere(pos, PLANET_RADIUS);
}

float sdAsteroids(vec3 pos) {
  vec3 asteroidPosition = vec3(15.0, (PLANET_RADIUS + 5.0), 0.0);
  vec4 i = vec4((pos - asteroidPosition), 1.0);
  vec3 inversePos = i.xyz;

  float box = sdBox(inversePos, vec3(5.0, 5.0, 5.0));
  float displacement = 2.0 * sin(1.0 * inversePos.x) + cos(1.0 * inversePos.y) + sin(1.0 * inversePos.z);
  box += displacement;
  float sphere = sdEllipsoid(inversePos, vec3(7.0, 6.0, 5.0));
  float result = max(box, sphere);

  for (int j = 1; j < 20; j++) {
    float sineVal = sin(radians(float(j) * 18.0));
    float cosineVal = cos(radians(float(j) * 18.0));
    vec3 asteroidPosition = vec3(15.0, (PLANET_RADIUS + 5.0) * cosineVal, (PLANET_RADIUS + 5.0) * sineVal);
    if (j % 2 == 1) {
      asteroidPosition.x *= -1.0;
    }
    i = vec4((pos - asteroidPosition), 1.0);
    inversePos = i.xyz;

    box = sdBox(inversePos, vec3(5.0, 5.0, 5.0));
    displacement = 2.0 * sin(1.0 * inversePos.x) + cos(1.0 * inversePos.y) + sin(1.0 * inversePos.z);
    box += displacement;
    sphere = sdEllipsoid(inversePos, vec3(7.0, 6.0, 5.0));
    result = min(result, max(box, sphere));
  }
  return result;
}

float sdRocket(vec3 pos) {
  float sineVal = sin(u_Time * M_PI * 0.001);
  float cosineVal = cos(u_Time * M_PI * 0.001);
  vec3 rocketPosition = vec3(0.0, (PLANET_RADIUS + 10.0) * cosineVal, (PLANET_RADIUS + 10.0) * sineVal);
  vec3 planetOffset = vec3(0.0, PLANET_RADIUS + 10.0, 0.0);
  vec4 planetOffset4 = vec4(planetOffset, 1.0);
  vec3 rocketAxis = vec3(0.0, 0.0, 1.0);
  vec4 i = rotationMatrix(vec3(0.0, 0.0, 1.0), -u_Time * M_PI * 0.1) * rotationMatrix(vec3(1.0, 0.0, 0.0), -u_Time * M_PI * 0.001) * vec4((pos - rocketPosition), 1.0);
  vec3 inversePos = i.xyz;

  float rocketLength = 7.0;
  float noseLength = 2.0;
  vec3 bodyBottom = vec3(0.0, 0.0, 0.0);
  vec3 bodyTop = vec3(0.0, 0.0, rocketLength);
  vec3 noseTop = vec3(0.0, 0.0, rocketLength + noseLength);
  float bodyRadius = 1.0;
  float noseRadius = 0.2;
  float finThickness = 0.1;
  float finHeight = 2.0;
  float nozzleTopRadius = 0.3;
  float nozzleBottomRadius = 0.6;
  float nozzleHeight = 1.0;
  float nozzleThickness = 0.1;
  vec3 nozzleBottom = vec3(0.0, -2.0 * nozzleThickness, 0.0);

  float body = sdCylinder(inversePos, bodyBottom, bodyTop, bodyRadius);
  float nose = sdRoundCone(inversePos, bodyTop, noseTop, bodyRadius, noseRadius);
  float bodyNose = smin(body, nose, 0.1);

  vec4 j = rotationMatrix(vec3(1.0, 0.0, 0.0), -M_PI / 2.0) * i;
  vec3 nozzleInversePos = j.xyz;

  float outerNozzle = sdCappedCone(nozzleInversePos, nozzleHeight, nozzleBottomRadius, nozzleTopRadius);
  float innerNozzle = sdCappedCone(nozzleInversePos - nozzleBottom, nozzleHeight - nozzleThickness, nozzleBottomRadius - nozzleThickness, nozzleTopRadius - nozzleThickness);
  float nozzle = max(-innerNozzle, outerNozzle);
  float bodyNoseNozzle = smin(bodyNose, nozzle, 0.4);

  vec3 fina1 = bodyBottom + vec3(bodyRadius, -finThickness * 0.5, finHeight * 0.1);
  vec3 finb1 = bodyBottom + vec3(bodyRadius, -finThickness * 0.5, finHeight);
  vec3 finc1 = bodyBottom + vec3(bodyRadius * 2.5, -finThickness * 0.5, finHeight * -0.25);
  vec3 fina2 = bodyBottom + vec3(bodyRadius, finThickness * 0.5, finHeight * 0.1);
  vec3 finb2 = bodyBottom + vec3(bodyRadius, finThickness * 0.5, finHeight);
  vec3 finc2 = bodyBottom + vec3(bodyRadius * 2.5, finThickness * 0.5, finHeight * -0.25);
  float fin = sdTriangularPrism(inversePos, fina1, finb1, finc1, fina2, finb2, finc2);
  float withFin = smin(bodyNoseNozzle, fin, 0.3);

  for (int k = 1; k < int(u_Fins); k++) {
    mat4 matrix = rotationMatrix(rocketAxis, radians(float(k) * 360.0 / u_Fins));
    vec3 ofina1 = (matrix * vec4(fina1, 1.0)).xyz;
    vec3 ofinb1 = (matrix * vec4(finb1, 1.0)).xyz;
    vec3 ofinc1 = (matrix * vec4(finc1, 1.0)).xyz;;
    vec3 ofina2 = (matrix * vec4(fina2, 1.0)).xyz;;
    vec3 ofinb2 = (matrix * vec4(finb2, 1.0)).xyz;;
    vec3 ofinc2 = (matrix * vec4(finc2, 1.0)).xyz;;
    float ofin = sdTriangularPrism(inversePos, ofina1, ofinb1, ofinc1, ofina2, ofinb2, ofinc2);
    withFin = smin(withFin, ofin, 0.3);
  }
  return withFin;
}

float sdExhaust(vec3 pos) {
  float sineVal = sin(u_Time * M_PI * 0.001);
  float cosineVal = cos(u_Time * M_PI * 0.001);
  vec3 rocketPosition = vec3(0.0, (PLANET_RADIUS + 10.0) * cosineVal, (PLANET_RADIUS + 10.0) * sineVal);
  vec4 i = rotationMatrix(vec3(1.0, 0.0, 0.0), -u_Time * M_PI * 0.001) * vec4((pos - rocketPosition), 1.0);
  vec3 inversePos = i.xyz;

  vec3 bigExhaustSize = vec3(0.6, 0.6, 1.0 + 0.04 * (u_Thrust - 1.0));
  vec3 bigExhaustPos = vec3(0.0, 0.0, 1.5 + 0.04 * (u_Thrust - 1.0));
  int exhaustCount = 6;
  float result = sdEllipsoid(inversePos + bigExhaustPos, bigExhaustSize);
  for (int i = 1; i < exhaustCount; i++) {
    float zLength = 0.8 + 0.04 * (u_Thrust - 1.0) + triangle_wave(u_Time, 1.0, 0.08 * u_Thrust);
    vec3 newSize = vec3(pow(0.8, float(i)) * bigExhaustSize.x, pow(0.8, float(i)) * bigExhaustSize.y, pow(zLength, float(i)) * bigExhaustSize.z);
    vec3 newPos = vec3(0.0, 0.0, zLength * (1.0 - pow(zLength, float(i))) / (1.0 - zLength)) + bigExhaustPos;
    result = smin(result, sdEllipsoid(inversePos + newPos, newSize), 0.5);
  }
  return result;
}

vec3 estimateNormal(vec3 p) {
  float x1 = min(sdRocket(vec3(p.x + NORMAL_EPSILON, p.y, p.z)), min(sdPlanet(vec3(p.x + NORMAL_EPSILON, p.y, p.z)), sdExhaust(vec3(p.x + NORMAL_EPSILON, p.y, p.z))));
  float x2 = min(sdRocket(vec3(p.x - NORMAL_EPSILON, p.y, p.z)), min(sdPlanet(vec3(p.x - NORMAL_EPSILON, p.y, p.z)), sdExhaust(vec3(p.x - NORMAL_EPSILON, p.y, p.z))));
  float y1 = min(sdRocket(vec3(p.x, p.y + NORMAL_EPSILON, p.z)), min(sdPlanet(vec3(p.x, p.y + NORMAL_EPSILON, p.z)), sdExhaust(vec3(p.x, p.y + NORMAL_EPSILON, p.z))));
  float y2 = min(sdRocket(vec3(p.x, p.y - NORMAL_EPSILON, p.z)), min(sdPlanet(vec3(p.x, p.y - NORMAL_EPSILON, p.z)), sdExhaust(vec3(p.x, p.y - NORMAL_EPSILON, p.z))));
  float z1 = min(sdRocket(vec3(p.x, p.y, p.z + NORMAL_EPSILON)), min(sdPlanet(vec3(p.x, p.y, p.z + NORMAL_EPSILON)), sdExhaust(vec3(p.x, p.y, p.z + NORMAL_EPSILON))));
  float z2 = min(sdRocket(vec3(p.x, p.y, p.z - NORMAL_EPSILON)), min(sdPlanet(vec3(p.x, p.y, p.z - NORMAL_EPSILON)), sdExhaust(vec3(p.x, p.y, p.z - NORMAL_EPSILON))));
  x1 = min(x1, sdAsteroids(vec3(p.x + NORMAL_EPSILON, p.y, p.z)));
  x2 = min(x2, sdAsteroids(vec3(p.x - NORMAL_EPSILON, p.y, p.z)));
  y1 = min(y1, sdAsteroids(vec3(p.x, p.y + NORMAL_EPSILON, p.z)));
  y2 = min(y2, sdAsteroids(vec3(p.x, p.y - NORMAL_EPSILON, p.z)));
  z1 = min(z1, sdAsteroids(vec3(p.x, p.y, p.z + NORMAL_EPSILON)));
  z2 = min(z2, sdAsteroids(vec3(p.x, p.y, p.z - NORMAL_EPSILON)));
  return normalize(vec3(x1 - x2, y1 - y2, z1 - z2));
}

bool intersectsBox(float minX, float maxX, float minY, float maxY, float minZ, float maxZ, vec3 eye, vec3 rayDirec) {
  float tnear = -MAX_DIST * 10.0;
  float tfar = MAX_DIST * 10.0;

  if (abs(rayDirec.x) < EPSILON) {
    if (eye.x < minX || eye.x > maxX) {
      return false;
    }
  }
  float t0 = (minX - eye.x) / rayDirec.x;
  float t1 = (maxX - eye.x) / rayDirec.x;
  if (t0 > t1) {
    float s = t0;
    t0 = t1;
    t1 = s;
  }
  if (t0 > tnear) {
    tnear = t0;
  }
  if (t1 < tfar) {
    tfar = t1;
  }

  if (abs(rayDirec.y) < EPSILON) {
    if (eye.y < minX || eye.y > maxX) {
      return false;
    }
  }
  t0 = (minX - eye.y) / rayDirec.y;
  t1 = (maxX - eye.y) / rayDirec.y;
  if (t0 > t1) {
    float s = t0;
    t0 = t1;
    t1 = s;
  }
  if (t0 > tnear) {
    tnear = t0;
  }
  if (t1 < tfar) {
    tfar = t1;
  }

  if (abs(rayDirec.z) < EPSILON) {
    if (eye.z < minX || eye.z > maxX) {
      return false;
    }
  }
  t0 = (minX - eye.z) / rayDirec.z;
  t1 = (maxX - eye.z) / rayDirec.z;
  if (t0 > t1) {
    float s = t0;
    t0 = t1;
    t1 = s;
  }
  if (t0 > tnear) {
    tnear = t0;
  }
  if (t1 < tfar) {
    tfar = t1;
  }
  if (tnear > tfar) {
    return false;
  }
   return true;
}

bool intersectsRocketOrExhaust(vec3 rayDirec) {
  float sineVal = sin(u_Time * M_PI * 0.001);
  float cosineVal = cos(u_Time * M_PI * 0.001);
  vec3 rocketPosition = vec3(0.0, (PLANET_RADIUS + 10.0) * cosineVal, (PLANET_RADIUS + 10.0) * sineVal);
  mat4 matrix = rotationMatrix(vec3(1.0, 0.0, 0.0), -u_Time * M_PI * 0.001);
  vec4 inversePos = matrix * vec4((u_Eye - rocketPosition), 1.0);
  vec4 inverseRay = matrix * vec4(rayDirec, 1.0);

  float minZ = -2.5;
  float maxZ = 2.5;
  float minX = -20.0;
  float maxX = 9.0;
  float minY = -2.5;
  float maxY = 2.5;

  return intersectsBox(minX, maxX, minY, maxY, minZ, maxZ, inversePos.xyz, inverseRay.xyz);
}

bool intersectsPlanetOrAsteroids(vec3 rayDirec) {
  float minZ = -PLANET_RADIUS - 7.0;
  float maxZ = PLANET_RADIUS + 7.0;
  float minX = -PLANET_RADIUS - 7.0;
  float maxX = PLANET_RADIUS + 7.0;
  float minY = -PLANET_RADIUS - 7.0;
  float maxY = PLANET_RADIUS + 7.0;
  return intersectsBox(minX, maxX, minY, maxY, minZ, maxZ, u_Eye, rayDirec);
}

bool intersectsObject(vec3 rayDirec) {
  float sineVal = sin(u_Time * M_PI * 0.001);
  float cosineVal = cos(u_Time * M_PI * 0.001);
  vec3 rocketPosition = vec3(0.0, (PLANET_RADIUS + 10.0) * cosineVal, (PLANET_RADIUS + 10.0) * sineVal);
  mat4 matrix = rotationMatrix(vec3(1.0, 0.0, 0.0), -u_Time * M_PI * 0.001);
  vec4 inversePos = matrix * vec4((u_Eye - rocketPosition), 1.0);
  vec4 inverseRay = matrix * vec4(rayDirec, 1.0);

  float minZ = -PLANET_RADIUS - 7.0;
  float maxZ = PLANET_RADIUS + 7.0;
  float minX = -PLANET_RADIUS - 7.0;
  float maxX = PLANET_RADIUS + 7.0;
  float minY = -PLANET_RADIUS - 7.0;
  float maxY = PLANET_RADIUS + 10.0 + 2.5;

  return intersectsBox(minX, maxX, minY, maxY, minZ, maxZ, inversePos.xyz, inverseRay.xyz);
}

vec4 rayMarch(vec3 rayDirec) {
  if (!intersectsObject(rayDirec)) {
    return vec4(0.0, 0.0, 0.0, -1.0);
  } else {
    bool intersectsRocket = intersectsRocketOrExhaust(rayDirec);
    bool intersectsPlanetOrAsteroids = intersectsPlanetOrAsteroids(rayDirec);
    float depth = 0.0;
    for (int i = 0; i < MAX_MARCHING_STEPS; i++) {
      vec3 current = u_Eye + depth * rayDirec;
      float rocketDist = sdRocket(current);
      float exhaustDist = sdExhaust(current);
      float asteroidDist = sdAsteroids(current);
      float planetDist = sdPlanet(current);
      float dist = MAX_DIST;
      if (intersectsRocket) {
        dist = min(dist, min(rocketDist, exhaustDist));
      }
      if (intersectsPlanetOrAsteroids) {
        dist = min(dist, min(planetDist, asteroidDist));
      }
      if (dist < EPSILON) {
        if (abs(rocketDist - dist) < EPSILON) {
          return vec4(current, 0.0);
        }
        if (abs(exhaustDist - dist) < EPSILON) {
          return vec4(current, 2.0);
        }
        if (abs(asteroidDist - dist) < EPSILON) {
          return vec4(current, 3.0);
        }
        if (abs(planetDist - dist) < EPSILON) {
          return vec4(current, 1.0);
        }
      }
      depth += dist;
      if (depth >= MAX_DIST) {
        return vec4(current, -1.0);
      }
    }
    return vec4(0.0, 0.0, 0.0, -2.0);
  }
}

float bias(float b, float t) {
  return pow(t, log(b) / log(0.5));
}

void main() {
  float fovy = 90.0;
  vec3 look = u_Ref - u_Eye;
  vec3 right = normalize(cross(look, u_Up));
  float aspect = float(u_Dimensions.x) / float(u_Dimensions.y);
  float tan_fovy2 = tan(fovy / 2.0);
  vec3 h = right * length(look) * aspect * tan_fovy2;
  vec3 v = u_Up * length(look) * tan_fovy2;
  vec3 p = u_Ref + fs_Pos.x * h + fs_Pos.y * v;

  vec3 sunPos = vec3(0.0, 0.0, 20000.0);

  vec4 target = rayMarch(normalize(p - u_Eye));
  if (target.w > -0.5) {
    if (target.w < 0.5) {
      float sineVal = sin(u_Time * M_PI * 0.001);
      float cosineVal = cos(u_Time * M_PI * 0.001);
      vec3 rocketPosition = vec3(0.0, (PLANET_RADIUS + 10.0) * cosineVal, (PLANET_RADIUS + 10.0) * sineVal);
      vec4 i = rotationMatrix(vec3(0.0, 0.0, 1.0), -u_Time * M_PI * 0.1) * rotationMatrix(vec3(1.0, 0.0, 0.0), -u_Time * M_PI * 0.001) * vec4((target.xyz - rocketPosition), 1.0);
      vec3 inversePos = i.xyz;
      if (inversePos.x * inversePos.x + inversePos.y * inversePos.y > 1.1) {
        float orangeMix = rocketTexture(inversePos.x * inversePos.x + inversePos.y * inversePos.y, inversePos.x * inversePos.x + inversePos.y * inversePos.y, inversePos.x * inversePos.x + inversePos.y * inversePos.y);
        out_Col = vec4(mix(vec3(1.0, 1.0, 1.0), vec3(1.0, 0.5, 0.0), orangeMix), 1.0);
      } else if (inversePos.z < -0.01) {
        float blackMix = rocketTexture(inversePos.x, inversePos.y, inversePos.z);
        out_Col = vec4(mix(vec3(0.0, 0.0, 0.0), vec3(0.333, 0.29, 0.235), blackMix), 1.0);
      } else if (inversePos.z > 7.1) {
        float redMix = rocketTexture(inversePos.z, inversePos.z, inversePos.z) / 2.0;
        out_Col = vec4(mix(vec3(0.0, 0.0, 0.0), vec3(1.0, 0.0, 0.0), redMix), 1.0);
      } else {
        float purpleMix = rocketTexture(inversePos.x, inversePos.y, inversePos.z);
        out_Col = vec4(mix(vec3(0.0, 0.0, 0.0), vec3(1.0, 0.0, 1.0), purpleMix), 1.0);
      }
    } else if (target.w < 1.5) {
      float height = planetTexture(target.x, target.y, target.z);
      if (height > 0.5) {
        out_Col = vec4(0.0, 0.0, 1.0, 1.0);
      } else {
        out_Col = vec4(0.0, 0.5, 0.0, 1.0);
      }
    } else if (target.w < 2.5) {
      float flameMix = rocketTexture(target.x * 3.0, target.y * 3.0, target.z * 3.0);
      out_Col = vec4(mix(vec3(1.0, 1.0, 1.0), vec3(0.902, 0.655, 0.122), flameMix), 1.0);
    } else if (target.w < 3.5) {
      float elevation = planetTexture(target.x * 100.0, target.y * 100.0, target.z * 100.0);
      out_Col = vec4(mix(vec3(0.463, 0.463, 0.463), vec3(0.353, 0.333, 0.298), bias(0.7, elevation)), 1.0);
    }
    if (target.w < 1.5 || target.w > 2.5) {
      vec3 lightVec = sunPos - target.xyz;
      vec3 normal = estimateNormal(target.xyz);
      float diffuseTerm = dot(normalize(normal), normalize(lightVec));
      diffuseTerm = clamp(diffuseTerm, 0.0, 1.0);
      float ambientTerm = 0.1;
      float lightIntensity = diffuseTerm + ambientTerm;

      vec3 cameraPosition = vec3(0.0, (PLANET_RADIUS + 20.0) * cos((u_Time - 0.25) * M_PI * 0.001), (PLANET_RADIUS + 20.0) * sin((u_Time - 0.25) * M_PI * 0.001));
      vec3 view = normalize(target.xyz - cameraPosition);
      vec3 light = normalize(target.xyz - sunPos);
      float specularIntensity = dot(reflect(-light, normal), view);
      specularIntensity = pow(max(specularIntensity, 0.0), 50.0);

      out_Col = vec4(out_Col.xyz * lightIntensity + vec3(0.0, 0.0, 0.0) * specularIntensity, 1.0);
    }
  } else {
    out_Col = vec4(0.0, 0.0, 0.0, 1.0);
  }
}
