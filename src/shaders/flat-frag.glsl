#version 300 es
precision highp float;

uniform vec3 u_Eye, u_Ref, u_Up;
uniform vec2 u_Dimensions;
uniform float u_Time;

in vec2 fs_Pos;
out vec4 out_Col;

#define EPSILON 0.000001
#define MAX_DIST 100.0
#define MAX_MARCHING_STEPS 500

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

float dot2(in vec3 v) {
  return dot(v,v);
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

float smin(float a, float b, float k) {
  float h = clamp(0.5 + 0.5 * (b - a) / k, 0.0, 1.0);
  return mix(b, a, h) - k * h * (1.0 - h);
}

float sdRocket(vec3 pos) {
  vec3 rocketAxis = vec3(0.0, 1.0, 0.0);
  vec3 bodyBottom = vec3(0.0, 0.0, 0.0);
  vec3 bodyTop = vec3(0.0, 5.0, 0.0);
  vec3 noseTop = vec3(0.0, 7.5, 0.0);
  float bodyRadius = 1.0;
  float noseRadius = 0.2;
  float finThickness = 0.1;
  float finHeight = 2.0;

  float body = sdCylinder(pos, bodyBottom, bodyTop, bodyRadius);
  float nose = sdRoundCone(pos, bodyTop, noseTop, bodyRadius, noseRadius);
  float bodyNose = smin(body, nose, 0.1);

  vec3 fin1a1 = bodyBottom + vec3(bodyRadius, finHeight * 0.1, -finThickness * 0.5);
  vec3 fin1b1 = bodyBottom + vec3(bodyRadius, finHeight, -finThickness * 0.5);
  vec3 fin1c1 = bodyBottom + vec3(bodyRadius * 2.5, finHeight * -0.25, -finThickness * 0.5);
  vec3 fin1a2 = bodyBottom + vec3(bodyRadius, finHeight * 0.1, finThickness * 0.5);
  vec3 fin1b2 = bodyBottom + vec3(bodyRadius, finHeight, finThickness * 0.5);
  vec3 fin1c2 = bodyBottom + vec3(bodyRadius * 2.5, finHeight * -0.25, finThickness * 0.5);
  float fin1 = sdTriangularPrism(pos, fin1a1, fin1b1, fin1c1, fin1a2, fin1b2, fin1c2);

  vec4 fin2a1 = rotationMatrix(rocketAxis, radians(120.0)) * vec4(fin1a1, 1.0);
  vec4 fin2b1 = rotationMatrix(rocketAxis, radians(120.0)) * vec4(fin1b1, 1.0);
  vec4 fin2c1 = rotationMatrix(rocketAxis, radians(120.0)) * vec4(fin1c1, 1.0);
  vec4 fin2a2 = rotationMatrix(rocketAxis, radians(120.0)) * vec4(fin1a2, 1.0);
  vec4 fin2b2 = rotationMatrix(rocketAxis, radians(120.0)) * vec4(fin1b2, 1.0);
  vec4 fin2c2 = rotationMatrix(rocketAxis, radians(120.0)) * vec4(fin1c2, 1.0);
  float fin2 = sdTriangularPrism(pos, fin2a1.xyz, fin2b1.xyz, fin2c1.xyz, fin2a2.xyz, fin2b2.xyz, fin2c2.xyz);

  vec4 fin3a1 = rotationMatrix(rocketAxis, radians(240.0)) * vec4(fin1a1, 1.0);
  vec4 fin3b1 = rotationMatrix(rocketAxis, radians(240.0)) * vec4(fin1b1, 1.0);
  vec4 fin3c1 = rotationMatrix(rocketAxis, radians(240.0)) * vec4(fin1c1, 1.0);
  vec4 fin3a2 = rotationMatrix(rocketAxis, radians(240.0)) * vec4(fin1a2, 1.0);
  vec4 fin3b2 = rotationMatrix(rocketAxis, radians(240.0)) * vec4(fin1b2, 1.0);
  vec4 fin3c2 = rotationMatrix(rocketAxis, radians(240.0)) * vec4(fin1c2, 1.0);
  float fin3 = sdTriangularPrism(pos, fin3a1.xyz, fin3b1.xyz, fin3c1.xyz, fin3a2.xyz, fin3b2.xyz, fin3c2.xyz);

  float withFin1 = smin(bodyNose, fin1, 0.2);
  float withFin2 = smin(withFin1, fin2, 0.2);
  float withFin3 = smin(withFin2, fin3, 0.2);
  return withFin3;
}

float rayMarch(vec3 rayDirec) {
  float depth = 0.0;
  for (int i = 0; i < MAX_MARCHING_STEPS; i++) {
    float dist = sdRocket(u_Eye + depth * rayDirec);
    if (dist < EPSILON) {
      return depth;
    }
    depth += dist;
    if (depth >= MAX_DIST) {
      return MAX_DIST;
    }
  }
  return MAX_DIST;
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

  float dist = rayMarch(normalize(p - u_Eye));
  if (dist < MAX_DIST) {
    out_Col = vec4(0.0, 0.0, 1.0, 1.0);
  } else {
    out_Col = vec4(0.5, 0.5, 0.5, 1.0);
  }
  // out_Col = vec4(0.5 * (fs_Pos + vec2(1.0)), 0.5 * (sin(u_Time * 3.14159 * 0.01) + 1.0), 1.0);
}
