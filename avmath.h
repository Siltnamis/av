#ifndef AV_MATH_HEADER
#define AV_MATH_HEADER

//TODO: Add a SIMD support
//TODO: Add quaternion

#include <math.h>

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

#define PI                      3.14159265358979323846264338327950288f
#define AV_MATH_PI              3.14159265358979323846264338327950288f
#define AV_MATH_EPSILON         1.19209290e-7f
#define AV_MATH_ZERO            0.0f
#define AV_MATH_ONE             1.0f
#define AV_MATH_TWO_THIRDS      0.666666666666666666666666666666666666667f

#define AV_MATH_TAU             6.28318530717958647692528676655900576f
#define AV_MATH_PI              3.14159265358979323846264338327950288f
#define AV_MATH_ONE_OVER_TAU    0.636619772367581343075535053490057448f
#define AV_MATH_ONE_OVER_PI     0.159154943091895335768883763372514362f

#define AV_MATH_TAU_OVER_2      3.14159265358979323846264338327950288f
#define AV_MATH_TAU_OVER_4      1.570796326794896619231321691639751442f
#define AV_MATH_TAU_OVER_8      0.785398163397448309615660845819875721f

#define AV_MATH_E               2.71828182845904523536f
#define AV_MATH_SQRT_TWO        1.41421356237309504880168872420969808f
#define AV_MATH_SQRT_THREE      1.73205080756887729352744634150587236f
#define AV_MATH_SQRT_FIVE       2.23606797749978969640917366873127623f

#define AV_MATH_LOG_TWO         0.693147180559945309417232121458176568f
#define AV_MATH_LOG_TEN         2.30258509299404568401799145468436421f
union vec2
{
    struct { float x, y; };
    struct { float s, t; };
    float e[2];
};

union vec3
{
    struct { float x, y, z; };
    struct { float r, g, b; };
    struct { float s, t, p; };
    vec2 xy;
    float e[3];
};

union vec4
{
    struct { float x, y, z, w; };
    struct { float r, g, b, a; };
    struct { float s, t, p, q; };
    struct { vec2 xy, zw; };
    vec3 xyz;
    vec3 rgb;
    float e[4];
};

/*inline vec4 vec4(vec3 v, float w)
{
    vec4 result;
    result.xyz = v;
    result.w = w;
    return result;
}
*/

union quat
{
    struct { float x, y, z, w; };
    vec4 xyzw;
    vec3 xyz;
    float e[4];
};

union mat2
{
    struct { vec2 x, y; };
    struct { float m00, m01, m10, m11; };
    vec2 col[2];
    float e[4];
};
inline mat2 mat2_identity()
{
    return {1.f, 0.f, 0.f, 1.f};
}

union mat3
{
    struct {vec3 x, y, z; };
    struct {
            float m00, m01, m02; // 1st column
            float m10, m11, m12; // 2nd column
            float m20, m21, m22; // 3rd column
            };
    vec3 col[3];
    float e[9];
};
inline mat3 mat3_identity()
{
    return {1, 0, 0, 0, 1, 0, 0, 0, 1};
}

union mat4
{
    struct { vec4 x, y, z, w; };
    struct {
            float m00, m01, m02, m03; // 1st column
            float m10, m11, m12, m13; // 2nd column
            float m20, m21, m22, m23; // 3rd column
            float m30, m31, m32, m33; // 4th column 
           };

    vec4 col[4];
    float e[16];
};
inline mat4 mat4_identity()
{
    mat4 result = {};
    result.m00 = 1;
    result.m11 = 1;
    result.m22 = 1;
    result.m33 = 1;
    return result;
}

//==========util==========
inline float radians(float angle)
{
    return angle*PI/180.f;
}

inline float degrees(float angle)
{
    return angle*180.f/PI;
}

//==========Vec2 stuff==========
inline vec2 operator+(vec2 A, vec2 B)
{
    vec2 result;
    result.x = A.x + B.x;
    result.y = A.y + B.y;
    return result;
}
inline vec2& operator+=(vec2& A, vec2 B)
{
    return (A = A + B);
}
inline vec2 operator-(vec2 A, vec2 B)
{
    vec2 result;
    result.x = A.x - B.x;
    result.y = A.y - B.y;
    return result;
}
inline vec2& operator-=(vec2& A, vec2 B)
{
    return (A = A - B);
}
inline vec2 operator-(vec2 A)
{
    vec2 result;
    result.x = -A.x;
    result.y = -A.y;
    return result;
}
inline vec2 operator*(vec2 v, float scalar)
{
    vec2 result;
    result.x = v.x*scalar;
    result.y = v.y*scalar;
    return result;
}
inline vec2 operator*(float scalar, vec2 v)
{
    vec2 result;
    result.x = v.x*scalar;
    result.y = v.y*scalar;
    return result;
}
inline vec2& operator*=(vec2& v, float scalar)
{
    return ( v = v*scalar);
}
inline vec2 operator/(vec2 v, float scalar)
{
    vec2 result;
    result.x = v.x/scalar;
    result.y = v.y/scalar;
    return result;
}
inline vec2& operator/=(vec2& v, float scalar)
{
    return (v = v/scalar);
}

inline float dot(vec2 A, vec2 B)
{
    return A.x*B.x + A.y*B.y;
}

inline float length(vec2 v)
{
    return sqrtf(v.x*v.x + v.y*v.y);
}

inline vec2 normalize(vec2 v)
{
    float l = length(v);
    if(l > 0)
        return v/l;
    else
        return {0.f, 0.f};
}

//==========Vec3 stuff==========
inline vec3 operator+(vec3 A, vec3 B)
{
    vec3 result;
    result.x = A.x + B.x;
    result.y = A.y + B.y;
    result.z = A.z + B.z;
    return result;
}
inline vec3& operator+=(vec3& A, vec3 B)
{
    return (A = A + B);
}
inline vec3 operator-(vec3 A, vec3 B)
{
    vec3 result;
    result.x = A.x - B.x;
    result.y = A.y - B.y;
    result.z = A.z - B.z;
    return result;
}
inline vec3& operator-=(vec3& A, vec3 B)
{
    return (A = A - B);
}
inline vec3 operator-(vec3 v)
{
    vec3 result;
    result.x = -v.x;
    result.y = -v.y;
    result.z = -v.z;
    return result;
}
inline vec3 operator*(vec3 v, float scalar)
{
    vec3 result;
    result.x = v.x * scalar;
    result.y = v.y * scalar;
    result.z = v.z * scalar;
    return result;
}
inline vec3 operator*(float scalar, vec3 v)
{
    vec3 result;
    result.x = v.x * scalar;
    result.y = v.y * scalar;
    result.z = v.z * scalar;
    return result;
}
inline vec3& operator*=(vec3& v, float scalar)
{
    return (v = v*scalar);
}

inline vec3 operator*(vec3 A, vec3 B)
{
    vec3 result;
    result.x = A.x*B.x;
    result.y = A.y*B.y;
    result.z = A.z*B.z;
    return result;
}

inline vec3 operator/(vec3 v, float scalar)
{
    vec3 result;
    result.x = v.x / scalar;
    result.y = v.y / scalar;
    result.z = v.z / scalar;
    return result;
}
inline vec3& operator/=(vec3& v, float scalar)
{
    return (v = v/scalar);
}

inline float dot(vec3 A, vec3 B)
{
    return A.x*B.x + A.y*B.y + A.z*B.z; 
}

inline vec3 cross(vec3 A, vec3 B)
{
    vec3 result;
    result.x = A.y*B.z - A.z*B.y;
    result.y = A.z*B.x - A.x*B.z;
    result.z = A.x*B.y - A.y*B.x;
    return result;
}

inline float length(vec3 v)
{
    return sqrtf(v.x*v.x + v.y*v.y + v.z*v.z);
}

inline vec3 normalize(vec3 v)
{
    float l = length(v);
    if(l > 0)
        return v/l;
    else 
        return {0.f, 0.f, 0.f};
}

inline vec3 reflect(vec3 v, vec3 n)
{
    vec3 r;
    r = v - 2*(dot(v, n))*n; 
    return r;
}

//==========Vec4 stuff==========
inline vec4 operator+(vec4 A, vec4 B)
{
    vec4 result;
    result.x = A.x + B.x;
    result.y = A.y + B.y;
    result.z = A.z + B.z;
    result.w = A.w + B.w;
    return result;
}
inline vec4& operator+=(vec4& A, vec4 B)
{
    return (A = A + B);
}
inline vec4 operator-(vec4 A, vec4 B)
{
    vec4 result;
    result.x = A.x - B.x;
    result.y = A.y - B.y;
    result.z = A.z - B.z;
    result.w = A.w - B.w;
    return result;
}
inline vec4& operator-=(vec4& A, vec4 B)
{
    return (A = A - B);
}
inline vec4 operator-(vec4 v)
{
    vec4 result;
    result.x = -v.x;
    result.y = -v.y;
    result.z = -v.z;
    result.w = -v.w;
    return result;
}
inline vec4 operator*(vec4 v, float scalar)
{
    vec4 result;
    result.x = v.x * scalar;
    result.y = v.y * scalar;
    result.z = v.z * scalar;
    result.w = v.w * scalar;
    return result;
}
inline vec4 operator*(float scalar, vec4 v)
{
    vec4 result;
    result.x = v.x * scalar;
    result.y = v.y * scalar;
    result.z = v.z * scalar;
    result.w = v.w * scalar;
    return result;
}
inline vec4& operator*=(vec4& v, float scalar)
{
    return (v = v*scalar);
}
inline vec4 operator/(vec4 v, float scalar)
{
    vec4 result;
    result.x = v.x / scalar;
    result.y = v.y / scalar;
    result.z = v.z / scalar;
    result.w = v.w / scalar;
    return result;
}
inline vec4& operator/=(vec4& v, float scalar)
{
    return (v = v/scalar);
}

inline float dot(vec4 A, vec4 B)
{
    return A.x*B.x + A.y*B.y + A.z*B.z + A.w*B.w;
}

inline float length(vec4 v)
{
    return sqrtf(v.x*v.x + v.y*v.y + v.z*v.z + v.w*v.w);
}

inline vec4 normalize(vec4 v)
{
    float l = length(v);
    if(l > 0)
        return v/l;
    else 
        return {0.f, 0.f, 0.f, 0.f};
}

//==========mat2stuff==========
inline mat2 operator+(mat2 ml, mat2 mr)
{
    return  {
                ml.m00 + mr.m00, ml.m01 + mr.m01,
                ml.m10 + mr.m10, ml.m11 + mr.m11
            };
            
}
inline mat2& operator+=(mat2& ml, mat2 mr)
{
    return (ml = ml + mr);
}
inline mat2 operator-(mat2 ml, mat2 mr)
{
    return  {
                ml.m00 - mr.m00, ml.m01 - mr.m01,
                ml.m10 - mr.m10, ml.m11 - mr.m11
            };
}
inline mat2& operator-=(mat2& ml, mat2 mr)
{
    return (ml = ml - mr);
}

inline mat2 operator*(mat2 ml, mat2 mr)
{
#if 1
    mat2 result;
    for(int j = 0; j < 2; ++j){
        for(int i = 0; i < 2; ++i){
            result.e[j*2 + i] = ml.e[2*0 + i]*mr.e[j*2 + 0] + ml.e[2*1 +i]*mr.e[j*2 + 1];
        }
    }
    return result;
#else
    
#endif
}

inline mat2& operator*=(mat2& ml, mat2 mr)
{
    return (ml = ml * mr);
}

inline mat2 operator*(mat2 m, float scalar)
{
    return {m.e[0]*scalar, m.e[1]*scalar, m.e[2]*scalar, m.e[3]*scalar};
}
inline mat2& operator*=(mat2& m, float scalar)
{
    return (m = m*scalar);
}

inline vec2 operator*(mat2 m, vec2 v)
{
    return {m.m00*v.x + m.m10*v.y, m.m01*v.x + m.m11*v.y};
}

//==========mat3stuff==========
inline mat3 operator*(mat3 ml, mat3 mr)
{
    mat3 result;
    for(int j = 0; j < 3; ++j){
        for(int i = 0; i < 3; ++i){
            result.e[j*3 + i] =     ml.e[3*0 + i] * mr.e[j*3 + 0]
                                +   ml.e[3*1 + i] * mr.e[j*3 + 1]  
                                +   ml.e[3*2 + i] * mr.e[j*3 + 2];
        }
    }
    return result;
}

inline mat3& operator*=(mat3& ml, mat3 mr)
{
    return (ml = ml * mr);
}

inline vec3 operator*(mat3 m, vec3 v)
{
    return {
            m.m00*v.x + m.m10*v.y + m.m20*v.z, 
            m.m01*v.x + m.m11*v.y + m.m21*v.z,
            m.m02*v.x + m.m12*v.y + m.m22*v.z
            };
}

//==========mat4stuff==========
inline mat4 operator*(mat4 ml, mat4 mr)
{
    mat4 result;
    for(int j = 0; j < 4; ++j){
        for(int i = 0; i < 4; ++i){
            result.e[j*4 + i] =     ml.e[4*0 + i] * mr.e[j*4 + 0]
                                +   ml.e[4*1 + i] * mr.e[j*4 + 1]
                                +   ml.e[4*2 + i] * mr.e[j*4 + 2]
                                +   ml.e[4*3 + i] * mr.e[j*4 + 3];
        }
    }
    return result;
}
inline mat4& operator*=(mat4& ml, mat4 mr)
{
    return (ml = ml * mr);
}
inline vec4 operator*(mat4 m, vec4 v)
{
    vec4 result;
    result.x = m.m00*v.x + m.m10*v.y + m.m20*v.z + m.m30*v.w;
    result.y = m.m01*v.x + m.m11*v.y + m.m21*v.z + m.m31*v.w;
    result.z = m.m02*v.x + m.m12*v.y + m.m22*v.z + m.m32*v.w;
    result.w = m.m03*v.x + m.m13*v.y + m.m23*v.z + m.m33*v.w;
    return result;
}

inline mat4 mat4_inverse(mat4 m)
{
    mat4 o;

    float ood;
    float tmp;

    float sf00 = m.m22 * m.m33 - m.m32 * m.m23;
    float sf01 = m.m21 * m.m33 - m.m31 * m.m23;
    float sf02 = m.m21 * m.m32 - m.m31 * m.m22;
    float sf03 = m.m20 * m.m33 - m.m30 * m.m23;
    float sf04 = m.m20 * m.m32 - m.m30 * m.m22;
    float sf05 = m.m20 * m.m31 - m.m30 * m.m21;
    float sf06 = m.m12 * m.m33 - m.m32 * m.m13;
    float sf07 = m.m11 * m.m33 - m.m31 * m.m13;
    float sf08 = m.m11 * m.m32 - m.m31 * m.m12;
    float sf09 = m.m10 * m.m33 - m.m30 * m.m13;
    float sf10 = m.m10 * m.m32 - m.m30 * m.m12;
    float sf11 = m.m11 * m.m33 - m.m31 * m.m13;
    float sf12 = m.m10 * m.m31 - m.m30 * m.m11;
    float sf13 = m.m12 * m.m23 - m.m22 * m.m13;
    float sf14 = m.m11 * m.m23 - m.m21 * m.m13;
    float sf15 = m.m11 * m.m22 - m.m21 * m.m12;
    float sf16 = m.m10 * m.m23 - m.m20 * m.m13;
    float sf17 = m.m10 * m.m22 - m.m20 * m.m12;
    float sf18 = m.m10 * m.m21 - m.m20 * m.m11;

    o.m00 = +(m.m11 * sf00 - m.m12 * sf01 + m.m13 * sf02);
    o.m10 = -(m.m10 * sf00 - m.m12 * sf03 + m.m13 * sf04);
    o.m20 = +(m.m10 * sf01 - m.m11 * sf03 + m.m13 * sf05);
    o.m30 = -(m.m10 * sf02 - m.m11 * sf04 + m.m12 * sf05);

    o.m01 = -(m.m01 * sf00 - m.m02 * sf01 + m.m03 * sf02);
    o.m11 = +(m.m00 * sf00 - m.m02 * sf03 + m.m03 * sf04);
    o.m21 = -(m.m00 * sf01 - m.m01 * sf03 + m.m03 * sf05);
    o.m31 = +(m.m00 * sf02 - m.m01 * sf04 + m.m02 * sf05);

    o.m02 = +(m.m01 * sf06 - m.m02 * sf07 + m.m03 * sf08);
    o.m12 = -(m.m00 * sf06 - m.m02 * sf09 + m.m03 * sf10);
    o.m22 = +(m.m00 * sf11 - m.m01 * sf09 + m.m03 * sf12);
    o.m32 = -(m.m00 * sf08 - m.m01 * sf10 + m.m02 * sf12);

    o.m03 = -(m.m01 * sf13 - m.m02 * sf14 + m.m03 * sf15);
    o.m13 = +(m.m00 * sf13 - m.m02 * sf16 + m.m03 * sf17);
    o.m23 = -(m.m00 * sf14 - m.m01 * sf16 + m.m03 * sf18);
    o.m33 = +(m.m00 * sf15 - m.m01 * sf17 + m.m02 * sf18);

    ood = 1.0f / (m.m00 * o.m00 +
                  m.m01 * o.m10 +
                  m.m02 * o.m20 +
                  m.m03 * o.m30);

    o.m00 *= ood;
    o.m01 *= ood;
    o.m02 *= ood;
    o.m03 *= ood;
    o.m10 *= ood;
    o.m11 *= ood;
    o.m12 *= ood;
    o.m13 *= ood;
    o.m20 *= ood;
    o.m21 *= ood;
    o.m22 *= ood;
    o.m23 *= ood;
    o.m30 *= ood;
    o.m31 *= ood;
    o.m32 *= ood;
    o.m33 *= ood;

    return o;
}


//==========quaternion stuff==========
inline quat operator+(quat A, quat B)
{
    quat result;
    result.xyzw = A.xyzw + B.xyzw; 
    return result;
}
inline quat& operator+=(quat& A, quat B)
{
    return (A = A + B);
}
inline quat operator-(quat A, quat B)
{
    quat result;
    result.xyzw = A.xyzw - B.xyzw;
    return result;
}
inline quat& operator-=(quat& A, quat B)
{
    return (A = A - B);
}
inline quat operator-(quat q)
{
    quat result;
    result.xyzw = -q.xyzw;
    return result;
}

inline quat operator*(quat A, quat B)
{
    quat result;
    result.x = A.w * B.x + A.x * B.w + A.y * B.z - A.z * B.y;
	result.y = A.w * B.y - A.x * B.z + A.y * B.w + A.z * B.x;
	result.z = A.w * B.z + A.x * B.y - A.y * B.x + A.z * B.w;
	result.w = A.w * B.w - A.x * B.x - A.y * B.y - A.z * B.z;
    return result;
}

inline quat& operator*=(quat& A, quat B)
{
    return (A = A*B);
}

inline quat conjugate(quat q)
{
    q.xyz = -q.xyz;
    return q;
}               
                
inline quat operator/(quat q, float scalar)
{
    quat result;
    result.xyzw = q.xyzw / scalar;
    return result;
}

inline quat& operator/=(quat& q, float scalar)
{
    return (q = q/scalar);
}

inline float dot(quat A, quat B)
{
    return dot(A.xyzw, B.xyzw);
}

inline quat quat_inverse(quat q)
{
    return conjugate(q)/dot(q, q);
}

inline quat operator/(quat A, quat B)
{
    return(quat_inverse(B)*A);
}

inline quat& operator/=(quat& A, quat B)
{
    return (A = A/B);
}
               
inline quat operator*(quat q, float scalar)
{
    quat result;
    result.xyzw = q.xyzw * scalar;
    return result;
}
inline quat operator*(float scalar, quat q)
{
    quat result;
    result.xyzw = q.xyzw * scalar;
    return result;
}

inline quat& operator*=(quat& q, float scalar)
{
    return (q = q*scalar);
}


inline float length(quat q)
{
    return length(q.xyzw);
}

inline quat normalize(quat q)
{
    float l = length(q);
    if(l > 0)
        return q/l;
    else 
        return {0.f, 0.f, 0.f, 0.f};
}

inline quat quat_identity()
{
    return quat{0, 0, 0, 1};
}


inline quat quat_axis_angle(vec3 axis, float angle)
{
    quat q;
    q.xyz = normalize(axis);
    q.xyz *= sinf(angle/2.f);
    q.w = cosf(angle/2.f);
    return q;
}

/*
inline quat quat_euler_angles(float yaw, float pitch, float roll)
{
    //assert(false);
    quat result;
#if 1
#if 1
    float t0 = cos(yaw/2.f); 
    float t1 = sin(yaw/2.f); 
    float t2 = cos(roll/2.f); 
    float t3 = sin(roll/2.f); 
    float t4 = cos(pitch/2.f); 
    float t5 = sin(pitch/2.f); 

    result.w = t0*t2*t4 + t1*t3*t5;
    result.x = t0*t3*t4 - t1*t2*t5;
    result.y = t0*t2*t5 + t1*t3*t4;
    result.z = t1*t2*t4 - t0*t3*t5;
#else
    vec3 c;
    c.x = cosf(pitch/2.f);
    c.y = cosf(yaw/2.f);
    c.z = cosf(roll/2.f);

    vec3 s;
    s.x = sinf(pitch/2.f);
    s.y = sinf(yaw/2.f);
    s.z = sinf(roll/2.f);

    result.w = c.x * c.y * c.z + s.x * s.y * s.z;
    result.x = s.x * c.y * c.z - c.x * s.y * s.z;
    result.y = c.x * s.y * c.z + s.x * c.y * s.z;
    result.z = c.x * c.y * s.z - s.x * s.y * c.z;


#endif
#else
    quat p = quat_axis_angle({1, 0, 0}, pitch);
    quat y = quat_axis_angle({0, 1, 0}, yaw);
    quat r = quat_axis_angle({0, 0, 1}, roll);
    result = y*p*r;
    result = normalize(result);
#endif
    result = normalize(result);
    return result; 
}
*/

inline quat slerp(quat q, quat r, float t)
{
 float cos_half_theta = dot (q, r);
	// as found here http://stackoverflow.com/questions/2886606/flipping-issue-when-interpolating-rotations-using-quaternions
	// if dot product is negative then one quaternion should be negated, to make
	// it take the short way around, rather than the long way
	// yeah! and furthermore Susan, I had to recalculate the d.p. after this
	if (cos_half_theta < 0.0f) {
		for (int i = 0; i < 4; i++) {
			q.e[i] *= -1.0f;
		}
		cos_half_theta = dot (q, r);
	}
	// if qa=qb or qa=-qb then theta = 0 and we can return qa
	if (fabs (cos_half_theta) >= 1.0f) {
		return q;
	}
	// Calculate temporary values
	float sin_half_theta = sqrt (1.0f - cos_half_theta * cos_half_theta);
	// if theta = 180 degrees then result is not fully defined
	// we could rotate around any axis normal to qa or qb
	quat result;
	if (fabs (sin_half_theta) < 0.001f) {
		for (int i = 0; i < 4; i++) {
			result.e[i] = (1.0f - t) * q.e[i] + t * r.e[i];
		}
		return result;
	}
	float half_theta = acos (cos_half_theta);
	float a = sin ((1.0f - t) * half_theta) / sin_half_theta;
	float b = sin (t * half_theta) / sin_half_theta;
	for (int i = 0; i < 4; i++) {
		result.e[i] = q.e[i] * a + r.e[i] * b;
	}
	return result;   
}

inline mat4 quat_to_mat4(quat q)
{
    mat4 result = mat4_identity();
    q = normalize(q);
    float xx, yy, zz, xy, xz, yz, wx, wy, wz;
    xx = q.x*q.x; yy = q.y*q.y; zz = q.z*q.z;
    xy = q.x*q.y; xz = q.x*q.z; yz = q.y*q.z;
    wx = q.w*q.x; wy = q.w*q.y; wz = q.w*q.z;

    result.m00 = 1.f - 2.f*(yy + zz);
    result.m01 = 2.f*(xy + wz);
    result.m02 = 2.f*(xz - wy);

    result.m10 = 2.f*(xy - wz);
    result.m11 = 1.f - 2.f*(xx + zz);
    result.m12 = 2.f*(yz + wx);

    result.m20 = 2.f*(xz + wy);
    result.m21 = 2.f*(yz - wx);
    result.m22 = 1.f - 2.f*(xx + yy);
    
    return result;
};




//==========gl related mat stuff==========
inline mat4 look_at(vec3 eye, vec3 center, vec3 up)
{
    vec3 f, s, u;
    f = center - eye;
    f = normalize(f);
    s = cross(f, up);
    s = normalize(s);
    u = cross(s, f);

    mat4 result = mat4_identity();

    result.m00 = +s.x;
    result.m10 = +s.y;
    result.m20 = +s.z;

    result.m01 = +u.x;
    result.m11 = +u.y;
    result.m21 = +u.z;

    result.m02 = -f.x;
    result.m12 = -f.y;
    result.m22 = -f.z;

    result.m30 = -dot(s, eye);
    result.m31 = -dot(u, eye);
    result.m32 =  dot(f, eye);

    return result;
}

inline mat4 perspective(float fov, float aspect_ratio, float z_near, float z_far)
{
    mat4 result = {};
    float tan_half_fov = tanf(fov/2.f);
    
    result.m00 = 1.f / (aspect_ratio*tan_half_fov);
    result.m11 = 1.f / (tan_half_fov);
    result.m22 = -(z_far + z_near)/(z_far - z_near);
    result.m23 = -1.f;
    result.m32 = -2.f*z_far*z_near/(z_far - z_near);

    return result;
}

inline mat3 mat3_ortho(float left, float right, float bottom, float top)
{
    mat3 result = {};

    result.m00 = 2.f/(right - left);
    result.m11 = 2.f/(top - bottom);
    result.m22 = 1.f;

    result.m20 = -(right + left)/(right - left);
    result.m21 = -(top + bottom)/(top - bottom);
    //result.m20 = -0.5;

    return result;
}

inline mat4 mat4_ortho(float left, float right, float bottom, float top, float nearz, float farz)
{
    mat4 result = {};
    result.m00 = 2.f/(right - left);
    result.m11 = 2.f/(top - bottom);
    result.m22 = -2.f/(farz - nearz);
    result.m30 = -(right + left)/(right - left);
    result.m31 = -(top + bottom)/(top - bottom);
    result.m32 = -(farz + nearz)/(farz - nearz);
    result.m33 = 1.f;
    return result;
}

//==========gl transformation stuff==========

inline mat4 mat4_translation(vec3 pos)
{
    mat4 result = mat4_identity();
    result.col[3].xyz = pos;
    return result;
}

inline mat4 mat4_scale(vec3 size)
{
    mat4 result = {};
    result.m00 = size.x;
    result.m11 = size.y;
    result.m22 = size.z;
    result.m33 = 1;
    return result;
}

inline mat4 mat4_scale(float size)
{
    return mat4_scale(vec3{size, size, size});
}

inline mat4 mat4_rotationx(float angle)
{
    mat4 result = {};
    result.m33 = 1;

    result.m00 = 1;

    result.col[1].y = cos(angle);
    result.col[1].z = sin(angle);

    result.col[2].y = -sin(angle);
    result.col[2].z = cos(angle);

    return result;
}

inline mat4 mat4_rotationy(float angle)
{
    mat4 result = {};
    result.m33 = 1;

    result.col[0].x = cos(angle);
    result.col[0].z = -sin(angle);

    result.m11 = 1;

    result.col[2].x = sin(angle);
    result.col[2].z = cos(angle);

    return result;
}

inline mat4 mat4_rotationz(float angle)
{
    mat4 result = {};
    result.m33 = 1;

    result.col[0].x = cos(angle);
    result.col[0].y = sin(angle);

    result.col[1].x = -sin(angle);
    result.col[1].y = cos(angle);

    result.m22 = 1;

    return result;
}

inline mat4 mat4_orient(vec3 look, vec3 up)
{
    mat4 rot = mat4_identity();

    vec3 right = cross(look, up);
    up = cross(right, look);
    rot.x.xyz = normalize(right);
    rot.y.xyz = normalize(look);
    rot.z.xyz = normalize(up);

    return rot;
}



#endif //AV_MATH_HEADER