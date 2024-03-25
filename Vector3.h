#include <math.h>
#include <algorithm>

//A class that represents a Vector, Point or Normal in 3D Space
template<typename T>
class Vector3
{
  public:
    Vector3() : x(T(0)), y(T(0)), z(T(0)) {}
    Vector3(const T &vv) : x(T(vv)), y(T(vv)), z(T(vv)) {}
    Vector3(T xx, T yy, T zz) : x(T(xx)), y(T(yy)), z(T(zz)) {}
    T x, y, z;

    //Returns the length of our vector
    T length()
    {
      return sqrt(x*x + y*y + z*z);
    } 

    //Returns the dot product of this and another vector
    T dot(const Vector3<T> &v) const
    {
      return x*v.x + y*v.y + z*v.z;
    }

    //Returns the cross product of this and another vector
    Vector3<T> cross(const Vector3<T> &v) const
    {
      return Vector3<T>(
          y*v.z - z*v.y,
          z*v.x - x*v.z,
          x*v.y - y*v.x);

    }

    //Overloaded arithmitic operators
    Vector3<T> operator + (const Vector3<T> &v) const
    { return Vector3<T>(x+v.x, y+v.y, z+v.z); }

    Vector3<T> operator - (const Vector3<T> &v) const
    { return Vector3<T>(x-v.x, y-v.y, z-v.z); }
    
    //Scale vector operator
    Vector3<T> operator * (const T &r) const
    { return Vector3<T>(x*r, y*r, z*r); }

    //[] operators
    const T& operator [] (uint8_t i) const { return (&x)[i]; }
    T& operator [] (uint8_t i) { return (&x)[i]; }
};

template<typename T>
std::ostream& operator<<(std::ostream& os, const Vector3<T> &v)
{
  os << "X:"<< v.x << ", Y: "<< v.y << ", Z:"<< v.z;
  return os;
}

typedef Vector3<float> Vector3f;

//normalize a vector in place
template<typename T>
void normalize(Vector3<T> &v)
{
  T len = v.length();
  if(len > 0)
  {
    T inverseLen = 1/len;
    v.x *= inverseLen, v.y *= inverseLen, v.z *= inverseLen;
  }
}

template<typename T>
T dot(const Vector3<T> &a, const Vector3<T> &b)
{ return a.x*b.x + a.y*b.y + a.z*b.z; }

template<typename T>
Vector3<T> cross(const Vector3<T> &a, const Vector3<T> &b)
{
  return Vector3<T>(
      a.y*b.z - a.z*b.y,
      a.z*b.x - a.x*b.z,
      a.x*b.y - a.y*b.x);
}

//Convert's a spherical coordinate encoding to cartesian coordiantes
template<typename T>
Vector3<T> sphericalToCartesian(const T &theta, const T &phi)
{
  return Vector3<T>(cos(phi)*sin(theta), sin(phi)*sin(theta), cos(theta));
};

//Get theta from cartesian coordinates
template<typename T>
inline T sphericalTheta(const Vector3<T> &v)
{
  return acos(clamp<T>(v[2], -1, 1));
}

template<typename T>
inline T sphericalPhi(const Vector3<T> &v)
{
  T phi = atan2(v[1], v[0]);
  return (phi < 0) ? phi+2*M_PI : phi;
}

template<typename T> inline T cosTheta(const Vector3<T> &w) { return cos(w[2]); }

template<typename T>
inline T sinTheta2(const Vector3<T> &w)
{
  return std::max(T(0), 1 - cosTheta(w) * cosTheta(w));
}

template<typename T>
inline T sinTheta(const Vector3<T> &w)
{
  return sqrt(sinTheta2(w));
}

template<typename T>
inline T cosPhi(const Vector3<T> &w)
{
  T sinTheta = sinTheta(w);
  if(sinTheta) { return clamp<T>(w[0]/sinTheta, -1, 1); }
  else return 1;
}

template<typename T>
inline T sinPhi(const Vector3<T> &w)
{
  T sinTheta = sinTheta(w);
  if(sinTheta) { return clamp<T>(w[1]/sinTheta, -1, 1); }
  else return 1;
}
