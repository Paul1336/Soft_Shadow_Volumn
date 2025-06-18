#ifndef VEC2H
#define VEC2H

#include <cmath>
#include <iostream>

class vec2
{
public:
    vec2() : x_(0), y_(0) {}
    vec2(float x, float y) : x_(x), y_(y) {}

    inline float x() const { return x_; }
    inline float y() const { return y_; }

    inline float operator[](int i) const { return i == 0 ? x_ : y_; }
    inline float &operator[](int i) { return i == 0 ? x_ : y_; }

    inline vec2 operator-() const { return vec2(-x_, -y_); }

    inline vec2 &operator+=(const vec2 &v)
    {
        x_ += v.x_;
        y_ += v.y_;
        return *this;
    }

    inline vec2 &operator-=(const vec2 &v)
    {
        x_ -= v.x_;
        y_ -= v.y_;
        return *this;
    }

    inline vec2 &operator*=(float s)
    {
        x_ *= s;
        y_ *= s;
        return *this;
    }

    inline vec2 &operator/=(float s)
    {
        float inv = 1.0f / s;
        x_ *= inv;
        y_ *= inv;
        return *this;
    }

    inline float length() const
    {
        return std::sqrt(x_ * x_ + y_ * y_);
    }

    inline float squared_length() const
    {
        return x_ * x_ + y_ * y_;
    }

private:
    float x_, y_;
};

// Operators
inline vec2 operator+(const vec2 &a, const vec2 &b)
{
    return vec2(a.x() + b.x(), a.y() + b.y());
}

inline vec2 operator-(const vec2 &a, const vec2 &b)
{
    return vec2(a.x() - b.x(), a.y() - b.y());
}

inline vec2 operator*(const vec2 &a, float s)
{
    return vec2(a.x() * s, a.y() * s);
}

inline vec2 operator*(float s, const vec2 &a)
{
    return vec2(a.x() * s, a.y() * s);
}

inline vec2 operator/(const vec2 &a, float s)
{
    return vec2(a.x() / s, a.y() / s);
}

inline float dot(const vec2 &a, const vec2 &b)
{
    return a.x() * b.x() + a.y() * b.y();
}

inline std::ostream &operator<<(std::ostream &os, const vec2 &v)
{
    os << "(" << v.x() << ", " << v.y() << ")";
    return os;
}
inline bool operator<(const vec2 &a, const vec2 &b)
{
    if (a.x() != b.x())
        return a.x() < b.x();
    return a.y() < b.y();
}

inline bool operator>(const vec2 &a, const vec2 &b)
{
    return b < a;
}

inline bool operator<=(const vec2 &a, const vec2 &b)
{
    return !(b < a);
}

inline bool operator>=(const vec2 &a, const vec2 &b)
{
    return !(a < b);
}
#endif // VEC2H
