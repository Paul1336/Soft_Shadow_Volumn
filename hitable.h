#ifndef HITABLEH
#define HITABLEH

#include "vec3.h"
#include "ray.h"

constexpr float EPSILON = 1e-6f;
struct hit_record
{
    float t;
    vec3 p;
    vec3 normal;
    float w_r, w_t;
    float n_t;
    vec3 Kd;
};

class hitable
{
public:
    virtual ~hitable() = default;
    virtual bool hit(const ray &r, float t_min, float t_max, hit_record &rec) const = 0;
};

class sphere : public hitable
{
public:
    vec3 center;
    float radius;
    vec3 Kd;
    float w_r, w_t, n_t;

    sphere() {}
    sphere(vec3 c, float r, vec3 kd, float wr = 0, float wt = 0, float nt = 1.0)
        : center(c), radius(r), Kd(kd), w_r(wr), w_t(wt), n_t(nt) {}

    virtual bool hit(const ray &r, float t_min, float t_max, hit_record &rec) const override
    {
        vec3 oc = r.origin() - center;
        float a = dot(r.direction(), r.direction());
        float b = dot(oc, r.direction());
        float c = dot(oc, oc) - radius * radius;
        float discriminant = b * b - a * c;

        if (discriminant > 0)
        {
            float temp = (-b - sqrt(discriminant)) / a;
            if (temp < t_max && temp > t_min)
            {
                rec.t = temp;
                rec.p = r.point_at_parameter(rec.t);
                rec.normal = (rec.p - center) / radius;
                rec.Kd = Kd;
                rec.w_r = w_r;
                rec.w_t = w_t;
                rec.n_t = n_t;
                return true;
            }
            temp = (-b + sqrt(discriminant)) / a;
            if (temp < t_max && temp > t_min)
            {
                rec.t = temp;
                rec.p = r.point_at_parameter(rec.t);
                rec.normal = (rec.p - center) / radius;
                rec.Kd = Kd;
                rec.w_r = w_r;
                rec.w_t = w_t;
                rec.n_t = n_t;
                return true;
            }
        }
        return false;
    }
};

class triangle : public hitable
{
public:
    vec3 v0, v1, v2;
    vec3 normal;
    vec3 Kd;

    triangle(const vec3 &_v0, const vec3 &_v1, const vec3 &_v2, const vec3 &_kd)
        : v0(_v0), v1(_v1), v2(_v2), Kd(_kd)
    {
        // 計算法向量，確保方向正確
        vec3 edge1 = v1 - v0;
        vec3 edge2 = v2 - v0;
        normal = unit_vector(cross(edge1, edge2));

        // 檢查法向量是否為零向量
        if (normal.length() < EPSILON)
        {
            normal = vec3(0, 1, 0); // 默認向上
        }
    }

    virtual bool hit(const ray &r, float t_min, float t_max, hit_record &rec) const override
    {
        vec3 edge1 = v1 - v0;
        vec3 edge2 = v2 - v0;
        vec3 h = cross(r.direction(), edge2);
        float a = dot(edge1, h);

        if (abs(a) < EPSILON)
            return false;

        float f = 1.0f / a;
        vec3 s = r.origin() - v0;
        float u = f * dot(s, h);
        if (u < 0.0f || u > 1.0f)
            return false;

        vec3 q = cross(s, edge1);
        float v = f * dot(r.direction(), q);
        if (v < 0.0f || u + v > 1.0f)
            return false;

        float t = f * dot(edge2, q);
        if (t > t_min && t < t_max)
        {
            rec.t = t;
            rec.p = r.point_at_parameter(t);

            // 確保法向量朝向光線來源方向
            rec.normal = normal;
            if (dot(rec.normal, r.direction()) > 0)
            {
                rec.normal = -rec.normal;
            }

            rec.Kd = Kd;
            rec.w_r = rec.w_t = 0.0f;
            rec.n_t = 1.0f;
            return true;
        }
        return false;
    }
};
class plane : public hitable
{
public:
    vec3 p;
    vec3 normal;
    vec3 Kd;
    virtual bool hit(const ray &r, float t_min, float t_max, hit_record &rec) const;
    plane(const vec3 &p, const vec3 &n, const vec3 &kd)
        : p(p), normal(unit_vector(n)), Kd(kd) {}
};

bool plane::hit(const ray &r, float t_min, float t_max, hit_record &rec) const
{
    float denom = dot(normal, r.direction());
    if (fabs(denom) > EPSILON)
    {
        float t = dot(p - r.origin(), normal) / denom;
        if (t >= t_min && t <= t_max)
        {
            rec.t = t;
            rec.p = r.point_at_parameter(t);
            rec.normal = denom < 0 ? normal : -normal; // 保證法向量朝外
            rec.Kd = Kd;
            rec.w_r = rec.w_t = 0.0f;
            rec.n_t = 1.0f;
            return true;
        }
    }
    return false;
}
#endif