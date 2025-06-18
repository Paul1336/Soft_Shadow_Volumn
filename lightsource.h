#ifndef LIGHTSOURCEH
#define LIGHTSOURCEH
#include "vec3.h"
#include "vec2.h"
#include <vector>
#include <limits.h>
#include "edge.h"
using namespace std;

struct QueryResponse
{
    float percentage;
    vec3 pos;
    long long time;
};

class LightSource
{
public:
    vec3 light_position;
    vec3 light_u;
    vec3 light_v;
    vec3 light_intensity;
    vec3 normal;
    int resolution;
    std::vector<std::vector<int>> buffer;
    std::vector<std::vector<bool>> visited;
    LightSource(vec3 position = vec3(10, 40, 20), vec3 u = vec3(10.0f, 0, 0), vec3 v = vec3(0, 0, 10.0f), vec3 intensity = vec3(1, 1, 1), int res = 64)
        : light_position(position), light_u(u), light_v(v), light_intensity(intensity), resolution(res)
    {
        buffer.assign(resolution, vector<int>(resolution, 0));
        normal = unit_vector(cross(light_u, light_v));
        visited.assign(resolution, std::vector<bool>(resolution, true));
    }
    ~LightSource() = default;
    void initBuffer()
    {
        buffer.assign(resolution, std::vector<int>(resolution, 0));
    }
    vec2 worldToUV(const vec3 &q, const vec3 &p) const
    {
        vec3 n = unit_vector(normal);
        vec3 dir = q - p;
        float denom = dot(dir, n);
        float t = dot(light_position - p, n) / denom;
        vec3 p_proj = p + (t * dir);

        vec3 d = p_proj - light_position;
        float uu = (d.x() * light_u.x() + d.y() * light_u.y() + d.z() * light_u.z()) / (light_u.x() * light_u.x() + light_u.y() * light_u.y() + light_u.z() * light_u.z());

        uu = uu * resolution;
        float vv = (d.x() * light_v.x() + d.y() * light_v.y() + d.z() * light_v.z()) / (light_v.x() * light_v.x() + light_v.y() * light_v.y() + light_v.z() * light_v.z());

        vv = vv * resolution;

        return vec2(uu, vv);
    }
    vector<vector<vec3>> getSamplePositions() const
    {
        vector<vector<vec3>> samples(
            resolution,
            vector<vec3>(resolution));
        // (i, j) 映射到 u_offset ∈ [-0.5,0.5], v_offset ∈ [-0.5,0.5]
        for (int j = 0; j < resolution; ++j)
        {
            float v_off = (j / float(resolution - 1)) - 0.5f;
            for (int i = 0; i < resolution; ++i)
            {
                float u_off = (i / float(resolution - 1)) - 0.5f;
                samples[j][i] = light_position + light_u * u_off + light_v * v_off;
            }
        }
        return samples;
    }
    long long update(const Edge &e, const vec3 &ref, const vec3 &p)
    {
        auto update_start = chrono::high_resolution_clock::now();
        vec2 ref_uv = worldToUV(ref, p);
        vec2 h = worldToUV(e.v1, p);
        vec2 l = worldToUV(e.v2, p);
        if (l.y() > h.y())
        {
            vec2 tmp = h;
            h = l;
            l = tmp;
        }
        float yhat = (h.y() - l.y()) / (h.x() - l.x()) * (ref.x() - l.x());
        float diff = yhat - ref.y();
        int sign = 1;
        if (diff > 1e-6f)
        { // 左側(往右-)
            sign = -1;
        }

        int lv = max(0, (int)floor(l.y() + 0.5));
        int x0 = round(l.x()), y0 = round(l.y());
        int x1 = round(h.x()), y1 = round(h.y());
        int dx = abs(x1 - x0), dy = -abs(y1 - y0);
        int sx = (x0 < x1) ? 1 : -1;
        int sy = (y0 < y1) ? 1 : -1;
        int err = dx + dy;
        while (y0 < resolution)
        {
            int y0_eps = y0;
            if (x1 == x0 && y1 == y0)
                break;
            int e2 = 2 * err;
            if (e2 >= dy)
            {
                err += dy;
                x0 += sx;
            }
            if (e2 <= dx)
            {
                err += dx;
                y0 += sy;
                if (y0 >= 0 && y0 < resolution)
                    for (int i = max(x0, 0); i < resolution; ++i)
                        buffer[y0][i] += sign;
            }
        }
        auto update_end = chrono::high_resolution_clock::now();
        return chrono::duration_cast<chrono::nanoseconds>(update_end - update_start).count();
    }

    QueryResponse integrate()
    {
        auto integrate_start = chrono::high_resolution_clock::now();
        int d_min = INT_MAX;
        vec2 min_pos;
        int min_n = 0;
        for (int y = 0; y < resolution; ++y)
            for (int x = 0; x < resolution; ++x)
            {
                if (buffer[y][x] < d_min)
                {
                    d_min = buffer[y][x];
                    min_pos = vec2(x, y);
                    min_n = 1;
                }
                else if (buffer[y][x] == d_min)
                {
                    min_n++;
                }
            }
        visited.assign(resolution, std::vector<bool>(resolution, true));
        vector<vec2> next;
        for (int y = 0; y < resolution; ++y)
            for (int x = 0; x < resolution; ++x)
            {
                if (buffer[y][x] == d_min)
                {
                    visited[y][x] = false;
                }
                else
                {
                    visited[y][x] = true;
                }
            }
        int max_size = INT_MIN;
        vec2 mid_point;
        for (int y = 0; y < resolution; ++y)
            for (int x = 0; x < resolution; ++x)
            {
                next.push_back(vec2(x, y));
                int size = 0;
                vec2 mid = vec2(0, 0);
                while (next.size() > 0)
                {
                    vec2 cur = next.back();
                    next.pop_back();
                    if (!visited[(int)cur.y()][(int)cur.x()])
                    {
                        visited[(int)cur.y()][(int)cur.x()] = 1;
                        mid += cur;
                        size++;
                        int y_p = cur.y() + 1;
                        int x_p = cur.x();
                        if (y_p >= 0 && y_p < resolution && x_p >= 0 && x_p < resolution)
                        {
                            next.push_back(vec2(x_p, y_p));
                        }
                        y_p = cur.y();
                        x_p = cur.x() + 1;
                        if (y_p >= 0 && y_p < resolution && x_p >= 0 && x_p < resolution)
                        {
                            next.push_back(vec2(x_p, y_p));
                        }
                        y_p = cur.y() - 1;
                        x_p = cur.x();
                        if (y_p >= 0 && y_p < resolution && x_p >= 0 && x_p < resolution)
                        {
                            next.push_back(vec2(x_p, y_p));
                        }
                        y_p = cur.y();
                        x_p = cur.x() - 1;
                        if (y_p >= 0 && y_p < resolution && x_p >= 0 && x_p < resolution)
                        {
                            next.push_back(vec2(x_p, y_p));
                        }
                    }
                }
                if (size > max_size)
                {
                    mid_point = mid / size;
                    max_size = size;
                }
            }
        if (((int)mid_point.x() >= 0 && (int)mid_point.x() < resolution &&
             (int)mid_point.y() >= 0 && (int)mid_point.y() < resolution) &&
            buffer[(int)mid_point.x()][(int)mid_point.y()] == d_min)
            min_pos = vec2(mid_point.x(), mid_point.y());

        float percentage = float(min_n) / float(resolution * resolution);
        float u = (min_pos.x() / float(resolution)) - 0.5f;
        float v = (min_pos.y() / float(resolution)) - 0.5f;
        vec3 pos = light_position + light_u * u + light_v * v;
        auto integrate_end = chrono::high_resolution_clock::now();

        return {percentage, pos, chrono::duration_cast<chrono::nanoseconds>(integrate_end - integrate_start).count()};
    }
};

#endif
