#ifndef HEMICUBEH
#define HEMICUBEH

#include "vec3.h"
#include "vec2.h"
#include <vector>
#include <cmath>
#include "edge.h"
#include "lightsource.h"
#include <iostream>
#include <algorithm>
#include <limits>
using namespace std;

class Hemicube
{
public:
    vec3 position;
    vec3 normal;
    vec3 up, right;
    int resolution = 256;
    float size = 50.0f;
    std::vector<std::vector<std::vector<int>>> buffer;

    Hemicube(const vec3 corner[4], LightSource lightSource, vec3 pos, const vec3 norm, float s = 2.0f, int res = 64)
        : position(pos), normal(unit_vector(norm)), size(s), resolution(res)
    {
        vec3 light_center = lightSource.light_position + 0.5 * (lightSource.light_u + lightSource.light_v);
        vec3 border[4];
        position = vec3(0, 0, 0);
        float _max[2] = {FLT_MIN, FLT_MIN};
        float _min[2] = {FLT_MAX, FLT_MAX};
        for (int i = 0; i < 4; i++)
        {
            vec3 p = corner[i];
            vec3 dir = p - light_center;
            float denom = dot(dir, normal);
            if (fabs(denom) < 1e-6f)
            {
                // 射線與平面平行，不相交或在平面上
                border[i] = vec3(INFINITY, INFINITY, INFINITY);
            }
            else
            {
                float t = dot(position - p, normal) / denom;
                border[i] = p + dir * t;
            }
            if (border[i].x() < _min[0])
            {
                _min[0] = border[i].x();
            }
            if (border[i].x() > _max[0])
            {
                _max[0] = border[i].x();
            }
            if (border[i].z() < _min[1])
            {
                _min[1] = border[i].z();
            }
            if (border[i].z() > _max[1])
            {
                _max[1] = border[i].z();
            }
            position += border[i];
        }
        position /= 4;
        cout << "hemicube center: " << position << endl;
        size = max((_max[0] - _min[0]), (_max[1] - _min[1])) * 10;
        cout << "hemicube size: " << size << endl;
        right = unit_vector(cross(normal, vec3(0, 1, 0)));
        up = unit_vector(cross(right, normal));
        buffer.resize(resolution, std::vector<std::vector<int>>(resolution));
    }
    ~Hemicube() = default;
    bool is_silhouette(const LightSource &lightSource, const vec3 &v1, const vec3 &v2, const vec3 &n1, const vec3 &n2)
    {
        vec3 corner[4] = {lightSource.light_position, lightSource.light_position + lightSource.light_u, lightSource.light_position + lightSource.light_v,
                          lightSource.light_position + lightSource.light_u + lightSource.light_v};
        for (int i = 0; i < 4; i++)
        {
            vec3 lightDir = unit_vector((v1 + v2) * 0.5f - corner[i]);
            if (dot(lightDir, n1) * dot(lightDir, n2) < 0)
            {
                return true;
            }
        }
        return false; // 一正一負，表示方向相反
    }

    vec2 world_to_buffer(const vec3 &p) const
    {
        float px = p.x(), pz = p.z();
        if (p.y() != 0)
        {
            px = (p.x() * right.x() + p.y() * right.y() + p.z() * right.z()) / (right.x() * right.x() + right.y() * right.y() + right.z() * right.z());
            pz = (p.x() * up.x() + p.y() * up.y() + p.z() * up.z()) / (up.x() * up.x() + up.y() * up.y() + up.z() * up.z());
        }

        float half_size = size / 2.0f;
        float x_norm = (px - position.x() + half_size) / size; // 映射到 0~1
        float z_norm = (pz - position.z() + half_size) / size;
        float x_pixel = x_norm * resolution;
        float z_pixel = z_norm * resolution;

        return vec2(x_pixel, z_pixel);
    }
    vec3 buffer_to_world(const vec2 &pixel) const
    {
        // 1. 把像素座標映射回 [0,1]
        float x_norm = pixel.x() / float(resolution);
        float z_norm = pixel.y() / float(resolution);

        // 2. 回到 [-half_size, +half_size]
        float half_size = size / 2.0f;
        float x_local = x_norm * size - half_size;
        float z_local = z_norm * size - half_size;

        // 3. 加上平面中心位置
        float x_world = x_local + position.x();
        float z_world = z_local + position.z();

        // 4. y 座標保持在平面上（跟 position.y 一樣）
        return vec3(x_world, position.y(), z_world);
    }
    bool simple_query(const vec3 &p, LightSource &lightSource, map<Edge, pair<int, int>> &edgeMap, map<int, Edge> &edgeIdMap, vector<unique_ptr<hitable>> &obj)
    {
        vec2 id = world_to_buffer(p);
        vector<int> edge_ids = {};
        if (id.y() < resolution && id.y() >= 0 && id.x() < resolution && id.x() >= 0)
            edge_ids = buffer[id.y()][id.x()];

        return (edge_ids.size() > 0);
    }
    QueryResponse query(const vec3 &p, LightSource &lightSource, map<Edge, pair<int, int>> &edgeMap,
                        map<int, Edge> &edgeIdMap, vector<unique_ptr<hitable>> &obj, long long &time)
    {

        lightSource.initBuffer();
        vec2 id = world_to_buffer(p);
        vector<int> edge_ids = {};
        if (id.y() < resolution && id.y() >= 0 && id.x() < resolution && id.x() >= 0)
        {
            edge_ids = buffer[id.y()][id.x()];
        }

        sort(edge_ids.begin(), edge_ids.end());
        int cur = -1;
        for (int i = 0; i < edge_ids.size(); ++i)
        {
            if (edge_ids[i] != cur)
            {

                cur = edge_ids[i];
                Edge e = edgeIdMap[cur];
                auto itF = edgeMap.find(e);
                if (itF == edgeMap.end())
                    continue;
                hitable *h1 = obj[itF->second.first].get();
                hitable *h2 = obj[itF->second.second].get();
                triangle *tri1 = dynamic_cast<triangle *>(h1);
                triangle *tri2 = dynamic_cast<triangle *>(h2);
                float d = dot(lightSource.normal, tri1->normal);
                if (d < 1e-6f)
                    time += lightSource.update(e, ((tri1->v0 + tri1->v1) / 2 + tri1->v2) / 2, p);
                else
                    time += lightSource.update(e, ((tri2->v0 + tri2->v1) / 2 + tri2->v2) / 2, p);
            }
        }

        return {lightSource.integrate()};
    }
    void Bresenhamline(float _x0, float _y0, float _x1, float _y1,
                       const int &edge_id,
                       std::vector<std::vector<std::vector<int>>> &buffer)
    {
        int x0 = round(_x0), y0 = round(_y0);
        int x1 = round(_x1), y1 = round(_y1);
        int dx = abs(x1 - x0), dy = -abs(y1 - y0);
        int sx = (x0 < x1) ? 1 : -1;
        int sy = (y0 < y1) ? 1 : -1;
        int err = dx + dy;

        while (true)
        {
            if (x0 >= 0 && y0 >= 0 && x0 < resolution && y0 << resolution)
            {
                buffer[y0][x0].push_back(edge_id);
            }
            int y0_eps = y0;
            if (((abs(_y1 - _y0) / abs(_x1 - _x0)) * abs((float)x0) - _x0) < y0)
                y0_eps -= 1;
            else
                y0_eps += 1;
            if (x0 >= 0 && y0_eps >= 0 && x0 < resolution && y0_eps << resolution)
                buffer[y0_eps][x0].push_back(edge_id);
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
            }
        }
    }
    void rasterize_triangle(
        const vec2 &v0, const vec2 &v1, const vec2 &v2,
        const int &edge_id,
        std::vector<std::vector<std::vector<int>>> &buffer)
    {
        float x0 = v0.x(), x1 = v1.x(), x2 = v2.x();
        float y0 = v0.y(), y1 = v1.y(), y2 = v2.y();
        Bresenhamline(x0, y0, x1, y1, edge_id, buffer);
        Bresenhamline(x1, y1, x2, y2, edge_id, buffer);
        Bresenhamline(x0, y0, x2, y2, edge_id, buffer);

        float mnx = std::min(std::min(x0, x1), x2);
        float mxx = std::max(std::max(x0, x1), x2);
        float mny = std::min(std::min(y0, y1), y2);
        float mxy = std::max(std::max(y0, y1), y2);

        int minX = std::max(0, int(std::floor(mnx) - 1));
        int maxX = std::min(int(buffer[0].size()) - 1, int(std::ceil(mxx) + 1));
        int minY = std::max(0, int(std::floor(mny) - 1));
        int maxY = std::min(int(buffer.size()) - 1, int(std::ceil(mxy) + 1));
        float denom = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0);

        for (int y = minY; y <= maxY; ++y)
        {
            for (int x = minX; x <= maxX; ++x)
            {
                float w0 = (x1 - x) * (y2 - y) - (x2 - x) * (y1 - y);
                float w1 = (x2 - x) * (y0 - y) - (x0 - x) * (y2 - y);
                float w2 = 1 - w0 - w1;

                if (w0 >= 0 && w1 >= 0 && w2 >= 0)
                {
                    buffer[y][x].push_back(edge_id);
                }
            }
        }
    }
    static float cross2D(const vec2 &a, const vec2 &b, const vec2 &c)
    {
        return (b.x() - a.x()) * (c.y() - a.y()) - (b.y() - a.y()) * (c.x() - a.x());
    }

    // Monotone Chain Convex Hull: returns points in CCW order
    std::vector<vec2> convexHull(std::vector<vec2> pts)
    {
        int n = int(pts.size());
        if (n < 3)
            return pts;
        std::sort(pts.begin(), pts.end());

        std::vector<vec2> lower;
        for (const auto &p : pts)
        {
            while (lower.size() >= 2 && cross2D(lower[lower.size() - 2], lower.back(), p) <= 0)
                lower.pop_back();
            lower.push_back(p);
        }

        std::vector<vec2> upper;
        for (int i = n - 1; i >= 0; --i)
        {
            const auto &p = pts[i];
            while (upper.size() >= 2 && cross2D(upper[upper.size() - 2], upper.back(), p) <= 0)
                upper.pop_back();
            upper.push_back(p);
        }

        lower.pop_back();
        upper.pop_back();
        lower.insert(lower.end(), upper.begin(), upper.end());
        return lower;
    }
    bool isPointInConvexPolygon(const vec2 &p, const std::vector<vec2> &poly)
    {
        int n = int(poly.size());
        if (n < 3)
            return false;
        bool isCCW = cross2D(poly[0], poly[1], poly[2]) > 0;
        for (int i = 0; i < n; ++i)
        {
            const vec2 &a = poly[i];
            const vec2 &b = poly[(i + 1) % n];
            float cval = cross2D(a, b, p);
            if (isCCW)
            {
                if (cval < 0)
                    return false;
            }
            else
            {
                if (cval > 0)
                    return false;
            }
        }
        return true;
    }

    std::vector<vec2> get_polygon(const vec2 input_pts[8])
    {
        std::vector<vec2> vertices = {};
        for (int i = 0; i < 8; i++)
        {
            vertices.push_back(input_pts[i]);
        }
        return convexHull(vertices);
    }
    void project_and_store_edge(const LightSource &lightSource, const vec3 &v1, const vec3 &v2, const int &edge_id)
    {
        vec3 light_corner[4] = {lightSource.light_position, lightSource.light_position + lightSource.light_u, lightSource.light_position + lightSource.light_v,
                                lightSource.light_position + lightSource.light_u + lightSource.light_v};
        vec3 project_corner[8];
        for (int i = 0; i < 4; i++)
        {
            vec3 p = light_corner[i];
            vec3 dir[2] = {v1 - light_corner[i], v2 - light_corner[i]};
            for (int j = 0; j < 2; j++)
            {
                float denom = dot(dir[j], normal);
                if (fabs(denom) < 1e-6f)
                {
                    project_corner[j * 4 + i] = vec3(INFINITY, INFINITY, INFINITY);
                }
                else
                {
                    float t = dot(position - p, normal) / denom;
                    project_corner[j * 4 + i] = p + dir[j] * t;
                }
            }
        }

        vec2 project_corner_2D[8];
        for (int i = 0; i < 8; ++i)
        {
            project_corner_2D[i] = world_to_buffer(project_corner[i]);
        }

        std::vector<vec2> poly = get_polygon(project_corner_2D);

        if (poly.size() == 6)
        {
            rasterize_triangle(poly[0], poly[1], poly[3], edge_id, buffer);
            rasterize_triangle(poly[0], poly[3], poly[4], edge_id, buffer);
            rasterize_triangle(poly[1], poly[2], poly[3], edge_id, buffer);
            rasterize_triangle(poly[0], poly[4], poly[5], edge_id, buffer);
        }
        else if (poly.size() == 4)
        {
            rasterize_triangle(poly[0], poly[1], poly[2], edge_id, buffer);
            rasterize_triangle(poly[0], poly[2], poly[3], edge_id, buffer);
        }
        else
        {
            cout << "abnormal polygon vertice: " << poly.size() << endl;
            // while (poly.size() > 2)
            // {
            //     rasterize_triangle(poly[0], poly[1], poly.back(), edge_id, buffer);
            //     poly.pop_back();
            // }
        }
    }
    void printBuffer() const
    {
        // return;
        for (int y = 0; y < resolution; ++y)
        {
            cout << "line: " << y;
            for (int x = 0; x < resolution; ++x)
            {
                const auto &cell = buffer[y][x];
                if (!cell.empty())
                {
                    std::cout << "(" << x << "," << y << "):";
                    for (int id : cell)
                    {
                        std::cout << id << ",";
                    }
                    std::cout << "\t";
                }
                else
                {
                    cout << 0;
                }
            }
            std::cout << "\n";
        }
    }
};
#endif