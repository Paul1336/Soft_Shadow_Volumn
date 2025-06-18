#ifndef SCENEH
#define SCENEH

#include <vector>
#include <memory>
#include "vec3.h"
#include "ray.h"
#include "hitable.h"
#include "edge.h"
#include "hemicube.h"
#include <map>
#include "lightsource.h"

class Scene
{
public:
    vector<unique_ptr<hitable>> objects;
    map<Edge, pair<int, int>> edgeMap;
    std::map<int, Edge> edgeIdMap;
    LightSource lightSource;
    Scene(vec3 position = vec3(10, 40, 20), vec3 u = vec3(10.0f, 0, 0), vec3 v = vec3(0, 0, 10.0f), vec3 intensity = vec3(1, 1, 1)) : lightSource(position, u, v, intensity)
    {
    }
    void add_object(unique_ptr<hitable> obj)
    {
        objects.push_back(move(obj));
    }
    void add_triangle_with_edge(const vec3 &a, const vec3 &b, const vec3 &c, const vec3 &color)
    {
        int faceIndex = objects.size(); // 對應當前 triangle 的 index

        add_object(make_unique<triangle>(a, b, c, color));

        // 建立三條邊並註冊
        Edge(a, b, faceIndex).register_edge(edgeMap);
        Edge(b, c, faceIndex).register_edge(edgeMap);
        Edge(c, a, faceIndex).register_edge(edgeMap);
    }
    bool hit(const ray &r, float t_min, float t_max, hit_record &closest_rec) const
    {
        hit_record temp_rec;
        bool hit_anything = false;
        float closest_so_far = t_max;

        for (const auto &obj : objects)
        {
            if (obj->hit(r, t_min, closest_so_far, temp_rec))
            {
                hit_anything = true;
                closest_so_far = temp_rec.t;
                closest_rec = temp_rec;
            }
        }
        return hit_anything;
    }
};

void add_cube(Scene &scene, vec3 center, float size, vec3 color)
{
    float s = size / 2.0f;
    vec3 v[8] = {
        center + vec3(-s, -s, -s), center + vec3(s, -s, -s),
        center + vec3(s, s, -s), center + vec3(-s, s, -s),
        center + vec3(-s, -s, s), center + vec3(s, -s, s),
        center + vec3(s, s, s), center + vec3(-s, s, s)};

    int faces[12][3] = {
        {0, 1, 2}, {2, 3, 0}, // back
        {4, 5, 6},
        {6, 7, 4}, // front
        {0, 4, 7},
        {7, 3, 0}, // left
        {1, 5, 6},
        {6, 2, 1}, // right
        {3, 2, 6},
        {6, 7, 3}, // top
        {0, 1, 5},
        {5, 4, 0} // bottom
    };

    for (int i = 0; i < 12; ++i)
    {
        scene.add_triangle_with_edge(v[faces[i][0]], v[faces[i][1]], v[faces[i][2]], color);
    }
}

void add_tetrahedron(Scene &scene, vec3 center, float size, vec3 color)
{
    float s = size;
    vec3 v0 = center + vec3(0, s, 0);
    vec3 v1 = center + vec3(-s, -s, s);
    vec3 v2 = center + vec3(s, -s, s);
    vec3 v3 = center + vec3(0, -s, -s);

    scene.add_triangle_with_edge(v0, v1, v2, color);
    scene.add_triangle_with_edge(v0, v2, v3, color);
    scene.add_triangle_with_edge(v0, v3, v1, color);
    scene.add_triangle_with_edge(v1, v3, v2, color);
}

void add_pyramid(Scene &scene, vec3 center, float size, vec3 color)
{
    float s = size / 2.0f;
    vec3 top = center + vec3(0, s, 0);
    vec3 base1 = center + vec3(-s, -s, -s);
    vec3 base2 = center + vec3(s, -s, -s);
    vec3 base3 = center + vec3(s, -s, s);
    vec3 base4 = center + vec3(-s, -s, s);

    // base (split into two triangles)
    scene.add_triangle_with_edge(base1, base2, base3, color);
    scene.add_triangle_with_edge(base3, base4, base1, color);

    // sides
    scene.add_triangle_with_edge(top, base1, base2, color);
    scene.add_triangle_with_edge(top, base2, base3, color);
    scene.add_triangle_with_edge(top, base3, base4, color);
    scene.add_triangle_with_edge(top, base4, base1, color);
}

#endif