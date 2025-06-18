// util
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <float.h>
#include <memory>
#include <sstream>
#include <random>
#include <algorithm>
#include <chrono>

// classes
#include "vec3.h"
#include "hitable.h"
#include "scene.h"

using namespace std;

// constants
const float PI = 3.14159265359f;
const int MAX_DEPTH = 5;
const int SHADOW_SAMPLES = 200;
const int IMAGE_WEIGHT = 1600;
const int IMAGE_HEIGHT = 800;
// const int SHADOW_SAMPLES = 200;
long long SHADOW_RENDERING_TIME = 0;
long long TRACING_TIME = 0;
long long HEMICUBE_CONSTRUCTION_TIME = 0;
long long QUERY_TIME = 0;
long long UPDATE_TIME = 0;
long long INTEGRATE_TIME = 0;
bool USE_SHADOW_VOLUME = true;
bool USE_SIMPLIFY = false;
bool SHOW_PENUMBRA = false;
const char *FILE_NAME = "volumn_complex_200_h64l64.ppm";

// Random number generator
random_device rd;
mt19937 gen(rd());
uniform_real_distribution<float> dis(0.0f, 1.0f);
float random_float()
{
    return dis(gen);
}
float random_float(float min, float max)
{
    return min + (max - min) * random_float();
}

// rendering functions
vec3 reflect(const vec3 &v, const vec3 &n)
{
    return v - 2 * dot(v, n) * n;
}
bool refract(const vec3 &v, const vec3 &n, float ni_over_nt, vec3 &refracted)
{
    vec3 uv = unit_vector(v);
    float dt = dot(uv, n);
    float discriminant = 1.0f - ni_over_nt * ni_over_nt * (1 - dt * dt);
    if (discriminant > 0)
    {
        refracted = ni_over_nt * (uv - n * dt) - n * sqrt(discriminant);
        return true;
    }
    return false;
}
vec3 shading_single(const hit_record &rec, Scene &scene, Hemicube &hemicube)
{
    vec3 color(0, 0, 0);
    int hits = 0;

    vec3 light_point = scene.lightSource.light_position;
    vec3 light_dir = unit_vector(light_point - rec.p);

    float bias = max(0.01f, abs(rec.t) * 1e-3f);
    vec3 shadow_origin = rec.p + bias * rec.normal;
    ray shadow_ray(shadow_origin, light_dir);

    hit_record temp;
    float light_distance = (light_point - rec.p).length();
    if (!scene.hit(shadow_ray, bias, light_distance - bias, temp))
    {
        float diff = max(0.0f, dot(rec.normal, light_dir));
        color += diff * scene.lightSource.light_intensity * rec.Kd;
        hits++;
    }

    return color;
}
vec3 shading_stochastic(const hit_record &rec, Scene &scene, Hemicube &hemicube)
{
    vec3 color(0, 0, 0);
    int hits = 0;
    for (int i = 0; i < SHADOW_SAMPLES; ++i)
    {

        float u = random_float(-0.5f, 0.5f);
        float v = random_float(-0.5f, 0.5f);
        vec3 light_point = scene.lightSource.light_position +
                           u * scene.lightSource.light_u + v * scene.lightSource.light_v;
        vec3 light_dir = unit_vector(light_point - rec.p);

        // 增加陰影光線的偏差值
        float bias = max(0.01f, abs(rec.t) * 1e-3f); // 增加偏差值
        vec3 shadow_origin = rec.p + bias * rec.normal;
        ray shadow_ray(shadow_origin, light_dir);

        hit_record temp;
        float light_distance = (light_point - rec.p).length();

        if (!scene.hit(shadow_ray, bias, light_distance - bias, temp))
        {
            float diff = max(0.0f, dot(rec.normal, light_dir));
            color += diff * scene.lightSource.light_intensity * rec.Kd;
            hits++;
        }
    }
    color /= float(SHADOW_SAMPLES);
    return color;
}
vec3 shading_shadow_volumn(const hit_record &rec, Scene &scene, Hemicube &hemicube, const QueryResponse &lowest_complexity)
{
    vec3 color(0, 0, 0);
    int hits = 0;
    vec3 light_point = lowest_complexity.pos;

    float percentage = lowest_complexity.percentage;
    vec3 light_dir = unit_vector(light_point - rec.p);

    float bias = max(0.01f, abs(rec.t) * 1e-3f);
    vec3 shadow_origin = rec.p + bias * rec.normal;
    ray shadow_ray(shadow_origin, light_dir);

    hit_record temp;
    float light_distance = (light_point - rec.p).length();

    if (!scene.hit(shadow_ray, bias, light_distance - bias, temp))
    {
        float diff = max(0.0f, dot(rec.normal, light_dir));
        color += diff * scene.lightSource.light_intensity * rec.Kd;
    }

    color *= percentage;
    return color;
}
vec3 shadowing(const hit_record &rec, Scene &scene, Hemicube &hemicube)
{
    auto shadow_start = chrono::high_resolution_clock::now();
    vec3 color;
    if (USE_SHADOW_VOLUME)
    {
        if (USE_SIMPLIFY)
        {
            if (hemicube.simple_query(rec.p, scene.lightSource, scene.edgeMap, scene.edgeIdMap, scene.objects))
                color = shading_stochastic(rec, scene, hemicube);
            else
                color = shading_single(rec, scene, hemicube);
        }
        else
        { // shadow volume
            auto query_start = chrono::high_resolution_clock::now();
            QueryResponse q = hemicube.query(rec.p, scene.lightSource, scene.edgeMap, scene.edgeIdMap, scene.objects, UPDATE_TIME);
            INTEGRATE_TIME += q.time;
            auto query_end = chrono::high_resolution_clock::now();
            QUERY_TIME += chrono::duration_cast<chrono::nanoseconds>(query_end - query_start).count();
            color = shading_shadow_volumn(rec, scene, hemicube, q);
        }
    }
    else // stochastic
        color = shading_stochastic(rec, scene, hemicube);

    color += 0.1f * rec.Kd; // ambient
    auto shadow_end = chrono::high_resolution_clock::now();
    SHADOW_RENDERING_TIME += chrono::duration_cast<chrono::nanoseconds>(shadow_end - shadow_start).count();

    return color;
}
vec3 skybox(const ray &r)
{
    vec3 unit_direction = unit_vector(r.direction());
    float t = 0.5 * (unit_direction.y() + 1.0);
    return (1.0 - t) * vec3(1.0, 1.0, 1.0) + t * vec3(0.5, 0.7, 1.0); // gradient sky
}
vec3 trace(const ray &r, Scene &scene, int depth, Hemicube &hemicube)
{
    auto trace_start = chrono::high_resolution_clock::now();
    hit_record rec;
    if (scene.hit(r, 0.001f, FLT_MAX, rec))
    {
        vec3 local = shadowing(rec, scene, hemicube);

        // 如果達到最大深度，只返回局部著色
        if (depth >= MAX_DEPTH)
        {
            auto trace_end = chrono::high_resolution_clock::now();
            TRACING_TIME += chrono::duration_cast<chrono::nanoseconds>(trace_end - trace_start).count();
            return local;
        }

        vec3 reflected_color(0, 0, 0), transmitted_color(0, 0, 0);

        if (rec.w_r > 0.0f && depth < MAX_DEPTH)
        {
            vec3 reflected = reflect(r.direction(), rec.normal);
            float bias = max(0.001f, abs(rec.t) * 1e-4f);
            reflected_color = trace(ray(rec.p + bias * rec.normal, reflected), scene, depth + 1, hemicube);
        }

        if (rec.w_t > 0.0f && depth < MAX_DEPTH)
        {
            vec3 refracted;
            bool entering = dot(r.direction(), rec.normal) < 0;
            float ni_over_nt = entering ? 1.0f / rec.n_t : rec.n_t;
            vec3 n = entering ? rec.normal : -rec.normal;

            if (refract(r.direction(), n, ni_over_nt, refracted))
            {
                float bias = max(0.001f, abs(rec.t) * 1e-4f);
                vec3 refract_origin = entering ? rec.p - bias * n : rec.p + bias * n;
                transmitted_color = trace(ray(refract_origin, refracted), scene, depth + 1, hemicube);
            }
            else
            {
                // 全內反射時使用天空盒顏色
                transmitted_color = skybox(r);
            }
        }
        auto trace_end = chrono::high_resolution_clock::now();
        TRACING_TIME += chrono::duration_cast<chrono::nanoseconds>(trace_end - trace_start).count();
        return (1.0f - rec.w_t) * ((1.0f - rec.w_r) * local + rec.w_r * reflected_color) + rec.w_t * transmitted_color;
    }
    auto trace_end = chrono::high_resolution_clock::now();
    TRACING_TIME += chrono::duration_cast<chrono::nanoseconds>(trace_end - trace_start).count();
    return skybox(r); // 沒有擊中物體時返回天空盒顏色
}

int main()
{
    auto start = chrono::high_resolution_clock::now();
    // Image settings
    const int nx = IMAGE_WEIGHT;
    const int ny = IMAGE_HEIGHT;
    // const int nx = 400;
    // const int ny = 200;
    cout << "Image width: " << nx << ", height: " << ny << endl;

    // Camera settings - adjusted for better view
    vec3 lookfrom(0, 2, 5);
    vec3 lookat(0, 1, -1);
    vec3 vup(0, 1, 0);
    float vfov = 60;
    float aspect = float(nx) / float(ny);

    float theta = vfov * PI / 180.0f;
    float half_height = tan(theta / 2);
    float half_width = aspect * half_height;

    vec3 w = unit_vector(lookfrom - lookat);
    vec3 u = unit_vector(cross(vup, w));
    vec3 v = cross(w, u);

    vec3 lower_left_corner = lookfrom - half_width * u - half_height * v - w;
    vec3 lower_right_corner = lookfrom + half_width * u - half_height * v - w;
    vec3 top_left_corner = lookfrom - half_width * u + half_height * v - w;
    vec3 top_right_corner = lookfrom + half_width * u + half_height * v - w;
    cout << "x: " << lower_left_corner[0] << "y: " << lower_left_corner[1] << "z: " << lower_left_corner[2] << endl;
    cout << "x: " << lower_right_corner[0] << "y: " << lower_right_corner[1] << "z: " << lower_right_corner[2] << endl;
    cout << "x: " << top_left_corner[0] << "y: " << top_left_corner[1] << "z: " << top_left_corner[2] << endl;
    cout << "x: " << top_right_corner[0] << "y: " << top_right_corner[1] << "z: " << top_right_corner[2] << endl;
    vec3 horizontal = 2 * half_width * u;
    vec3 vertical = 2 * half_height * v;
    auto ckp_setting = chrono::high_resolution_clock::now();
    chrono::duration<double, nano> duration = ckp_setting - start;
    cout << "Setting time:" << chrono::duration_cast<chrono::nanoseconds>(duration).count() << " ns" << endl;

    // Create scene
    Scene scene(vec3(10, 20, 20), vec3(10.0f, -5.0, 0), vec3(0, -5.0, 10.0f), vec3(1, 1, 1)); // correct scene
    scene.add_object(make_unique<plane>(vec3(0, 0, 0), vec3(0, 1, 0), vec3(0.8f, 0.8f, 0.0f)));
    // 調整物體位置讓它們更容易看到，並確保沒有反射/折射屬性
    add_cube(scene, vec3(-2.0f, 2.5f, -1), 1.0f, vec3(0.8f, 0.3f, 0.3f)); // 紅色立方體correct scene
    add_pyramid(
        scene,
        vec3(-3.0f, 2.0f, 0.0f), // 位置稍微往左
        0.5f,                    // 邊長 0.5
        vec3(0.9f, 0.6f, 0.1f)   // 黃色系
    );
    add_tetrahedron(scene, vec3(0.0f, 2.0f, -0.2), 0.4f, vec3(0.2f, 0.8f, 0.3f)); // 綠色四面體
    add_cube(
        scene,
        vec3(3.0f, 2.5f, 0.5f), // 位置稍微往右、稍高
        0.7f,                   // 邊長 0.7
        vec3(0.3f, 0.7f, 0.9f)  // 青藍色系
    );
    add_pyramid(scene, vec3(2.0f, 1.5f, -0.5), 1.0f, vec3(0.2f, 0.4f, 0.9f)); // 藍色金字塔
    add_tetrahedron(
        scene,
        vec3(0.5f, 1.0f, -0.5f), // 位置往上、往後
        0.4f,                    // 邊長 0.4
        vec3(0.8f, 0.2f, 0.7f)   // 粉紫色系
    );

    //  Open output file
    auto ckp_scene = chrono::high_resolution_clock::now();
    duration = ckp_scene - ckp_setting;
    cout << "Costruction scene time:" << chrono::duration_cast<chrono::nanoseconds>(duration).count() << " ns" << endl;
    int i, j = 0;

    // constructing hemicube
    vec3 border[4] = {lower_left_corner, lower_right_corner, top_left_corner, top_right_corner};
    Hemicube hemicube(border, scene.lightSource, vec3(0, 0, 0), vec3(0, 1, 0));
    for (map<Edge, pair<int, int>>::iterator it = scene.edgeMap.begin(); it != scene.edgeMap.end(); ++it, ++i)
    {
        Edge e = it->first;
        pair<int, int> faces = it->second;

        hitable *hit_ptr1 = scene.objects[faces.first].get();
        triangle *tri_ptr1 = dynamic_cast<triangle *>(hit_ptr1);
        hitable *hit_ptr2 = scene.objects[faces.second].get();
        triangle *tri_ptr2 = dynamic_cast<triangle *>(hit_ptr2);
        if (tri_ptr1 && tri_ptr2)
        {

            if (tri_ptr1->normal != tri_ptr2->normal)
            {
                j++;
                scene.edgeIdMap.emplace(j, e);
                if (hemicube.is_silhouette(scene.lightSource, e.v1, e.v2, tri_ptr1->normal, tri_ptr2->normal))
                {
                    hemicube.project_and_store_edge(scene.lightSource, e.v1, e.v2, j);
                }
            }
        }
    }
    for (int i = 0; i < hemicube.resolution; ++i)
    {
        for (int j = 0; j < hemicube.resolution; ++j)
        {
            if (hemicube.buffer[i][j].size() > 0)
            {
                if (SHOW_PENUMBRA)
                {
                    vec3 tmp = hemicube.buffer_to_world(vec2(j, i));
                    add_pyramid(scene, vec3(tmp.x(), 0.25f, tmp.z()), 0.5f, vec3(0.2f, 0.4f, 0.9f));
                }
            }
        }
    }
    auto ckp_hemicube = chrono::high_resolution_clock::now();
    duration = ckp_hemicube - ckp_scene;
    cout << "Costruction hemicube time:" << chrono::duration_cast<chrono::nanoseconds>(duration).count() << " ns" << endl;
    ofstream file(FILE_NAME);
    if (!file)
    {
        cerr << "Error: Cannot create output file ray.ppm" << endl;
        return 1;
    }

    file << "P3\n"
         << nx << " " << ny << "\n255\n";

    // Render
    long long total_time = 0;
    long long last_100_time = 0;
    cout << "Rendering..." << endl;
    SHADOW_RENDERING_TIME = 0;
    TRACING_TIME = 0;
    for (int j = ny - 1; j >= 0; --j)
    {
        auto render_start = chrono::high_resolution_clock::now();
        if (j % 100 == 0)
        {
            cout << "Scanlines remaining: " << j << endl;
            cout << "Total render time per line: " << total_time / 1000000 << " ms" << endl;
            cout << "Avg render time per line(last 100): " << last_100_time / 100 / 1000000 << " ms" << endl;

            last_100_time = 0;
        }

        for (int i = 0; i < nx; ++i)
        {
            float u_coord = float(i) / float(nx);
            float v_coord = float(j) / float(ny);

            // Fixed ray direction calculation
            vec3 direction = lower_left_corner + u_coord * horizontal + v_coord * vertical - lookfrom;
            ray r(lookfrom, unit_vector(direction));
            vec3 col = trace(r, scene, 0, hemicube);

            // Clamp color values
            int ir = int(255.99f * min(max(col.r(), 0.0f), 1.0f));
            int ig = int(255.99f * min(max(col.g(), 0.0f), 1.0f));
            int ib = int(255.99f * min(max(col.b(), 0.0f), 1.0f));

            file << ir << " " << ig << " " << ib << "\n";
        }
        auto render_end = chrono::high_resolution_clock::now();
        total_time += chrono::duration_cast<chrono::nanoseconds>(render_end - render_start).count();
        last_100_time += chrono::duration_cast<chrono::nanoseconds>(render_end - render_start).count();
    }
    cout << "Toral shadow rendering time: " << SHADOW_RENDERING_TIME / 1000000 << " ms" << endl;
    cout << "Toral tracing time: " << TRACING_TIME / 1000000 << " ms" << endl;
    cout << "Toral query time: " << QUERY_TIME / 1000000 << " ms" << endl;
    cout << "Toral update time: " << UPDATE_TIME / 1000000 << " ms" << endl;
    cout << "Toral integrate time: " << INTEGRATE_TIME / 1000000 << " ms" << endl;
    cout << "Avg render time per line: " << total_time / ny / 1000000 << " ms" << endl;
    file.close();
    cout << "Rendering complete! Output saved to ray.ppm" << endl;

    return 0;
}