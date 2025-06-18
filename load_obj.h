// load obj
bool load_obj_with_filter(const string &filename, Scene &scene, const vec3 &kd = vec3(0.7f, 0.7f, 0.7f),
                          float scale = 1.0f, const vec3 &translation = vec3(0, 0, 0), const vec3 &rotation = vec3(0, 0, 0),
                          bool filter_large_faces = true, float max_face_size = 10.0f)
{
    ifstream file(filename);
    if (!file)
    {
        cerr << "警告：無法開啟 OBJ 檔案：" << filename << endl;
        return false;
    }

    vector<vec3> vertices;
    string line;

    vec3 min_bound(FLT_MAX, FLT_MAX, FLT_MAX);
    vec3 max_bound(-FLT_MAX, -FLT_MAX, -FLT_MAX);

    // 第一遍：收集頂點
    while (getline(file, line))
    {
        istringstream iss(line);
        string type;
        iss >> type;

        if (type == "v")
        {
            float x, y, z;
            iss >> x >> y >> z;
            vec3 vertex(x, y, z);
            vertices.push_back(vertex);

            min_bound = vec3(min(min_bound.x(), x), min(min_bound.y(), y), min(min_bound.z(), z));
            max_bound = vec3(max(max_bound.x(), x), max(max_bound.y(), y), max(max_bound.z(), z));
        }
    }

    if (vertices.empty())
    {
        cerr << "警告：在 OBJ 檔案中找不到頂點：" << filename << endl;
        return false;
    }

    // 輸出模型資訊
    vec3 center = (min_bound + max_bound) * 0.5f;
    vec3 size = max_bound - min_bound;
    float max_size = max(max(size.x(), size.y()), size.z());

    cout << "\n=== OBJ 檔案分析 ===" << endl;
    cout << "檔案：" << filename << endl;
    cout << "頂點數量：" << vertices.size() << endl;
    cout << "邊界框：(" << min_bound.x() << "," << min_bound.y() << "," << min_bound.z() << ") 到 ("
         << max_bound.x() << "," << max_bound.y() << "," << max_bound.z() << ")" << endl;
    cout << "尺寸：" << size.x() << " x " << size.y() << " x " << size.z() << endl;
    cout << "最大維度：" << max_size << endl;

    // 變換頂點
    for (auto &vertex : vertices)
    {
        vertex = vertex - center;
        vertex = vertex * scale;

        // 旋轉變換
        if (rotation.x() != 0)
        {
            float cos_x = cos(rotation.x() * PI / 180.0f);
            float sin_x = sin(rotation.x() * PI / 180.0f);
            float y = vertex.y() * cos_x - vertex.z() * sin_x;
            float z = vertex.y() * sin_x + vertex.z() * cos_x;
            vertex = vec3(vertex.x(), y, z);
        }
        if (rotation.y() != 0)
        {
            float cos_y = cos(rotation.y() * PI / 180.0f);
            float sin_y = sin(rotation.y() * PI / 180.0f);
            float x = vertex.x() * cos_y + vertex.z() * sin_y;
            float z = -vertex.x() * sin_y + vertex.z() * cos_y;
            vertex = vec3(x, vertex.y(), z);
        }
        if (rotation.z() != 0)
        {
            float cos_z = cos(rotation.z() * PI / 180.0f);
            float sin_z = sin(rotation.z() * PI / 180.0f);
            float x = vertex.x() * cos_z - vertex.y() * sin_z;
            float y = vertex.x() * sin_z + vertex.y() * cos_z;
            vertex = vec3(x, y, vertex.z());
        }

        vertex = vertex + translation;
    }

    // 第二遍：處理面，並進行過濾
    file.clear();
    file.seekg(0, ios::beg);

    int total_faces = 0;
    int added_triangles = 0;
    int filtered_faces = 0;

    while (getline(file, line))
    {
        istringstream iss(line);
        string type;
        iss >> type;

        if (type == "f")
        {
            total_faces++;
            vector<int> idx;
            string vertex_str;
            while (iss >> vertex_str)
            {
                istringstream vss(vertex_str);
                string i_str;
                getline(vss, i_str, '/');
                int i = stoi(i_str) - 1;
                if (i >= 0 && i < vertices.size())
                {
                    idx.push_back(i);
                }
            }

            if (idx.size() >= 3)
            {
                // 檢查面的大小，過濾掉過大的面（可能是背景）
                bool should_add = true;

                if (filter_large_faces)
                {
                    for (size_t i = 1; i + 1 < idx.size(); ++i)
                    {
                        vec3 v0 = vertices[idx[0]];
                        vec3 v1 = vertices[idx[i]];
                        vec3 v2 = vertices[idx[i + 1]];

                        // 計算三角形的邊長
                        float edge1_len = (v1 - v0).length();
                        float edge2_len = (v2 - v0).length();
                        float edge3_len = (v2 - v1).length();

                        float max_edge = max(max(edge1_len, edge2_len), edge3_len);

                        // 如果任何一條邊太長，跳過這個三角形
                        if (max_edge > max_face_size)
                        {
                            should_add = false;
                            filtered_faces++;
                            cout << "過濾大面：最大邊長 = " << max_edge << endl;
                            break;
                        }

                        // 檢查三角形是否接近平行於某個坐標軸平面（可能是背景牆）
                        vec3 normal = cross(v1 - v0, v2 - v0);
                        if (normal.length() > EPSILON)
                        {
                            normal = unit_vector(normal);

                            // 如果法向量幾乎平行於X、Y或Z軸，且三角形很大，可能是背景
                            if ((abs(abs(normal.x()) - 1.0f) < 0.1f ||
                                 abs(abs(normal.y()) - 1.0f) < 0.1f ||
                                 abs(abs(normal.z()) - 1.0f) < 0.1f) &&
                                max_edge > 2.0f)
                            {
                                should_add = false;
                                filtered_faces++;
                                cout << "過濾背景面：法向量 = (" << normal.x() << "," << normal.y() << "," << normal.z()
                                     << "), 最大邊長 = " << max_edge << endl;
                                break;
                            }
                        }
                    }
                }

                if (should_add)
                {
                    // 三角化並加入場景
                    for (size_t i = 1; i + 1 < idx.size(); ++i)
                    {
                        vec3 v0 = vertices[idx[0]];
                        vec3 v1 = vertices[idx[i]];
                        vec3 v2 = vertices[idx[i + 1]];

                        vec3 edge1 = v1 - v0;
                        vec3 edge2 = v2 - v0;
                        vec3 cross_product = cross(edge1, edge2);

                        if (cross_product.length() > EPSILON)
                        {
                            scene.add_object(make_unique<triangle>(v0, v1, v2, kd));
                            added_triangles++;
                        }
                    }
                }
            }
        }
    }
    cout << "=== 載入結果 ===" << endl;
    cout << "總面數：" << total_faces << endl;
    cout << "過濾掉的面：" << filtered_faces << endl;
    cout << "成功載入的三角形：" << added_triangles << endl;
    cout << "==================" << endl;
    return true;
}