#ifndef EDGEH
#define EDGEH
using namespace std;
class Edge
{
public:
    vec3 v1, v2;
    int faceIndex1 = -1;
    int faceIndex2 = -1;
    Edge()
        : v1(0, 0, 0),
          v2(0, 0, 0),
          faceIndex1(-1),
          faceIndex2(-1)
    {
    }
    Edge(const vec3 &a, const vec3 &b, int faceIndex)
        : v1(min(a, b)), v2(max(a, b)), faceIndex1(faceIndex) {}

    // 排除 faceIndex1, faceIndex2 比較
    bool operator<(const Edge &other) const
    {
        if (v1 != other.v1)
            return v1 < other.v1;
        return v2 < other.v2;
    }

    void register_edge(map<Edge, pair<int, int>> &edgeMap)
    {
        Edge key(v1, v2, -1); // 只保留頂點資訊，避免 faceIndex 影響比較
        auto it = edgeMap.find(key);
        if (it == edgeMap.end())
        {
            edgeMap[key] = {faceIndex1, -1};
        }
        else
        {
            it->second.second = faceIndex1;
        }
    }
};
#endif