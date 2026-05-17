#include <cstdio>
#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <set>
#include <numeric>

using namespace std;

mt19937 rng(4257179);

vector<int> Dsatur(int n, int m, const vector<pair<int, int>>& edges) {
    vector<vector<int>> graph(n);
    for (int i = 0; i < m; ++i) {
        auto [u, v] = edges[i];
        graph[u].push_back(v);
        graph[v].push_back(u);
    }
    vector<int> coloring(n, -1);
    vector<set<int>> satur(n);

    for (int taken = 0; taken < n; ++taken) {
        int ver = -1;
        pair<int, int> pr = {-1, -1};
        for (int v = 0; v < n; ++v) {
            if (coloring[v] != -1) continue;
            pair<int, int> cur = make_pair(satur[v].size(), graph[v].size());
            if (cur > pr) {
                pr = cur;
                ver = v;
            }
        }
        set<int> bad_colors;
        for (int adj : graph[ver]) {
            if (coloring[adj] != -1) {
                bad_colors.insert(coloring[adj]);
            }
        }
        int c = 1;
        while (bad_colors.find(c) != bad_colors.end()) {
            ++c;
        }
        coloring[ver] = c;
        for (int adj : graph[ver]) {
            satur[adj].insert(c);
        }
    }
    return coloring;
}

int main(int argc, char** argv) {
    freopen(argv[1], "r", stdin);
    freopen("tmp/ans.txt", "w", stdout);

    int n, m;
    cin >> n >> m;
    vector<pair<int, int>> edges(m);
    for (int i = 0; i < m; ++i) {
        int u, v;
        cin >> u >> v;
        edges[i] = make_pair(u, v);
    }

    vector<int> coloring = Dsatur(n, m, edges);
    for (int c : coloring) {
        cout << c << ' ';
    }
    cout << endl;
    
    fclose(stdin);
    fclose(stdout);
    return 0;
}
