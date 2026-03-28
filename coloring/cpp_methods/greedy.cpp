#include <cstdio>
#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <set>
#include <numeric>

using namespace std;

mt19937 rng(4257179);

vector<int> Greedy(int n, int m, const vector<pair<int, int>>& edges) {
    vector<vector<int>> graph(n);
    for (int i = 0; i < m; ++i) {
        auto [u, v] = edges[i];
        graph[u].push_back(v);
        graph[v].push_back(u);
    }
    vector<int> best_coloring(n, -1);
    int best_val = n + 1;
    
    for (int iter = 0; iter < 1000; ++iter) {
        vector<int> order(n);
        iota(order.begin(), order.end(), 0);
        shuffle(order.begin(), order.end(), rng);
        
        vector<int> cur_coloring(n, -1);
        int cur_val = 0;

        for (int v : order) {
            set<int> colors;
            for (int adj : graph[v]) {
                if (cur_coloring[adj] != -1) {
                    colors.insert(cur_coloring[adj]);
                }
            }
            int c = 1;
            while (colors.find(c) != colors.end()) {
                ++c;
            }
            cur_val = max(cur_val, c);
            cur_coloring[v] = c;
            
        }

        if (cur_val < best_val) {
            best_val = cur_val;
            best_coloring = cur_coloring;
        }

    }
    
    return best_coloring;
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

    vector<int> coloring = Greedy(n, m, edges);
    for (int c : coloring) {
        cout << c << ' ';
    }
    cout << endl;
    
    fclose(stdin);
    fclose(stdout);
    return 0;
}
