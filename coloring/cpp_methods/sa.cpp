#include <cstdio>
#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <set>
#include <numeric>

using namespace std;

mt19937 rng(4257179);
uniform_real_distribution<double> dist(0.0, 1.0);

double prob() {
    return dist(rng);
}


vector<int> SimulatedAnnealing(int n, int m, const vector<pair<int, int>>& edges) {
    vector<vector<int>> graph(n);
    for (int i = 0; i < m; ++i) {
        auto [u, v] = edges[i];
        graph[u].push_back(v);
        graph[v].push_back(u);
    }
    vector<int> best_coloring(n, -1);
    int best_val = n + 1;

    vector<int> used(n + 2);
    
    auto CostOrder = [&](const vector<int>& order) {
        vector<int> cur_coloring(n, -1);
        int cur_val = 0;

        for (int v : order) {
            for (int adj : graph[v]) {
                if (cur_coloring[adj] != -1) {
                    used[cur_coloring[adj]] = 1;
                }
            }
            int c = 1;
            while (used[c]) {
                ++c;
            }
            cur_val = max(cur_val, c);
            cur_coloring[v] = c;

             for (int adj : graph[v]) {
                if (cur_coloring[adj] != -1) {
                    used[cur_coloring[adj]] = 0;
                }
            }
        }

        if (cur_val < best_val) {
            best_val = cur_val;
            best_coloring = cur_coloring;
        }
        
        return cur_val;
    };

    const int restart = 1000 / n;
    const int iter = 9e8 / (restart * m);
    
    for (int epoch = 0; epoch < restart; ++epoch) {
        vector<int> order(n);
        iota(order.begin(), order.end(), 0);
        shuffle(order.begin(), order.end(), rng);

        int cur_val = CostOrder(order);
        double temperature = 300;
        double step = 0.9;
        
        for (int it = 0; it < iter; ++it) {
            int snd = rng() % (n - 1) + 1;
            int fst = rng() % snd;
            swap(order[fst], order[snd]);
            int new_val = CostOrder(order);
            if (new_val < cur_val || prob() < exp((cur_val - new_val) / temperature)) {
                cur_val = new_val;
            } else {
                swap(order[fst], order[snd]);
            }
            temperature *= step;
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

    vector<int> coloring = SimulatedAnnealing(n, m, edges);
    for (int c : coloring) {
        cout << c << ' ';
    }
    cout << endl;
    
    fclose(stdin);
    fclose(stdout);
    return 0;
}
