#include <cstdio>
#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <set>
#include <numeric>
#include <cmath>
#include <chrono>

using namespace std;

mt19937 rng(4257179);

double Dist(pair<double, double> a, pair<double, double> b) {
    double dx = a.first - b.first;
    double dy = a.second - b.second;
    return sqrt(dx * dx + dy * dy);
}

vector<int> Greedy(int n, const vector<pair<double, double>>& points) {
    vector<int> order = {0};
    vector<bool> used(n, false);
    used[0] = true;
    while (order.size() != n) {
        int end = order.back();
        int opt_idx = -1;
        double opt_val;
        for (int j = 0; j < n; ++j) {
            if (used[j]) continue;
            double d = Dist(points[end], points[j]);
            
            if (opt_idx == -1 || d < opt_val) {
                opt_idx = j;
                opt_val = d;
            }
        }
        order.emplace_back(opt_idx);
        used[opt_idx] = true;
    }
    return order;
}


void LocalOpt(int n, const vector<pair<double, double>>& points, vector<int>& order) {
    double path_sm = 0;
    for (int i = 0; i < n; ++i) {
        int nxt = (i + 1 == n ? 0 : i + 1);
        path_sm += Dist(points[order[i]], points[order[nxt]]);
    }
    auto prv = [&](int i) {
        if (i == 0) return n - 1;
        return i - 1;
    };

    auto nxt = [&](int i) {
        if (i == n - 1) return 0;
        return i + 1;
    };

    bool found = true;
    auto start_time = chrono::steady_clock::now();
    while (found) {
        auto now = chrono::steady_clock::now();
        auto elapsed = chrono::duration_cast<chrono::seconds>(now - start_time).count();
        if (elapsed >= 60) break;
        
        found = false;
        
        for (int s = 1; s < n; ++s) {
            if (found) break;
            // [... s, ..., t, ...]
            // [..., t, ..., s, ...]
        
            for (int t = s + 1; t < n; ++t) {
                if (found) break;
                
                double was_val = Dist(points[order[s]], points[order[prv(s)]]) + 
                                 Dist(points[order[t]], points[order[nxt(t)]]);
                
                double next_val = Dist(points[order[s]], points[order[nxt(t)]]) + 
                                  Dist(points[order[t]], points[order[prv(s)]]);
    
                if (next_val < was_val) {
                    found = true;
                    reverse(order.begin() + s, order.begin() + t + 1);
                    path_sm += next_val - was_val;
                    break;
                }
            }
        }
    }
}

int main(int argc, char** argv) {
    freopen(argv[1], "r", stdin);
    freopen("tmp/ans.txt", "w", stdout);

    int n;
    cin >> n;
    vector<pair<double, double>> points(n);
    for (int i = 0; i < n; ++i) {
        double x, y;
        cin >> x >> y;
        points[i] = make_pair(x, y);
    }

    vector<int> ordering = Greedy(n, points);
    LocalOpt(n, points, ordering);
    for (int i : ordering) {
        cout << i << ' ';
    }
    cout << endl;
    
    fclose(stdin);
    fclose(stdout);
    return 0;
}
