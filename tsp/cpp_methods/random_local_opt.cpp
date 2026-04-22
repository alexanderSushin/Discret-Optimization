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


const int RESTART = 10000;
const int K = 20;
const double INF = 1e18;

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
        if (elapsed > 30) break;
        
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


vector<int> Annealing(int n, const vector<pair<double, double>>& points) {
    vector<int> best_order(n);
    iota(best_order.begin(), best_order.end(), 0);
    double best = INF;
    auto prv = [&](int i) {
        if (i == 0) return n - 1;
        return i - 1;
    };

    auto nxt = [&](int i) {
        if (i == n - 1) return 0;
        return i + 1;
    };

    auto relax = [&](double path_sm, const vector<int>& order) {
        if (best > path_sm) {
            best = path_sm;
            best_order = order;
        }
    };

    auto cost = [&](const vector<int>& order) {
        double path_sm = 0;
        for (int i = 0; i < n; ++i) {
            int nxt = (i + 1 == n ? 0 : i + 1);
            path_sm += Dist(points[order[i]], points[order[nxt]]);
        }
        return path_sm;
    };

    auto init = Greedy(n, points);
    LocalOpt(n, points, init);
    relax(cost(init), init);

    auto start_time = chrono::steady_clock::now();
    while (1) {
        auto now = chrono::steady_clock::now();
        auto elapsed = chrono::duration_cast<chrono::seconds>(now - start_time).count();
        if (elapsed > 100) break;
        
        vector<int> order = best_order;
        for (int i = 0; i < K; ++i) {
            int s = rng() % n;
            int t = rng() % n;
            if (s > t) swap(s, t);
            reverse(order.begin() + s, order.begin() + t);
        }
        LocalOpt(n, points, order);
        relax(cost(order), order);
    }

    LocalOpt(n, points, best_order);
    return best_order;
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

    vector<int> ordering = Annealing(n, points);
    for (int i : ordering) {
        cout << i << ' ';
    }
    cout << endl;
    
    fclose(stdin);
    fclose(stdout);
    return 0;
}
