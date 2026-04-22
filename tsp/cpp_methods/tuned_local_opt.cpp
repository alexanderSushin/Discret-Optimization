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

struct PointDist {
    vector<pair<double, double>> pts;
    vector<vector<double>> dist_matrix;

    PointDist(const vector<pair<double, double>>& pts) : pts(pts) {
        if (pts.size() <= 5000) {
            int n = pts.size();
            dist_matrix.resize(n, vector<double>(n));
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    dist_matrix[i][j] = Dist(pts[i], pts[j]);
                }
            }
        }
    }

    PointDist() = default;

    double GetDist(int i, int j) {
        if (pts.size() != dist_matrix.size()) {
            return Dist(pts[i], pts[j]);
        }
        return dist_matrix[i][j];
    }
};

PointDist dists;

const int K = 20;
const double INF = 1e18;

vector<int> Greedy(int n) {
    vector<int> order = {0};
    vector<bool> used(n, false);
    used[0] = true;
    while (order.size() != n) {
        int end = order.back();
        int opt_idx = -1;
        double opt_val;
        for (int j = 0; j < n; ++j) {
            if (used[j]) continue;
            double d = dists.GetDist(end, j);
            
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

void LocalOpt(int n, vector<int>& order) {
    double path_sm = 0;
    for (int i = 0; i < n; ++i) {
        int nxt = (i + 1 == n ? 0 : i + 1);
        path_sm += dists.GetDist(order[i], order[nxt]);
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
                
                double was_val = dists.GetDist(order[s], order[prv(s)]) + 
                                 dists.GetDist(order[t], order[nxt(t)]);
                
                double next_val = dists.GetDist(order[s], order[nxt(t)]) + 
                                  dists.GetDist(order[t], order[prv(s)]);
    
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


void double_bridge(vector<int>& order) {
    int n = order.size();

    int a = rng() % n;
    int b = rng() % n;
    int c = rng() % n;
    int d = rng() % n;

    vector<int> cuts = {a, b, c, d};
    sort(cuts.begin(), cuts.end());

    a = cuts[0];
    b = cuts[1];
    c = cuts[2];
    d = cuts[3];

    vector<int> new_order;
    
    new_order.insert(new_order.end(), order.begin(), order.begin() + a);
    new_order.insert(new_order.end(), order.begin() + c, order.begin() + d);
    new_order.insert(new_order.end(), order.begin() + b, order.begin() + c);
    new_order.insert(new_order.end(), order.begin() + a, order.begin() + b);
    new_order.insert(new_order.end(), order.begin() + d, order.end());

    order = new_order;
}

vector<int> Annealing(int n) {
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

    auto relax = [&](double path_sm, const vector<int>& order) -> double {
        if (best > path_sm) {
            double diff = best - path_sm;
            best = path_sm;
            best_order = order;

            return diff;
        }
        return 0;
    };

    auto cost = [&](const vector<int>& order) {
        double path_sm = 0;
        for (int i = 0; i < n; ++i) {
            int nxt = (i + 1 == n ? 0 : i + 1);
            path_sm += dists.GetDist(order[i], order[nxt]);
        }
        return path_sm;
    };

    auto init = Greedy(n);
    LocalOpt(n, init);
    relax(cost(init), init);
    
    auto start_time = chrono::steady_clock::now();
    while (1) {
        auto now = chrono::steady_clock::now();
        auto elapsed = chrono::duration_cast<chrono::seconds>(now - start_time).count();
        if (elapsed > 250) break;
        
        vector<int> order = best_order;
        double_bridge(order);
        LocalOpt(n, order);
        relax(cost(order), order);
    }

    LocalOpt(n, best_order);
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

    dists = PointDist(points);

    vector<int> ordering = Annealing(n);
    for (int i : ordering) {
        cout << i << ' ';
    }
    cout << endl;
    
    fclose(stdin);
    fclose(stdout);
    return 0;
}
