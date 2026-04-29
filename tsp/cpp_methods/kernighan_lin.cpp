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
const double EPS = 1e-9;

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

void cyclic_reverse(vector<int>& order, int l, int r) {
    int n = order.size();

    while (true) {
        swap(order[l], order[r]);

        if (l == r || (l + 1) % n == r)
            break;

        l = (l + 1) % n;
        r = (r - 1 + n) % n;
    }
}

void KernighanLin(int n, vector<int>& order) {
    auto cost = [&](const vector<int>& cur_order) {
        double path_sm = 0;
        for (int i = 0; i < n; ++i) {
            int nxt = (i + 1 == n ? 0 : i + 1);
            path_sm += dists.GetDist(cur_order[i], cur_order[nxt]);
        }
        return path_sm;
    };

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
    vector<vector<bool>> used(n, vector<bool>(n));

    while (found) {
        found = false;
        double cycle_cost = cost(order);

        vector<int> b_order(n);
        iota(b_order.begin(), b_order.end(), 0);
        shuffle(b_order.begin(), b_order.end(), rng);

        for (int init_b : b_order) {;

            vector<pair<int, int>> deleted;

            auto delete_edge = [&](int x, int y) {
              deleted.push_back({x, y});
              used[x][y] = used[y][x] = true;
            };

            if (found) break;
            vector<int> now_order = order;

            // [..., a - b, ...]
            int pos_b = init_b;
            int pos_a = prv(init_b);

            delete_edge(now_order[pos_a], now_order[pos_b]);

            double ref_cost = cycle_cost - dists.GetDist(now_order[pos_a], now_order[pos_b]);

            while (1) {
                auto now = chrono::steady_clock::now();
                auto elapsed = chrono::duration_cast<chrono::seconds>(now - start_time).count();
                if (elapsed > 30) return;
                int opt_d = -1;
                double opt_val = INF;

                for (int pos_d = 0; pos_d < n; ++pos_d) {
                    // [..., a - b, ... d]
                    if (pos_d == pos_a || pos_d == pos_b) continue;
                    int pos_c = prv(pos_d);
                    if (pos_c == pos_b) continue;

                    if (used[now_order[pos_c]][now_order[pos_d]]) continue;

                    double add_edge = dists.GetDist(now_order[pos_b], now_order[pos_d]);
                    if (add_edge < opt_val) {
                        opt_val = add_edge;
                        opt_d = pos_d;
                    }
                }

                ref_cost += opt_val;

                if (opt_d == -1 || ref_cost >= cycle_cost) {
                    break;
                }
                int pos_d = opt_d;
                int pos_c = prv(pos_d);

                ref_cost -= dists.GetDist(now_order[pos_c], now_order[pos_d]);
                delete_edge(now_order[pos_c], now_order[pos_d]);

                double estimation = ref_cost + dists.GetDist(now_order[pos_a], now_order[pos_c]);
                if (estimation < cycle_cost - EPS) {
                    found = true;
                }

                cyclic_reverse(now_order, pos_b, pos_c);
                if (found) {
                    order = now_order;
                    break;
                }

                // [b ... c] -> [c ... b]
                // reverse(pos_b, pos_c);
            }


            for (auto [x, y] : deleted) {
                used[x][y] = used[y][x] = false;
            }
        }


    }
}


vector<int> LocalSearch(int n) {
    vector<int> best_order(n);
    iota(best_order.begin(), best_order.end(), 0);
    double best = INF;

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
    auto start_time = chrono::steady_clock::now();
    while (1) {
        auto now = chrono::steady_clock::now();
        auto elapsed = chrono::duration_cast<chrono::seconds>(now - start_time).count();
        if (elapsed > 600) break;

        KernighanLin(n, init);
        relax(cost(init), init);
        shuffle(init.begin(), init.end(), rng);
    }

    // cout << cost(best_order) << endl;
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

    vector<int> ordering = LocalSearch(n);
    for (int i : ordering) {
        cout << i << ' ';
    }
    cout << endl;
    
    fclose(stdin);
    fclose(stdout);
    return 0;
}
