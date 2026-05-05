#include <cstdio>
#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <set>
#include <numeric>
#include <cmath>
#include <chrono>
#include <cassert>

using namespace std;

mt19937 rng(4257179);
const double EPS = 1e-9;

double Dist(pair<double, double> a, pair<double, double> b) {
    double dx = a.first - b.first;
    double dy = a.second - b.second;
    return sqrt(dx * dx + dy * dy);
}

struct PointDist {
    vector<pair<double, double> > pts;
    vector<vector<double> > dist_matrix;

    PointDist(const vector<pair<double, double> > &pts) : pts(pts) {
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

    int GetPointsSize() {
        return pts.size();
    }
};

PointDist dists;

const double INF = 1e18;

double Cost(const vector<int> &cur_order) {
    int n = cur_order.size();
    double path_sm = 0;
    for (int i = 0; i < n; ++i) {
        int nxt = (i + 1 == n ? 0 : i + 1);
        path_sm += dists.GetDist(cur_order[i], cur_order[nxt]);
    }
    return path_sm;
}


struct Cycle {
    vector<int> order;
    double taken;

    Cycle(vector<int> order, double taken) : order(order), taken(taken) {}

    int Size() {
        return order.size();
    }

    int GetPoint(int i) {
        return order[i];
    }

    double GetTaken() {
        return taken;
    }
};

vector<vector<int>> Greedy(int n, int v, double cap, const vector<double> &req) {
    vector<Cycle> cycles;
    for (int i = 1; i <= n; ++i) {
        cycles.emplace_back(vector<int>{i}, req[i]);
    }
    bool found = true;
    while (cycles.size() > v || found) {
        found = false;
        int cycle_id1, cycle_id2;
        int cycle_pos1, cycle_pos2;
        double best_add = INF;
        for (int i = 0; i < cycles.size(); ++i) {
            for (int j = i + 1; j < cycles.size(); ++j) {
                if (cycles[i].GetTaken() + cycles[j].GetTaken() > cap) continue;

                for (int p1 : {0, cycles[i].Size() - 1}) {
                    for (int p2 : {0, cycles[j].Size() - 1}) {
                        double add = 0;
                        add -= dists.GetDist(cycles[i].GetPoint(p1), 0);
                        add -= dists.GetDist(cycles[j].GetPoint(p2), 0);
                        add += dists.GetDist(cycles[i].GetPoint(p1), cycles[j].GetPoint(p2));
                        if (add < best_add) {
                            best_add = add;
                            cycle_id1 = i;
                            cycle_id2 = j;
                            cycle_pos1 = p1;
                            cycle_pos2 = p2;
                        }
                    }
                }
            }
        }

        if (best_add < 0) {
            found = true;
        }

        if (best_add == INF) {
            break;
        }

        Cycle lhs = cycles[cycle_id1];
        Cycle rhs = cycles[cycle_id2];
        cycles.erase(cycles.begin() + cycle_id2);
        cycles.erase(cycles.begin() + cycle_id1);

        if (cycle_pos1 + 1 != lhs.Size()) {
            reverse(lhs.order.begin(), lhs.order.end());
        }
        if (cycle_pos2 != 0) {
            reverse(rhs.order.begin(), rhs.order.end());
        }
        lhs.order.insert(lhs.order.end(), rhs.order.begin(), rhs.order.end());
        lhs.taken += rhs.taken;
        cycles.push_back(lhs);
    }

    vector<vector<int>> decomp;
    for (int i = 0; i < cycles.size(); ++i) {
        decomp.push_back(cycles[i].order);
        decomp.back().push_back(0);
    }
    while (decomp.size() < v) {
        decomp.push_back({0});
    }

    return decomp;
}



int main(int argc, char **argv) {
    freopen(argv[1], "r", stdin);
    freopen("tmp/ans.txt", "w", stdout);

    int n, v;
    cin >> n >> v;
    --n;
    double cap;
    cin >> cap;

    vector<double> req(n + 1);
    vector<pair<double, double> > points(n + 1);

    for (int i = 0; i <= n; ++i) {
        cin >> req[i];

        double x, y;
        cin >> x >> y;
        points[i] = make_pair(x, y);
    }

    dists = PointDist(points);

    vector<vector<int>> decomp = Greedy(n, v, cap, req);
    for (int i = 0; i < v; ++i) {
        auto it = find(decomp[i].begin(), decomp[i].end(), 0);
        rotate(decomp[i].begin(), it, decomp[i].end());
        decomp[i].erase(decomp[i].begin());
        for (int x: decomp[i]) {
            cout << x << ' ';
        }
        cout << endl;
    }

    fclose(stdin);
    fclose(stdout);
    return 0;
}
