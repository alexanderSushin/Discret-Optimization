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
#include <functional>

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


const double Lambda = 1000;

double OptInsVertex(vector<int> &path, int ver) {
    double best_add = INF;
    int best_id = -1;

    for (int i = 0; i < path.size(); ++i) {
        int nxt = i;
        int prv = (i == 0 ? (int) path.size() - 1 : i - 1);
        double add = -dists.GetDist(path[prv], path[nxt])
                     + dists.GetDist(path[prv], ver)
                     + dists.GetDist(ver, path[nxt]);

        if (add < best_add) {
            best_add = add;
            best_id = i;
        }
    }

    path.insert(path.begin() + best_id, ver);
    double result = Cost(path);
    return result;
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

    auto GreedySolve = [&](function<double(double, double)> decide) {
        bool found = true;
        while (cycles.size() > v || found) {
            found = false;
            int cycle_id1, cycle_id2;
            int cycle_pos1, cycle_pos2;
            double best_add = INF;
            for (int i = 0; i < cycles.size(); ++i) {
                for (int j = i + 1; j < cycles.size(); ++j) {
                    // if (cycles[i].GetTaken() + cycles[j].GetTaken() > cap) continue;
                    for (int p1 : {0, cycles[i].Size() - 1}) {
                        for (int p2 : {0, cycles[j].Size() - 1}) {
                            double add = decide(cycles[i].GetTaken(), cycles[j].GetTaken());

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

            if (best_add >= INF) {
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
    };


    auto decide1 = [cap](double lhs, double rhs) {
        if (lhs + rhs > cap) return 2 * INF;
        return 0.0;
    };

    GreedySolve(decide1);

    // cout << cycles.size() << endl;
    if (cycles.size() > v) {
        auto decide2 = [&](double lhs, double rhs) {
          return max(0.0, lhs + rhs - cap) * Lambda;
        };
        GreedySolve(decide2);
    }


    vector<vector<int>> decomp;
    for (int i = 0; i < cycles.size(); ++i) {
        decomp.push_back(cycles[i].order);
    }

    for (int i = 0; i < decomp.size(); ++i) {
        decomp[i].push_back(0);
    }

    while (decomp.size() < v) {
        decomp.push_back({0});
    }


    return decomp;
}

vector<vector<int>> LocalSearch
(
    int n,
    int v,
    double cap,
    const vector<double>& req,
    vector<vector<int>>& decomp
) {
    vector<double> taken(v);
    vector<int> mapping(n + 1);
    double result = 0;
    for (int i = 0; i < v; ++i) {
        result += Cost(decomp[i]);
        for (auto x : decomp[i]) {
            if (x != 0) {
                taken[i] += req[x];
                mapping[x] = i;
            }
        }
    }
    for (int i = 0; i < v; ++i) {
        result += max(0.0, taken[i] - cap) * Lambda;
    }
    auto start_time = chrono::steady_clock::now();
    while (1) {
        // cout << result << endl;
        auto now = chrono::steady_clock::now();
        auto elapsed = chrono::duration_cast<chrono::seconds>(now - start_time).count();
        if (elapsed > 60) {
            break;
        }

        double best_result = result;
        int best_id = -1;
        int best_change = -1;

        for (int id = 1; id <= n; ++id) {
            for (int path_id = 0; path_id < v; ++path_id) {
                if (mapping[id] == path_id) continue;
                auto path1 = decomp[mapping[id]];
                auto path2 = decomp[path_id];
                double cur_result = result;
                cur_result -= max(0.0, taken[mapping[id]] - cap) * Lambda;
                cur_result += max(0.0, taken[mapping[id]] - req[id] - cap) * Lambda;

                cur_result -= max(0.0, taken[path_id] - cap) * Lambda;
                cur_result += max(0.0, taken[path_id] + req[id] - cap) * Lambda;

                cur_result -= Cost(path1);
                cur_result -= Cost(path2);

                path1.erase(find(path1.begin(), path1.end(), id));
                cur_result += Cost(path1);
                cur_result += OptInsVertex(path2, id);

                if (cur_result < best_result) {
                    best_result = cur_result;
                    best_id = id;
                    best_change = path_id;
                }

            }
        }

        if (best_result < result) {
            result = best_result;

            auto path1 = decomp[mapping[best_id]];
            auto path2 = decomp[best_change];

            path1.erase(find(path1.begin(), path1.end(), best_id));
            OptInsVertex(path2, best_id);
            decomp[mapping[best_id]] = path1;
            taken[mapping[best_id]] -= req[best_id];
            decomp[best_change] = path2;
            taken[best_change] += req[best_id];
            mapping[best_id] = best_change;
        } else {
            break;
        }
    }
    return decomp;
}

int main(int argc, char **argv) {
    freopen(argv[1], "r", stdin);
    freopen("tmp/ans.txt", "w", stdout);

    int n, v;
    cin >> n >> v;
    double cap;
    cin >> cap;

    vector<double> req(n + 1);
    vector<pair<double, double> > points(n + 1);
    points[0] = make_pair(0, 0);

    for (int i = 1; i <= n; ++i) {
        cin >> req[i];

        double x, y;
        cin >> x >> y;
        points[i] = make_pair(x, y);
    }

    dists = PointDist(points);

    vector<vector<int>> decomp = Greedy(n, v, cap, req);
    LocalSearch(n, v, cap, req, decomp);
    // double result = 0;
    // for (int i = 0; i < v; ++i) {
    //     result += Cost(decomp[i]);
    // }
    // cout << result << endl;

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
