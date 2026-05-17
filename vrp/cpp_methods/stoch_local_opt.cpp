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
uniform_real_distribution<double> dst(0.0, 1.0);

double prob() {
    return dst(rng);
}

const double EPS = 1e-9;
double THRESHOLD;

using TimePoint = chrono::steady_clock::time_point;

bool TimeIsOver(const TimePoint& deadline) {
    return chrono::steady_clock::now() >= deadline;
}

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

void GreedySolve(
    int n,
    int v,
    double cap,
    const vector<double>& req,
    vector<Cycle>& cycles,
    function<double(double, double)> decide,
    TimePoint deadline
) {
    bool found = true;
    while ((cycles.size() > v || found) && !TimeIsOver(deadline)) {
        found = false;
        int cycle_id1, cycle_id2;
        int cycle_pos1, cycle_pos2;
        double best_add = INF;

        for (int i = 0; i < cycles.size(); ++i) {
            if (TimeIsOver(deadline)) return;

            for (int j = i + 1; j < cycles.size(); ++j) {
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
}

vector<vector<int>> Greedy(
    int n,
    int v,
    double cap,
    const vector<double> &req,
    vector<Cycle>& cycles,
    TimePoint deadline
) {
    auto decide1 = [cap](double lhs, double rhs) {
        if (lhs + rhs > cap) return 2 * INF;
        return 0.0;
    };

    GreedySolve(n, v, cap, req, cycles, decide1, deadline);

    if (cycles.size() > v && !TimeIsOver(deadline)) {
        auto decide2 = [&](double lhs, double rhs) {
          return max(0.0, lhs + rhs - cap) * Lambda;
        };
        GreedySolve(n, v, cap, req, cycles, decide2, deadline);
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

double GetDecompValue(int n, int v, double cap, const vector<double>& req, const vector<vector<int>>& decomp) {
    double result = 0;
    vector<double> taken(v);
    for (int i = 0; i < v; ++i) {
        result += Cost(decomp[i]);
        for (auto x : decomp[i]) {
            if (x != 0) {
                taken[i] += req[x];
            }
        }
    }
    for (int i = 0; i < v; ++i) {
        result += max(0.0, taken[i] - cap) * Lambda;
    }
    return result;
}

vector<vector<int>> LocalSearch
(
    int n,
    int v,
    double cap,
    const vector<double>& req,
    vector<vector<int>>& decomp,
    TimePoint deadline
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

    while (!TimeIsOver(deadline)) {

        double best_result = result;
        int best_id = -1;
        int best_change = -1;

        for (int id = 1; id <= n; ++id) {
            if (TimeIsOver(deadline)) return decomp;

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
            continue;
        }

        int best_id1, best_id2;
        for (int id1 = 1; id1 <= n; ++id1) {
            if (TimeIsOver(deadline)) return decomp;

            for (int id2 = id1 + 1; id2 <= n; ++id2) {
                if (mapping[id1] == mapping[id2]) continue;
                double cur_result = result;
                cur_result -= max(0.0, taken[mapping[id1]] - cap) * Lambda;
                cur_result -= max(0.0, taken[mapping[id2]] - cap) * Lambda;

                auto path1 = decomp[mapping[id1]];
                auto path2 = decomp[mapping[id2]];

                cur_result -= Cost(path1);
                cur_result -= Cost(path2);

                path1.erase(find(path1.begin(), path1.end(), id1));
                cur_result += OptInsVertex(path1, id2);

                path2.erase(find(path2.begin(), path2.end(), id2));
                cur_result += OptInsVertex(path2, id1);

                cur_result += max(0.0, taken[mapping[id1]] - req[id1] + req[id2] - cap) * Lambda;
                cur_result += max(0.0, taken[mapping[id2]] - req[id2] + req[id1] - cap) * Lambda;

                if (cur_result < best_result) {
                    best_result = cur_result;
                    best_id1 = id1;
                    best_id2 = id2;
                }
            }
        }

        if (best_result < result) {
            result = best_result;
            auto path1 = decomp[mapping[best_id1]];
            auto path2 = decomp[mapping[best_id2]];

            path1.erase(find(path1.begin(), path1.end(), best_id1));
            path2.erase(find(path2.begin(), path2.end(), best_id2));
            OptInsVertex(path1, best_id2);
            OptInsVertex(path2, best_id1);

            taken[mapping[best_id1]] += req[best_id2] - req[best_id1];
            taken[mapping[best_id2]] += req[best_id1] - req[best_id2];
            swap(mapping[best_id1], mapping[best_id2]);

            decomp[mapping[best_id1]] = path2;
            decomp[mapping[best_id2]] = path1;
        } else {
            break;
        }
    }

    return decomp;
}

vector<vector<int>> RestartLocalSearch(int n, int v, double cap, const vector<double>& req) {
    auto deadline = chrono::steady_clock::now() + chrono::seconds(6);

    vector<Cycle> cycles;
    for (int i = 1; i <= n; ++i) {
        cycles.emplace_back(vector<int>{i}, req[i]);
    }

    auto best_decomp = Greedy(n, v, cap, req, cycles, deadline);
    double best_result = GetDecompValue(n, v, cap, req, best_decomp);

    while (!TimeIsOver(deadline)) {
        auto sol = best_decomp;

        vector<Cycle> new_cycles;
        for (int i = 0; i < v; ++i) {
            if (TimeIsOver(deadline)) break;

            vector<int> ord;
            for (int j : sol[i]) {
                if (j != 0) {
                    ord.push_back(j);
                }
            }
            vector<int> path;
            double taken = 0;

            for (int x : ord) {
                if (prob() <= THRESHOLD) {
                    if (!path.empty()) {
                        new_cycles.emplace_back(path, taken);
                    }
                    path.clear();
                    taken = 0;
                }
                path.push_back(x);
                taken += req[x];
            }

            if (!path.empty()) {
                new_cycles.emplace_back(path, taken);
            }
        }

        if (TimeIsOver(deadline)) {
            break;
        }

        sol = Greedy(n, v, cap, req, new_cycles, deadline);

        if (TimeIsOver(deadline)) {
            break;
        }

        sol = LocalSearch(n, v, cap, req, sol, deadline);

        if (TimeIsOver(deadline)) {
            break;
        }

        auto cur_result = GetDecompValue(n, v, cap, req, sol);
        if (cur_result < best_result) {
            best_result = cur_result;
            best_decomp = sol;
        }
    }

    return best_decomp;
}

int main(int argc, char **argv) {
    freopen(argv[1], "r", stdin);
    freopen("tmp/ans.txt", "w", stdout);
    THRESHOLD = stod(argv[2]);
    
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
    auto decomp = RestartLocalSearch(n, v, cap, req);

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