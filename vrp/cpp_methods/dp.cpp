
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
#define route vector<int>

vector<vector<int>> BruteForce(int n, int v, double cap, const vector<double>& req) {
    vector<vector<double>> dp(1 << n, vector<double>(n, INF));
    vector<vector<route>> res(1 << n, vector<route>(n));

    for (int i = 0; i < n; ++i) {
        dp[1 << i][i] = dists.GetDist(0, i + 1);
        res[1 << i][i] = {i};
    }
    for (int msk = 1; msk < (1 << n); ++msk) {
        for (int last = 0; last < n; ++last) {
            if ((msk & (1 << last)) == 0) continue;
            for (int adj = 0; adj < n; ++adj) {
                if (msk & (1 << adj)) continue;
                int nmsk = msk ^ (1 << adj);
                double nres = dp[msk][last] + dists.GetDist(last + 1, adj + 1);
                if (dp[nmsk][adj] > nres) {
                    dp[nmsk][adj] = nres;
                    res[nmsk][adj] = res[msk][last];
                    res[nmsk][adj].push_back(adj);
                }
            }
        }
    }

    vector<double> opt_val(1 << n, INF);
    vector<route> opt_route(1 << n);
    for (int msk = 1; msk < (1 << n); ++msk) {
        double taken = 0;
        for (int last = 0; last < n; ++last) {
            if ((msk & (1 << last)) == 0) continue;
            taken += req[last + 1];
        }

        if (taken > cap) continue;

        for (int last = 0; last < n; ++last) {
            if ((msk & (1 << last)) == 0) continue;
            dp[msk][last] += dists.GetDist(last + 1, 0);

            if (opt_val[msk] > dp[msk][last]) {
                opt_val[msk] = dp[msk][last];
                opt_route[msk] = res[msk][last];
            }
        }

    }

    vector<vector<double>> covered(v + 1, vector<double>(1 << n, INF));
    vector<vector<int>> prev(v + 1, vector<int>(1 << n, 0));
    covered[0][0] = 0;

    for (int layer = 0; layer < v; ++layer) {
        for (int msk = 0; msk < (1 << n); ++msk) {
            covered[layer + 1][msk] = covered[layer][msk];
            for (int s = msk; s; s = (s - 1) & msk) {
                if (covered[layer + 1][msk] > covered[layer][msk ^ s] + opt_val[s]) {
                    covered[layer + 1][msk] = covered[layer][msk ^ s] + opt_val[s];
                    prev[layer + 1][msk] = s;
                }
            }
        }
    }

    vector<vector<int>> ans(v);
    int msk = (1 << n) - 1;
    for (int i = v - 1; i >= 0; --i) {
        int taken_here = prev[i + 1][msk];
        for (int x : opt_route[taken_here]) {
            ans[i].push_back(x + 1);
        }
        msk ^= taken_here;
    }
    return ans;
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
    vector<vector<int>> decomp(v);
    if (n <= 22) {
        decomp = BruteForce(n, v, cap, req);
    }
    for (int i = 0; i < v; ++i) {
        for (int x : decomp[i]) {
            cout << x << ' ';
        }
        cout << '\n';
    }
    fclose(stdin);
    fclose(stdout);
    return 0;
}

