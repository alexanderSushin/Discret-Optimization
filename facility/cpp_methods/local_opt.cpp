#include <cstdio>
#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <set>
#include <numeric>
#include <chrono>
#include <cmath>

using namespace std;

mt19937 rng(4257179);
uniform_real_distribution<double> dst(0.0, 1.0);

double prob() {
    return dst(rng);
}

double Lambda = 1e6;

double Dist(pair<double, double> a, pair<double, double> b) {
    double dx = a.first - b.first;
    double dy = a.second - b.second;
    return sqrt(dx * dx + dy * dy);
}

struct PointDist {
    vector<vector<double>> dist_matrix;

    PointDist(vector<pair<double, double>> lhs, vector<pair<double, double>> rhs) {
        dist_matrix.resize(lhs.size(), vector<double>(rhs.size()));

        for (int i = 0; i < (int)lhs.size(); ++i) {
            for (int j = 0; j < (int)rhs.size(); ++j) {
                dist_matrix[i][j] = Dist(lhs[i], rhs[j]);
            }
        }
    }

    PointDist() = default;

    double GetDist(int i, int j) {
        return dist_matrix[i][j];
    }
};

PointDist dists;

const double INF = 1e18;

double FacilityCost(
    int sid,
    int used_count,
    double load,
    const vector<double>& price,
    const vector<double>& cap
) {
    double result = 0.0;

    if (used_count > 0) {
        result += price[sid];
    }

    result += max(0.0, load - cap[sid]) * Lambda;
    return result;
}

double GetCost(
    int n,
    int m,
    const vector<int>& choice,
    const vector<double>& price,
    const vector<double>& cap,
    const vector<double>& calls
) {
    double result = 0.0;

    vector<int> used(n, 0);
    vector<double> ncap(n, 0.0);

    for (int i = 0; i < m; ++i) {
        int sid = choice[i];

        if (sid < 0 || sid >= n) {
            return INF;
        }

        used[sid]++;
        ncap[sid] += calls[i];
        result += dists.GetDist(i, sid);
    }

    for (int sid = 0; sid < n; ++sid) {
        result += FacilityCost(sid, used[sid], ncap[sid], price, cap);
    }

    return result;
}

vector<int>
LocalOpt(
    int n,
    int m,
    const vector<double>& price,
    const vector<double>& cap,
    const vector<pair<double, double>>& shops,
    const vector<double>& calls,
    const vector<pair<double, double>>& customers
) {
    vector<int> match(m, -1);
    vector<double> ncap(n, 0.0);
    vector<int> used(n, 0);

    double value = 0.0;

    for (int cid = 0; cid < m; ++cid) {
        int sid = rng() % n;
        match[cid] = sid;
        used[sid]++;
        ncap[sid] += calls[cid];
        value += dists.GetDist(cid, sid);
    }

    for (int sid = 0; sid < n; ++sid) {
        value += FacilityCost(sid, used[sid], ncap[sid], price, cap);
    }

    bool found = true;
    const double EPS = 1e-9;

    while (found) {
        found = false;

        double best_delta = INF;
        int best_cid = -1;
        int best_sid = -1;

        for (int cid = 0; cid < m; ++cid) {
            int old_sid = match[cid];

            for (int new_sid = 0; new_sid < n; ++new_sid) {
                if (old_sid == new_sid) continue;

                double delta = 0.0;

                delta += dists.GetDist(cid, new_sid);
                delta -= dists.GetDist(cid, old_sid);

                double before =
                    FacilityCost(old_sid, used[old_sid], ncap[old_sid], price, cap) +
                    FacilityCost(new_sid, used[new_sid], ncap[new_sid], price, cap);

                double after =
                    FacilityCost(
                        old_sid,
                        used[old_sid] - 1,
                        ncap[old_sid] - calls[cid],
                        price,
                        cap
                    ) +
                    FacilityCost(
                        new_sid,
                        used[new_sid] + 1,
                        ncap[new_sid] + calls[cid],
                        price,
                        cap
                    );

                delta += after - before;

                if (delta < best_delta) {
                    best_delta = delta;
                    best_cid = cid;
                    best_sid = new_sid;
                }
            }
        }

        if (best_delta < -EPS) {
            int old_sid = match[best_cid];
            int new_sid = best_sid;

            value += best_delta;

            used[old_sid]--;
            ncap[old_sid] -= calls[best_cid];

            used[new_sid]++;
            ncap[new_sid] += calls[best_cid];

            match[best_cid] = new_sid;

            found = true;
            continue;
        }

        best_delta = 0.0;

        int best_c1 = -1;
        int best_c2 = -1;

        for (int c1 = 0; c1 < m; ++c1) {
            int s1 = match[c1];

            for (int c2 = c1 + 1; c2 < m; ++c2) {
                int s2 = match[c2];

                if (s1 == s2) continue;

                double delta = 0.0;

                delta += dists.GetDist(c1, s2);
                delta += dists.GetDist(c2, s1);
                delta -= dists.GetDist(c1, s1);
                delta -= dists.GetDist(c2, s2);

                double before =
                    FacilityCost(s1, used[s1], ncap[s1], price, cap) +
                    FacilityCost(s2, used[s2], ncap[s2], price, cap);

                double new_load_s1 = ncap[s1] - calls[c1] + calls[c2];
                double new_load_s2 = ncap[s2] - calls[c2] + calls[c1];

                double after =
                    FacilityCost(s1, used[s1], new_load_s1, price, cap) +
                    FacilityCost(s2, used[s2], new_load_s2, price, cap);

                delta += after - before;

                if (delta < best_delta) {
                    best_delta = delta;
                    best_c1 = c1;
                    best_c2 = c2;
                }
            }
        }

        if (best_delta < -EPS) {
            int s1 = match[best_c1];
            int s2 = match[best_c2];

            value += best_delta;

            ncap[s1] = ncap[s1] - calls[best_c1] + calls[best_c2];
            ncap[s2] = ncap[s2] - calls[best_c2] + calls[best_c1];

            swap(match[best_c1], match[best_c2]);

            found = true;
        }
    }

    return match;
}

vector<int>
RepeatedLocalOpt(
    int n,
    int m,
    const vector<double>& price,
    const vector<double>& cap,
    const vector<pair<double, double>>& shops,
    const vector<double>& calls,
    const vector<pair<double, double>>& customers
) {
    double best_value = INF;
    vector<int> best_choice(m, -1);

    auto relax = [&](const vector<int>& choice) {
        double value = GetCost(n, m, choice, price, cap, calls);

        if (value < best_value) {
            best_value = value;
            best_choice = choice;
        }

        return value;
    };

    auto start = chrono::steady_clock::now();

    while (true) {
        auto now = chrono::steady_clock::now();
        double elapsed = chrono::duration<double>(now - start).count();

        if (elapsed > 60.0) break;

        auto choice = LocalOpt(n, m, price, cap, shops, calls, customers);
        relax(choice);
    }

    return best_choice;
}

int main(int argc, char** argv) {
    freopen(argv[1], "r", stdin);
    freopen("tmp/ans.txt", "w", stdout);

    int n, m;
    cin >> n >> m;

    vector<double> price(n);
    vector<double> cap(n);
    vector<pair<double, double>> shops(n);

    for (int i = 0; i < n; ++i) {
        cin >> price[i] >> cap[i];

        double x, y;
        cin >> x >> y;

        shops[i] = make_pair(x, y);
    }

    vector<double> calls(m);
    vector<pair<double, double>> customers(m);

    for (int i = 0; i < m; ++i) {
        cin >> calls[i];

        double x, y;
        cin >> x >> y;

        customers[i] = make_pair(x, y);
    }

    dists = PointDist(customers, shops);

    auto match = RepeatedLocalOpt(n, m, price, cap, shops, calls, customers);

    vector<int> shop_list;

    for (int i = 0; i < m; ++i) {
        shop_list.push_back(match[i]);
    }

    sort(shop_list.begin(), shop_list.end());
    shop_list.resize(unique(shop_list.begin(), shop_list.end()) - shop_list.begin());

    vector<int> f(m);

    for (int i = 0; i < m; ++i) {
        f[i] = lower_bound(shop_list.begin(), shop_list.end(), match[i]) - shop_list.begin();
    }

    for (auto i : shop_list) {
        cout << i << ' ';
    }

    cout << endl;

    for (auto i : f) {
        cout << i << ' ';
    }

    cout << endl;

    fclose(stdin);
    fclose(stdout);

    return 0;
}