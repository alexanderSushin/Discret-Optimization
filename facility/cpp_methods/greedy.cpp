#include <cstdio>
#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <set>
#include <numeric>
#include <chrono>

using namespace std;

mt19937 rng(4257179);

double Dist(pair<double, double> a, pair<double, double> b) {
    double dx = a.first - b.first;
    double dy = a.second - b.second;
    return sqrt(dx * dx + dy * dy);
}

const double INF = 1e18;
double RATIO = 0.3;

struct MatchResult {
    double cost;
    vector<int> f;
};

MatchResult GetMatching(
    const vector<int>& choice,
    const vector<double>& cap,
    const vector<double>& calls,
    const vector<pair<double, double>>& shops,
    const vector<pair<double, double>>& customers
) {
    int m = (int)calls.size();
    int k = (int)choice.size();

    vector<double> rem(k);
    for (int j = 0; j < k; ++j) rem[j] = cap[choice[j]];

    vector<int> ord(m);
    iota(ord.begin(), ord.end(), 0);

    sort(ord.begin(), ord.end(), [&](int a, int b) {
        if (calls[a] != calls[b]) return calls[a] > calls[b];
        return a < b;
    });

    vector<int> f(m, -1);
    double total_dist = 0.0;

    for (int id : ord) {
        int best = -1;
        double best_d = INF;

        for (int j = 0; j < k; ++j) {
            if (rem[j] < calls[id]) continue;
            double d = Dist(customers[id], shops[choice[j]]);
            if (d < best_d) {
                best_d = d;
                best = j;
            }
        }

        if (best == -1) {
            return {INF, {}};
        }

        f[id] = best;
        rem[best] -= calls[id];
        total_dist += best_d;
    }

    return {total_dist, f};
}

pair<vector<int>, vector<int>>
Greedy(
    int n,
    int m,
    const vector<double>& price,
    const vector<double>& cap,
    const vector<pair<double, double>>& shops,
    const vector<double>& calls,
    const vector<pair<double, double>>& customers
) {
    vector<int> order(n);
    iota(order.begin(), order.end(), 0);

    sort(order.begin(), order.end(), [&](int a, int b) {
        if (price[a] != price[b]) return price[a] < price[b];
        if (cap[a] != cap[b]) return cap[a] > cap[b];
        return a < b;
    });

    double best_value = INF;
    vector<int> best_choice, best_f;

    auto try_order = [&](const vector<int>& ord) {
        vector<int> choice;
        double shop_cost = 0.0;

        for (int pos = 0; pos < n; ++pos) {
            int shop_id = ord[pos];
            choice.push_back(shop_id);
            shop_cost += price[shop_id];

            auto res = GetMatching(choice, cap, calls, shops, customers);
            if (res.cost < INF) {
                double value = shop_cost + res.cost;
                if (value < best_value) {
                    best_value = value;
                    best_choice = choice;
                    best_f = res.f;
                }
                break; 
            }
        }
    };

    auto eval_choice = [&](const vector<int>& choice) -> pair<double, vector<int>> {
        double shop_cost = 0.0;
        for (int x : choice) shop_cost += price[x];
    
        auto res = GetMatching(choice, cap, calls, shops, customers);
        if (res.cost >= INF) return {INF, {}};
        return {shop_cost + res.cost, res.f};
    };

    try_order(order);
    auto start = chrono::steady_clock::now();

    while (true) {
        auto now = chrono::steady_clock::now();
        double elapsed = chrono::duration<double>(now - start).count();
        if (elapsed > 30.0) break;
    
        shuffle(order.begin(), order.end(), rng);
        try_order(order);
    }

    start = chrono::steady_clock::now();

    while (true) {
        auto now = chrono::steady_clock::now();
        double elapsed = chrono::duration<double>(now - start).count();
        if (elapsed > 30.0 || best_choice.empty()) break;
        if (best_choice.empty()) break;
    
        vector<int> cand = best_choice;
        vector<bool> in(n, false);
        for (int x : cand) in[x] = true;
    
        int max_remove = max(1, (int)ceil(RATIO * cand.size()));
        uniform_int_distribution<int> dist(1, max_remove);
        int k_remove = dist(rng);
    
        for (int t = 0; t < k_remove; ++t) {
            int pos = uniform_int_distribution<int>(0, (int)cand.size() - 1)(rng);
            in[cand[pos]] = false;
            cand.erase(cand.begin() + pos);
        }
    
        vector<int> new_sets;
        for (int j = 0; j < n; ++j) {
            if (!in[j]) new_sets.push_back(j);
        }
        shuffle(new_sets.begin(), new_sets.end(), rng);
    
        for (int j : new_sets) {
            cand.push_back(j);
            
            auto [val, f] = eval_choice(cand);
        
            if (val < best_value) {
                best_value = val;
                best_choice = cand;
                best_f = f;
            }
        }
        
    }
    return {best_choice, best_f};
}

int main(int argc, char** argv) {
    freopen(argv[1], "r", stdin);
    freopen("tmp/ans.txt", "w", stdout);
    RATIO = stod(argv[2]);
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

    auto [choice, f] = Greedy(n, m, price, cap, shops, calls, customers);
    
    for (auto i : choice) {
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
