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

        for (int i = 0; i < lhs.size(); ++i) {
            for (int j = 0; j < rhs.size(); ++j) {
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

    double avg_call = 0;
    for (int i = 0; i < m; ++i) {
        avg_call += calls[i];
    }
    avg_call /= m;

    vector<int> f(m, -1);
    double total_dist = 0.0;
    double bad = 0;
    
    for (int id : ord) {
        int best = -1;
        double best_d = INF;

        for (int j = 0; j < k; ++j) {
            if (rem[j] < calls[id]) continue;
            double d = dists.GetDist(id, choice[j]);
            if (d < best_d) {
                best_d = d;
                best = j;
            }
        }

        if (best == -1) {
            bad += calls[id] / avg_call;
        } else {
            f[id] = best;
            rem[best] -= calls[id];
            total_dist += best_d;
        }

    }

    total_dist += Lambda * bad;

    return {total_dist, f};
}

const int EPOCH_COUNT = 1e4;

double STEP = 0.999;

pair<vector<int>, vector<int>>
Annealing(
    int n,
    int m,
    const vector<double>& price,
    const vector<double>& cap,
    const vector<pair<double, double>>& shops,
    const vector<double>& calls,
    const vector<pair<double, double>>& customers
) {

    double best_value = INF;
    vector<int> best_choice(n);
    vector<int> best_f(m, -1);
    

    auto eval_choice = [&](const vector<int>& choice) -> pair<double, vector<int>> {
        double shop_cost = 0.0;
        vector<int> shop_list;
        for (int i = 0; i < n; ++i) {
            if (choice[i]) {
                shop_cost += price[i];
                shop_list.push_back(i);
            }
        }
    
        auto res = GetMatching(shop_list, cap, calls, shops, customers);
        return {shop_cost + res.cost, res.f};
    };

    auto relax = [&](const vector<int>& choice) {
        auto [value, f] = eval_choice(choice);
        if (value < best_value) {
            best_value = value;
            best_choice = choice;
            best_f = f;
        }
        return value;
    };


    auto start = chrono::steady_clock::now();
    
    while (true) {
        auto now = chrono::steady_clock::now();
        double elapsed = chrono::duration<double>(now - start).count();
        if (elapsed > 60.0) break;
        
        vector<int> choice(n);
        double temperature = 100.0;
        
        for (int i = 0; i < n; ++i) {
            if (rng() % 2) choice[i] = 1;
        }
        double was_value = relax(choice);

        for (int epoch = 0; epoch < EPOCH_COUNT; ++epoch) {
            int bit = rng() % n;
            choice[bit] = (choice[bit] ^ true);
            double new_value = relax(choice);

            if (new_value < was_value || prob() < exp((was_value - new_value) / temperature)) {
                was_value = new_value;
            } else {
                choice[bit] = (choice[bit] ^ true);
            }

            temperature *= STEP;
        }
        
    }


    vector<int> best_shops;
    for (int i = 0; i < n; ++i) {
        if (best_choice[i]) best_shops.push_back(i);
    }
    
    return make_pair(best_shops, best_f);
}

int main(int argc, char** argv) {
    freopen(argv[1], "r", stdin);
    freopen("tmp/ans.txt", "w", stdout);
    STEP = stod(argv[2]);
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
    
    auto [choice, f] = Annealing(n, m, price, cap, shops, calls, customers);
    
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
