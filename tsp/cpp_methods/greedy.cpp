#include <cstdio>
#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <set>
#include <numeric>

using namespace std;

mt19937 rng(4257179);

double Dist(pair<double, double> a, pair<double, double> b) {
    double dx = a.first - b.first;
    double dy = a.second - b.second;
    return dx * dx + dy * dy;
}

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

    vector<int> ordering = Greedy(n, points);
    for (int i : ordering) {
        cout << i << ' ';
    }
    cout << endl;
    
    fclose(stdin);
    fclose(stdout);
    return 0;
}
