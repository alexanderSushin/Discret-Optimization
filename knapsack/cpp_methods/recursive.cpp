#include <cstdio>
#include <iostream>
#include <vector>
#include <utility>

using namespace std;

int n, W;
vector<int> cost, weight;

int opt = 0;
vector<int> answer;

vector<int> order;

int GetLowerBound(int ptr, int total_cost, int total_weight) {
    for (int p = ptr; p < n; ++p) {
        int i = order[p];
        if (total_weight + weight[i] <= W) {
            total_cost += cost[i];
            total_weight += weight[i];
        }
    }
    return total_cost;
}

long double GetUpperBound(int ptr, int total_cost, int total_weight) {
    long double cur_cost = total_cost;
    int cur_weight = total_weight;
    for (int p = ptr; p < n; ++p) {
        int i = order[p];
        if (cur_weight + weight[i] <= W) {
            cur_cost += cost[i];
            cur_weight += weight[i];
        } else {
            int can = W - cur_weight;
            cur_cost += can * (long double)cost[i] / weight[i];
            return cur_cost;
        }
    }
    return cur_cost;
}

pair<int, long double> GetBounds(int ptr, int total_cost, int total_weight) {
    return make_pair
                (
                    GetLowerBound(ptr, total_cost, total_weight),
                    GetUpperBound(ptr, total_cost, total_weight)
                );
}


void Rec(int ptr, vector<int> ids, int total_cost, int total_weight) {
    if (ptr == n || total_weight == W) {
        if (total_cost > opt) {
            opt = total_cost;
            answer = ids;
        }
        return;
    }
    auto [lower_bound, upper_bound] = GetBounds(ptr, total_cost, total_weight);
    if (upper_bound <= opt) {
        return;
    }

    if (lower_bound > opt) {
        answer = ids;
        opt = lower_bound;
        int copy_total_weight = total_weight;
        for (int p = ptr; p < n; ++p) {
            int i = order[p];
            if (copy_total_weight + weight[i] <= W) {
                answer.emplace_back(i);
                copy_total_weight += weight[i];
            }
        }

    }
    int i = order[ptr];
    if (total_weight + weight[i] <= W) {
        vector<int> to_rec = ids;
        to_rec.emplace_back(i);
        Rec(ptr + 1, to_rec, total_cost + cost[i], total_weight + weight[i]);
    }
    Rec(ptr + 1, ids, total_cost, total_weight);
}

int main(int argc, char** argv) {
    freopen(argv[1], "r", stdin);
    freopen("tmp/ans.txt", "w", stdout);

    cin >> n >> W;
    cost.resize(n);
    weight.resize(n);

    for (int i = 0; i < n; ++i) {
        cin >> cost[i] >> weight[i];
    }

    order.resize(n);
    for (int i = 0; i < n; ++i) {
        order[i] = i;
    }
    sort(order.begin(), order.end(), [&](int i, int j) {
        return cost[i] / (long double)weight[i] > cost[j] / (long double)weight[j];
    });

    // initial solution
    int taken = 0;
    for (int i : order) {
        if (taken + weight[i] <= W) {
            taken += weight[i];
            opt += cost[i];
            answer.emplace_back(i);
        }
    }

    Rec(0, vector<int>{}, 0, 0);

    for (int i : answer) {
        cout << i << ' ';
    }
    cout << endl;
    
    fclose(stdin);
    fclose(stdout);
    return 0;
}
