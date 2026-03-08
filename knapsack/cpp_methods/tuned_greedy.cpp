#include <cstdio>
#include <iostream>
#include <vector>

using namespace std;

int n, W;
vector<int> cost, weight;
vector<int> order;

vector<int> Greedy(int ptr, int& score) {
    vector<int> answer;
    int taken = 0;
    
    for (int p = ptr; p < n; ++p) {
        int i = order[p];
        if (taken + weight[i] <= W) {
            score += cost[i];
            taken += weight[i];
            answer.emplace_back(i);
        }
    }
    return answer;
}

int main(int argc, char** argv) {
    freopen(argv[1], "r", stdin);
    freopen("tmp/ans.txt", "w", stdout);

    cin >> n >> W;
    cost.resize(n), weight.resize(n);
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
    
    vector<int> chosen;
    int opt = 0;
    for (int ptr = 0; ptr < n; ++ptr) {
        int score = 0;
        auto now = Greedy(ptr, score);
        if (score > opt) {
            chosen = now;
            opt = score;
        }
    }
    
    for (int i : chosen) {
        cout << i << ' ';
    }
    cout << endl;
    
    fclose(stdin);
    fclose(stdout);
    return 0;
}
