#include <cstdio>
#include <iostream>
#include <vector>

using namespace std;

vector<int> Greedy(int n, int W, const vector<int>& cost, const vector<int>& weight) {
    vector<int> order(n);
    for (int i = 0; i < n; ++i) {
        order[i] = i;
    }
    sort(order.begin(), order.end(), [&](int i, int j) {
        return cost[i] / (long double)weight[i] > cost[j] / (long double)weight[j];
    });
    vector<int> answer;
    int taken = 0;
    for (int i : order) {
        if (taken + weight[i] <= W) {
            taken += weight[i];
            answer.emplace_back(i);
        }
    }
    return answer;
}

int main(int argc, char** argv) {
    freopen(argv[1], "r", stdin);
    freopen("tmp/ans.txt", "w", stdout);

    int n, W;
    cin >> n >> W;
    vector<int> cost(n), weight(n);
    for (int i = 0; i < n; ++i) {
        cin >> cost[i] >> weight[i];
    }

    vector<int> chosen = Greedy(n, W, cost, weight);
    for (int i : chosen) {
        cout << i << ' ';
    }
    cout << endl;
    
    fclose(stdin);
    fclose(stdout);
    return 0;
}
