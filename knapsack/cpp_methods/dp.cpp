#include <cstdio>
#include <iostream>
#include <vector>

using namespace std;

#include <cstdio>
#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

vector<int> DP(int n, int W, const vector<int>& cost, const vector<int>& weight) {
    vector<int> curr(W + 1, 0);
    vector<vector<bool>> chosen(n, vector<bool>(W + 1, false));
    
    for (int i = 0; i < n; ++i) {
        for (int j = W; j >= weight[i]; --j) {
            int new_val = curr[j - weight[i]] + cost[i];
            if (new_val > curr[j]) {
                curr[j] = new_val;
                chosen[i][j] = true;
            }
        }
    }
    
    int opt = 0;
    for (int w = 0; w <= W; ++w) {
        if (curr[w] > curr[opt]) {
            opt = w;
        }
    }
    
    vector<int> answer;
    int remaining = opt;
    for (int i = n - 1; i >= 0; --i) {
        if (chosen[i][remaining]) {
            answer.push_back(i);
            remaining -= weight[i];
        }
    }
    
    reverse(answer.begin(), answer.end());
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

    vector<int> chosen = DP(n, W, cost, weight);
    for (int i : chosen) {
        cout << i << ' ';
    }
    cout << endl;
    
    fclose(stdin);
    fclose(stdout);
    return 0;
}
