#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#include <queue>
#include <set>

using namespace std;

template <typename T = double>
class TSP2OptFastSym {
public:
    TSP2OptFastSym(const vector<T>& distMatrix, int n, int k_nn = 20)
        : d(distMatrix), n(n), k(k_nn) {
        build_candidate_list();
    }

    void solve(vector<int>& tour) {
        bool improved = true;
        T best_cost = compute_cost(tour);

        while (improved) {
            improved = false;

            for (int i = 0; i < n; ++i) {
                int a = tour[i];
                int b = tour[(i + 1) % n];

                for (int neighbor : candidate[a]) {
                    int j = find_index(tour, neighbor);
                    int c = tour[j];
                    int d = tour[(j + 1) % n];

                    if (i == j || (i + 1) % n == j || (j + 1) % n == i)
                        continue;  // skip adjacent or same

                    T delta = -dist(a, b) - dist(c, d) + dist(a, c) + dist(b, d);
                    if (delta < -1e-9) {
                        if ((i + 1) % n < j)
                            reverse_segment(tour, i + 1, j);
                        else
                            reverse_segment(tour, j + 1, i);
                        best_cost += delta;
                        improved = true;
                        goto next_iter;
                    }
                }
            }
        next_iter:;
        }
    }

    T compute_cost(const vector<int>& tour) const {
        T cost = 0;
        for (int i = 0; i < n; ++i)
            cost += dist(tour[i], tour[(i + 1) % n]);
        return cost;
    }

private:
    const vector<T>& d;
    int n;
    int k;
    vector<vector<int>> candidate;

    inline T dist(int i, int j) const {
        return d[i * n + j];
    }

    void build_candidate_list() {
        candidate.resize(n);
        for (int i = 0; i < n; ++i) {
            priority_queue<pair<T, int>> pq;
            for (int j = 0; j < n; ++j) {
                if (i == j) continue;
                T dij = dist(i, j);
                if ((int)pq.size() < k)
                    pq.push({dij, j});
                else if (dij < pq.top().first) {
                    pq.pop();
                    pq.push({dij, j});
                }
            }
            while (!pq.empty()) {
                candidate[i].push_back(pq.top().second);
                pq.pop();
            }
        }
    }

    int find_index(const vector<int>& tour, int city) const {
        for (int i = 0; i < n; ++i)
            if (tour[i] == city) return i;
        return -1;
    }

    void reverse_segment(vector<int>& tour, int i, int j) const {
        while (i < j) {
            swap(tour[i % n], tour[j % n]);
            ++i;
            --j;
        }
    }
};


template <typename T = double>
class TSP3OptFastSym {
public:
    TSP3OptFastSym(const vector<T>& distMatrix, int n, int k_nn = 20)
        : d(distMatrix), n(n), k(k_nn) {
        build_candidate_list();
    }

    void solve(vector<int>& tour) {
        bool improved = true;
        T best_cost = compute_cost(tour);

        while (improved) {
            improved = false;

            for (int i = 0; i < n - 2; ++i) {
                int a = tour[i];
                int b = tour[(i + 1) % n];

                for (int x : candidate[a]) {
                    int j = find_index(tour, x);
                    if (j <= i + 1 || j >= n - 1) continue;

                    int c = tour[j];
                    int d = tour[(j + 1) % n];

                    for (int y : candidate[c]) {
                        int k = find_index(tour, y);
                        if (k <= j + 1 || k >= n) continue;

                        int e = tour[k];
                        int f = tour[(k + 1) % n];

                        // Try a single reconnection scheme (most effective one)
                        T delta = -dist(a, b) - dist(c, d) - dist(e, f)
                                  + dist(a, c) + dist(b, e) + dist(d, f);

                        if (delta < -1e-6) {
                            reverse_segment(tour, i + 1, j);
                            reverse_segment(tour, j + 1, k);
                            improved = true;
                            best_cost += delta;
                            goto next_iter;
                        }
                    }
                }
            }
        next_iter:;
        }
    }

    T compute_cost(const vector<int>& tour) const {
        T cost = 0;
        for (int i = 0; i < n; ++i)
            cost += dist(tour[i], tour[(i + 1) % n]);
        return cost;
    }

private:
    const vector<T>& d;
    int n;
    int k;
    vector<vector<int>> candidate;

    T dist(int i, int j) const {
        return d[i * n + j];
    }

    void build_candidate_list() {
        candidate.resize(n);
        for (int i = 0; i < n; ++i) {
            priority_queue<pair<T, int>> pq;
            for (int j = 0; j < n; ++j) {
                if (i == j) continue;
                T dij = dist(i, j);
                if ((int)pq.size() < k) pq.push({dij, j});
                else if (dij < pq.top().first) {
                    pq.pop();
                    pq.push({dij, j});
                }
            }
            while (!pq.empty()) {
                candidate[i].push_back(pq.top().second);
                pq.pop();
            }
        }
    }

    int find_index(const vector<int>& tour, int city) const {
        for (int i = 0; i < n; ++i)
            if (tour[i] == city) return i;
        return -1;
    }

    void reverse_segment(vector<int>& tour, int i, int j) const {
        while (i < j) {
            swap(tour[i], tour[j]);
            ++i;
            --j;
        }
    }
};


// Access dist[i][j] from flattened array
template <typename T> 
inline T getDist(const T* distMatrix, int n, int i, int j) {
    return distMatrix[i * n + j];
}

template <typename T> 
T calculateTourCost(const vector<int>& tour, const T* distMatrix, int n) {
    T cost = 0.0;
    for (int i = 0; i < n; ++i) {
        int from = tour[i];
        int to = tour[(i + 1) % n];
        cost += getDist(distMatrix, n, from, to);
    }
    return cost;
}

template <typename T> 
bool twoOptMove(std::vector<int>& tour, const T* distMatrix, int n, bool sym=false) {
    bool improved = false;
    T best_gain = 0.0;
    int best_i = -1, best_j = -1;

    for (int i = 0; i < n - 1; ++i) {
        for (int j = i + 2; j < n; ++j) {
            if (i == 0 && j == n - 1) continue;  // avoid full reversal

            int a = tour[i];
            int b = tour[i + 1];
            int c = tour[j];
            int d = tour[(j + 1) % n];

            // Old edges
            T oldCost = getDist(distMatrix, n, a, b)
                           + getDist(distMatrix, n, c, d);

            // New edges
            T newCost = getDist(distMatrix, n, a, c)
                           + getDist(distMatrix, n, b, d);

            // Internal edges (reversed)
            if(!sym){
                for (int k = i + 1; k < j; ++k) {
                    int u = tour[k];
                    int v = tour[k + 1];
                    oldCost += getDist<T>(distMatrix, n, u, v);
                    newCost += getDist<T>(distMatrix, n, v, u);  // reversed direction
                }
            }

            T gain = oldCost - newCost;
            if (gain > 0) {                
                reverse(tour.begin() + i + 1, tour.begin() + j + 1);
                return true;
            }
        }
    }


    return false;
}
