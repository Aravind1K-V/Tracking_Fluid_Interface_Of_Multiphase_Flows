#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>
#include <queue>
#include <map>
#include <set>

using namespace std;

struct Point {
    double x, y;
    Point(double x = 0, double y = 0) : x(x), y(y) {}

    bool operator<(const Point& other) const {
        if (x != other.x) {
            return x < other.x;
        }
        return y < other.y;
    }
};

double distance(const Point& a, const Point& b) {
    return sqrt(pow(a.x - b.x, 2) + pow(a.y - b.y, 2));
}

double angleWithCentroid(const Point& centroid, const Point& p) {
    return atan2(p.y - centroid.y, p.x - centroid.x);
}

Point calculateCentroid(const vector<Point>& points) {
    double x_sum = 0, y_sum = 0;
    for (const auto& p : points) {
        x_sum += p.x;
        y_sum += p.y;
    }
    return Point(x_sum / points.size(), y_sum / points.size());
}

double trapeziumArea(const Point& a, const Point& b, const Point& c, const Point& d) {
    vector<Point> points = {a, b, c, d};
    Point centroid = calculateCentroid(points);

    sort(points.begin(), points.end(), [&centroid](const Point& a, const Point& b) {
        return angleWithCentroid(centroid, a) < angleWithCentroid(centroid, b);
    });

    double area = 0.0;
    int n = points.size();
    for (int i = 0; i < n; ++i) {
        int j = (i + 1) % n;
        area += points[i].x * points[j].y - points[j].x * points[i].y;
    }

    return fabs(area) / 2.0;
}

double area(const Point& a, const Point& b, const Point& c, const Point& int_i, const Point& ext_i, const Point& int_i_plus_1, const Point& ext_i_plus_1) {
    return trapeziumArea(a, b, int_i, int_i_plus_1);
}

Point interpolate(const Point& external, const Point& internal, double x) {
    return Point(x * internal.x + (1 - x) * external.x,
                 x * internal.y + (1 - x) * external.y);
}

double calculate_cost(const Point& M_i, const Point& M_i_plus_1, const vector<vector<float>>& grid, double lambda, int p, const Point& ext_i, const Point& int_i, const Point& ext_i_plus_1, const Point& int_i_plus_1) {
    int x1 = (ext_i.x + int_i.x + ext_i_plus_1.x + int_i_plus_1.x) / 4;
    int y1 = (ext_i.y + int_i.y + ext_i_plus_1.y + int_i_plus_1.y) / 4;
    double calculated_volume_fraction = area(M_i, M_i_plus_1, Point(M_i.x, M_i_plus_1.y), int_i, ext_i, int_i_plus_1, ext_i_plus_1);
    double vol_target = grid[x1][y1];
    double vol_diff = fabs(calculated_volume_fraction - vol_target);
    double length = distance(M_i, M_i_plus_1);
    return pow(vol_diff, p) + lambda * length;
}

vector<Point> dp_algorithm(const vector<pair<pair<int, int>, pair<int, int>>>& pairs, 
                           const vector<vector<float>>& grid, 
                           double lambda, int p, int num_interpolations, vector<Point>& internalarr, vector<Point>& externalarr) {
    int N = pairs.size();
    vector<vector<double>> dp(N + 1, vector<double>(num_interpolations, numeric_limits<double>::max()));
    
    //Stores the optimal previous step and is used to check for the current index.
    vector<vector<pair<int, int>>> prev(N + 1, vector<pair<int, int>>(num_interpolations, {-1, -1}));
    
    
    //Setting the cost as 0 for all interpolations for the segment 0.
    for (int j = 0; j < num_interpolations; ++j) {
        dp[0][j] = 0;
    }

    for (int i = 1; i <= N; i++) {
        
        //Pick the external and internal points from the pairs.
        Point external(pairs[i % N].first.first, pairs[i % N].first.second);
        Point internal(pairs[i % N].second.first, pairs[i % N].second.second);
        
        //Running for all interpolation points.
        for (int j = 0; j < num_interpolations; ++j) {
            double x_i = static_cast<double>(j) / (num_interpolations - 1);
            Point M_i = interpolate(external, internal, x_i);
            
            //Storing the previous best.
            for (int k = 0; k < num_interpolations; ++k) {
                Point prev_external(pairs[i - 1].first.first, pairs[i - 1].first.second);
                Point prev_internal(pairs[i - 1].second.first, pairs[i - 1].second.second);
                double x_i_minus_1 = static_cast<double>(k) / (num_interpolations - 1);
                Point M_i_minus_1 = interpolate(prev_external, prev_internal, x_i_minus_1);
                
                //Cost to reach from the previous best to current index.
                double cost = calculate_cost(M_i_minus_1, M_i, grid, lambda, p, prev_external, prev_internal, external, internal);
                if (dp[i - 1][k] + cost < dp[i][j]) {
                    dp[i][j] = dp[i - 1][k] + cost;
                    prev[i][j] = {i - 1, k};
                }
            }
        }
    }
    
    
    //Finding the minimum interpolation value possible.
    int min_j = 0;
    for (int j = 0; j < num_interpolations; ++j) {
        if (dp[N][j] < dp[N][min_j]) {
            min_j = j;
        }
    }
    
    //This retrieves the previous indices prev_i and prev_j from the prev table. 
    // The prev table was built during the DP process to remember the previous state that led to the current optimal state. 
    //Essentially, it stores the traceback pointers.
    
    vector<Point> optimal_points;
    int i = N, j = min_j;
    while (i > 0) {
        Point external(pairs[i % N].first.first, pairs[i % N].first.second);
        Point internal(pairs[i % N].second.first, pairs[i % N].second.second);
        double x = static_cast<double>(j) / num_interpolations;
        optimal_points.push_back(interpolate(external, internal, x));
        internalarr.push_back(internal);
        externalarr.push_back(external);
        auto [prev_i, prev_j] = prev[i][j];
        i = prev_i;
        j = prev_j;
    }

    return optimal_points;
}

void check_volume_fractions(const vector<Point>& optimal_points, const vector<vector<float>>& grid, vector<Point>& internalarr, vector<Point>& externalarr) {
    cout << "\nChecking Volume Fractions:\n";


    for (int i = 0; i < optimal_points.size(); i++) {
        int j = (i + 1) % optimal_points.size();
        Point M_i = optimal_points[i];
        Point M_i_plus_1 = optimal_points[j];
        int x1 = (externalarr[i].x + internalarr[i].x + externalarr[j].x + internalarr[j].x) / 4;
        int y1 = (externalarr[i].y + internalarr[i].y + externalarr[j].y + internalarr[j].y) / 4;
        double calculated_volume_fraction = area(M_i, M_i_plus_1, Point(M_i.x, M_i_plus_1.y), internalarr[i], externalarr[i], internalarr[j], externalarr[j]);
        double expected_volume_fraction = grid[x1][y1];

        cout << "From Point: (" << M_i.x << ", " << M_i.y << ") To Point: (" << M_i_plus_1.x << ", " << M_i_plus_1.y << ") => Expected Volume Fraction is: " << expected_volume_fraction << ", Calculated Volume Fraction is: " << calculated_volume_fraction << endl;
    }
}

int main() {
    vector<vector<float>> grid = {
        {0, 0, 0, 0, 0},
        {0, 0.175, 0.5, 0.05, 0},
        {0, 0.7, 1, 0.2, 0},
        {0, 0.21, 0.6, 0.06, 0},
        {0, 0, 0, 0, 0}
    };

    vector<vector<int>> neighbours(grid.size(), vector<int>(grid[0].size(), 0));
    for (int i = 1; i <= grid.size() - 2; i++) {
        for (int j = 1; j <= grid[0].size() - 2; j++) {
            int cnt = 0;
            if (grid[i][j] > 0) {
                if (grid[i - 1][j] > 0) cnt++;
                if (grid[i + 1][j] > 0) cnt++;
                if (grid[i][j - 1] > 0) cnt++;
                if (grid[i][j + 1] > 0) cnt++;
                neighbours[i][j] = cnt;
            }
        }
    }

    for (int i = 0; i < neighbours.size(); i++) {
        for (int j = 0; j < neighbours[0].size(); j++) {
            cout << neighbours[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;

    vector<vector<string>> tempgrid(grid.size(), vector<string>(grid[0].size(), " "));
    map<pair<pair<int, int>, pair<int, int>>, bool> mp1;

    for (int i = 1; i <= grid.size() - 2; i++) {
        for (int j = 1; j <= grid[0].size() - 2; j++) {
            if (neighbours[i][j] != 0 && grid[i][j] <= 0.5) {
                if (neighbours[i][j] <= 2 && neighbours[i - 1][j] == 0 && neighbours[i + 1][j] != 0 && neighbours[i][j - 1] == 0 && neighbours[i][j + 1] != 0) {
                    tempgrid[i][j] = "/";
                    mp1[{{i + 1, j}, {i, j + 1}}] = true;
                } else if (neighbours[i][j] == 3 && neighbours[i - 1][j] == 0 && neighbours[i][j - 1] != 0 && neighbours[i][j + 1] != 0 && neighbours[i + 1][j] != 0) {
                    tempgrid[i][j] = "-";
                    mp1[{{i, j}, {i, j + 1}}] = true;
                } else if (neighbours[i][j] <= 2 && neighbours[i - 1][j] == 0 && neighbours[i][j - 1] != 0 && neighbours[i][j + 1] == 0) {
                    tempgrid[i][j] = "\\";
                    mp1[{{i, j}, {i + 1, j + 1}}] = true;
                } else if (neighbours[i][j] == 3 && neighbours[i - 1][j] != 0 && neighbours[i][j - 1] == 0 && neighbours[i + 1][j] != 0 && neighbours[i][j + 1] != 0) {
                    tempgrid[i][j] = "|";
                    mp1[{{i, j}, {i + 1, j}}] = true;
                } else if (neighbours[i][j] == 3 && neighbours[i - 1][j] != 0 && neighbours[i][j - 1] != 0 && neighbours[i][j + 1] == 0 && neighbours[i + 1][j] != 0) {
                    tempgrid[i][j] = "|";
                    mp1[{{i, j + 1}, {i + 1, j + 1}}] = true;
                } else if (neighbours[i][j] <= 2 && neighbours[i - 1][j] != 0 && neighbours[i][j - 1] == 0 && neighbours[i][j + 1] != 0 && neighbours[i + 1][j] == 0) {
                    tempgrid[i][j] = "\\";
                    mp1[{{i, j}, {i + 1, j + 1}}] = true;
                } else if (neighbours[i][j] == 3 && neighbours[i - 1][j] != 0 && neighbours[i][j - 1] != 0 && neighbours[i][j + 1] != 0 && neighbours[i + 1][j] == 0) {
                    tempgrid[i][j] = "-";
                    mp1[{{i + 1, j}, {i + 1, j + 1}}] = true;
                } else if (neighbours[i][j] <= 2 && neighbours[i - 1][j] != 0 && neighbours[i + 1][j] == 0 && neighbours[i][j - 1] != 0 && neighbours[i][j + 1] == 0) {
                    tempgrid[i][j] = "/";
                    mp1[{{i + 1, j}, {i, j + 1}}] = true;
                } else {
                    tempgrid[i][j] = "$";
                }
            } else if (neighbours[i][j] != 0 && grid[i][j] > 0.5) {
                if (neighbours[i][j] <= 2 && neighbours[i - 1][j] == 0 && neighbours[i + 1][j] != 0 && neighbours[i][j - 1] == 0 && neighbours[i][j + 1] != 0) {
                    mp1[{{i + 1, j}, {i, j}}] = true;
                    mp1[{{i, j}, {i, j + 1}}] = true;
                } else if (neighbours[i][j] == 3 && neighbours[i - 1][j] == 0 && neighbours[i][j - 1] != 0 && neighbours[i][j + 1] != 0 && neighbours[i + 1][j] != 0) {
                    tempgrid[i][j] = "-";
                    mp1[{{i, j}, {i, j + 1}}] = true;
                } else if (neighbours[i][j] <= 2 && neighbours[i - 1][j] == 0 && neighbours[i][j - 1] != 0 && neighbours[i][j + 1] == 0) {
                    tempgrid[i][j] = "\\";
                    mp1[{{i, j}, {i, j + 1}}] = true;
                    mp1[{{i, j + 1}, {i + 1, j + 1}}] = true;
                } else if (neighbours[i][j] == 3 && neighbours[i - 1][j] != 0 && neighbours[i][j - 1] == 0 && neighbours[i + 1][j] != 0 && neighbours[i][j + 1] != 0) {
                    tempgrid[i][j] = "|";
                    mp1[{{i, j}, {i + 1, j}}] = true;
                } else if (neighbours[i][j] == 3 && neighbours[i - 1][j] != 0 && neighbours[i][j - 1] != 0 && neighbours[i][j + 1] == 0 && neighbours[i + 1][j] != 0) {
                    tempgrid[i][j] = "|";
                    mp1[{{i, j + 1}, {i + 1, j + 1}}] = true;
                } else if (neighbours[i][j] <= 2 && neighbours[i - 1][j] != 0 && neighbours[i][j - 1] == 0 && neighbours[i][j + 1] != 0 && neighbours[i + 1][j] == 0) {
                    tempgrid[i][j] = "\\";
                    mp1[{{i, j}, {i + 1, j}}] = true;
                    mp1[{{i + 1, j}, {i + 1, j + 1}}] = true;
                } else if (neighbours[i][j] == 3 && neighbours[i - 1][j] != 0 && neighbours[i][j - 1] != 0 && neighbours[i][j + 1] != 0 && neighbours[i + 1][j] == 0) {
                    tempgrid[i][j] = "-";
                    mp1[{{i + 1, j}, {i + 1, j + 1}}] = true;
                } else if (neighbours[i][j] <= 2 && neighbours[i - 1][j] != 0 && neighbours[i + 1][j] == 0 && neighbours[i][j - 1] != 0 && neighbours[i][j + 1] == 0) {
                    tempgrid[i][j] = "/";
                    mp1[{{i + 1, j}, {i + 1, j + 1}}] = true;
                    mp1[{{i + 1, j + 1}, {i, j + 1}}] = true;
                } else {
                    tempgrid[i][j] = "$";
                }
            }
        }
    }

    for (int i = 0; i < tempgrid.size(); i++) {
        for (int j = 0; j < tempgrid[0].size(); j++) {
            cout << tempgrid[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;

    set<Point> internalCurve, externalCurve;
    map<pair<int, int>, bool> mp2;
    int count = 0;
    vector<pair<int, int>> orderedpairs;
    queue<pair<int, int>> q1;
    map<pair<int, int>, bool> vis;

    cout << "External" << endl;
    for (auto& i : mp1) {
        if (count == 0) {
            count++;
            q1.push({i.first.first.first, i.first.first.second});
        }
        cout << "edge from:{" << i.first.first.first << "," << i.first.first.second << "} to :{" << i.first.second.first << "," << i.first.second.second << "}" << endl;

        externalCurve.insert(Point(i.first.first.first, i.first.first.second));
        externalCurve.insert(Point(i.first.second.first, i.first.second.second));
        mp2[i.first.first] = true;
        mp2[i.first.second] = true;
    }
    cout << endl;

    while (!q1.empty()) {
        pair<int, int> fr = q1.front();
        q1.pop();
        orderedpairs.push_back(fr);
        vis[fr] = true;
        bool isfound = false;

        for (auto hell : mp1) {
            pair<pair<int, int>, pair<int, int>> t4 = hell.first;
            if ((!isfound) && fr.first == t4.first.first && fr.second == t4.first.second) {
                if (!vis[t4.second]) {
                    q1.push(t4.second);
                    isfound = true;
                }
            }
            if ((!isfound) && fr.first == t4.second.first && fr.second == t4.second.second) {
                if (!vis[t4.first]) {
                    q1.push(t4.first);
                    isfound = true;
                }
            }

            if (isfound) {
                break;
            }
        }
        if (!isfound) {
            break;
        }
    }

    vector<vector<int>> arr(tempgrid.size(), vector<int>(tempgrid[0].size(), 0));
    for (int i = 1; i < arr.size() - 1; i++) {
        for (int j = 1; j < arr[0].size() - 1; j++) {
            if (tempgrid[i - 1][j - 1] != " " && tempgrid[i][j - 1] != " " && tempgrid[i - 1][j] != " " && tempgrid[i][j] != " ") {
                if (!mp2[{i, j}]) {
                    arr[i][j] = 1;
                }
            }
        }
    }

    for (int i = 0; i < arr.size(); i++) {
        for (int j = 0; j < arr[0].size(); j++) {
            cout << arr[i][j] << " ";
        }
        cout << endl;
    }

    map<pair<int, int>, vector<pair<int, int>>> mp;
    vector<vector<int>> dir = {{-1, 0}, {1, 0}, {0, 1}, {0, -1}};
    set<pair<pair<int, int>, pair<int, int>>> uniqueEdges;

    for (int i = 0; i < arr.size(); i++) {
        for (int j = 0; j < arr[0].size(); j++) {
            if (arr[i][j] == 1) {
                queue<vector<int>> q;
                q.push({i, j});
                while (!q.empty()) {
                    int x1 = q.front()[0];
                    int y1 = q.front()[1];
                    q.pop();
                    arr[x1][y1] = 2;
                    for (const auto& k : dir) {
                        int x = x1 + k[0];
                        int y = y1 + k[1];
                        if (x >= 1 && x < arr.size() - 1) {
                            if (y >= 1 && y < arr[0].size() - 1) {
                                if (arr[x][y] == 1) {
                                    q.push({x, y});
                                    pair<int, int> p1 = {x1, y1};
                                    pair<int, int> p2 = {x, y};
                                    if (p1 > p2) swap(p1, p2);
                                    uniqueEdges.insert({p1, p2});
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    cout << endl;
    cout << "Internal" << endl;
    for (const auto& edge : uniqueEdges) {
        cout << "Edge from: {" << edge.first.first << "," << edge.first.second << "} to: {"
             << edge.second.first << "," << edge.second.second << "}" << endl;
        internalCurve.insert(Point(edge.first.first, edge.first.second));
        internalCurve.insert(Point(edge.second.first, edge.second.second));
    }

    vector<pair<int, int>> dir1 = {{-1, 0}, {1, 0}, {0, 1}, {0, -1}};
    vector<pair<pair<int, int>, pair<int, int>>> pairs;
    for (auto i : orderedpairs) {
        vector<pair<int, int>> temp12;
        for (auto j : dir1) {
            int x = i.first + j.first;
            int y = i.second + j.second;
            if (internalCurve.find(Point(x, y)) != internalCurve.end()) {
                pairs.push_back({{i.first, i.second}, {x, y}});
                temp12.push_back({x, y});
            }
        }
        if (temp12.size() == 2) {
            mp[temp12[0]].push_back(temp12[1]);
        }
    }

    cout << "Internal" << endl;
    for (auto i : internalCurve) {
        cout << "(" << i.x << "," << i.y << ") ";
    }
    cout << endl;

    cout << "External" << endl;
    for (auto i : externalCurve) {
        cout << "(" << i.x << "," << i.y << ") ";
    }
    cout << endl;
    cout << endl;

    cout << "Pairs for DP of {Ei,Ii}" << endl;
    for (auto i : pairs) {
        cout << "(" << i.first.first << "," << i.first.second << ") and (" << i.second.first << "," << i.second.second << ")" << endl;
    }
    cout << endl;

    double lambda = 0.001;
    int p = 2;
    int num_interpolations = 1000;

    vector<Point> internalarr;
    vector<Point> externalarr;
    vector<Point> optimal_points = dp_algorithm(pairs, grid, lambda, p, num_interpolations, internalarr, externalarr);

    cout << "Optimal Interpolated Points:" << endl;
    for (const auto& p : optimal_points) {
        cout << "(" << p.x << ", " << p.y << ")" << endl;
    }

    double total_cost = 0;
    for (size_t i = 0; i < optimal_points.size(); ++i) {
        total_cost += calculate_cost(optimal_points[i], optimal_points[(i + 1) % optimal_points.size()], grid, lambda, p, externalarr[i], internalarr[i], externalarr[(i + 1) % optimal_points.size()], internalarr[(i + 1) % optimal_points.size()]);
    }

    cout << "Total Cost: " << total_cost << endl;

    check_volume_fractions(optimal_points, grid, internalarr, externalarr);

    return 0;
}