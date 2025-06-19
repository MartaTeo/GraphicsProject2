#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

using namespace std;

struct Vec3 {
    double x, y, z;
};

struct Vec2 {
    double x, y;
};

vector<Vec3> V;
vector<vector<int>> F;
vector<vector<int>> A;
vector<bool> B;
vector<Vec2> UV;

void read_obj(string file) {
    ifstream in(file);
    string s;
    while (in >> s) {
        if (s == "v") {
            Vec3 v;
            in >> v.x >> v.y >> v.z;
            V.push_back(v);
        } else if (s == "f") {
            vector<int> f(3);
            in >> f[0] >> f[1] >> f[2];
            for (int& x : f) x -= 1;
            F.push_back(f);
        }
    }
}

void build_graph() {
    int n = V.size();
    A.assign(n, {});
    B.assign(n, false);
    for (auto& f : F) {
        for (int i = 0; i < 3; i++) {
            int a = f[i];
            int b = f[(i + 1) % 3];
            if (find(A[a].begin(), A[a].end(), b) == A[a].end()) A[a].push_back(b);
            if (find(A[b].begin(), A[b].end(), a) == A[b].end()) A[b].push_back(a);
        }
    }

    for (int i = 0; i < n; i++) {
        for (int j : A[i]) {
            int count = 0;
            for (auto& f : F) {
                for (int k = 0; k < 3; k++) {
                    int u = f[k], v = f[(k + 1) % 3];
                    if ((u == i && v == j) || (u == j && v == i))
                        count++;
                }
            }
            if (count == 1)
                B[i] = true;
        }
    }
}

void map_boundary() {
    int n = V.size();
    UV.assign(n, {0, 0});
    vector<int> boundary;
    for (int i = 0; i < n; i++)
        if (B[i]) boundary.push_back(i);
    int m = boundary.size();
    for (int i = 0; i < m; i++) {
        double t = 2 * M_PI * i / m;
        UV[boundary[i]].x = cos(t);
        UV[boundary[i]].y = sin(t);
    }
}

void solve_linear(vector<vector<double>>& M, vector<double>& b, vector<double>& x) {
    int n = M.size();
    x = vector<double>(n, 0);
    for (int i = 0; i < n; i++) {
        int pivot = i;
        for (int j = i + 1; j < n; j++)
            if (abs(M[j][i]) > abs(M[pivot][i]))
                pivot = j;
        swap(M[i], M[pivot]);
        swap(b[i], b[pivot]);

        double diag = M[i][i];
        for (int j = i; j < n; j++)
            M[i][j] /= diag;
        b[i] /= diag;

        for (int j = i + 1; j < n; j++) {
            double factor = M[j][i];
            for (int k = i; k < n; k++)
                M[j][k] -= factor * M[i][k];
            b[j] -= factor * b[i];
        }
    }

    for (int i = n - 1; i >= 0; i--) {
        x[i] = b[i];
        for (int j = i + 1; j < n; j++)
            x[i] -= M[i][j] * x[j];
    }
}

void embed() {
    int n = V.size();
    vector<vector<double>> M(n, vector<double>(n, 0));
    vector<double> bx(n, 0), by(n, 0);

    for (int i = 0; i < n; i++) {
        if (B[i]) {
            M[i][i] = 1;
            bx[i] = UV[i].x;
            by[i] = UV[i].y;
        } else {
            int d = A[i].size();
            M[i][i] = d;
            for (int j : A[i])
                M[i][j] = -1;
        }
    }

    vector<double> x, y;
    solve_linear(M, bx, x);
    solve_linear(M, by, y);

    for (int i = 0; i < n; i++)
        UV[i] = {x[i], y[i]};
}

void write_obj(string file) {
    ofstream out(file);
    for (auto& v : V)
        out << "v " << v.x << " " << v.y << " " << v.z << "\n";
    for (auto& uv : UV)
        out << "vt " << uv.x << " " << uv.y << "\n";
    for (auto& f : F)
        out << "f " << f[0]+1 << "/" << f[0]+1 << " " << f[1]+1 << "/" << f[1]+1 << " " << f[2]+1 << "/" << f[2]+1 << "\n";
}

int main() {
    read_obj("goethe.obj");
    build_graph();
    map_boundary();
    embed();
    write_obj("goethe_uv.obj");
    return 0;
}
