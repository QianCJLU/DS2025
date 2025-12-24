#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <ctime>
#include <cmath>

using namespace std;

// 复数类
class Complex {
public:
    double real, imag;
    Complex(double r = 0, double i = 0) : real(r), imag(i) {}
    double modulus() const { return sqrt(real * real + imag * imag); }
    bool operator==(const Complex& other) const {
        return real == other.real && imag == other.imag;
    }
    bool operator<(const Complex& other) const {
        double m1 = modulus(), m2 = other.modulus();
        if (m1 != m2) return m1 < m2;
        return real < other.real;
    }
    friend ostream& operator<<(ostream& os, const Complex& c) {
        os << "(" << c.real << "," << c.imag << ")";
        return os;
    }
};

// 随机生成无序复数向量
vector<Complex> randomComplexVector(int n, int minVal = -10, int maxVal = 10) {
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dis(minVal, maxVal);
    vector<Complex> vec;
    for (int i = 0; i < n; ++i) {
        vec.emplace_back(dis(gen), dis(gen));
    }
    return vec;
}

// 置乱
void shuffleVector(vector<Complex>& vec) {
    random_device rd;
    mt19937 gen(rd());
    shuffle(vec.begin(), vec.end(), gen);
}

// 查找
int findComplex(const vector<Complex>& vec, const Complex& target) {
    for (size_t i = 0; i < vec.size(); ++i) {
        if (vec[i] == target) return i;
    }
    return -1;
}

// 插入
void insertComplex(vector<Complex>& vec, const Complex& c) {
    vec.push_back(c);
}

// 删除
bool deleteComplex(vector<Complex>& vec, const Complex& c) {
    auto it = find(vec.begin(), vec.end(), c);
    if (it != vec.end()) {
        vec.erase(it);
        return true;
    }
    return false;
}

// 唯一化
void uniqueComplex(vector<Complex>& vec) {
    sort(vec.begin(), vec.end(), [](const Complex& a, const Complex& b) {
        if (a.real != b.real) return a.real < b.real;
        return a.imag < b.imag;
    });
    vec.erase(unique(vec.begin(), vec.end()), vec.end());
}

// 冒泡排序
void bubbleSort(vector<Complex>& vec) {
    for (size_t i = 0; i < vec.size(); ++i) {
        for (size_t j = 0; j < vec.size() - i - 1; ++j) {
            if (vec[j + 1] < vec[j]) {
                swap(vec[j], vec[j + 1]);
            }
        }
    }
}

// 归并排序
void merge(vector<Complex>& vec, int l, int m, int r) {
    int n1 = m - l + 1, n2 = r - m;
    vector<Complex> L(vec.begin() + l, vec.begin() + m + 1);
    vector<Complex> R(vec.begin() + m + 1, vec.begin() + r + 1);
    int i = 0, j = 0, k = l;
    while (i < n1 && j < n2) {
        if (L[i] < R[j]) vec[k++] = L[i++];
        else vec[k++] = R[j++];
    }
    while (i < n1) vec[k++] = L[i++];
    while (j < n2) vec[k++] = R[j++];
}
void mergeSort(vector<Complex>& vec, int l, int r) {
    if (l < r) {
        int m = l + (r - l) / 2;
        mergeSort(vec, l, m);
        mergeSort(vec, m + 1, r);
        merge(vec, l, m, r);
    }
}

// 区间查找
vector<Complex> intervalSearch(const vector<Complex>& vec, double m1, double m2) {
    vector<Complex> res;
    for (const auto& c : vec) {
        double mod = c.modulus();
        if (mod >= m1 && mod < m2) res.push_back(c);
    }
    return res;
}

int main() {
    int n = 20;
    vector<Complex> vec = randomComplexVector(n);

    cout << "原始无序向量:\n";
    for (auto& c : vec) cout << c << " ";
    cout << endl;

    // 置乱
    shuffleVector(vec);
    cout << "置乱后:\n";
    for (auto& c : vec) cout << c << " ";
    cout << endl;

    // 查找
    Complex target = vec[5];
    int idx = findComplex(vec, target);
    cout << "查找 " << target << " 的位置: " << idx << endl;

    // 插入
    insertComplex(vec, Complex(1, 2));
    cout << "插入 (1,2) 后:\n";
    for (auto& c : vec) cout << c << " ";
    cout << endl;

    // 删除
    deleteComplex(vec, Complex(1, 2));
    cout << "删除 (1,2) 后:\n";
    for (auto& c : vec) cout << c << " ";
    cout << endl;

    // 唯一化
    uniqueComplex(vec);
    cout << "唯一化后:\n";
    for (auto& c : vec) cout << c << " ";
    cout << endl;

    // 排序效率比较
    vector<Complex> vec1 = vec, vec2 = vec, vec3 = vec;
    shuffleVector(vec1); // 乱序
    sort(vec2.begin(), vec2.end()); // 顺序
    reverse(vec3.begin(), vec3.end()); // 逆序

    clock_t t1 = clock();
    bubbleSort(vec1);
    clock_t t2 = clock();
    cout << "冒泡排序（乱序）耗时: " << t2 - t1 << " ms\n";

    t1 = clock();
    bubbleSort(vec2);
    t2 = clock();
    cout << "冒泡排序（顺序）耗时: " << t2 - t1 << " ms\n";

    t1 = clock();
    bubbleSort(vec3);
    t2 = clock();
    cout << "冒泡排序（逆序）耗时: " << t2 - t1 << " ms\n";

    // 归并排序
    shuffleVector(vec1);
    sort(vec2.begin(), vec2.end());
    reverse(vec3.begin(), vec3.end());

    t1 = clock();
    mergeSort(vec1, 0, vec1.size() - 1);
    t2 = clock();
    cout << "归并排序（乱序）耗时: " << t2 - t1 << " ms\n";

    t1 = clock();
    mergeSort(vec2, 0, vec2.size() - 1);
    t2 = clock();
    cout << "归并排序（顺序）耗时: " << t2 - t1 << " ms\n";

    t1 = clock();
    mergeSort(vec3, 0, vec3.size() - 1);
    t2 = clock();
    cout << "归并排序（逆序）耗时: " << t2 - t1 << " ms\n";

    // 区间查找
    double m1 = 5, m2 = 10;
    vector<Complex> subvec = intervalSearch(vec2, m1, m2);
    cout << "模在[" << m1 << "," << m2 << ")区间的元素:\n";
    for (auto& c : subvec) cout << c << " ";
    cout << endl;

    return 0;
}