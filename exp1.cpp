#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <stack>
#include <string>
#include <cctype>
#include <stdexcept>

using namespace std;

// 一、复数类及相关操作
class Complex {
private:
    double real;   // 实部
    double imag;   // 虚部

public:
    // 构造函数
    Complex(double r = 0, double i = 0) : real(r), imag(i) {}

    // 获取模长
    double modulus() const {
        return sqrt(real * real + imag * imag);
    }

    // 重载比较运算符(用于排序)
    bool operator<(const Complex& other) const {
        if (modulus() != other.modulus()) {
            return modulus() < other.modulus();
        }
        return real < other.real;  // 模相等时比较实部
    }

    // 判断两个复数是否相等(实部和虚部均相同)
    bool operator==(const Complex& other) const {
        return (fabs(real - other.real) < 1e-9) && (fabs(imag - other.imag) < 1e-9);
    }

    // 输出复数
    friend ostream& operator<<(ostream& os, const Complex& c) {
        os << "(" << c.real << ", " << c.imag << ")";
        return os;
    }
};

// 复数向量生成
vector<Complex> generateRandomComplexVector(int size, double minVal, double maxVal) {
    vector<Complex> vec;
    for (int i = 0; i < size; ++i) {
        double real = minVal + (maxVal - minVal) * rand() / RAND_MAX;
        double imag = minVal + (maxVal - minVal) * rand() / RAND_MAX;
        vec.push_back(Complex(real, imag));
    }
    return vec;
}

// 复数向量置乱
void shuffleComplexVector(vector<Complex>& vec) {
    for (int i = vec.size() - 1; i > 0; --i) {
        int j = rand() % (i + 1);
        swap(vec[i], vec[j]);
    }
}

// 复数查找
int findComplex(const vector<Complex>& vec, const Complex& target) {
    for (int i = 0; i < vec.size(); ++i) {
        if (vec[i] == target) {
            return i;
        }
    }
    return -1;
}

// 复数插入
void insertComplex(vector<Complex>& vec, int pos, const Complex& elem) {
    if (pos >= 0 && pos <= vec.size()) {
        vec.insert(vec.begin() + pos, elem);
    }
}

// 复数删除
void deleteComplex(vector<Complex>& vec, int pos) {
    if (pos >= 0 && pos < vec.size()) {
        vec.erase(vec.begin() + pos);
    }
}

// 复数向量唯一化
void uniqueComplexVector(vector<Complex>& vec) {
    sort(vec.begin(), vec.end());
    auto last = unique(vec.begin(), vec.end());
    vec.erase(last, vec.end());
}

// 复数起泡排序
void bubbleSortComplex(vector<Complex>& vec) {
    int n = vec.size();
    for (int i = 0; i < n - 1; ++i) {
        for (int j = 0; j < n - i - 1; ++j) {
            if (!(vec[j] < vec[j + 1])) {
                swap(vec[j], vec[j + 1]);
            }
        }
    }
}

// 复数归并排序辅助函数
void mergeComplex(vector<Complex>& vec, int left, int mid, int right) {
    int n1 = mid - left + 1;
    int n2 = right - mid;

    vector<Complex> L(n1), R(n2);
    for (int i = 0; i < n1; ++i) L[i] = vec[left + i];
    for (int i = 0; i < n2; ++i) R[i] = vec[mid + 1 + i];

    int i = 0, j = 0, k = left;
    while (i < n1 && j < n2) {
        if (L[i] < R[j]) {
            vec[k++] = L[i++];
        }
        else {
            vec[k++] = R[j++];
        }
    }

    while (i < n1) vec[k++] = L[i++];
    while (j < n2) vec[k++] = R[j++];
}

// 复数归并排序
void mergeSortComplex(vector<Complex>& vec, int left, int right) {
    if (left < right) {
        int mid = left + (right - left) / 2;
        mergeSortComplex(vec, left, mid);
        mergeSortComplex(vec, mid + 1, right);
        mergeComplex(vec, left, mid, right);
    }
}

// 复数区间查找
vector<Complex> rangeSearchComplex(const vector<Complex>& vec, double m1, double m2) {
    vector<Complex> result;
    for (const auto& c : vec) {
        double mod = c.modulus();
        if (mod >= m1 && mod < m2) {
            result.push_back(c);
        }
    }
    return result;
}

// 二、基于栈的字符串计算器
template <typename T>
class Stack {
private:
    T* data;
    int topIdx;
    int capacity;
public:
    Stack(int size = 1000) : capacity(size), topIdx(-1) {
        data = new T[capacity];
    }
    ~Stack() {
        delete[] data;
    }
    bool isEmpty() const {
        return topIdx == -1;
    }
    void push(const T& item) {
        if (topIdx == capacity - 1) {
            throw runtime_error("Stack overflow");
        }
        data[++topIdx] = item;
    }
    T pop() {
        if (isEmpty()) {
            throw runtime_error("Stack underflow");
        }
        return data[topIdx--];
    }
    T top() const {
        if (isEmpty()) {
            throw runtime_error("Stack is empty");
        }
        return data[topIdx];
    }
};

enum Operator { ADD, SUB, MUL, DIV, POW, FAC, L_P, R_P, EOE };
const int N_OPTR = 9;

const char pri[N_OPTR][N_OPTR] = {
    { '>','>','<','<','<','<','<','>','>' },
    { '>','>','<','<','<','<','<','>','>' },
    { '>','>','>','>','<','<','<','>','>' },
    { '>','>','>','>','<','<','<','>','>' },
    { '>','>','>','>','>','<','<','>','>' },
    { '>','>','>','>','>','>',' ','>','>' },
    { '<','<','<','<','<','<','<','=',' ' },
    { ' ',' ',' ',' ',' ',' ',' ',' ',' ' },
    { '<','<','<','<','<','<','<',' ','=' }
};

int opToIdx(char op) {
    switch (op) {
    case '+': return ADD;
    case '-': return SUB;
    case '*': return MUL;
    case '/': return DIV;
    case '^': return POW;
    case '!': return FAC;
    case '(': return L_P;
    case ')': return R_P;
    case '\0': return EOE;
    default: throw runtime_error("Invalid operator");
    }
}

char getPriority(char opTop, char opCur) {
    int idxTop = opToIdx(opTop);
    int idxCur = opToIdx(opCur);
    return pri[idxTop][idxCur];
}

double calculate(double a, double b, char op) {
    switch (op) {
    case '+': return a + b;
    case '-': return a - b;
    case '*': return a * b;
    case '/':
        if (fabs(b) < 1e-9) throw runtime_error("Division by zero");
        return a / b;
    case '^': {
        double res = 1;
        int exp = (int)b;
        if (exp < 0) throw runtime_error("Negative exponent not supported");
        for (int i = 0; i < exp; i++) res *= a;
        return res;
    }
    case '!': {
        if (a < 0 || fabs(a - (int)a) > 1e-9)
            throw runtime_error("Factorial only supports non-negative integers");
        int n = (int)a;
        double res = 1;
        for (int i = 1; i <= n; i++) res *= i;
        return res;
    }
    default: throw runtime_error("Invalid operator in calculation");
    }
}

double stringCalculator(const string& expr) {
    Stack<double> numStack;
    Stack<char> opStack;
    opStack.push('\0');

    int i = 0;
    while (i < expr.size() || !opStack.isEmpty()) {
        if (i < expr.size() && (isdigit(expr[i]) || expr[i] == '.')) {
            double num = 0;
            int decimal = 0;
            bool isDecimal = false;
            while (i < expr.size() && (isdigit(expr[i]) || expr[i] == '.')) {
                if (expr[i] == '.') {
                    if (isDecimal) throw runtime_error("Invalid number format");
                    isDecimal = true;
                }
                else {
                    if (isDecimal) {
                        decimal = decimal * 10 + (expr[i] - '0');
                    }
                    else {
                        num = num * 10 + (expr[i] - '0');
                    }
                }
                i++;
            }
            while (decimal > 0) {
                num /= 10;
                decimal /= 10;
            }
            numStack.push(num);
        }
        else if (i < expr.size()) {
            char opCur = expr[i];
            while (true) {
                char opTop = opStack.top();
                char p = getPriority(opTop, opCur);
                if (p == '<') {
                    opStack.push(opCur);
                    i++;
                    break;
                }
                else if (p == '=') {
                    opStack.pop();
                    i++;
                    break;
                }
                else if (p == '>') {
                    char op = opStack.pop();
                    if (op == '!') {
                        double a = numStack.pop();
                        numStack.push(calculate(a, 0, op));
                    }
                    else {
                        if (numStack.isEmpty()) throw runtime_error("Invalid expression");
                        double b = numStack.pop();
                        if (numStack.isEmpty()) throw runtime_error("Invalid expression");
                        double a = numStack.pop();
                        numStack.push(calculate(a, b, op));
                    }
                }
                else {
                    throw runtime_error("Invalid expression");
                }
            }
        }
        else {
            char op = opStack.pop();
            if (op == '\0') break;
            if (op == '!') {
                double a = numStack.pop();
                numStack.push(calculate(a, 0, op));
            }
            else {
                if (numStack.isEmpty()) throw runtime_error("Invalid expression");
                double b = numStack.pop();
                if (numStack.isEmpty()) throw runtime_error("Invalid expression");
                double a = numStack.pop();
                numStack.push(calculate(a, b, op));
            }
        }
    }

    if (numStack.isEmpty()) throw runtime_error("Invalid expression");
    double result = numStack.pop();
    if (!numStack.isEmpty()) throw runtime_error("Invalid expression");
    return result;
}

// 三、柱状图最大矩形面积
int largestRectangleArea(vector<int>& heights) {
    stack<int> stk;
    heights.push_back(0);
    int maxArea = 0;

    for (int i = 0; i < heights.size(); ++i) {
        while (!stk.empty() && heights[i] < heights[stk.top()]) {
            int h = heights[stk.top()];
            stk.pop();
            int w = stk.empty() ? i : i - stk.top() - 1;
            maxArea = max(maxArea, h * w);
        }
        stk.push(i);
    }

    heights.pop_back();
    return maxArea;
}

vector<int> generateRandomHeights(int size) {
    vector<int> heights(size);
    for (int i = 0; i < size; ++i) {
        heights[i] = rand() % 105;
    }
    return heights;
}

void printHeights(const vector<int>& heights) {
    int n = heights.size();
    if (n <= 20) {
        for (int h : heights) cout << h << " ";
    }
    else {
        for (int i = 0; i < 10; ++i) cout << heights[i] << " ";
        cout << "... ";
        for (int i = n - 10; i < n; ++i) cout << heights[i] << " ";
    }
    cout << endl;
}

// 打印复数向量（简化输出）
void printComplexVector(const vector<Complex>& vec, const string& msg = "") {
    if (!msg.empty()) cout << msg << endl;
    int limit = min(5, (int)vec.size());
    for (int i = 0; i < limit; ++i) {
        cout << vec[i] << " ";
    }
    if (vec.size() > 10) cout << "... ";
    for (int i = max(0, (int)vec.size() - limit); i < vec.size(); ++i) {
        cout << vec[i] << " ";
    }
    cout << endl << "Size: " << vec.size() << endl << endl;
}

int main() {
    srand(time(0));

    // 测试1：复数类操作
    cout << "=== 复数类操作测试 ===" << endl;
    vector<Complex> compVec = generateRandomComplexVector(10, 0, 10);
    printComplexVector(compVec, "初始随机复数向量:");

    shuffleComplexVector(compVec);
    printComplexVector(compVec, "置乱后的复数向量:");

    if (!compVec.empty()) {
        Complex target = compVec[2];
        int pos = findComplex(compVec, target);
        cout << "查找 " << target << ": " << (pos != -1 ? "找到，位置 " + to_string(pos) : "未找到") << endl << endl;
    }

    insertComplex(compVec, 3, Complex(100, 200));
    printComplexVector(compVec, "插入元素后的向量:");

    deleteComplex(compVec, 5);
    printComplexVector(compVec, "删除元素后的向量:");

    compVec.push_back(Complex(1, 2));
    compVec.push_back(Complex(1, 2));
    printComplexVector(compVec, "添加重复元素后:");
    uniqueComplexVector(compVec);
    printComplexVector(compVec, "唯一化后:");

    sort(compVec.begin(), compVec.end());
    printComplexVector(compVec, "排序后的向量:");

    vector<Complex> rangeResult = rangeSearchComplex(compVec, 2, 5);
    printComplexVector(rangeResult, "模介于[2,5)的元素:");

    // 排序效率测试
    cout << "\n=== 复数排序效率测试 ===" << endl;
    vector<Complex> largeCompVec = generateRandomComplexVector(5000, 0, 100);
    vector<Complex> sortedComp = largeCompVec;
    sort(sortedComp.begin(), sortedComp.end());
    vector<Complex> reversedComp = sortedComp;
    reverse(reversedComp.begin(), reversedComp.end());

    clock_t start, end;

    start = clock();
    bubbleSortComplex(sortedComp);
    end = clock();
    cout << "起泡排序(顺序): " << (double)(end - start) / CLOCKS_PER_SEC << "s" << endl;

    start = clock();
    bubbleSortComplex(largeCompVec);
    end = clock();
    cout << "起泡排序(乱序): " << (double)(end - start) / CLOCKS_PER_SEC << "s" << endl;

    start = clock();
    bubbleSortComplex(reversedComp);
    end = clock();
    cout << "起泡排序(逆序): " << (double)(end - start) / CLOCKS_PER_SEC << "s" << endl;

    start = clock();
    mergeSortComplex(sortedComp, 0, sortedComp.size() - 1);
    end = clock();
    cout << "归并排序(顺序): " << (double)(end - start) / CLOCKS_PER_SEC << "s" << endl;

    start = clock();
    mergeSortComplex(largeCompVec, 0, largeCompVec.size() - 1);
    end = clock();
    cout << "归并排序(乱序): " << (double)(end - start) / CLOCKS_PER_SEC << "s" << endl;

    start = clock();
    mergeSortComplex(reversedComp, 0, reversedComp.size() - 1);
    end = clock();
    cout << "归并排序(逆序): " << (double)(end - start) / CLOCKS_PER_SEC << "s" << endl << endl;

    // 测试2：字符串计算器
    cout << "=== 字符串计算器测试 ===" << endl;
    vector<string> exprs = {
        "3+4*2",
        "(3+4)^2",
        "5!",
        "3*2+5!-10/2",
        "10/(2+3)",
        "3^2+4^2",
        "((5+3)*2)!",
        "10/0",        // 错误案例
        "3+(4*2",      // 错误案例
        "7!-2*3"
    };

    for (const string& expr : exprs) {
        try {
            cout << expr << " = " << stringCalculator(expr) << endl;
        }
        catch (const exception& e) {
            cout << expr << " = 无效表达式 (" << e.what() << ")" << endl;
        }
    }
    cout << endl;

    // 测试3：柱状图最大矩形面积
    cout << "=== 柱状图最大矩形最大矩形面积测试 ===" << endl;
    for (int i = 0; i < 10; ++i) {
        int size = rand() % 105 + 1;  // 长度范围：1~105
        vector<int> heights = generateRandomHeights(size);

        cout << "第" << i + 1 << "组 (长度: " << size << "):" << endl;
        printHeights(heights);
        cout << "最大面积: " << largestRectangleArea(heights) << "\n" << endl;
    }
    return 0;
}