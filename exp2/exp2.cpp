#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <cstring>
#include <fstream>
#include <queue>
#include <vector>
#include <string>
#include <cctype>
#include <limits> // 用于忽略文件读取后的换行符
using namespace std;

// 定义 Rank 类型（与 size_t 匹配，避免类型转换丢失数据）
typedef size_t Rank;

// 位图类 Bitmap（修复浅拷贝问题，实现深拷贝）
class Bitmap {
private:
    unsigned char* M;
    Rank N, _sz; // N：字节数，_sz：有效位个数
    void init(Rank n) {
        M = new unsigned char[N = (n + 7) / 8]; // 向上取整到字节
        memset(M, 0, N);
        _sz = 0;
    }
    void expand(Rank k) { // 确保能容纳第k位（k从0开始）
        if (k < 8 * N) return;
        Rank newN = (k + 7) / 8;
        unsigned char* oldM = M;
        init(newN);
        memcpy(M, oldM, N); // 复制原数据
        delete[] oldM;
    }
public:
    // 构造函数
    Bitmap(Rank n = 8) { init(n); }
    // 拷贝构造函数（深拷贝）
    Bitmap(const Bitmap& other) {
        N = other.N;
        _sz = other._sz;
        M = new unsigned char[N];
        memcpy(M, other.M, N);
    }
    // 赋值运算符（深拷贝，防止自赋值和重复释放）
    Bitmap& operator=(const Bitmap& other) {
        if (this == &other) return *this; // 避免自赋值
        delete[] M; // 释放当前对象的内存
        N = other.N;
        _sz = other._sz;
        M = new unsigned char[N];
        memcpy(M, other.M, N);
        return *this;
    }
    // 文件读取构造函数
    Bitmap(char* file, Rank n = 8) {
        init(n);
        FILE* fp = fopen(file, "r");
        if (fp) {
            fread(M, sizeof(char), N, fp);
            fclose(fp);
        }
        _sz = 0;
        for (Rank k = 0; k < n; k++) _sz += test(k);
    }
    // 析构函数
    ~Bitmap() { delete[] M; M = NULL; _sz = 0; }

    Rank size() { return _sz; }
    void set(Rank k) {
        expand(k);
        if (!test(k)) _sz++; // 仅未设置时计数+1
        M[k >> 3] |= (0x80 >> (k & 0x07));
    }
    void clear(Rank k) {
        expand(k);
        if (test(k)) _sz--; // 仅已设置时计数-1
        M[k >> 3] &= ~(0x80 >> (k & 0x07));
    }
    bool test(Rank k) {
        expand(k);
        return (M[k >> 3] & (0x80 >> (k & 0x07))) != 0;
    }
    void dump(char* file) {
        FILE* fp = fopen(file, "w");
        if (fp) {
            fwrite(M, sizeof(char), N, fp);
            fclose(fp);
        }
    }
    char* bits2string(Rank n) {
        if (n == 0) return new char[1]{ '\0' };
        expand(n - 1);
        char* s = new char[n + 1];
        s[n] = '\0';
        for (Rank i = 0; i < n; i++) s[i] = test(i) ? '1' : '0';
        return s;
    }
};

// 二叉树节点结构（存储频率+字母索引）
template <typename T>
struct BinNode {
    T freq; // 频率
    int charIdx; // 字母索引（0-25，-1表示非叶子节点）
    BinNode<T>* left;
    BinNode<T>* right;
    // 构造函数
    BinNode(T f, int idx = -1) : freq(f), charIdx(idx), left(NULL), right(NULL) {}
};

// 二叉树类 BinTree（统一管理节点内存，避免重复释放）
template <typename T>
class BinTree {
protected:
    BinNode<T>* root;
public:
    BinTree() : root(NULL) {}
    ~BinTree() { destroy(root); } // 析构时释放所有节点
    BinNode<T>* getRoot() { return root; }
    void setRoot(BinNode<T>* node) { root = node; }
private:
    // 后序遍历销毁节点（确保每个节点仅释放一次）
    void destroy(BinNode<T>* node) {
        if (!node) return;
        destroy(node->left);
        destroy(node->right);
        delete node;
        node = NULL; // 避免野指针
    }
};

// Huffman树类（修复节点内存重复释放问题）
class HuffTree : public BinTree<int> {
public:
    void buildHuffTree(int freq[]) {
        // 优先队列：小顶堆（频率小的节点优先）
        priority_queue<BinNode<int>*, vector<BinNode<int>*>, Compare> pq;
        // 初始化叶子节点（仅频率>0的字母）
        for (int i = 0; i < 26; i++) {
            if (freq[i] > 0) {
                pq.push(new BinNode<int>(freq[i], i)); // 绑定字母索引
            }
        }
        // 处理只有1个叶子节点的边界情况
        if (pq.size() == 1) {
            BinNode<int>* leaf = pq.top();
            pq.pop();
            BinNode<int>* rootNode = new BinNode<int>(leaf->freq);
            rootNode->left = leaf;
            setRoot(rootNode);
            return;
        }
        // 构建Huffman树核心逻辑
        while (pq.size() > 1) {
            BinNode<int>* left = pq.top(); pq.pop();
            BinNode<int>* right = pq.top(); pq.pop();
            // 创建父节点（频率求和，非叶子节点索引为-1）
            BinNode<int>* parent = new BinNode<int>(left->freq + right->freq);
            parent->left = left;
            parent->right = right;
            pq.push(parent);
        }
        // 设置根节点（取出后队列清空，避免重复访问）
        if (!pq.empty()) {
            setRoot(pq.top());
            pq.pop();
        }
    }
private:
    // 比较器：用于优先队列排序
    struct Compare {
        bool operator()(BinNode<int>* a, BinNode<int>* b) {
            return a->freq > b->freq; // 小顶堆：频率小的优先
        }
    };
};

// Huffman编码类（基于Bitmap，优化赋值逻辑）
class HuffCode {
private:
    Bitmap codeBits;
    Rank length;
public:
    HuffCode() : length(0) {}
    // 从字符串构建编码
    void fromString(const string& codeStr) {
        length = codeStr.size();
        codeBits = Bitmap(length); // 现在支持深拷贝，无内存问题
        for (Rank i = 0; i < length; i++) {
            if (codeStr[i] == '1') codeBits.set(i);
            else codeBits.clear(i);
        }
    }
    // 编码转换为字符串
    string toString() {
        if (length == 0) return "";
        char* cstr = codeBits.bits2string(length);
        string str(cstr);
        delete[] cstr; // 释放临时内存
        return str;
    }
    Rank getLength() { return length; }
};

// Huffman编码工具类（优化索引合法性检查）
class HuffEncoder {
private:
    HuffTree huffTree;
    HuffCode codeTable[26]; // 0-25对应a-z
    int freq[26]; // 缓存字母频率
public:
    // 构造函数：初始化频率数组
    HuffEncoder() {
        memset(freq, 0, 26 * sizeof(int));
    }
    // 统计频率（输出到外部数组）
    void countFrequency(const string& text, int freqOut[]) {
        memset(freqOut, 0, 26 * sizeof(int));
        for (Rank i = 0; i < text.size(); i++) {
            char c = text[i];
            if (isalpha(c)) {
                c = tolower(c);
                int idx = c - 'a';
                if (idx >= 0 && idx < 26) { // 确保索引合法
                    freqOut[idx]++;
                }
            }
        }
    }
    // 统计频率（缓存到内部成员）
    void countFrequency(const string& text) {
        memset(freq, 0, 26 * sizeof(int));
        for (Rank i = 0; i < text.size(); i++) {
            char c = text[i];
            if (isalpha(c)) {
                c = tolower(c);
                int idx = c - 'a';
                if (idx >= 0 && idx < 26) { // 确保索引合法
                    freq[idx]++;
                }
            }
        }
    }
    // 构建Huffman树和编码表
    void buildCodeTable(int freq[]) {
        // 初始化编码表为空
        for (int i = 0; i < 26; i++) {
            codeTable[i].fromString("");
        }
        huffTree.buildHuffTree(freq);
        generateCodes(huffTree.getRoot(), "");
    }
    // 对单个单词编码
    HuffCode encodeWord(const string& word) {
        HuffCode result;
        string fullCode;
        for (Rank i = 0; i < word.size(); i++) {
            char c = word[i];
            if (!isalpha(c)) continue; // 过滤非字母字符
            c = tolower(c);
            int idx = c - 'a';
            if (idx < 0 || idx >= 26) continue; // 过滤非法索引
            // 获取编码（单个字母时默认编码"0"）
            string code = codeTable[idx].toString();
            if (code.empty() && freq[idx] > 0) {
                code = "0";
            }
            fullCode += code;
        }
        result.fromString(fullCode);
        return result;
    }
private:
    // 递归生成编码（左子树0，右子树1）
    void generateCodes(BinNode<int>* node, string currentCode) {
        if (!node) return;
        // 叶子节点：绑定编码到对应字母
        if (!node->left && !node->right) {
            int idx = node->charIdx;
            if (idx >= 0 && idx < 26) { // 确保索引合法
                codeTable[idx].fromString(currentCode);
            }
            return;
        }
        // 递归遍历左右子树
        generateCodes(node->left, currentCode + "0");
        generateCodes(node->right, currentCode + "1");
    }
};

// 主函数（修改为从TXT文档读取文本）
int main() {
    // -------------------------- 硬编码文件路径（核心修改）--------------------------
    string speech;       // 存储从文件读取的文本
    // 直接写入你的绝对路径（注意：Windows系统中路径分隔符用双反斜杠 \\，或单斜杠 /）
    string filePath = "D:\\code\\数据结构实验\\实验2\\exp2\\exp2\\x64\\Debug\\i have a dream.txt";
    ifstream file;       // 文件输入流

    // 打开文件并检查是否成功
    file.open(filePath, ios::in); // ios::in 表示只读模式
    if (!file.is_open()) {
        cerr << "错误：无法打开文件！请检查路径是否正确。" << endl;
        system("pause");
        return 1; // 异常退出
    }

    // 读取文件全部内容到 speech 中（逐行读取，保留完整文本）
    string line;
    while (getline(file, line)) {
        speech += line + " "; // 每行末尾加空格，避免单词拼接
    }
    file.close(); // 读取完成后关闭文件，释放资源
    // -------------------------- 硬编码结束 --------------------------

    // 初始化编码器
    HuffEncoder encoder;
    int freq[26];

    // 统计频率（使用从文件读取的 speech）
    encoder.countFrequency(speech, freq);
    encoder.countFrequency(speech); // 缓存到内部

    // 构建编码表
    encoder.buildCodeTable(freq);

    // 测试编码（可修改 testWords 为需要测试的单词）
    vector<string> testWords = { "dream", "equal", "nation", "freedom", "a", "abc" };
    for (const string& word : testWords) {
        HuffCode code = encoder.encodeWord(word);
        cout << "单词 '" << word << "' 的Huffman编码：" << code.toString() << endl;
    }

    // 防止控制台一闪而过（Visual Studio专用）
    system("pause");
    return 0;
}