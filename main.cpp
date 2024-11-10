#include <algorithm>
#include <chrono>
#include <cstdint>  // For uint16_t
#include <fstream>
#include <iostream>
#include <iterator>
#include <numeric>
#include <queue>
#include <string>

using std::cout;
using std::endl;
using std::max;
using std::min;
using std::pair;
using std::queue;

// 方向数组，用于移动到上下左右四个相邻的格子
constexpr int directions[4][2] = {{1, 0}, {-1, 0}, {0, 1}, {0, -1}};

// 将二维坐标转换为一维索引的辅助函数
inline int to1D(const int x, const int y, const int cols) { return x * cols + y; }

// 模板化的 BFS 搜索来找到数堆的行列范围
template <typename T>
pair<pair<int, int>, pair<int, int>> bfs(T *grid, int rows, int cols,
                                         int startX, int startY) {
    int minRow = startX, maxRow = startX, minCol = startY, maxCol = startY;
    queue<pair<int, int>> q;
    q.emplace(startX, startY);
    grid[to1D(startX, startY, cols)] = 0;  // 标记为已访问

    while (!q.empty()) {
        auto [x, y] = q.front();
        q.pop();

        // 更新当前数堆的边界
        minRow = min(minRow, x);
        maxRow = max(maxRow, x);
        minCol = min(minCol, y);
        maxCol = max(maxCol, y);

        // 遍历四个方向
        for (const auto &direction : directions) {
            int nx = x + direction[0], ny = y + direction[1];
            // 检查新位置是否在边界内且是非零元素
            if (nx >= 0 && nx < rows && ny >= 0 && ny < cols &&
                grid[to1D(nx, ny, cols)] != 0) {
                grid[to1D(nx, ny, cols)] = 0;  // 标记为已访问
                q.emplace(nx, ny);
            }
        }
    }

    return {{minRow, maxRow}, {minCol, maxCol}};
}

// 模板化的函数来寻找数堆的行列范围
template <typename T>
void findClusterRanges(T *grid, int rows, int cols, float *result) {
    const bool *visited = new bool[rows * cols]();  // 用于标记已经访问过的格子

    // // 遍历整个数组
    // for (int i = 0; i < rows; i++) {
    //     for (int j = 0; j < cols; j++) {
    //         if (grid[to1D(i, j, cols)] != 0) {             // 找到一个新的数堆
    //             auto range = bfs(grid, rows, cols, i, j);  // 使用模板化的 BFS
    //             计算数堆的行列范围 cout << "行范围: (" << range.first.first <<
    //             ", " << range.first.second
    //                  << "), \t列范围: (" << range.second.first << ", "
    //                  << range.second.second << ")";
    //             // 输出行列宽度
    //             cout << "    \t行宽度: " << range.first.second -
    //             range.first.first + 1
    //                  << ", 列宽度: " << range.second.second -
    //                  range.second.first + 1;
    //             // 输出中心
    //             cout << "  \t中心: (" << (range.first.first +
    //             range.first.second) / 2 << ", "
    //                  << (range.second.first + range.second.second) / 2 <<
    //                  ")\n";
    //         }
    //     }
    // }

    int numClusters = 0;  // 记录数堆的数量
    // 传参形式返回 包括行列范围，行列宽度，中心坐标，指针传递
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            if (grid[to1D(i, j, cols)] != 0) {  // 找到一个新的数堆
                auto range =
                    bfs(grid, rows, cols, i, j);  // 使用模板化的 BFS 计算数堆的行列范围
                result[numClusters * 8 + 0] = range.first.first;
                result[numClusters * 8 + 1] = range.first.second;
                result[numClusters * 8 + 2] = range.second.first;
                result[numClusters * 8 + 3] = range.second.second;
                result[numClusters * 8 + 4] =
                    range.first.second - range.first.first + 1;
                result[numClusters * 8 + 5] =
                    range.second.second - range.second.first + 1;
                result[numClusters * 8 + 6] =
                    (range.first.first + range.first.second) / 2;
                result[numClusters * 8 + 7] =
                    (range.second.first + range.second.second) / 2;
                numClusters++;
            }
        }
    }

    // 释放内存

    delete[] visited;
    visited = nullptr;
}

template <typename T>
void cropImage(const T *inputImage, int rows, int cols, int centerX,
               int centerY, int rowWidth, int colWidth, T *outputImage) {
    // 计算裁剪区域的起始和结束行列
    const int minRow = std::max(centerX - rowWidth, 0);
    const int maxRow = std::min(centerX + rowWidth, rows - 1);
    const int minCol = std::max(centerY - colWidth, 0);
    const int maxCol = std::min(centerY + colWidth, cols - 1);

    // 计算裁剪区域的尺寸
    const int croppedRows = maxRow - minRow + 1;
    const int croppedCols = maxCol - minCol + 1;

    // 将裁剪区域的像素复制到输出图像
    for (int i = 0; i < croppedRows; ++i) {
        for (int j = 0; j < croppedCols; ++j) {
            outputImage[i * croppedCols + j] =
                inputImage[(minRow + i) * cols + (minCol + j)];
        }
    }

    // 输出裁剪区域的尺寸
    // cout << "裁剪区域的尺寸: " << croppedRows << " x " << croppedCols << endl;

    // 存储裁剪后的图像 txt 文件 uint16_t数据
    std::ofstream file("croppedImage.txt");
    if (!file) {
        std::cerr << "无法打开文件: croppedImage.txt" << std::endl;
        return;
    }

    // 将裁剪后的图像数据写入文件
    for (int i = 0; i < croppedRows; ++i) {
        for (int j = 0; j < croppedCols; ++j) {
            file << outputImage[i * croppedCols + j] << " ";
        }
        file << std::endl;
    }

    file.close();
}

template <typename T>
bool readRawImage(const std::string &filename, T *imageData, int imgHeight,
                  int imgWidth, long long skipBytes = 0) {
    // 打开文件，二进制读模式
    std::ifstream file(filename, std::ios::binary);
    if (!file) {
        std::cerr << "无法打开文件: " << filename << std::endl;
        return false;
    }

    // 计算每张图像的字节数（16位图像，每像素2字节）
    const int imageSize = imgHeight * imgWidth;

    // 跳过指定的字节
    file.seekg(skipBytes, std::ios::beg);

    // 读取图像数据
    file.read(reinterpret_cast<char *>(imageData), imageSize * sizeof(T));

    // 检查是否成功读取到图像数据
    if (!file) {
        std::cerr << "读取图像数据失败，可能是文件大小不够。" << std::endl;
        return false;
    }

    file.close();
    return true;
}

template <typename T>
void rotateImage90CounterClockwise(T *imageData, const int imgHeight,
                                   const int imgWidth) {
    // 创建新的数组来存储旋转后的图像
    T *rotatedData = new T[imgHeight * imgWidth];

    // 对每个像素进行重新排列
    for (int row = 0; row < imgHeight; ++row) {
        for (int col = 0; col < imgWidth; ++col) {
            // 计算旋转后的新位置
            const int newRow = imgWidth - col - 1;
            const int newCol = row;
            rotatedData[newRow * imgHeight + newCol] =
                imageData[row * imgWidth + col];
        }
    }

    // 将旋转后的数据复制回原始数组
    std::copy(rotatedData, rotatedData + (imgHeight * imgWidth), imageData);

    // 释放旋转数据的临时数组
    delete[] rotatedData;
    rotatedData = nullptr;
}

int main0() {
    const auto start = std::chrono::high_resolution_clock::now();

    const int rows = 2176;
    // 设定数组的行数和列数
    const int cols = 2176;

    // 读取原始图像数据
    auto *imageData = new unsigned short[rows * cols];
    unsigned short *croppedImage = nullptr;
    float *result = nullptr;
    const std::string filename = "2176x2176x599_skip2_normal.raw";
    for (int i = 0; i < 597; i++) {
        // 读取原始图像数据
        readRawImage(filename, imageData, rows, cols, 1LL * (2 + i) * rows * cols * 2);

        rotateImage90CounterClockwise(imageData, rows, cols);

        // 计算最小值的索引
        unsigned short *minValIt =
            std::min_element(imageData, imageData + rows * cols);

        unsigned short minVal = *minValIt;
        const int minIndex = std::distance(imageData, minValIt);

        const int minRowIndex = minIndex / cols;
        const int minColIndex = minIndex % cols;

        constexpr int rowWidth = 9999;
        // cout << "最小值的行列索引: (" << minRowIndex << ", " << minColIndex <<
        // ")" << endl;
        constexpr int colWidth = 100;

        // 裁剪图像
        // 计算裁剪区域的起始和结束行列
        const int minRow = std::max(minRowIndex - rowWidth, 0);
        const int maxRow = std::min(minRowIndex + rowWidth, rows - 1);
        const int minCol = std::max(minColIndex - colWidth, 0);
        const int maxCol = std::min(minColIndex + colWidth, cols - 1);

        // 计算裁剪区域的尺寸
        const int croppedRows = maxRow - minRow + 1;
        const int croppedCols = maxCol - minCol + 1;
        // 裁剪图像
        if (croppedImage == nullptr) {
            croppedImage = new unsigned short[croppedRows * croppedCols];
        }

        // cout << "裁剪图像的尺寸: " << croppedRows << " x " << croppedCols <<
        // endl;
        cropImage(imageData, rows, cols, minRowIndex, minColIndex, rowWidth,
                  colWidth, croppedImage);

        // 计算平均值
        const double avgVal = std::accumulate(imageData, imageData + rows * cols, 0.0) /
                              (rows * cols);

        // cout << "图像最小值: " << minVal << endl;
        // cout << "图像平均值: " << avgVal << endl;

        // 阈值
        const double threshold = minVal + 0.2 * (avgVal - minVal);
        // 二值化图像 低于阈值的像素设为1，高于阈值的像素设为0

        for (int i = 0; i < croppedRows; i++) {
            for (int j = 0; j < croppedCols; j++) {
                croppedImage[to1D(i, j, croppedCols)] =
                    croppedImage[to1D(i, j, croppedCols)] < threshold ? 1 : 0;
            }
        }

        // std::fstream file("croppedImageBin.txt", std::ios::out);
        // if (!file) {
        //     std::cerr << "无法打开文件: croppedImageBin.txt" << std::endl;
        //     return 1;
        // }

        // // 将二值化图像数据写入文件
        // for (int i = 0; i < croppedRows; ++i) {
        //     for (int j = 0; j < croppedCols; ++j) {
        //         file << croppedImage[i * croppedCols + j] << " ";
        //     }
        //     file << std::endl;
        // }

        // file.close();

        result = new float[30 * 8];
        findClusterRanges(croppedImage, croppedRows, croppedCols, result);

        // for (int i = 0; i < 13; i++) {
        //     cout << "行范围: (" << result[i * 8 + 0] + minRow << ", " << result[i * 8 + 1] + minRow
        //          << "), \t列范围: (" << result[i * 8 + 2] + minCol << ", "
        //          << result[i * 8 + 3] + minCol << ")";
        //     // 输出行列宽度
        //     cout << "    \t行宽度: " << result[i * 8 + 4] << ", 列宽度: " << result[i * 8 + 5];
        //     // 输出中心
        //     cout << "  \t中心: (" << result[i * 8 + 6] << ", " << result[i * 8 + 7] << ")\n";
        // }

        // for (int i = 0; i < rows; i++) {
        //     for (int j = 0; j < cols; j++) {
        //         imageData[to1D(i, j, cols)] = imageData[to1D(i, j, cols)] <
        //         threshold ? 1 : 0;
        //     }
        // }
        // findClusterRanges(imageData, rows, cols);
        // 找到数堆的行列范围

        // 双排分割线
        std::cout << "-------------------------------------------\n";
        std::cout << "-------------------------------------------\n";
    }

    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed = end - start;
    std::cout << "程序运行时间: " << elapsed.count() << " 秒" << std::endl;
    // 释放内存
    delete[] imageData;
    delete[] croppedImage;
    delete[] result;

    imageData = nullptr;
    croppedImage = nullptr;

    return 0;
}
int main1() {
    // cut10有10张图像将cut10.raw复制60份在后面追加写入，文件名为cut10_60.raw
    std::ifstream file("cut10.raw", std::ios::binary);
    if (!file) {
        std::cerr << "无法打开文件: cut10.raw" << std::endl;
        return 1;
    }
    std::ofstream file_out("cut10_60.raw", std::ios::binary);
    if (!file_out) {
        std::cerr << "无法打开文件: cut10_60.raw" << std::endl;
        return 1;
    }
    // 读取原始图像数据
    unsigned short *imageData = new unsigned short[2176 * 2176 * 10];
    file.read(reinterpret_cast<char *>(imageData),
              2176 * 2176 * 10 * sizeof(unsigned short));
    // 写入60份
    for (int i = 0; i < 60; i++) {
        file_out.write(reinterpret_cast<char *>(imageData),
                       2176 * 2176 * 10 * sizeof(unsigned short));
    }
    file.close();
    file_out.close();
    delete[] imageData;
    imageData = nullptr;
    return 0;
}

int main() {
    main0();
    return 0;
}
