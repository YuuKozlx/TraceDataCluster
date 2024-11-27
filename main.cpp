#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>  // For uint16_t
#include <cstring>
#include <fstream>
#include <iostream>
#include <iterator>
#include <numeric>
#include <queue>
#include <string>
#include <vector>

using std::cout;
using std::endl;
using std::max;
using std::min;
using std::pair;
using std::queue;
using std::sqrt;
using std::string;
using std::vector;

// 方向数组，用于移动到上下左右四个相邻的格子
constexpr int directions[4][2] = {{1, 0}, {-1, 0}, {0, 1}, {0, -1}};

// 将二维坐标转换为一维索引的辅助函数
inline int to1D(const int x, const int y, const int cols) { return x * cols + y; }

template <typename T>
pair<pair<int, int>, pair<int, int>> bfs(T* grid, int rows, int cols,
                                         int startX, int startY,
                                         bool* visited) {
    int minRow = startX, maxRow = startX, minCol = startY, maxCol = startY;
    queue<pair<int, int>> q;
    q.emplace(startX, startY);
    visited[to1D(startX, startY, cols)] = true;  // 标记为已访问

    while (!q.empty()) {
        auto [x, y] = q.front();
        q.pop();

        // 更新当前数堆的边界
        minRow = min(minRow, x);
        maxRow = max(maxRow, x);
        minCol = min(minCol, y);
        maxCol = max(maxCol, y);

        // 遍历四个方向
        for (const auto& direction : directions) {
            int nx = x + direction[0], ny = y + direction[1];
            // 检查新位置是否在边界内且未访问且是非零元素
            if (nx >= 0 && nx < rows && ny >= 0 && ny < cols &&
                !visited[to1D(nx, ny, cols)] && grid[to1D(nx, ny, cols)] != 0) {
                visited[to1D(nx, ny, cols)] = true;  // 标记为已访问
                q.emplace(nx, ny);
            }
        }
    }

    return {{minRow, maxRow}, {minCol, maxCol}};
}

template <typename T>
void findClusterRanges(T* grid, int rows, int cols, float* result) {
    bool* visited = new bool[rows * cols]();  // 初始化访问标记数组，全为 false

    int numClusters = 0;  // 记录数堆的数量

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            if (!visited[to1D(i, j, cols)] && grid[to1D(i, j, cols)] != 0) {  // 找到一个新的数堆
                auto range = bfs(grid, rows, cols, i, j, visited);            // 计算数堆的行列范围

                result[numClusters * 8 + 0] = range.first.first;
                result[numClusters * 8 + 1] = range.first.second;
                result[numClusters * 8 + 2] = range.second.first;
                result[numClusters * 8 + 3] = range.second.second;
                result[numClusters * 8 + 4] =
                    range.first.second - range.first.first + 1;
                result[numClusters * 8 + 5] =
                    range.second.second - range.second.first + 1;

                // 计算质心
                for (int i = 0; i < numClusters; i++) {
                    float rowWeightedSum = 0, rowWeightTotal = 0;
                    float colWeightedSum = 0, colWeightTotal = 0;

                    int rowStart = max(0, (int)result[i * 8 + 0] - (int)result[i * 8 + 4] / 2);
                    int rowEnd = min(rows - 1, (int)result[i * 8 + 1] + (int)result[i * 8 + 4] / 2);
                    int colStart = max(0, (int)result[i * 8 + 2] - (int)result[i * 8 + 5] / 2);
                    int colEnd = min(cols - 1, (int)result[i * 8 + 3] + (int)result[i * 8 + 5] / 2);

                    // 遍历裁剪范围内的所有像素
                    for (int r = rowStart; r <= rowEnd; ++r) {
                        for (int c = colStart; c <= colEnd; ++c) {
                            float value = grid[to1D(r, c, cols)];
                            if (value > 1e-6) {
                                rowWeightedSum += r * value;
                                rowWeightTotal += value;

                                colWeightedSum += c * value;
                                colWeightTotal += value;
                            }
                        }
                    }

                    // 计算行质心和列质心，避免除以零的情况
                    result[i * 8 + 6] = rowWeightTotal > 1e-6 ? rowWeightedSum / rowWeightTotal : -1.0f;
                    result[i * 8 + 7] = colWeightTotal > 1e-6 ? colWeightedSum / colWeightTotal : -1.0f;
                }

                numClusters++;  // 更新数堆计数
            }
        }
    }

    cout << "清理前数堆的数量: " << numClusters << endl;
    // 清理异常数堆，若数堆的行宽度或列宽度小于10，则认为是异常数堆;若数堆最小行 < 10 或最大行 > rows - 10 或最小列 < 10 或最大列 > cols - 10 也认为是异常数堆

    for (int i = 0; i < numClusters; i++) {
        if (result[i * 8 + 4] < 10 || result[i * 8 + 5] < 10 || result[i * 8 + 0] < 10 || result[i * 8 + 1] > rows - 10 || result[i * 8 + 2] < 10 || result[i * 8 + 3] > cols - 10) {
            result[i * 8 + 0] = -1;
            result[i * 8 + 1] = -1;
            result[i * 8 + 2] = -1;
            result[i * 8 + 3] = -1;
            result[i * 8 + 4] = -1;
            result[i * 8 + 5] = -1;
            result[i * 8 + 6] = -1;
            result[i * 8 + 7] = -1;
        }
    }

    // 将有效簇移到前面
    int j = 0;
    for (int i = 0; i < numClusters; i++) {
        if (result[i * 8 + 0] != -1) {
            if (i != j) {
                for (int k = 0; k < 8; k++) {
                    result[j * 8 + k] = result[i * 8 + k];
                }
            }
            j++;
        }
    }

    cout << "清理后数堆的数量: " << j << endl;

    // 释放内存
    delete[] visited;
}

// // dataCluster 只在竖直方向有重合，水平方向没有重合，找出所有的dataCluster
// template <typename T>
// void findClusterRanges(T* grid, int rows, int cols, float* result) {
//     int numClusters = 0;  // 记录簇的数量

//     // 逐行求和
//     int* sumRow = new int[rows]();
//     for (int i = 0; i < rows; i++) {
//         for (int j = 0; j < cols; j++) {
//             sumRow[i] += grid[to1D(i, j, cols)];
//         }
//     }

//     // 计算簇的行范围 由于数据是浮点数，所以需要设置一个阈值 1e-6
//     for (int i = 0; i < rows; i++) {
//         if (sumRow[i] > 1e-6) {
//             int j = i;
//             while (j < rows && sumRow[j] > 1e-6) {
//                 j++;
//             }
//             result[numClusters * 8 + 0] = i;
//             result[numClusters * 8 + 1] = j - 1;
//             numClusters++;
//             i = j;
//         }
//     }

//     // 求簇的列范围
//     for (int i = 0; i < numClusters; i++) {
//         int minCol = cols, maxCol = 0;
//         for (int j = 0; j < cols; j++) {
//             for (int k = result[i * 8 + 0]; k <= result[i * 8 + 1]; k++) {
//                 if (grid[to1D(k, j, cols)] > 1e-6) {
//                     minCol = min(minCol, j);
//                     maxCol = max(maxCol, j);
//                 }
//             }
//         }
//         result[i * 8 + 2] = minCol;
//         result[i * 8 + 3] = maxCol;
//         result[i * 8 + 4] = result[i * 8 + 1] - result[i * 8 + 0] + 1;
//         result[i * 8 + 5] = result[i * 8 + 3] - result[i * 8 + 2] + 1;
//     }

//     // 求行质心
//     for (int i = 0; i < numClusters; i++) {
//         float rowWeightedSum = 0, rowWeightTotal = 0;
//         float colWeightedSum = 0, colWeightTotal = 0;

//         int rowStart = max(0, (int)result[i * 8 + 0] - (int)result[i * 8 + 4] / 2);
//         int rowEnd = min(rows - 1, (int)result[i * 8 + 1] + (int)result[i * 8 + 4] / 2);
//         int colStart = max(0, (int)result[i * 8 + 2] - (int)result[i * 8 + 5] / 2);
//         int colEnd = min(cols - 1, (int)result[i * 8 + 3] + (int)result[i * 8 + 5] / 2);

//         // 遍历裁剪范围内的所有像素
//         for (int r = rowStart; r <= rowEnd; ++r) {
//             for (int c = colStart; c <= colEnd; ++c) {
//                 float value = grid[to1D(r, c, cols)];
//                 if (value > 1e-6) {
//                     rowWeightedSum += r * value;
//                     rowWeightTotal += value;

//                     colWeightedSum += c * value;
//                     colWeightTotal += value;
//                 }
//             }
//         }

//         // 计算行质心和列质心，避免除以零的情况
//         result[i * 8 + 6] = rowWeightTotal > 1e-6 ? rowWeightedSum / rowWeightTotal : -1.0f;
//         result[i * 8 + 7] = colWeightTotal > 1e-6 ? colWeightedSum / colWeightTotal : -1.0f;
//     }

//     // 清理异常簇，若簇的行宽度或列宽度小于10，则认为是异常簇;若簇最小行 < 10 或最大行 > rows - 10 或最小列 < 10 或最大列 > cols - 10 也认为是异常簇

//     for (int i = 0; i < numClusters; i++) {
//         if (result[i * 8 + 4] < 10 || result[i * 8 + 5] < 10 || result[i * 8 + 0] < 10 || result[i * 8 + 1] > rows - 10 || result[i * 8 + 2] < 10 || result[i * 8 + 3] > cols - 10) {
//             result[i * 8 + 0] = -1;
//             result[i * 8 + 1] = -1;
//             result[i * 8 + 2] = -1;
//             result[i * 8 + 3] = -1;
//             result[i * 8 + 4] = -1;
//             result[i * 8 + 5] = -1;
//             result[i * 8 + 6] = -1;
//             result[i * 8 + 7] = -1;
//         }
//     }

//     // 将有效簇移到前面
//     int j = 0;
//     for (int i = 0; i < numClusters; i++) {
//         if (result[i * 8 + 0] != -1) {
//             if (i != j) {
//                 for (int k = 0; k < 8; k++) {
//                     result[j * 8 + k] = result[i * 8 + k];
//                 }
//             }
//             j++;
//         }
//     }

//     // 释放内存
//     delete[] sumRow;
//     sumRow = nullptr;

//     // // 输出簇的中心
//     // for (int i = 0; i < numClusters; i++) {
//     //     cout << "簇的中心: (" << result[i * 8 + 6] << ", " << result[i * 8 + 7] << ")\n";
//     // }
// }

template <typename T>
void cropImage(const T* inputImage, int rows, int cols, int centerX,
               int centerY, int rowWidth, int colWidth, T* outputImage) {
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
    cout << "裁剪区域的尺寸: " << croppedRows << " x " << croppedCols << endl;

    // // 保存为raw float32 位图像 ，方便matlab读取
    // std::fstream file("croppedImage.raw", std::ios::out | std::ios::binary);
    // if (!file) {
    //     std::cerr << "无法打开文件: croppedImage.raw" << std::endl;
    //     return;
    // }

    // // 将裁剪区域的像素写入文件
    // for (int i = 0; i < croppedRows; ++i) {
    //     for (int j = 0; j < croppedCols; ++j) {
    //         file.write(reinterpret_cast<const char*>(&outputImage[i * croppedCols + j]), sizeof(T));
    //     }
    // }

    // file.close();
}

template <typename T>
bool readRawImage(const std::string& filename, T* imageData, int imgHeight,
                  int imgWidth, int imgViews = 1, long long skipBytes = 0) {
    // 打开文件，二进制读模式
    std::ifstream file(filename, std::ios::binary);
    if (!file) {
        std::cerr << "无法打开文件: " << filename << std::endl;
        return false;
    }

    // 计算每张图像的字节数（16位图像，每像素2字节）
    const int imageSize = imgHeight * imgWidth * imgViews;

    // 跳过指定的字节
    file.seekg(skipBytes, std::ios::beg);

    // 读取图像数据
    file.read(reinterpret_cast<char*>(imageData), imageSize * sizeof(T));

    // 检查是否成功读取到图像数据
    if (!file) {
        std::cerr << "读取图像数据失败，可能是文件大小不够。" << std::endl;
        return false;
    }

    file.close();
    return true;
}

template <typename T>
void rotateImage90CounterClockwise(T* imageData, const int imgHeight,
                                   const int imgWidth) {
    // 创建新的数组来存储旋转后的图像
    T* rotatedData = new T[imgHeight * imgWidth];

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

// 计数排序，用于选出第k大 或 第m小的数  排序对象为 unsigned short类型的数组
int countingSort(unsigned short* arr, int n, int k, int m, float& kthMax, float& mthMin) {
    // 找出数组中的最大值
    int maxVal = 65535;

    // 创建计数数组并初始化为0
    int* count = new int[maxVal + 1]();
    // 统计每个元素的出现次数
    for (int i = 0; i < n; ++i) {
        count[arr[i]]++;
    }

    // 从计数数组中找出第k大的数
    int sum = 0;
    for (int i = maxVal; i >= 0; --i) {
        sum += count[i];
        if (sum >= k) {
            kthMax = i;
            break;
        }
    }

    // 从计数数组中找出第m小的数

    sum = 0;
    for (int i = 0; i <= maxVal; ++i) {
        sum += count[i];
        if (sum >= m) {
            mthMin = i;
            break;
        }
    }

    // 释放内存
    delete[] count;
    count = nullptr;

    return 0;
}

// 根据前k大的数的值，以及前m小的数的值，对图像进行平滑处理，平滑方式为5x5的高斯滤波
int smoothImage(unsigned short* arr, int rows, int cols, float kthMax, float mthMin) {
    // 创建新的数组来存储平滑后的图像
    unsigned short* smoothedData = new unsigned short[rows * cols];

    // 备份原始数据
    std::copy(arr, arr + (rows * cols), smoothedData);

    // 高斯滤波核
    float kernel[5][5] = {
        {0.003, 0.013, 0.022, 0.013, 0.003},
        {0.013, 0.059, 0.097, 0.059, 0.013},
        {0.022, 0.097, 0.159, 0.097, 0.022},
        {0.013, 0.059, 0.097, 0.059, 0.013},
        {0.003, 0.013, 0.022, 0.013, 0.003}};

    // 对每个像素进行重新排列
    for (int row = 2; row < rows - 2; ++row) {
        for (int col = 2; col < cols - 2; ++col) {
            // 计算新像素值
            if (arr[row * cols + col] > kthMax || arr[row * cols + col] < mthMin) {
                float sum = 0;
                for (int i = -2; i <= 2; ++i) {
                    for (int j = -2; j <= 2; ++j) {
                        sum += kernel[i + 2][j + 2] * smoothedData[(row + i) * cols + (col + j)];
                    }
                }
                arr[row * cols + col] = static_cast<unsigned short>(sum);
            }
        }
    }

    // 释放平滑数据的临时数组
    delete[] smoothedData;
    smoothedData = nullptr;

    return 0;
}

// 对数化处理，用前面的函数求出的第k大的值进行对数化处理
int logTransform(unsigned short* arr, float* logData, int rows, int cols, float kthMax) {
    // 对每个像素进行重新排列
    float ivs = 1.0 / kthMax;
    for (int row = 0; row < rows; ++row) {
        for (int col = 0; col < cols; ++col) {
            // 计算新像素值
            logData[row * cols + col] = -log((1 + arr[row * cols + col]) * ivs);
        }
    }

    return 0;
}
// 输入为当前帧的第k个数堆的行列坐标，以及第二帧的全部行列坐标，以及第二帧的有效数堆数目，以及第二帧的平均距离。
int matchNeighbourView(std::vector<float>& matchRow, std::vector<float>& matchCol, float* rowCoordAt_cur_View, float* colCoordAt_cur_View, int effitiveNumClusterAt_cur_View, int avgDistAt_cur_View) {
    float dr = 0;
    float dc = 0;
    float dist2 = 0;
    for (int i = 0; i < effitiveNumClusterAt_cur_View; i++) {
        dr = rowCoordAt_cur_View[i] - matchRow.back();
        dc = colCoordAt_cur_View[i] - matchCol.back();
        dist2 = pow(dr, 2) + pow(dc, 2);
        if (dist2 < 0.25 * pow(avgDistAt_cur_View, 2)) {
            // 匹配成功
            matchRow.push_back(rowCoordAt_cur_View[i]);
            matchCol.push_back(colCoordAt_cur_View[i]);

            return 0;
        }
    }

    return -1;
}

int evalCenterCoord(float* arr, int maxNumCluster, int numViews, float* matchRowCoord, float* matchColCoord) {
    // 从arr中提取中心坐标  第6个为行坐标 第7个为列坐标
    // 每个view的大小为 30*8

    float* centerRowCoord = new float[maxNumCluster * numViews];
    float* centerColCoord = new float[maxNumCluster * numViews];

    for (int i = 0; i < numViews; i++) {
        for (int j = 0; j < maxNumCluster; j++) {
            centerRowCoord[i * maxNumCluster + j] = arr[i * 30 * 8 + j * 8 + 6];
            centerColCoord[i * maxNumCluster + j] = arr[i * 30 * 8 + j * 8 + 7];
        }
    }

    // 图像旋转的过程中，会出现新的数堆以及已有数堆消失的情况，导致相邻帧图像的数堆发生错位，需要对数堆进行重新匹配
    // 可根据当前帧除去首尾的数堆的平均距离来判断，相邻帧的数堆距离若小于该值的0.8倍，则认为是同一数堆
    // 当前帧的数堆匹配好后，对匹配后的数堆进行排序，按照数堆的中心坐标进行排序，从而得到数堆的匹配关系

    // 逐帧计算有效数堆，行列均为0的数堆不要
    int* effectiveNumCluster = new int[numViews]();
    for (int i = 0; i < numViews; i++) {
        for (int j = 0; j < maxNumCluster; j++) {
            if (centerRowCoord[i * maxNumCluster + j] > 1e-2 || centerColCoord[i * maxNumCluster + j] > 1e-2) {
                effectiveNumCluster[i]++;
            }
        }
    }

    // 计算每帧的数堆的平均距离, 除去首尾的数堆
    float* avgDist = new float[numViews]();
    for (int i = 0; i < numViews; i++) {
        float sumDist = 0;
        for (int j = 1; j < effectiveNumCluster[i] - 1; j++) {
            sumDist += sqrt(pow(centerRowCoord[i * maxNumCluster + j + 1] - centerRowCoord[i * maxNumCluster + j], 2) +
                            pow(centerColCoord[i * maxNumCluster + j + 1] - centerColCoord[i * maxNumCluster + j], 2));
        }
        avgDist[i] = sumDist / (effectiveNumCluster[i] - 3);
    }

    std::vector<float> matchRowTemp;
    std::vector<float> matchColTemp;
    int matchCount = 0;  // 存储匹配计数

    // 先求第一个小球的匹配关系验证
    for (int j = 0; j < effectiveNumCluster[0]; j++) {
        matchRowTemp.push_back(centerRowCoord[j]);
        matchColTemp.push_back(centerColCoord[j]);
        bool matched = true;
        for (int i = 1; i < numViews; i++) {
            if (matchNeighbourView(matchRowTemp, matchColTemp, centerRowCoord + i * maxNumCluster, centerColCoord + i * maxNumCluster, effectiveNumCluster[i], avgDist[i]) != 0) {
                cout << "匹配失败" << endl;
                matched = false;
                matchColTemp.clear();
                matchRowTemp.clear();

                break;  // 如果匹配失败，退出内层循环
            }
        }

        // 只有在所有视图都匹配成功时才计数
        if (matched) {
            cout << "匹配成功" << endl;
            // 保存匹配结果
            for (int i = 0; i < numViews; i++) {
                matchRowCoord[i * maxNumCluster + matchCount] = matchRowTemp[i];
                matchColCoord[i * maxNumCluster + matchCount] = matchColTemp[i];
            }
            matchCount++;
            matchColTemp.clear();
            matchRowTemp.clear();
        }
    }

    // 将结果保存至文件 matchResult.txt 列宽为匹配的数堆数目 行宽为视图数目
    std::ofstream file("matchResult.txt");
    if (!file) {
        std::cerr << "无法打开文件: matchResult.txt" << std::endl;
        return 1;
    }

    // 将结果写入文件  对应行列坐标在一行  先写行坐标 再写列坐标
    for (int i = 0; i < numViews; i++) {
        for (int j = 0; j < matchCount; j++) {
            file << matchRowCoord[i * maxNumCluster + j] << " " << matchColCoord[i * maxNumCluster + j] << " ";
        }
        file << std::endl;
    }

    file.close();

    // 释放内存
    delete[] centerRowCoord;
    delete[] centerColCoord;
    delete[] effectiveNumCluster;
    delete[] avgDist;

    return 0;
}

int main0() {
    const auto start = std::chrono::high_resolution_clock::now();

    const int rows = 2176;
    // 设定数组的行数和列数
    const int cols = 2176;
    const int numViews = 599;

    // 读取原始图像数据
    auto* imageData = new unsigned short[rows * cols];
    auto* logData = new float[rows * cols];
    float* croppedImage = nullptr;
    float* result = nullptr;
    float* allResult = new float[30 * 8 * 597];
    const std::string filename = "G:/Code/TraceDataCluster/2176x2176x599_skip2_normal.raw";

    for (int i = 0; i < 597; i++) {
        // 读取原始图像数据
        readRawImage(filename, imageData, rows, cols, 1, 1LL * (2 + i) * rows * cols * 2);

        rotateImage90CounterClockwise(imageData, rows, cols);

        // 计数排序，选出第k大的数
        float kthMax = 0;
        float mthMin = 0;
        countingSort(imageData, rows * cols, 1000, 100, kthMax, mthMin);

        // 平滑处理
        smoothImage(imageData, rows, cols, kthMax, mthMin);

        // 对数化处理
        logTransform(imageData, logData, rows, cols, kthMax);

        // 计算最大值的索引
        auto* maxValIt =
            std::max_element(logData, logData + rows * cols);

        auto maxVal = *maxValIt;
        const int maxIndex = std::distance(logData, maxValIt);

        const int maxRowIndex = maxIndex / cols;
        const int maxColIndex = maxIndex % cols;

        constexpr int rowWidth = 9999;
        cout << "最大值的行列索引: (" << maxRowIndex << ", " << maxColIndex << ")" << endl;
        constexpr int colWidth = 9999;

        // 裁剪图像
        // 计算裁剪区域的起始和结束行列
        const int minRow = std::max(maxRowIndex - rowWidth, 0);
        const int maxRow = std::min(maxRowIndex + rowWidth, rows - 1);
        const int minCol = std::max(maxColIndex - colWidth, 0);
        const int maxCol = std::min(maxColIndex + colWidth, cols - 1);

        // 计算裁剪区域的尺寸
        const int croppedRows = maxRow - minRow + 1;
        const int croppedCols = maxCol - minCol + 1;
        // 裁剪图像
        if (croppedImage == nullptr) {
            croppedImage = new float[croppedRows * croppedCols];
        }

        // cout << "裁剪图像的尺寸: " << croppedRows << " x " << croppedCols <<
        // endl;
        cropImage(logData, rows, cols, maxRowIndex, maxColIndex, rowWidth,
                  colWidth, croppedImage);

        // 计算平均值
        const double avgVal = std::accumulate(logData, logData + rows * cols, 0.0) /
                              (rows * cols);
        // 阈值
        const double threshold = maxVal + 0.5 * (avgVal - maxVal);
        //  低于阈值的像素值设为0，高于阈值的像素值保持不变

        for (int i = 0; i < croppedRows; i++) {
            for (int j = 0; j < croppedCols; j++) {
                croppedImage[to1D(i, j, croppedCols)] =
                    croppedImage[to1D(i, j, croppedCols)] < threshold ? 0 : croppedImage[to1D(i, j, croppedCols)];
            }
        }

        // //  若croppedImage.raw存在，则删除  重新写入
        // //  追加写入 croppedImage.raw
        // if (i == 0) {
        //     std::remove("croppedImage.raw");
        // }

        // // 保存为raw float32 位图像 ，方便matlab读取
        // std::fstream file("croppedImage.raw", std::ios::out | std::ios::app | std::ios::binary);

        // // 将裁剪区域的像素写入文件
        // for (int i = 0; i < croppedRows; ++i) {
        //     for (int j = 0; j < croppedCols; ++j) {
        //         file.write(reinterpret_cast<const char*>(&croppedImage[i * croppedCols + j]), sizeof(float));
        //     }
        // }

        // file.close();

        result = new float[30 * 8];
        findClusterRanges(croppedImage, croppedRows, croppedCols, result);

        // 对数堆的行列坐标进行修正
        for (int i = 0; i < 30; i++) {
            result[i * 8 + 0] += minRow;
            result[i * 8 + 1] += minRow;
            result[i * 8 + 2] += minCol;
            result[i * 8 + 3] += minCol;
            result[i * 8 + 6] += minRow;
            result[i * 8 + 7] += minCol;
        }

        memcpy(allResult + i * 30 * 8, result, 30 * 8 * sizeof(float));

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
        // std::cout << "-------------------------------------------\n";
        // std::cout << "-------------------------------------------\n";
    }

    // 将结果写入文件 只写入行列坐标

    std::ofstream file("allResult.txt");
    if (!file) {
        std::cerr << "无法打开文件: allResult.txt" << std::endl;
        return 1;
    }

    // 将结果写入文件  对应行列坐标在一行  先写行坐标 再写列坐标
    for (int i = 0; i < numViews - 2; i++) {
        for (int j = 0; j < 30; j++) {
            for (int k = 6; k < 8; k++) {
                file << allResult[i * 30 * 8 + j * 8 + k] << " ";
            }
        }
        file << std::endl;
    }
    file.close();

    // 保存所有数堆的行宽度和列宽度
    std::ofstream file2("allResult-width.txt");
    if (!file2) {
        std::cerr << "无法打开文件: allResult-width.txt" << std::endl;
        return 1;
    }

    // 将结果写入文件  对应行列坐标在一行  先写行坐标 再写列坐标
    for (int i = 0; i < numViews - 2; i++) {
        for (int j = 0; j < 30; j++) {
            for (int k = 4; k < 6; k++) {
                file2 << allResult[i * 30 * 8 + j * 8 + k] << " ";
            }
        }
        file2 << std::endl;
    }

    file2.close();

    // 匹配相邻帧的数堆

    float* matchRowCoord = new float[30 * 597];
    float* matchColCoord = new float[30 * 597];

    evalCenterCoord(allResult, 30, 597, matchRowCoord, matchColCoord);

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
    unsigned short* imageData = new unsigned short[2176 * 2176 * 10];
    file.read(reinterpret_cast<char*>(imageData),
              2176 * 2176 * 10 * sizeof(unsigned short));
    // 写入60份
    for (int i = 0; i < 60; i++) {
        file_out.write(reinterpret_cast<char*>(imageData),
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
