#pragma once

#include <algorithm>
#include <fstream>
#include <iostream>
#include <numeric>
#include <string>

/**
 * brief 图像属性类(模板类)
 * imageBuffer: 图像数据缓冲区
 * rows: 图像行数
 * cols: 图像列数
 * views: 视图数
 * pixelSize: 图像像素大小
 *
 */

template <typename T> class ImageData {
  ~ImageData() {
    if (imageBuffer != nullptr) {
      delete[] imageBuffer;
      imageBuffer = nullptr;
    }
  }

private:
  T *imageBuffer = nullptr;
  T *imageData = nullptr;
  int rows = 0;
  int cols = 0;
  int views = 0;
  float pixelSize = 0.0;

public:
  void setImageBuffer(T *imageData) { this->imageData = imageData; }
  void setRows(int rows, int cols, int views) {
    this->rows = rows;
    this->cols = cols;
  }
  void setPixelSize(float pixelSize) { this->pixelSize = pixelSize; }

  T *getImageBuffer() { return imageBuffer; }
  int getRows() { return rows; }
  int getCols() { return cols; }
  float getPixelSize() { return pixelSize; }

public:
  T *LoadImage(int viewIndex);
};

template <typename T> T *ImageData<T>::LoadImage(int viewIndex) {
  // 读取原始图像数据
  unsigned short *imageData = new T[rows * cols];
  std::memcpy(imageBuffer, imageData + viewIndex * rows * cols,
              rows * cols * sizeof(T));
}
