#include "iuwtdecomposition.h"

#include "../system/threadpool.h"

void IUWTDecomposition::DecomposeMT(ThreadPool& threadPool, const float* input,
                                    float* scratch, bool includeLargest) {
  Image& i1(_scales.back().Coefficients());
  i1 = Image(_width, _height);

  // The first iteration of the loop, unrolled, so that we don't have to
  // copy the input into i0.
  Image& coefficients0 = _scales[0].Coefficients();
  coefficients0 = Image(_width, _height);
  convolveMT(threadPool, i1.data(), input, scratch, _width, _height, 1);
  convolveMT(threadPool, coefficients0.data(), i1.data(), scratch, _width,
             _height, 1);

  // coefficients = i0 - i2
  differenceMT(threadPool, coefficients0.data(), input, coefficients0.data(),
               _width, _height);

  // i0 = i1;
  Image i0(i1);

  for (int scale = 1; scale != int(_scaleCount); ++scale) {
    Image& coefficients = _scales[scale].Coefficients();
    coefficients = Image(_width, _height);
    convolveMT(threadPool, i1.data(), i0.data(), scratch, _width, _height,
               scale + 1);
    convolveMT(threadPool, coefficients.data(), i1.data(), scratch, _width,
               _height, scale + 1);

    // coefficients = i0 - i2
    differenceMT(threadPool, coefficients.data(), i0.data(),
                 coefficients.data(), _width, _height);

    // i0 = i1;
    if (scale + 1 != int(_scaleCount))
      memcpy(i0.data(), i1.data(), sizeof(float) * _width * _height);
  }

  // The largest (residual) scales are in i1, but since the
  // largest scale is aliased to i1, it's already stored there.
  // Hence we can skip this:
  // if(includeLargest)
  //	_scales.back().Coefficients() = i1;

  // Do free the memory of the largest scale if it is not necessary:
  if (!includeLargest) _scales.back().Coefficients().reset();
}

void IUWTDecomposition::convolveMT(ThreadPool& threadPool, float* output,
                                   const float* image, float* scratch,
                                   size_t width, size_t height, int scale) {
  ConvolveHorizontalPartialFunc hFunc;
  hFunc._output = scratch;
  hFunc._image = image;
  hFunc._width = width;
  hFunc._scale = scale;
  for (size_t t = 0; t != threadPool.size(); ++t) {
    hFunc._startY = (height * t) / threadPool.size();
    hFunc._endY = (height * (t + 1)) / threadPool.size();
    threadPool.queue(hFunc);
  }
  threadPool.wait_for_all_tasks();

  ConvolveVerticalPartialFunc vFunc;
  vFunc._output = output;
  vFunc._image = scratch;
  vFunc._width = width;
  vFunc._height = height;
  vFunc._scale = scale;
  for (size_t t = 0; t != threadPool.size(); ++t) {
    vFunc._startX = (width * t) / threadPool.size();
    vFunc._endX = (width * (t + 1)) / threadPool.size();
    threadPool.queue(vFunc);
  }
  threadPool.wait_for_all_tasks();
}

void IUWTDecomposition::differenceMT(class ThreadPool& threadPool, float* dest,
                                     const float* lhs, const float* rhs,
                                     size_t width, size_t height) {
  DifferencePartialFunc func;
  func._dest = dest;
  func._lhs = lhs;
  func._rhs = rhs;
  func._width = width;
  for (size_t t = 0; t != threadPool.size(); ++t) {
    func._startY = (height * t) / threadPool.size();
    func._endY = (height * (t + 1)) / threadPool.size();
    threadPool.queue(func);
  }
  threadPool.wait_for_all_tasks();
}

void IUWTDecomposition::convolveHorizontalFast(float* output,
                                               const float* image, size_t width,
                                               size_t height, int scale) {
  const size_t H_SIZE = 5;
  const float h[H_SIZE] = {1.0 / 16.0, 4.0 / 16.0, 6.0 / 16.0, 4.0 / 16.0,
                           1.0 / 16.0};
  int scaleDist = (1 << scale);
  int dist[H_SIZE];
  size_t minX[H_SIZE], maxX[H_SIZE];
  for (int hIndex = 0; hIndex != H_SIZE; ++hIndex) {
    int hShift = hIndex - H_SIZE / 2;
    dist[hIndex] = (scaleDist - 1) * hShift;
    minX[hIndex] = std::max<int>(0, -dist[hIndex]);
    maxX[hIndex] = std::min<int>(width, width - dist[hIndex]);
  }
  for (size_t y = 0; y != height; ++y) {
    float* outputPtr = &output[y * width];
    const float* inputPtr = &image[y * width];

    for (size_t x = 0; x != minX[1]; ++x) {
      outputPtr[x] = inputPtr[x + dist[2]] * h[2] +
                     inputPtr[x + dist[3]] * h[3] +
                     inputPtr[x + dist[4]] * h[4];
    }

    for (size_t x = minX[1]; x != minX[0]; ++x) {
      outputPtr[x] =
          inputPtr[x + dist[2]] * h[2] + inputPtr[x + dist[1]] * h[1] +
          inputPtr[x + dist[3]] * h[3] + inputPtr[x + dist[4]] * h[4];
    }

    for (size_t x = minX[0]; x != maxX[4]; ++x) {
      outputPtr[x] =
          inputPtr[x + dist[2]] * h[2] + inputPtr[x + dist[1]] * h[1] +
          inputPtr[x + dist[0]] * h[0] + inputPtr[x + dist[3]] * h[3] +
          inputPtr[x + dist[4]] * h[4];
    }

    for (size_t x = maxX[4]; x != maxX[3]; ++x) {
      outputPtr[x] =
          inputPtr[x + dist[2]] * h[2] + inputPtr[x + dist[1]] * h[1] +
          inputPtr[x + dist[0]] * h[0] + inputPtr[x + dist[3]] * h[3];
    }

    for (size_t x = maxX[3]; x != width; ++x) {
      outputPtr[x] = inputPtr[x + dist[2]] * h[2] +
                     inputPtr[x + dist[1]] * h[1] +
                     inputPtr[x + dist[0]] * h[0];
    }
  }
}

// This version is not as fast as the one below.
void IUWTDecomposition::convolveVerticalPartialFastFailed(
    float* output, const float* image, size_t width, size_t height,
    size_t startX, size_t endX, int scale) {
  const size_t H_SIZE = 5;
  const float h[H_SIZE] = {1.0 / 16.0, 4.0 / 16.0, 6.0 / 16.0, 4.0 / 16.0,
                           1.0 / 16.0};
  int scaleDist = (1 << scale);
  for (size_t y = 0; y < height; ++y) {
    float* outputPtr = &output[y * width];
    const float* inputPtr = &image[y * width];
    for (size_t x = startX; x != endX; ++x) outputPtr[x] = inputPtr[x] * h[2];
  }
  for (int hIndex = 0; hIndex != H_SIZE; ++hIndex) {
    if (hIndex != 2) {
      int hShift = hIndex - H_SIZE / 2;
      int dist = (scaleDist - 1) * hShift;
      size_t minY = std::max<int>(0, -dist),
             maxY = std::min<int>(height, height - dist);
      for (size_t y = minY; y < maxY; ++y) {
        float* outputPtr = &output[y * width];
        const float* inputPtr = &image[(y + dist) * width];
        for (size_t x = startX; x != endX; ++x)
          outputPtr[x] += inputPtr[x] * h[hIndex];
      }
    }
  }
}

void IUWTDecomposition::convolveVerticalPartialFast(float* output,
                                                    const float* image,
                                                    size_t width, size_t height,
                                                    size_t startX, size_t endX,
                                                    int scale) {
  const size_t H_SIZE = 5;
  const float h[H_SIZE] = {1.0 / 16.0, 4.0 / 16.0, 6.0 / 16.0, 4.0 / 16.0,
                           1.0 / 16.0};
  int scaleDist = (1 << scale);
  int dist[H_SIZE];
  size_t minY[H_SIZE], maxY[H_SIZE];
  for (int hIndex = 0; hIndex != H_SIZE; ++hIndex) {
    int hShift = hIndex - H_SIZE / 2;
    dist[hIndex] = (scaleDist - 1) * hShift;
    minY[hIndex] = std::max<int>(0, -dist[hIndex]);
    maxY[hIndex] = std::min<int>(height, height - dist[hIndex]);
  }

  for (size_t y = 0; y != minY[1]; ++y) {
    float* outputPtr = &output[y * width];
    const float* inputPtr2 = &image[(y + dist[2]) * width];
    const float* inputPtr3 = &image[(y + dist[3]) * width];
    const float* inputPtr4 = &image[(y + dist[4]) * width];
    for (size_t x = startX; x != endX; ++x) {
      outputPtr[x] =
          inputPtr2[x] * h[2] + inputPtr3[x] * h[3] + inputPtr4[x] * h[4];
    }
  }

  for (size_t y = minY[1]; y != minY[0]; ++y) {
    float* outputPtr = &output[y * width];
    const float* inputPtr1 = &image[(y + dist[1]) * width];
    const float* inputPtr2 = &image[(y + dist[2]) * width];
    const float* inputPtr3 = &image[(y + dist[3]) * width];
    const float* inputPtr4 = &image[(y + dist[4]) * width];
    for (size_t x = startX; x != endX; ++x) {
      outputPtr[x] = inputPtr1[x] * h[1] + inputPtr2[x] * h[2] +
                     inputPtr3[x] * h[3] + inputPtr4[x] * h[4];
    }
  }

  for (size_t y = minY[0]; y != maxY[4]; ++y) {
    float* outputPtr = &output[y * width];
    const float* inputPtr0 = &image[(y + dist[0]) * width];
    const float* inputPtr1 = &image[(y + dist[1]) * width];
    const float* inputPtr2 = &image[(y + dist[2]) * width];
    const float* inputPtr3 = &image[(y + dist[3]) * width];
    const float* inputPtr4 = &image[(y + dist[4]) * width];
    for (size_t x = startX; x != endX; ++x) {
      outputPtr[x] = inputPtr0[x] * h[0] + inputPtr1[x] * h[1] +
                     inputPtr2[x] * h[2] + inputPtr3[x] * h[3] +
                     inputPtr4[x] * h[4];
    }
  }

  for (size_t y = maxY[4]; y != maxY[3]; ++y) {
    float* outputPtr = &output[y * width];
    const float* inputPtr0 = &image[(y + dist[0]) * width];
    const float* inputPtr1 = &image[(y + dist[1]) * width];
    const float* inputPtr2 = &image[(y + dist[2]) * width];
    const float* inputPtr3 = &image[(y + dist[3]) * width];
    for (size_t x = startX; x != endX; ++x) {
      outputPtr[x] = inputPtr0[x] * h[0] + inputPtr1[x] * h[1] +
                     inputPtr2[x] * h[2] + inputPtr3[x] * h[3];
    }
  }

  for (size_t y = maxY[3]; y != height; ++y) {
    float* outputPtr = &output[y * width];
    const float* inputPtr0 = &image[(y + dist[0]) * width];
    const float* inputPtr1 = &image[(y + dist[1]) * width];
    const float* inputPtr2 = &image[(y + dist[2]) * width];
    for (size_t x = startX; x != endX; ++x) {
      outputPtr[x] =
          inputPtr0[x] * h[0] + inputPtr1[x] * h[1] + inputPtr2[x] * h[2];
    }
  }
}
