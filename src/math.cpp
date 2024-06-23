#include "math.hpp"
#include <cmath>
#include <vector>
#include <iostream>
#include <iostream>
#include <fftw3.h>
#include <complex>

fftw_complex* _kernel_fft;
fftw_complex* _grid_fft;
fftw_complex* _input_fft;
fftw_plan _gridPlan;
fftw_plan _inversePlan;

void setupConvolution(std::vector<std::vector<double>> kernel, int inputWidth, int inputHeight) {
    int kernelWidth = kernel.size();
    int kernelHeight = kernel[0].size();

    fftw_complex* kernelIn = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * inputWidth * inputHeight);
    _kernel_fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * inputWidth * inputHeight);
    fftw_plan kernelPlan = fftw_plan_dft_2d(inputWidth, inputHeight, kernelIn, _kernel_fft, FFTW_FORWARD, FFTW_PATIENT);

    for (int x = 0; x < inputWidth; x++) {
        for (int y = 0; y < inputHeight; y++) {
            kernelIn[x * inputHeight + y][0] = kernel[x][y];
            kernelIn[x * inputHeight + y][1] = 0;
        }
    }

    fftw_execute(kernelPlan);
    fftw_destroy_plan(kernelPlan);
    fftw_free(kernelIn);

    _grid_fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * inputWidth * inputHeight);
    _input_fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * inputWidth * inputHeight);

    _gridPlan = fftw_plan_dft_2d(inputWidth, inputHeight, _input_fft, _grid_fft, FFTW_FORWARD, FFTW_PATIENT);
    _inversePlan = fftw_plan_dft_2d(inputWidth, inputHeight, _grid_fft, _input_fft, FFTW_BACKWARD, FFTW_PATIENT);
}

void cleanupMath() {
    fftw_destroy_plan(_gridPlan);
    fftw_destroy_plan(_inversePlan);
    fftw_free(_grid_fft);
    fftw_free(_input_fft);
    fftw_free(_kernel_fft);
}

double bell(double x, double m, double s) {
    return std::exp(-std::pow((x-m) / s, 2) / 2);
}

int wrap(int x, int n) {
    return (x % n + n) % n;
}

double polynomial(double x, double m, double s, double a) {
    double r = std::abs(x - m);
    double r2 = r * r;
    double k2 = 9 * s*s;
    return r2 > k2 ? -1 : std::pow(1 - r2 / k2, a) * 2 - 1;
}

double wrapD(double x, double n) {
    return std::fmod(std::fmod(x, n) + n, n);
}

double clip(double x, double a, double b) {
    return x < a ? a : (x > b ? b : x);
}

std::vector<std::vector<double>> ogrid(int xStart, int xEnd, int yStart, int yEnd) {
    int xSize = std::abs(xEnd - xStart);
    int ySize = std::abs(yEnd - yStart);
    std::vector<std::vector<double>> grid(xSize, std::vector<double>(ySize));

    for (int i = xStart; i < xEnd; i++) {
        grid[0].push_back(i);
    }

    for (int i = yStart; i < yEnd; i++) {
        grid[1].push_back(i);
    }

    return grid;
}

std::vector<std::vector<double>> fftShift(std::vector<std::vector<double>> input) {
    std::vector<std::vector<double>> output(input.size(), std::vector<double>(input[0].size()));
    // numpy fftshift equivalent
    for (int i = 0; i < input.size(); i++) {
        for (int j = 0; j < input[0].size(); j++) {
            output[i][j] = input[wrap(i + input.size() / 2, input.size())][wrap(j + input[0].size() / 2, input[0].size())];
        }
    }
    return output;
}

std::vector<std::vector<double>> roll(std::vector<std::vector<double>> grid, int start_row, int start_col) {
    int rows = grid.size();
    int cols = grid[0].size();

    std::vector<std::vector<double>> newGrid(rows, std::vector<double>(cols));

    // Perform the roll operation
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            // Calculate the new indices after rolling
            int new_row = (i - start_row + rows) % rows;
            int new_col = (j - start_col + cols) % cols;

            // Assign the value to the new position
            newGrid[new_row][new_col] = grid[i][j];
        }
    }

    return newGrid;
}
std::vector<std::vector<double>> fastConv2D(std::vector<std::vector<double>> grid) {
    int inputWidth = grid.size();
    int inputHeight = grid[0].size();

    for (int x = 0; x < inputWidth; x++) {
        for (int y = 0; y < inputHeight; y++) {
            _input_fft[x * inputHeight + y][0] = grid[x][y];
            _input_fft[x * inputHeight + y][1] = 0;
        }
    }

    fftw_execute(_gridPlan);

    for (int x = 0; x < inputWidth; x++) {
        for (int y = 0; y < inputHeight; y++) {
            double real = _grid_fft[x * inputHeight + y][0] * _kernel_fft[x * inputHeight + y][0] - _grid_fft[x * inputHeight + y][1] * _kernel_fft[x * inputHeight + y][1];
            double imag = _grid_fft[x * inputHeight + y][0] * _kernel_fft[x * inputHeight + y][1] + _grid_fft[x * inputHeight + y][1] * _kernel_fft[x * inputHeight + y][0];
            _grid_fft[x * inputHeight + y][0] = real;
            _grid_fft[x * inputHeight + y][1] = imag;
        }
    }

    fftw_execute(_inversePlan);

    std::vector<std::vector<double>> output(inputWidth, std::vector<double>(inputHeight));
    for (int x = 0; x < inputWidth; x++) {
        for (int y = 0; y < inputHeight; y++) {
            output[x][y] = _input_fft[x * inputHeight + y][0] / (inputWidth * inputHeight);
        }
    }

    return roll(output, -1, -1);
}

double* convertToDoublePointer(std::vector<std::vector<double>>& vec, int subVecIndex) {
    double* arr = new double[vec.size()];
    for (int i = 0; i < vec.size(); i++) {
        arr[i] = vec[i][subVecIndex];
    }
    return arr;
}