#include <vector>
#include <cmath>
#include "math.hpp"
#include "kernel.hpp"
#include <iostream>

double _kernelCoreFunction(double r, double alpha, int coreId) {
    if (r > 1) {
        return 0;
    }
    return std::pow(4*r*(1-r), alpha);
    switch(coreId) {
        case 0: // bump4
            if (r > 1) {
                return 0;
            }
            return std::exp(alpha - alpha/4/r/(1-r));
        case 1: // quad4
            if (r > 1) {
                return 0;
            }
            return std::pow(4*r*(1-r), alpha);
    }

    return 0;
}

double _kernelFunction(double r, int layerId, std::vector<double> b, int bDiv, double alpha, int coreId) {
    r = std::abs(r);
    if (layerId == 0) {
        return _kernelCoreFunction(r, alpha, coreId);
    } else {
        if (r >= 1) {
            return 0;
        }

        double R1 = r * (layerId + 1);
        double B1 = std::floor(R1);
        return _kernelCoreFunction(std::fmod(R1, 1.0), alpha, coreId) * b[B1];
    }
}

std::vector<std::vector<double>> buildKernel(int size, int nSize, int coreId, int layerId, std::vector<double> b, double alpha) {
    std::vector<std::vector<double>> kernel(size, std::vector<double>(size));

    double weight = 0;

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            int ii = (i + size / 2) % size - size / 2;
            int jj = (j + size / 2) % size - size / 2;
            double r = std::sqrt(ii * ii + jj * jj) / nSize;
            double v = _kernelFunction(r, layerId, b, 12, alpha, coreId);
            kernel[i][j] = v;
            weight += v;
        }
    }

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            kernel[i][j] /= weight;
        }
    }
    
    return roll(kernel, 1, 1);
}
