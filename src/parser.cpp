#include "parser.hpp"
#include "grid.hpp"
#include <string>
#include <regex>
#include <iostream>
#include <map>
#include <cmath>

std::map<std::string, int> coreName = {
    {"bump4", 0},
    {"quad4", 1},
    {"trap1/5", 2},
    {"stpz1/4", 3},
    {"life", 4}
};

std::map<std::string, int> deltaName = {
    {"gaus", 0},
    {"quad4", 1},
    {"trap", 2},
    {"stpz", 3},
};  

std::vector<std::string> split(std::string st, std::regex re) {
    std::sregex_token_iterator first{st.begin(), st.end(), re, -1}, last;
    return {first, last};
}

bool endsWith(std::string st, std::string suffix) {
    if (st.size() < suffix.size()) {
        return false;
    }
    return st.compare(st.size() - suffix.size(), suffix.size(), suffix) == 0;
}

bool startsWith(std::string st, std::string prefix) {
    if (st.size() < prefix.size()) {
        return false;
    }
    return st.compare(0, prefix.size(), prefix) == 0;
}

double UnSimpFrac(std::string st) {
    std::regex re("\\/");
    std::vector<std::string> parts = split(st, re);
    if (parts.size() == 2) {
        return std::stod(parts[0]) / std::stod(parts[1]);
    }
    return std::stod(st);
}

GridTemplate ParseRule(std::string st) {
    GridTemplate gt;
    std::regex re("\\=|\\;|\\*");
    std::vector<std::string> parts = split(st, re);
    if (parts.size() > 1) {
        gt.kernelRadius = std::stoi(parts[1]);
    }
    if (parts.size() > 6) {
        gt.limitValue = !endsWith(parts[6], "+");
        gt.updateFrequency = std::round(1 / std::stod(parts[6]));
    }
    if (parts.size() > 3) {
        re = std::regex("\\(|\\)|\\,");
        std::vector<std::string> parts3 = split(parts[3], re);
        re = std::regex("/");
        std::vector<std::string> parts3p = split(parts3[0].substr(4), re);

        bool hasE = static_cast<int>(parts[3].find(")(")) >= 0;
        gt.mu = -1;
        gt.sigma = -1;
        gt.layerId = -1;
        gt.coreId = -1;
        if (parts3[0] == "bimo4") { parts3[0] = "quad4"; }
        else if (parts3[0] == "bist4") { parts3[0] = "stpz1/4"; }
        else if (parts3[0] == "trmo4") { parts3[0] = "quad4"; gt.layerId = 2; }
        if (parts3p.size() >= 3) { gt.kernel_B[0] = UnSimpFrac(parts3p[1] + "/" + parts3p[2]); gt.layerId = 1; }
        gt.coreId = coreName[parts3[0]];
        if (gt.layerId == -1) {
            if (hasE) gt.layerId = std::max((static_cast<int>(parts3.size()) - 4) / 2, 0);
            else gt.layerId = std::max(static_cast<int>(parts3.size()) - 2, 0);
        }
        gt.kernel_B = std::vector<double>(4);
        gt.kernel_E = std::vector<double>(4);
        for (int i = 0; i < 4; i++) {
            if (i <= gt.layerId) {
                gt.kernel_B[i] = UnSimpFrac(parts3[i+1]);
            }
            if (hasE) {
                if (i <= gt.layerId) {
                    gt.kernel_E[i] = UnSimpFrac(parts3[i+gt.layerId]);
                }
            }
        }
        if (gt.layerId <= 1 && gt.kernel_B[0] == 0) {
            gt.kernel_B[0] = 1;
        }

        if (parts.size() > 5) {
            re = std::regex("\\(|\\)|\\,");
            std::vector<std::string> parts5 = split(parts[5], re);
            gt.mu = std::stod(parts5[1]);
            gt.sigma = std::stod(parts5[2]);
        }
    }

    std::cout << "Kernel radius: " << gt.kernelRadius << "\n";
    std::cout << "Update frequency: " << gt.updateFrequency << "\n";
    std::cout << "Limit value: " << gt.limitValue << "\n";
    std::cout << "Mu: " << gt.mu << "\n";
    std::cout << "Sigma: " << gt.sigma << "\n";
    std::cout << "Layer ID: " << gt.layerId << "\n";
    std::cout << "Core ID: " << gt.coreId << "\n";
    return gt;
}