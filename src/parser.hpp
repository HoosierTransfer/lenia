#pragma once
#include "grid.hpp"
#include <string>
#include <regex>
#include <vector>

GridTemplate ParseRule(std::string st);

std::vector<std::string> split(std::string st, std::regex re);

bool endsWith(std::string st, std::string suffix);

bool startsWith(std::string st, std::string prefix);