#include "grid.hpp"
#include <SDL2/SDL.h>
#include <vector>
#include <iostream>
#include "math.hpp"
#include "kernel.hpp"
#include "implot.h"
#include "colormap.hpp"
#include <string>
#include <regex>
#include "parser.hpp"

const std::vector<std::string> statNames = {"Mass", "Growth", "Mass volume", "Growth volume", "Mass density", "Growth density", "Centroid speed", "Growth-centroid distance", "Centroid rotate speed", "Growth-centroid rotate speed", "Major axis rotate speed", "Moment of intertia", "Skewness", "Hu's 5", "Hu's 6", "Hu's 7", "Kurtosis", "Flusser's 8", "Flusser's 9", "Flusser's 10"};

const std::vector<std::vector<ImVec4>> colors = {
			{{13, 13, 13, 255}, {12, 12, 12, 255}, {11, 11, 11, 255}, {10, 10, 10, 255}, {9, 9, 9, 255}, {8, 8, 8, 255}, {7, 7, 7, 255}, {6, 6, 6, 255}, {5, 5, 5, 255}, {4, 4, 4, 255}, {1, 1, 1, 255}, /**/{0, 15, 15, 255}, {0, 9, 15, 255}, {7, 11, 15, 255}},  // B/W
			{{9, 15, 15, 255}, {6, 13, 15, 255}, {3, 11, 15, 255}, {0, 9, 15, 255}, {3, 8, 15, 255}, {6, 7, 15, 255}, {9, 6, 15, 255}, {11, 4, 15, 255}, {13, 2, 15, 255}, {15, 0, 15, 255}, {11, 0, 5, 255}, /**/{15, 15, 0, 255}, {15, 7, 0, 255}, {15, 11, 7, 255}},  // blue
			{{9, 15, 9, 255}, {7, 13, 7, 255}, {5, 11, 5, 255}, {3, 10, 3, 255}, {5, 10, 3, 255}, {8, 10, 1, 255}, {11, 8, 0, 255}, {13, 6, 3, 255}, {15, 4, 6, 255}, {15, 2, 9, 255}, {9, 0, 7, 255}, /**/{0, 15, 15, 255}, {0, 9, 15, 255}, {7, 11, 15, 255}}};


Grid::Grid(SDL_Renderer* renderer, int width, int height, GridTemplate template_)
{
    this->renderer = renderer;
    this->width = width;
    this->height = height;
    this->grid = std::vector<std::vector<double>>(width, std::vector<double>(height));
    this->delta = std::vector<std::vector<double>>(width, std::vector<double>(height));
    this->texture = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_RGBA8888, SDL_TEXTUREACCESS_STREAMING, width, height);
    this->textureBuffer = new uint32_t[width * height];
    this->kernel = buildKernel(width, template_.kernelRadius, template_.coreId, template_.layerId, template_.kernel_B, template_.alpha);
    this->updateFrequency = template_.updateFrequency;
    this->NS = template_.kernelRadius;
    this->stats = std::vector<std::vector<double>>(20, std::vector<double>(100));
    this->mu = template_.mu;
    this->sigma = template_.sigma;
    this->alpha = template_.alpha;
    this->statsI = 0;
    this->generation = 0;
    this->statsAverage = std::vector<double>(20);
    setupConvolution(kernel, width, height);
    resetCounters();
}

void Grid::zeroStats() {
    mX = 0;
    mY = 0;

    gX = 0;
    gY = 0;

    axisA = 0;
    axisB = 0;
    th = 0;
    ec = 0;
    cp = 0;

    volume = 0;
    growth = 0;
    massMomentOfInertia = 0;
    radiusOfGyration = 0;

    mass = 0;
    m10 = 0;
    m01 = 0;
    m11 = 0;
    m20 = 0;
    m02 = 0;
    m30 = 0;
    m03 = 0;
    m21 = 0;
    m12 = 0;
    m22 = 0;
    m31 = 0;
    m13 = 0;
    m40 = 0;
    m04 = 0;

    g00 = 0;
    g10 = 0;
    g01 = 0; 

    u11 = 0;
    u20 = 0;
    u02 = 0;
    u30 = 0;
    u03 = 0;
    u21 = 0; 
    u12 = 0;
    u22 = 0;
    u31 = 0;
    u13 = 0;
    u40 = 0;
    u04 = 0;

    n11 = 0;
    n20 = 0;
    n02 = 0;
    n30 = 0;
    n03 = 0;
    n21 = 0; 
    n12 = 0;
    n22 = 0;
    n31 = 0;
    n13 = 0;
    n40 = 0;
    n04 = 0;

    phi1 = 0;
    phi2 = 0;
    phi3 = 0;
    phi4 = 0;
    phi5 = 0;
    phi6 = 0;
    phi7 = 0;
    phi8 = 0;
    phi9 = 0;
    phi10 = 0;
    phi11 = 0;
    phi12 = 0;
    phi13 = 0;

    oldth = 0;
    dth = 0;
    th180 = 0;

    shiftX = 0;
    shiftY = 0;
    oldmX = 0;
    oldmY = 0;
    dM = 0;
    aM = 0;
    oldaM = 0;
    daM = 0;
    aMG = 0;
    dMG = 0;
    oldaMG = 0;
    daMG = 0;
}

void Grid::resetCounters() {
    zeroStats();
    oldth = 0;
    dth = 0;
    th180 = 0;
    oldaM = 0;
    shiftX = 0;
    shiftY = 0;
    mX = width/2;
    mY = width/2;
    oldmX = mX;
    oldmY = mY;
}

void Grid::clear() {
    resetCounters();
    statsI = 0;
    generation = 0;
    grid = std::vector<std::vector<double>>(width, std::vector<double>(height));
}

std::vector<double> TwoPoints(double x1, double y1, double x2, double y2, double olda) {
    double dX = x1 - x2;
    double dY = y1 - y2;
    double d = std::sqrt(dX*dX + dY*dY);
    double a = d < EPSILON && d > -EPSILON ? 0 : std::atan2(dY, dX) / M_PI * 180;
    double da = olda == 0 ? 0 : a - olda;
    da = da + 540 % 360 - 180;
    return {d, a, da};
}

std::vector<double> Grid::getStatValues() {
    return {mass/NS/NS, g00/NS/NS, volume/NS/NS, growth/NS/NS, mass/volume, g00/growth, dM/NS*updateFrequency, dMG/NS, dM/NS*updateFrequency<0.01 ? 0 : daM*updateFrequency, dMG/NS<0.01 ? 0 : daMG*updateFrequency, dth*updateFrequency, phi1, phi4, phi5, phi6, phi7, phi9, phi10, phi11, phi12};
}

void Sparkline(std::string id, const double* value, int count, double min_v, double max_v, int offset, const ImVec4& color, const ImVec2& size) {
    ImPlot::PushStyleVar(ImPlotStyleVar_PlotPadding, ImVec2(0,0));
    if (ImPlot::BeginPlot(id.c_str(), size, ImPlotFlags_CanvasOnly)) {
        ImPlot::SetupAxes(nullptr, nullptr, ImPlotAxisFlags_NoDecorations, ImPlotAxisFlags_NoDecorations);
        ImPlot::SetupAxesLimits(0, count - 1, min_v, max_v, ImGuiCond_Always);
        ImPlot::SetNextLineStyle(color);
        ImPlot::SetNextFillStyle(color, 0.25);
        ImPlot::PlotLine(id.c_str(), value, count, 1, 0, ImPlotLineFlags_Shaded, offset);
        ImPlot::EndPlot();
    }
    ImPlot::PopStyleVar();
}

ImVec4 StatRGB(int part) {
    int s = 2 - std::floor(part / 11);
    int c = (part < 11) ? part : part - 11;
    return (ImVec4){colors[s][c].x * 0x11 / 255, colors[s][c].y * 0x11 / 255, colors[s][c].z * 0x11 / 255, 1};
}

void Grid::plotStats() {
    ImGuiTableFlags flags = ImGuiTableFlags_BordersOuter | ImGuiTableFlags_BordersV | ImGuiTableFlags_RowBg | ImGuiTableFlags_Resizable | ImGuiTableFlags_Reorderable;
    if (ImGui::BeginTable("stats", 2, flags, ImVec2(-1, 0))) {
        ImGui::TableSetupColumn("Stat");
        ImGui::TableSetupColumn("Value");
        ImGui::TableHeadersRow();
        ImGui::TableSetColumnIndex(0);
        ImGui::Text("growth function");
        double j = 0;
        std::vector<double> x(1000);
        for (int i = 0; i < 1000; i++) {
            x[i] = polynomial(j, mu, sigma, alpha);
            j += 0.001;
        }
        ImGui::TableSetColumnIndex(1);
        if (ImPlot::BeginPlot("Growth function")) {
            ImPlot::PlotLine("Growth function", x.data(), 1000, 1.0/1000.0);
            ImPlot::EndPlot();
        }
        ImGui::TableSetColumnIndex(1);
        for (int i = 0; i < statNames.size(); i++) {
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0);
            ImGui::Text("%s", statNames[i].c_str());
            ImGui::TableSetColumnIndex(1);
            Sparkline(statNames[i], stats[i].data(), stats[i].size(), 0, 1, statsI, StatRGB(i), ImVec2(300, 60));
        }
        ImGui::EndTable();
    }
}   

void Grid::updateStats() {
    std::vector<double> statValues = getStatValues();
    for (int i = 0; i < statValues.size(); i++) {
        stats[i][wrap(statsI, 100)] = statValues[i];
        statsAverage[i] = statsAverage [i]+ statValues[i];
    }
    statsI++;
}

void Grid::calcMoments() {
    volume = growth = mass = m10 = m01 = m11 = m20 = m02 = m30 = m03 = m21 = m12 = m22 = m31 = m13 = m40 = m04 = 0;
    g00 = g10 = g01 = 0;

    oldmX = mX - shiftX;
    oldmY = mY - shiftY;
    oldaM = aM;
    oldaMG = aMG;
    for (int i = 0; i < width; i++) {
        for (int j = 0; j < height; j++) {
            int y = i<oldmY-width/2 ? i+width : i>oldmY+width/2 ? i-width : i;
			int x = j<oldmX-width/2 ? j+width : j>oldmX+width/2 ? j-width : j;
            double v = grid[i][j];
            double g = std::max(0.0, delta[i][j]);
            if (v > EPSILON) {
                volume++;
            }
            if (g > EPSILON) {
                growth++;
            }

            double vx = v*x;
            double vy = v*y;
            double vxx = vx*x;
            double vyy = vy*y;
            double vxy = vx*y;
            double vxxx = vxx*x;
            double vyyy = vyy*y;
            double vxxy = vxx*y;
            double vxyy = vyy*x;

            double gx = g*x;
            double gy = g*y;

            mass += v;
            m10 += vx;
            m01 += vy;
            m11 += vxy;
            m20 += vxx;
            m02 += vyy;
			m30 += vxxx;
            m03 += vyyy;
            m21 += vxxy;
            m12 += vxyy;
			m22 += vxxy*y;
            m31 += vxxx*y;
            m13 += vyyy*x;
            m40 += vxxx*x;
            m04 += vyyy*y;
			g00 += g;
			g10 += gx;
            g01 += gy;
        }
    }
    mX = mass == 0 ? width/2 : m10/mass;
    mY = mass == 0 ? height/2 : m01/mass;
    gX = g00 == 0 ? width/2 : g10/g00;
    gY = g00 == 0 ? height/2 : g01/g00;

    double X2 = mX * mX;
    double X3 = X2 * mX;
    double Y2 = mY * mY;
    double Y3 = Y2 * mY;
    double XY = mX * mY;

    u11 = m11 - mX*m01;
	u20 = m20 - mX*m10;
	u02 = m02 - mY*m01;
	u30 = m30 - 3*mX*m20 + 2*X2*m10;
	u03 = m03 - 3*mY*m02 + 2*Y2*m01;
	u21 = m21 - 2*mX*m11 - mY*m20 + 2*X2*m01;
	u12 = m12 - 2*mY*m11 - mX*m02 + 2*Y2*m10;
	u22 = m22 - 2*mY*m21 + Y2*m20 - 2*mX*m12 + 4*XY*m11 - 2*mX*Y2*m10 + X2*m02 - 2*X2*mY*m01 + X2*Y2*mass;
	u31 = m31 - mY*m30 + 3*mX*(mY*m20-m21) + 3*X2*(m11-mY*m10) + X3*(mY*mass-m01);
	u13 = m13 - mX*m03 + 3*mY*(mX*m02-m12) + 3*Y2*(m11-mX*m01) + Y3*(mX*mass-m10);
	u40 = m40 - 4*mX*m30 + 6*X2*m20 - 4*X3*m10 + X2*X2*mass;
	u04 = m04 - 4*mY*m03 + 6*Y2*m02 - 4*Y3*m01 + Y2*Y2*mass;
    if (mass < EPSILON && mass > -EPSILON) {
        n11 = 0;
        n20 = 0;
        n02 = 0;
        n30 = 0;
        n21 = 0; 
        n12 = 0;
        n22 = 0;
        n31 = 0;
        n13 = 0;
        n40 = 0;
        n04 = 0;
    } else {
        double m2 = mass*mass;
        double mA = mass * mass * std::sqrt(mass);
        double m3 = mass * mass * mass;
        n11 = u11 / m2;
        n20 = u20 / m2;
        n02 = u02 / m2;
        n30 = u30 / mA;
        n03 = u03 / mA;
        n21 = u21 / mA;
        n12 = u12 / mA;
        n22 = u22 / m3;
        n31 = u31 / m3;
        n13 = u13 / m3;
        n40 = u40 / m3;
        n04 = u04 / m3;
    }

    double Z = 2 * n11;

	double A = n20 + n02;
	double B = n20 - n02;
	double C = n30 + n12;
	double D = n30 - n12;
	double E = n03 + n21;
	double F = n03 - n21;
	double G = n30 - 3*n12;
	double H = 3*n21 - n03;
	double Y = 2*n22;
	double I = n40 + n04;
	double J = n40 - n04;
	double K = n31 + n13;
	double L = n31 - n13;
	double CC = C*C;
	double EE = E*E;
	double CC_EE = CC - EE;
	double CC_EE3 = CC - 3*EE;
	double CC3_EE = 3*CC - EE;
	double CE = C*E;
	double DF = D*F;

    phi1 = A;
	phi2 = B*B + Z*Z;  //dependent on others
	phi3 = G*G + H*H;  //dependent on others
	phi4 = CC + EE;
	phi5 = G*C*CC_EE3 + H*E*CC3_EE;
	phi7 = H*C*CC_EE3 - G*E*CC3_EE;
	phi6 = B*CC_EE + 2*Z*CE;
	phi8 = Z*CC_EE/2 - B*CE;
	phi9 = I + Y;
	phi10 =    J*CC_EE + 4*L*DF;
	phi11 = -2*K*CC_EE - 2*J*DF;
    double M = I - 3*Y;
	double t1 = CC_EE*CC_EE - 4*CE*CE;
	double t2 = 4*CE*CC_EE;
	phi12 =  4*L*t2 + M*t1;
	phi13 = -4*L*t1 + M*t2;

    double t3 = phi1 / 2 * mass;
    double t4 = std::sqrt(phi2) / 2 * mass;
    axisA = t3+t4;
	axisB = t3-t4;
    ec = std::sqrt(1 - axisB/axisA);

    cp = mass / (u20 + u02);


    th = std::atan2(Z, B) / 2 / M_PI*180;
    dth = 0;
    if (oldth > 30 && th < -30) {
        th180 = th180 + 180 % 360;
        dth = 180;
    } else if (oldth < -30 && th > 30) {
        th180 = th180 - 180 % 360;
        dth = -180;
    }

    dth += oldth == 0 ? 0 : th - oldth;
    dth = dth + 540 % 360 - 180;
    oldth = th;
    th += th180;

    std::vector<double> p2 = TwoPoints(mX, mY, oldmX, oldmY, oldaM);
    dM = p2[0];
    aM = p2[1];
    daM = p2[2];
    p2 = TwoPoints(mX, mY, gX, gY, oldaMG);
    dMG = p2[0];
    aMG = p2[1];
    daMG = p2[2];
    if (generation < 2) {
        daM = 0;
        daMG = 0;
    }
    mX = wrapD(mX, width);
    mY = wrapD(mY, width);
}

Grid::~Grid()
{
    // SDL_DestroyTexture(texture);
}

void Grid::set(int x, int y, double value)
{
    grid[wrap(x, width)][wrap(y, height)] = value;
}

double Grid::get(int x, int y)
{
    return grid[wrap(x, width)][wrap(y, height)];
}

void Grid::draw(int cellSize)
{
    for (int x = 0; x < width; x++)
    {
        for (int y = 0; y < height; y++)
        {   
            SDL_Color col = interpolateTurbo(grid[x][y]);
            uint32_t color = SDL_MapRGBA(SDL_AllocFormat(SDL_PIXELFORMAT_RGBA8888), col.r, col.g, col.b, 255);
            textureBuffer[x + y * width] = color;
        }
    }

    SDL_UpdateTexture(texture, NULL, textureBuffer, width * sizeof(uint32_t));
    SDL_Rect dest = {0, 0, width * cellSize, height * cellSize};
    SDL_RenderCopy(renderer, texture, NULL, &dest);
}


double Grid::growthFunc(double U)
{
    // return bell(U, mu, sigma) * 2 - 1;  
    return polynomial(U, mu, sigma, alpha);
}

void Grid::update()
{
    generation++;
    std::vector<std::vector<double>> newGrid = fastConv2D(grid);
    for (int x = 0; x < width; x++)
    {
        for (int y = 0; y < height; y++)
        {
            double d = growthFunc(newGrid[x][y]);
            delta[x][y] = d;
            grid[x][y] = clip(grid[x][y] + 1.0 / updateFrequency * d, 0, 1);
        }
    }
}

void Grid::setGridTemplate(GridTemplate gt) {
    template_ = gt;
    NS = gt.kernelRadius;
    updateFrequency = gt.updateFrequency;
    sigma = gt.sigma;
    mu = gt.mu;
    alpha = 4.0;
    kernel = buildKernel(width, gt.kernelRadius, gt.coreId, gt.layerId, gt.kernel_B, gt.alpha);
}

void Grid::setRule(std::string st) {
    GridTemplate gt = ParseRule(st);
    setGridTemplate(gt);
    clear();
}

std::string ToZip(int v) {
    return v == 0 ? "0" : v == 100 ? "1" : std::string(1, static_cast<char>(v + (ZIP_START - 1)));
}

int FromZip(std::string c) {
    return c == "0" ? 0 : c == "1" ? 100 : static_cast<int>(c[0]) - (ZIP_START - 1);
}

bool IsZip(std::string st) {
    return static_cast<int>(st[0]) >= ZIP_START;
}

int fromRepeatSt(std::string st) {
    return (st == "") ? 1 : IsZip(st) ? (st.size() == 1) ? FromZip(st) : FromZip(st.substr(0, 1)) * 100 + FromZip(st.substr(1, 1)) : std::stoi(st);
} 

// std::vector<std::vector<double>> Grid::ToCellVector(std::string st) {
//     bool isZip = startsWith(st, ZIP_HEADER);
//     bool isZip2 = startsWith(st, ZIP2_HEADER);
//     isZip = isZip || isZip2;
//     if (isZip) {
//         st = st.substr(isZip2 ? ZIP2_HEADER.size() : ZIP_HEADER.size());
//     }
//     std::string colSeperator = isZip ? "" : ",";

//     std::vector<std::string> cells = split(st, std::regex("\\/"));
//     int h = cells.size();
//     int w = 0;
//     for (int i = 0; i < h; i++) {
//         std::string row = std::regex_replace(cells[i], std::regex("^ +| +$"), "");
//         if (isZip) {
//             std::vector<std::string> row2 = split(row, std::regex("-"));
//             row = "";
//             for (int j = 0; j < row2.size(); j++) {
//                 std::vector<std::string> row3 = split(row2[j], std::regex("."));
//                 if (row3.size() == 1) {
//                     row += row2[j];
//                 } else {
//                     row += []{std::string s; for(int i = 0; i < 5; i++) s += ToZip(0); return s;}();
//                 }
//             }
//         }
//         // std::vector<double> rowVec = split(row, std::regex(colSeperator));
//     }
// }

void Grid::setRuleAndAddCells(std::string st) {
    std::regex re("\\;");

    std::vector<std::string> parts = split(st, re);
    std::string cellSt = "";
    if (parts.size() == 1) {
        cellSt = st;
    } else if (parts.size() >= 4) {
        re = std::regex("\\=");
        cellSt = split(parts[3], re)[1];
        std::string ruleSt = "";
        for (int i = 0; i < parts.size() - 1; i++) {
            ruleSt += parts[i] + ";";
        }
        setRule(ruleSt);
    }
    if (cellSt != "") {
    }
}