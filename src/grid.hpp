#pragma once

#include <vector>
#include <SDL2/SDL.h>
#include <string>

const std::string ZIP_HEADER = "(zip)";
const std::string ZIP2_HEADER = "(zip2)";
const int ZIP_START = 192;

const int B_DIV = 12;

struct GridTemplate {
    int kernelRadius;
    int updateFrequency;
    double mu;
    double sigma;
    std::vector<double> kernel_B;
    std::vector<double> kernel_E; // unused
    int coreId;
    int layerId;
    int deltaId; // TODO: implement
    double alpha;
    bool limitValue; // unused
};

class Grid {
public:
    Grid(SDL_Renderer* renderer, int width, int height, GridTemplate template_);
    ~Grid();

    double get(int x, int y);
    void calcMoments();
    void draw(int cellSize);
    void printStats();
    void plotStats();
    void set(int x, int y, double value);
    void update();
    void updateStats();
    void clear();

private:
    std::vector<std::vector<double>> convolve(std::vector<std::vector<double>> kernel);
    std::vector<double> getStatValues();
    void resetCounters();
    void zeroStats();
    void setRule(std::string st);

    double growthFunc(double U);

    void setGridTemplate(GridTemplate gt);
    void setRuleAndAddCells(std::string st);
    std::vector<std::vector<double>> ToCellVector(std::string st);

    GridTemplate template_;

    double mu;
    double sigma;
    double alpha;

    int height;
    int statsI;
    int updateFrequency;
    int width;
    int generation;
    uint32_t* textureBuffer;
    SDL_Renderer* renderer;
    SDL_Texture* texture;
    std::vector<std::vector<double>> delta;
    std::vector<std::vector<double>> grid;
    std::vector<std::vector<double>> kernel;
    std::vector<double> statsAverage;
    std::vector<std::vector<double>> stats;
    int NS = 0;

    double mX;
    double mY;
    
    double gX;
    double gY;

    double axisA;
    double axisB;
    double th;
    double ec;
    double cp;

    double volume;
    double growth;
    double massMomentOfInertia;
    double radiusOfGyration;

    double mass;
    // i dont know what i am doing i just copied this code
    double m10;
    double m01;
    double m11;
    double m20;
    double m02;
    double m30;
    double m03;
    double m21;
    double m12;
    double m22;
    double m31;
    double m13;
    double m40;
    double m04;

    double g00;
    double g10;
    double g01; 

    double u11;
    double u20;
    double u02;
    double u30;
    double u03;
    double u21; 
    double u12;
    double u22;
    double u31;
    double u13;
    double u40;
    double u04;

    double n11;
    double n20;
    double n02;
    double n30;
    double n03;
    double n21; 
    double n12;
    double n22;
    double n31;
    double n13;
    double n40;
    double n04;

    double phi1;
    double phi2;
    double phi3;
    double phi4;
    double phi5;
    double phi6;
    double phi7;
    double phi8;
    double phi9;
    double phi10;
    double phi11;
    double phi12;
    double phi13;

    double oldth;
    double dth;
    double th180;
    
    double shiftX;
    double shiftY;
    double oldmX;
    double oldmY;
    double dM;
    double aM;
    double oldaM;
    double daM;
    double dMG;
    double aMG;
    double oldaMG;
    double daMG;
};