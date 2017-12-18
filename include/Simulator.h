#include <iostream>
#include <vector>
#include <thread>
#include "cinder/gl/gl.h"

#define WALL
#define MATERIALS_COUNT 4

#define GRID_SIZE(w, h)     ((w) * (h))

#define GRID_INDEX(x, y, h) (((x) * (h)) + (y))
#define GRID_X(i, h)        ((i) / (h))
#define GRID_Y(i, h)        ((i) % (h))

#define LERP(v, a, b) (((b) - (a)) * (v))

#define OUTSIDE   0
#define INSIDE -1

class Simulator;

struct Material
{
    float
        mass           = 1.0f,
        restDensity    = 2.0f,
        stiffness      = 1.0f,
        bulkViscosity  = 1.0f,
        surfaceTension = 0.0f,
        kElastic       = 0.0f,
        maxDeformation = 0.0f,
        meltRate       = 0.0f,
        viscosity      = 0.02f,
        damping        = 0.001f,
        friction       = 0.0f,
        stickiness     = 0.0f,
        smoothing      = 0.02f,
        gravity        = 0.03f;

    ci::Color color = ci::Color::gray(0.5f);

    // Index of material relative to the material list.
    // Required for fast material index reverse lookup.
    int index = -1;
};

struct Particle
{
    ci::vec3   position;
    ci::vec3   trail;
    ci::ColorA color = ci::Color::white();

    const Material& material;

    // Particle physical properties.
    float
        // position
        x = 0, y = 0,
        // velocity
        u = 0, v = 0,
        // velocity gradient
        gu = 0, gv = 0,
        // per-particle stress tension
        T00 = 0, T01 = 0, T11 = 0,
        // last 3 values for quadratic interpolation
        px[3], py[3],
        gx[3], gy[3];

    // Particle grid position and index.
    int gridX = 0, gridY = 0, gridIndex = 0;

    // Reflection status.
    int status = OUTSIDE;

    Particle(const Simulator& simulator, const Material& material, float x = 0, float y = 0, float u = 0, float v = 0);

    // Quadratic interpolation kernel weights (magic happens here).
    void quadraticInterpolationKernel();
};

struct Node
{
    // Node physical properties.
    float
        mass = 0,
        particleDensity = 0,
        // force
        gx = 0, gy = 0,
        // velocity
        u = 0, v = 0,
        // previous velocity
        u2 = 0, v2 = 0,
        // acceleration
        ax = 0, ay = 0,
        // per-material force
        cgx[MATERIALS_COUNT], cgy[MATERIALS_COUNT];
    
    // Label for active/non-active nodes.
    bool active = false;
    
    Node();
};

class Simulator
{
    // Array of grids.
    Node* grid;

    // List of active grids.
    std::vector<Node*> activeGrids;

    // Density matrix multiplication.
    float uscip(
        float p00, float x00, float y00,
        float p01, float x01, float y01,
        float p10, float x10, float y10,
        float p11, float x11, float y11,
        float u, float v
    );
public:
    // Array of materials.
    Material               materials[MATERIALS_COUNT];

    // List of particles.
    std::vector<Particle>  particles;

    // Simulation grid size.
    const int              gridW, gridH;
    // Simulation grid scale.
    const float            scale;
    
    // Static polygon matrix (will directly be copied on update).
    bool* staticPolygonMatrix;
    // Dynamic polygon matrix.
    bool* polygonMatrix;
    // Normal matrix (result of marching squares).
    int* normalMatrix;

    Simulator(int gridWidth, int gridHeight, float scale = 1.0f);
    ~Simulator();
    
    // Main update function.
    void update(double deltaTime);
};
