#include <iostream>
#include <vector>
#include <thread>
#include "cinder/gl/gl.h"

//#define DISABLE_COLLISION

/// Special values for ambiguous marching squares values (TODO: rework on this)
#define NOT_AMBIGUOUS 0
#define AMBIGUOUS_LB 1
#define AMBIGUOUS_TR 2

#define GRID_SIZE(w, h)     ((w) * (h))

#define GRID_INDEX(x, y, h) (((x) * (h)) + (y))
#define GRID_X(i, h)        ((i) / (h))
#define GRID_Y(i, h)        ((i) % (h))

#define LERP(v, a, b) (((b) - (a)) * (v))

// Number of materials supported in simulation.
const int MATERIALS_COUNT = 4;

// Particle status relative to polygon.
const int
    // using negative values to avoid clashing with positive bitwise reflection index
    PARTICLE_OUTSIDE_POLYGON = -1,
    PARTICLE_INSIDE_POLYGON = -2;

/// Forward declaration of Simulator class.
class Simulator;

struct Material
{
    // Material properties.
    float
        mass = 1.0f,
        restDensity = 2.0f,
        stiffness = 1.0f,
        bulkViscosity = 1.0f,
        surfaceTension = 0.0f,
        kElastic = 0.0f,
        maxDeformation = 0.0f,
        meltRate = 0.0f,
        viscosity = 0.02f,
        damping = 0.001f,
        friction = 0.0f,
        stickiness = 0.0f,
        smoothing = 0.02f,
        gravity = 0.03f;

    // Material color.
    ci::Color color = ci::Color::gray(0.5f);

    // Index of material relative to the material list.
    // Required for fast material index reverse lookup.
    int index = -1;
};

struct Particle
{
    // Particle position vector (required for drawing).
    ci::vec3 position;
    // Particle trail vector (required for drawing).
    ci::vec3 trail;
    // Particle color (required for drawing).
    ci::ColorA color = ci::Color::white();

    // Particle material.
    const Material& material;

    // Particle physical properties.
    float
        x = 0, y = 0,
        u = 0, v = 0,
        gu = 0, gv = 0,
        T00 = 0, T01 = 0, T11 = 0,
        px[3], py[3],
        gx[3], gy[3];

    // Particle grid position and index.
    int gridX = 0, gridY = 0, gridIndex = 0;

    // Reflection status.
    int status = PARTICLE_OUTSIDE_POLYGON;
    // Ambiguous status.
    int ambiguous = NOT_AMBIGUOUS;

    Particle(const Simulator& simulator, const Material& material, float x = 0, float y = 0, float u = 0, float v = 0);

    // Quadratic interpolation kernel weights (magic happens here).
    inline void quadraticInterpolationKernelWeights();
};

struct Polygon
{
    // Polygon points.
    std::vector<ci::vec2> points;

    // Returns true if point is inside polygon, false otherwise.
    inline const bool isInside(const float x, const float y) const;
};

struct Node
{
    // Node physical properties.
    float
        mass = 0,
        particleDensity = 0,
        gx = 0, gy = 0,
        u = 0, v = 0,
        u2 = 0, v2 = 0,
        ax = 0, ay = 0,
        cgx[MATERIALS_COUNT], cgy[MATERIALS_COUNT];

    // Label for active/non-active nodes.
    bool active = false;

    Node();
};

class Simulator
{
    // Threads.
    std::vector<std::thread> threads;

    // Array of grids.
    Node* grid;

    // List of active grids.
    std::vector<Node*> activeGrids;

    // Dynamic polygon matrix.
    bool* solidMatrix;
    // Normal matrix (result of marching squares).
    int* normalMatrix;

    // Density matrix multiplication.
    inline const float uscip(
        const float p00, const float x00, const float y00,
        const float p01, const float x01, const float y01,
        const float p10, const float x10, const float y10,
        const float p11, const float x11, const float y11,
        const float u, const float v
    ) const;

public:
    // Array of materials.
    Material materials[MATERIALS_COUNT];

    // List of particles.
    std::vector<Particle> particles;

    // List of polygons.
    std::vector<Polygon> polygons;

    // Simulation result size.
    const int width, height;
    // Simulation grid size.
    const int gridWidth, gridHeight;
    // Simulation grid scale.
    const float gridScale;

    // Static polygon matrix (will directly be copied on update).
    bool* terrainMatrix;

    Simulator(int width, int height, float scale = 1.0f);
    ~Simulator();

    // Get number of threads.
    const int getThreadCount() const;
    // Set number of threads.
    void setThreadCount(const int threadCount);

    // Simulation update function.
    void update();
};
