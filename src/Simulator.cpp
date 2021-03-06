#include "Simulator.h"
#include <algorithm>

using namespace std;
using namespace ci;

#define SMALL 0.001f
#define TURN_R -1
#define COLLINEAR 0
#define TURN_L 1

/// ================================================================================================
/// Helper functions
/// ================================================================================================

inline const int turn(const vec2 p, const vec2 q, const vec2 r);

inline const int lineTest(
    const float xtest, const float ytest,
    const float x1, const float y1,
    const float x2, const float y2
);

inline const vec2 normalIndex(const int index, const int ambiguous);

/// ================================================================================================
/// Particle
/// ================================================================================================

Particle::Particle(
    const Simulator& simulator,
    const Material& material,
    float x, float y,
    float u, float v
) :
    material(material),
    position(x, y, 0),
    trail(x, y, 0),
    color(material.color),
    x(x), y(y),
    u(u), v(v),
    px(), py(),
    gx(), gy()
{
    // Set initial grid position and index.
    gridX = int(x - 0.5f);
    gridY = int(y - 0.5f);
    gridIndex = GRID_INDEX(gridX, gridY, simulator.gridHeight);

    quadraticInterpolationKernelWeights();
}

void Particle::quadraticInterpolationKernelWeights()
{
    float dx = gridX - x, dy = gridY - y;

    px[0] = 0.5f * dx * dx + 1.5f * dx + 1.125f;
    gx[0] = dx + 1.5f;
    dx++;
    px[1] = -dx * dx + 0.75f;
    gx[1] = -2.0f * dx;
    dx++;
    px[2] = 0.5f * dx * dx - 1.5f * dx + 1.125f;
    gx[2] = dx - 1.5f;

    py[0] = 0.5f * dy * dy + 1.5f * dy + 1.125f;
    gy[0] = dy + 1.5f;
    dy++;
    py[1] = -dy * dy + .75f;
    gy[1] = -2.0f * dy;
    dy++;
    py[2] = 0.5f * dy * dy - 1.5f * dy + 1.125f;
    gy[2] = dy - 1.5f;
}

/// ================================================================================================
/// Polygon
/// ================================================================================================

const bool Polygon::isInside(const float x, const float y) const
{
    const vec2 p(x, y);
    const int pointsCount = int(points.size());
    for (int i = 0; i < pointsCount; i++)
    {
        int j = i + 1;
        if (j >= pointsCount)
        {
            j = 0;
        }

        const vec2
            a = points[i],
            b = points[j];

        if (turn(a, b, p) != TURN_L)
        {
            return false;
        }
    }

    return true;
}

/// ================================================================================================
/// Node
/// ================================================================================================

Node::Node()
{
    memset(cgx, 0, sizeof(float) * 4);
    memset(cgy, 0, sizeof(float) * 4);
}

/// ================================================================================================
/// Simulator
/// ================================================================================================

Simulator::Simulator(int width, int height, float scale) :
    gridScale(scale),
    width(width),
    height(height),
    gridWidth(int(width / scale)),
    gridHeight(int(height / scale)),
    threads(thread::hardware_concurrency())
{
    // Resize active grid size.
    activeGrids.resize(gridWidth * gridHeight);

    // Initialize grid array.
    grid = new Node[gridWidth * gridHeight]();

    // Initialize static polygon matrix.
    terrainMatrix = new bool[gridWidth * gridHeight]();

    // Initialize polygon matrix.
    solidMatrix = new bool[gridWidth * gridHeight]();

    // Initialize normal matrix.
    normalMatrix = new int[gridWidth * gridHeight]();

    // Initialize material indexes.
    for (int i = 0; i < MATERIALS_COUNT; i++)
    {
        materials[i].index = i;
    }
}

Simulator::~Simulator()
{
    delete grid;
    delete terrainMatrix;
    delete solidMatrix;
    delete normalMatrix;
}

const float Simulator::uscip(
    const float p00, const float x00, const float y00,
    const float p01, const float x01, const float y01,
    const float p10, const float x10, const float y10,
    const float p11, const float x11, const float y11,
    const float u, const float v
) const
{
    const float
        dx = x00 - x01,
        dy = y00 - y10,
        a = p01 - p00,
        b = p11 - p10 - a,
        c = p10 - p00,
        d = y11 - y01;

    return ((((d - 2 * b - dy) * u - 2 * a + y00 + y01) * v + ((3 * b + 2 * dy - d) * u + 3 * a - 2 * y00 - y01)) * v + ((((2 * c - x00 - x10) * u + (3 * b + 2 * dx + x10 - x11)) * u - b - dy - dx) * u + y00)) * v + (((x11 - 2 * (p11 - p01 + c) + x10 + x00 + x01) * u + (3 * c - 2 * x00 - x10)) * u + x00) * u + p00;
}

const int Simulator::getThreadCount() const
{
    return int(threads.size());
}

void Simulator::setThreadCount(const int threadCount)
{
    threads = vector<thread>(threadCount);
}

void Simulator::update()
{
    const int
        particleCount = int(particles.size()),
        polygonCount = int(polygons.size()),
        threadCount = int(threads.size());

#ifndef DISABLE_COLLISION
    // ---------------------------------------------------------------------------------------------
    // Polygon processings
    // ---------------------------------------------------------------------------------------------

    // Reset solid matrix.
    memset(solidMatrix, false, sizeof(bool) * GRID_SIZE(gridWidth, gridHeight));

    // Copy static polygon matrix to polygon matrix.
    memcpy_s(
        solidMatrix, sizeof(bool) * GRID_SIZE(gridWidth, gridHeight),
        terrainMatrix, sizeof(bool) * GRID_SIZE(gridWidth, gridHeight)
    );

    // Set solid matrix values.
    for (int t = 0; t < threadCount; t++)
    {
        threads[t] = thread(bind([&](const int bi, const int ei, const int t) {
            for (int pi = bi; pi < ei; pi++)
            {
                const Polygon& polygon = polygons[pi];

                // Find minimum X and Y for each polygons.
                int
                    minX = gridWidth,
                    minY = gridHeight,
                    maxX = 0,
                    maxY = 0;
                for (auto jt = polygon.points.cbegin(); jt != polygon.points.cend(); jt++)
                {
                    const vec2& p = *jt;

                    minX = min(minX, int(p.x));
                    maxX = max(maxX, int(p.x));
                    minY = min(minY, int(p.y));
                    maxY = max(maxY, int(p.y));
                }

                // Loop [minX..maxX][minY..maxY] and set solid matrix if point is inside polygon.
                for (int x = minX; x <= maxX; x++)
                {
                    for (int y = minY; y <= maxY; y++)
                    {
                        int pos = GRID_INDEX(x, y, gridHeight);

                        if (!solidMatrix[pos])
                        {
                            auto result = polygon.isInside(float(x), float(y));
                            solidMatrix[pos] = result;
                        }
                    }
                }
            }
        }, t * polygonCount / threadCount, (t + 1) == threadCount ? polygonCount : (t + 1) * polygonCount / threadCount, t));
    }
    for_each(threads.begin(), threads.end(), [](thread& x) { x.join(); });

    // Calculate normal matrix.
    for (int x = 0; x < gridWidth - 1; x++)
    {
        for (int y = 0; y < gridHeight - 1; y++)
        {
            int
                tl0 = GRID_INDEX(x, y, gridHeight),
                tr1 = GRID_INDEX(x + 1, y, gridHeight),
                br2 = GRID_INDEX(x + 1, y + 1, gridHeight),
                bl3 = GRID_INDEX(x, y + 1, gridHeight);

            int index = 0b0000;
            if (solidMatrix[tl0]) index |= 0b1000;
            if (solidMatrix[tr1]) index |= 0b0100;
            if (solidMatrix[br2]) index |= 0b0010;
            if (solidMatrix[bl3]) index |= 0b0001;

            normalMatrix[GRID_INDEX(x, y, gridHeight)] = index;
        }
    }

    // Calculate reflection for each particle.
    for (int t = 0; t < threadCount; t++)
    {
        threads[t] = thread(bind([&](const int bi, const int ei, const int t) {
            for (int pi = bi; pi < ei; pi++)
            {
                Particle& P = particles[pi];

                float x = P.x, y = P.y;
                int
                    cx = P.gridX, cy = P.gridY,
                    index = normalMatrix[GRID_INDEX(cx, cy, gridHeight)];

                // Four cell midpoints.
                vec2
                    t(cx, cy - 0.5f),
                    b(cx, cy + 0.5f),
                    l(cx - 0.5f, cy),
                    r(cx + 0.5f, cy);

                // Test whether polygon reflects or not.
                bool reflect = false;
                int ambiguous = NOT_AMBIGUOUS;
                switch (index)
                {
                case 0b0001:
                    reflect = lineTest(x, y, l.x, l.y, b.x, b.y) < 0;
                    break;
                case 0b1110:
                    reflect = lineTest(x, y, l.x, l.y, b.x, b.y) > 0;
                    break;

                case 0b0010:
                    reflect = lineTest(x, y, b.x, b.y, r.x, r.y) < 0;
                    break;
                case 0b1101:
                    reflect = lineTest(x, y, b.x, b.y, r.x, r.y) > 0;
                    break;

                case 0b0100:
                    reflect = lineTest(x, y, t.x, t.y, r.x, r.y) > 0;
                    break;
                case 0b1011:
                    reflect = lineTest(x, y, t.x, t.y, r.x, r.y) < 0;
                    break;

                case 0b0111:
                    reflect = lineTest(x, y, l.x, l.y, t.x, t.y) < 0;
                    break;
                case 0b1000:
                    reflect = lineTest(x, y, l.x, l.y, t.x, t.y) < 0;
                    break;

                case 0b0011:
                    reflect = y > cy;
                    break;
                case 0b1100:
                    reflect = y < cy;
                    break;

                case 0b0110:
                    reflect = x > cx;
                    break;
                case 0b1001:
                    reflect = x < cx;
                    break;

                case 0b0101: {
                    const bool
                        lbtest = lineTest(x, y, l.x, l.y, b.x, b.y) < 0,
                        trtest = lineTest(x, y, t.x, t.y, r.x, r.y) > 0;
                    if (lbtest || trtest)
                    {
                        reflect = true;
                        ambiguous = trtest ? AMBIGUOUS_TR : AMBIGUOUS_LB;
                    }
                    break;
                }
                case 0b1010: {
                    const bool
                        lbtest = lineTest(x, y, l.x, l.y, b.x, b.y) > 0,
                        trtest = lineTest(x, y, t.x, t.y, r.x, r.y) < 0;
                    if (lbtest || trtest)
                    {
                        reflect = true;
                        ambiguous = lbtest ? AMBIGUOUS_LB : AMBIGUOUS_TR;
                    }
                    break;
                }

                case 0b0000:
                case 0b1111:
                    break;
                default:
                    throw Exception("Unknown index: " + index);
                }

                // Status:
                //   - If reflect: use index.
                //   - If not reflect: either inside polygon (index == 0b1111) or outside.
                P.status = reflect ? index : (index == 0b1111 ? PARTICLE_INSIDE_POLYGON : PARTICLE_OUTSIDE_POLYGON);
                P.ambiguous = ambiguous;
            }
        }, t * particleCount / threadCount, (t + 1) == threadCount ? particleCount : (t + 1) * particleCount / threadCount, t));
    }
    for_each(threads.begin(), threads.end(), [](thread& x) { x.join(); });
#endif

    // ---------------------------------------------------------------------------------------------
    // BEGIN MPM
    // ---------------------------------------------------------------------------------------------

    for (int t = 0; t < threadCount; t++)
    {
        threads[t] = thread(bind([&](const int bi, const int ei, const int t)
        {
            for (int pi = bi; pi < ei; pi++)
            {
                Particle& P = particles[pi];
                const Material& M = P.material;
                Node* n;

                // Calculate velocity gradient.
                float
                    gu = 0, gv = 0,
                    dudx = 0, dudy = 0,
                    dvdx = 0, dvdy = 0,
                    *ppx = P.px, *ppy = P.py,
                    *pgx = P.gx, *pgy = P.gy;

                n = &grid[P.gridIndex];
                for (int i = 0; i < 3; i++, n += gridHeight - 3)
                {
                    for (int j = 0; j < 3; j++, n++)
                    {
                        float
                            pxi = ppx[i], pyj = ppy[j],
                            gxi = pgx[i], gyj = pgy[j],
                            phi = pxi * pyj,
                            gx = gxi * pyj, gy = pxi * gyj;

                        gu += phi * n->u2;
                        gv += phi * n->v2;

                        dudx += n->u2 * gx;
                        dudy += n->u2 * gy;
                        dvdx += n->v2 * gx;
                        dvdy += n->v2 * gy;
                    }
                }

                // Apply stress tensor to particle.
                float
                    w1 = dudy - dvdx,
                    wT0 = 0.5f * w1 * (P.T01 + P.T01),
                    wT1 = 0.5f * w1 * (P.T00 - P.T11),
                    D00 = dudx,
                    D01 = 0.5f * (dudy + dvdx),
                    D11 = dvdy,
                    trace = 0.5f * (D00 + D11);
                P.T00 += 0.5f * (-wT0 + (D00 - trace) - M.meltRate * P.T00);
                P.T01 += 0.5f * (wT1 + D01 - M.meltRate * P.T01);
                P.T11 += 0.5f * (wT0 + (D11 - trace) - M.meltRate * P.T11);

                // Apply deformation.
                float normal = P.T00 * P.T00 + 2 * P.T01 * P.T01 + P.T11 * P.T11;
                if (normal > M.maxDeformation)
                {
                    P.T00 = P.T01 = P.T11 = 0;
                }

                // Apply velocity gradient.
                P.gu = gu;
                P.gv = gv;

                // Apply velocity.
                P.u += M.smoothing * (gu - P.u);
                P.v += M.smoothing * (gv - P.v);

                // Apply position.
                P.x += gu;
                P.y += gv;

#ifndef DISABLE_COLLISION
                // ---------------------------------------------------------------------------------
                // Calculate particle position and velocity due to reflection
                // ---------------------------------------------------------------------------------

                switch (P.status)
                {
                    // Do not do anything is particle is outside.
                case PARTICLE_OUTSIDE_POLYGON:
                    break;
                    // Push outward if particle is inside.
                case PARTICLE_INSIDE_POLYGON:
                    P.x -= P.gu * 2;
                    P.y -= P.gv * 2;
                    break;
                    // Recalculate particle position and velocity.
                default:
                    /*
                    // Flip particle if reflect on diagonals.
                    switch (P.status)
                    {
                    case 0b0001:
                    case 0b1110:
                    case 0b0010:
                    case 0b1101:
                    case 0b0100:
                    case 0b1011:
                    case 0b1000:
                    case 0b0111:
                    swap(P.u, P.v);
                    swap(P.gu, P.gv);
                    break;
                    default:
                    break;
                    }

                    // Get normal vector based on index.
                    vec2 normal = -normalIndex(P.status, P.ambiguous);

                    // Multiply velocities by normal sign.
                    P.gu *= normal.x;
                    P.gv *= normal.y;
                    P.u *= normal.x;
                    P.v *= normal.y;
                    */

                    // Get normal vector based on index.
                    vec2 _normal = normalIndex(P.status, P.ambiguous);
                    _normal /= _normal.length();

                    // Calculate reflected vector.
                    vec2
                        _gv(P.gu, P.gv),
                        _v(P.u, P.v);
                    _gv = -_gv + (2 * dot(_gv, _normal) * _normal);
                    _v = -_v + (2 * dot(_v, _normal) * _normal);

                    // Copy result to particle properties.
                    P.gu = _gv.x;
                    P.gv = _gv.y;
                    P.u = _v.x;
                    P.v = _v.y;

                    P.x += P.u + 0.01f * rand() / RAND_MAX;
                    P.y += P.v + 0.01f * rand() / RAND_MAX;

                    break;
                }
#endif

                // Hard boundary correction.
                if (P.x < 1)
                {
                    P.x = 1 + 0.01f * rand() / RAND_MAX;
                }
                else if (P.x >(gridWidth - 1) - 1)
                {
                    P.x = (gridWidth - 1) - 1 - 0.01f * rand() / RAND_MAX;
                }
                if (P.y < 1)
                {
                    P.y = 1 + 0.01f * rand() / RAND_MAX;
                }
                else if (P.y >(gridHeight - 1) - 1)
                {
                    P.y = (gridHeight - 1) - 1 - 0.01f * rand() / RAND_MAX;
                }

                // Assign final particle position, trail, and color.
                P.position = vec3(P.x, P.y, 0) * gridScale;
                P.trail = vec3(P.x - P.gu, P.y - P.gv, 0) * gridScale;
                P.color = M.color;

                // Update grid cell index and kernel weights.
                int
                    cx = P.gridX = (int)(P.x - 0.5f),
                    cy = P.gridY = (int)(P.y - 0.5f);
                P.gridIndex = GRID_INDEX(cx, cy, gridHeight);

                P.quadraticInterpolationKernelWeights();

                // Add particle mass, velocity and density gradient to grid.
                int materialIndex = P.material.index;
                float
                    m = P.material.mass,
                    mu = m * P.u, mv = m * P.v,
                    *px = P.px, *gx = P.gx,
                    *py = P.py, *gy = P.gy;

                n = &grid[P.gridIndex];
                for (int i = 0; i < 3; i++, n += gridHeight - 3)
                {
                    for (int j = 0; j < 3; j++, n++)
                    {
                        float
                            pxi = px[i], pyj = py[j],
                            gxi = gx[i], gyj = gy[j],
                            phi = pxi * pyj;

                        n->active = true;
                        n->mass += phi * m;
                        n->particleDensity += phi;
                        n->u += phi * mu;
                        n->v += phi * mv;
                        n->cgx[materialIndex] += gxi * pyj;
                        n->cgy[materialIndex] += pxi * gyj;
                    }
                }
            }
        }, t * particleCount / threadCount, (t + 1) == threadCount ? particleCount : (t + 1) * particleCount / threadCount, t));
    }
    for_each(threads.begin(), threads.end(), [](thread& x) { x.join(); });

    // Add active nodes to list.
    activeGrids.clear();
    for (int i = 0; i < gridWidth * gridHeight; i++)
    {
        Node& N = grid[i];

        if (N.active && N.mass > 0)
        {
            activeGrids.push_back(&N);

            N.active = false;
            N.ax = 0;
            N.ay = 0;
            N.gx = 0;
            N.gy = 0;
            N.u /= N.mass;
            N.v /= N.mass;

            for (int j = 0; j < MATERIALS_COUNT; j++)
            {
                N.gx += N.cgx[j];
                N.gy += N.cgy[j];
            }

            for (int j = 0; j < MATERIALS_COUNT; j++)
            {
                N.cgx[j] -= N.gx - N.cgx[j];
                N.cgy[j] -= N.gy - N.cgy[j];
            }
        }
    }
    int activeGridCount = (int)activeGrids.size();

    // Calculate pressure and add forces to grid.
    for (int t = 0; t < threadCount; t++)
    {
        threads[t] = thread(bind([&](const int bi, const int ei, const int t) {
            for (int pi = bi; pi < ei; pi++)
            {
                Particle& P = particles[pi];
                const Material& M = P.material;

                Node* n;

                float
                    // External force
                    fx = 0, fy = 0,
                    // Velocity gradient
                    dudx = 0, dudy = 0,
                    dvdx = 0, dvdy = 0,
                    // Surface tension
                    sx = 0, sy = 0,
                    *ppx = P.px, *pgx = P.gx,
                    *ppy = P.py, *pgy = P.gy;

                // Calculate total velocity gradient and surface tension.
                n = &grid[P.gridIndex];
                int materialIndex = M.index;
                for (int i = 0; i < 3; i++, n += gridHeight - 3)
                {
                    for (int j = 0; j < 3; j++, n++)
                    {
                        float
                            pxi = ppx[i], pyj = ppy[j],
                            gxi = pgx[i], gyj = pgy[j],
                            gx = gxi * pyj, gy = pxi * gyj,
                            phi = pxi * pyj;

                        // Velocity gradient
                        dudx += n->u * gx;
                        dudy += n->u * gy;
                        dvdx += n->v * gx;
                        dvdy += n->v * gy;

                        // Surface tension
                        sx += phi * n->cgx[materialIndex];
                        sy += phi * n->cgy[materialIndex];
                    }
                }

                // Get node for current particle.
                int cx = (int)P.x, cy = (int)P.y;
                const Node&
                    n1 = grid[GRID_INDEX(cx, cy, gridHeight)],
                    n2 = grid[GRID_INDEX(cx, cy + 1, gridHeight)],
                    n3 = grid[GRID_INDEX(cx + 1, cy, gridHeight)],
                    n4 = grid[GRID_INDEX(cx + 1, cy + 1, gridHeight)];

                // Calculate particle density and pressure.
                float
                    density = uscip(
                        n1.particleDensity, n1.gx, n1.gy,
                        n2.particleDensity, n2.gx, n2.gy,
                        n3.particleDensity, n3.gx, n3.gy,
                        n4.particleDensity, n4.gx, n4.gy,
                        P.x - cx, P.y - cy
                    ),
                    pressure = min(M.stiffness / M.restDensity * (density - M.restDensity), 2.0f);

                // Calculate stress tensor.
                float
                    w1 = dudy - dvdx,
                    wT0 = 0.5f * w1 * (P.T01 + P.T01),
                    wT1 = 0.5f * w1 * (P.T00 - P.T11),
                    D00 = dudx,
                    D01 = 0.5f * (dudy + dvdx),
                    D11 = dvdy,
                    trace = 0.5f * (D00 + D11);
                D00 -= trace;
                D11 -= trace;
                P.T00 += 0.5f * (-wT0 + D00 - M.meltRate * P.T00);
                P.T01 += 0.5f * (wT1 + D01 - M.meltRate * P.T01);
                P.T11 += 0.5f * (wT0 + D11 - M.meltRate * P.T11);

                // Calculate stress tensor fracture.
                float normal = P.T00 * P.T00 + 2 * P.T01 * P.T01 + P.T11 * P.T11;
                if (normal > M.maxDeformation)
                {
                    P.T00 = P.T01 = P.T11 = 0;
                }
                float
                    T00 = M.mass * (M.kElastic * P.T00 + M.viscosity * D00 + pressure + trace * M.bulkViscosity),
                    T01 = M.mass * (M.kElastic * P.T01 + M.viscosity * D01),
                    T11 = M.mass * (M.kElastic * P.T11 + M.viscosity * D11 + pressure + trace * M.bulkViscosity);

                // Calculate surface tension.
                float magnitude = sx * sx + sy * sy;
                if (magnitude > 0)
                {
                    float
                        length = sqrtf(magnitude),
                        a = M.mass * M.surfaceTension / length;

                    T00 -= a * (0.5f * magnitude - sx * sx);
                    T01 -= a * (-sx * sy);
                    T11 -= a * (0.5f * magnitude - sy * sy);
                }

                // ---------------------------------------------------------------------------------
                // Add wall forces
                // ---------------------------------------------------------------------------------

                if (P.x < 4)
                {
                    fx += (4 - P.x);
                }
                else if (P.x > gridWidth - 4)
                {
                    fx += (gridWidth - 4 - P.x);
                }

                if (P.y < 4)
                {
                    fy += (4 - P.y);
                }
                else if (P.y > gridHeight - 4)
                {
                    fy += (gridHeight - 4 - P.y);
                }

#ifndef DISABLE_COLLISION
                // ---------------------------------------------------------------------------------
                // Add force to push outward if particle is not outside polygon (inside/border)
                // ---------------------------------------------------------------------------------

                if (P.status != PARTICLE_OUTSIDE_POLYGON)
                {
                    vec2 f(P.gu, P.gv);
                    f *= (1.0f / f.length()) * 1.0f;

                    fx -= f.x;
                    fy -= f.y;
                }
                // Reset particle status.
                P.status = PARTICLE_OUTSIDE_POLYGON;

#endif // ! DISABLE_COLLISION

                // Add forces to grid.
                n = &grid[P.gridIndex];
                for (int i = 0; i < 3; i++, n += gridHeight - 3) {
                    for (int j = 0; j < 3; j++, n++)
                    {
                        float
                            pxi = ppx[i], pyj = ppy[j],
                            gxi = pgx[i], gyj = pgy[j],
                            gx = gxi * pyj, gy = pxi * gyj,
                            phi = pxi * pyj;

                        n->ax += -(gx * T00 + gy * T01) + fx * phi;
                        n->ay += -(gx * T01 + gy * T11) + fy * phi;
                    }
                }
            }
        }, t * particleCount / threadCount, (t + 1) == threadCount ? particleCount : (t + 1) * particleCount / threadCount, t));
    }
    for_each(threads.begin(), threads.end(), [](thread& x) { x.join(); });

    // Update active node acceleration vectors.
    for (int t = 0; t < threadCount; t++)
    {
        threads[t] = thread(bind([&](const int bi, const int ei, const int t) {
            for (int i = bi; i < ei; i++)
            {
                Node& N = *activeGrids[i];
                N.u2 = 0;
                N.v2 = 0;
                N.ax /= N.mass;
                N.ay /= N.mass;
            }
        }, t * activeGridCount / threadCount, (t + 1) == threadCount ? activeGridCount : (t + 1) * activeGridCount / threadCount, t));
    }
    for_each(threads.begin(), threads.end(), [](thread& x) {x.join(); });

    for (int t = 0; t < threadCount; t++)
    {
        threads[t] = thread(bind([&](const int bi, const int ei, const int t) {
            for (int pi = bi; pi < ei; pi++)
            {
                Particle& P = particles[pi];
                const Material& M = P.material;

                Node* n;

                float *px = P.px, *py = P.py;

                // Update particle velocities.
                n = &grid[P.gridIndex];
                for (int i = 0; i < 3; i++, n += gridHeight - 3)
                {
                    for (int j = 0; j < 3; j++, n++)
                    {
                        float
                            pxi = px[i], pyj = py[j],
                            phi = pxi * pyj;

                        P.u += phi * n->ax;
                        P.v += phi * n->ay;
                    }
                }

                // Add damping.
                P.u *= 1 - M.damping;
                P.v *= 1 - M.damping;

                // Add gravity.
                P.v += M.gravity;

                // Calculate particle momentum.
                float
                    m = P.material.mass,
                    mu = m * P.u, mv = m * P.v;

                // Add particle momentum back to the grid.
                n = &grid[P.gridIndex];
                for (int i = 0; i < 3; i++, n += gridHeight - 3)
                {
                    for (int j = 0; j < 3; j++, n++)
                    {
                        float
                            pxi = px[i], pyj = py[j],
                            phi = pxi * pyj;

                        n->u2 += phi * mu;
                        n->v2 += phi * mv;
                    }
                }
            }
        }, t * particleCount / threadCount, (t + 1) == threadCount ? particleCount : (t + 1) * particleCount / threadCount, t));
    }
    for_each(threads.begin(), threads.end(), [](thread& x) {x.join(); });

    // Reset node values for active grids.
    for (int t = 0; t < threadCount; t++)
    {
        threads[t] = thread(bind([&](const int bi, const int ei, const int t) {
            for (int gi = bi; gi < ei; gi++)
            {
                Node& N = *activeGrids[gi];

                N.u = 0;
                N.v = 0;
                N.u2 /= N.mass;
                N.v2 /= N.mass;
                N.mass = 0;
                N.particleDensity = 0;
                memset(N.cgx, 0, sizeof(float) * MATERIALS_COUNT);
                memset(N.cgy, 0, sizeof(float) * MATERIALS_COUNT);
            }
        }, t * activeGridCount / threadCount, (t + 1) == threadCount ? activeGridCount : (t + 1) * activeGridCount / threadCount, t));
    }
    for_each(threads.begin(), threads.end(), [](thread& x) { x.join(); });

    // ---------------------------------------------------------------------------------------------
    // END MPM
    // ---------------------------------------------------------------------------------------------
}

/// ================================================================================================
/// Helper functions
/// ================================================================================================

const int turn(const vec2 p, const vec2 q, const vec2 r)
{
    const float result = ((r.x - q.x) * (p.y - q.y)) - ((r.y - q.y) * (p.x - q.x));
    return (abs(result) < SMALL) ? COLLINEAR : (signbit(result) ? TURN_L : TURN_R);
}

const int lineTest(
    const float xtest, const float ytest,
    const float x1, const float y1,
    const float x2, const float y2
)
{
    float
        m = (y2 - y1) / (x2 - x1),
        c = y1 - m * x1,
        d = ytest - (m * xtest + c);

    return (abs(d) < SMALL) ? 0 : (signbit(d) ? 1 : -1);
}

const vec2 normalIndex(const int index, const int ambiguous)
{
    switch (index)
    {
    case 0b0001: return vec2(-1, 1);
    case 0b1110: return vec2(1, -1);

    case 0b0010: return vec2(-1, -1);
    case 0b1101: return vec2(1, 1);

    case 0b0100: return vec2(1, -1);
    case 0b1011: return vec2(-1, 1);

    case 0b1000: return vec2(1, 1);
    case 0b0111: return vec2(-1, -1);

    case 0b0011: return vec2(0, -1);
    case 0b1100: return vec2(0, 1);

    case 0b0110: return vec2(-1, 0);
    case 0b1001: return vec2(1, 0);

    case 0b0101:
        switch (ambiguous)
        {
        case AMBIGUOUS_LB:
            return vec2(1, -1);
        case AMBIGUOUS_TR:
            return vec2(-1, 1);
        default:
            throw Exception("Invalid ambiguous: " + ambiguous);
        }
    case 0b1010:
        switch (ambiguous)
        {
        case AMBIGUOUS_LB:
            return vec2(-1, -1);
        case AMBIGUOUS_TR:
            return vec2(1, 1);
        default:
            throw Exception("Invalid ambiguous: " + ambiguous);
        }

    case 0b0000:
    case 0b1111:
    default:
        throw Exception("Invalid index: " + index);
    }
}
