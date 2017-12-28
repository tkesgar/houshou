#include <iostream>
#include "boost/format.hpp"
#include "cinder/app/App.h"
#include "cinder/app/RendererGl.h"
#include "cinder/gl/gl.h"
#include "Simulator.h"

using namespace std;
using namespace boost;
using namespace ci;
using namespace ci::app;

// Color for terrain.
const Color
    TERRAIN_COLOR = Color(0.64f, 0.64f, 0.16f),
    POLYGON_COLOR = Color(0.0f, 0.0f, 0.75f);

// Struct to store benchmark data.
struct BenchmarkRecord
{
    double t;
    double u;
    double d;
    float fps;
};

/// ================================================================================================
/// Shaders for fluid rendering
/// ================================================================================================

const char* VERTEX_SHADER = CI_GLSL(150,
    uniform mat4 ciModelViewProjection;

    in vec4 ciPosition;
    in vec4 ciColor;
    in vec4 trailPosition;

    out vec4 vColor;
    out vec4 trailPos;

    void main(void) {
        gl_Position = ciModelViewProjection * ciPosition;
        trailPos = ciModelViewProjection * trailPosition;
        vColor = ciColor;
    }
);

const char* GEOMETRY_SHADER = CI_GLSL(150,
    layout(points) in;
    layout(line_strip, max_vertices = 3) out;

    in vec4 trailPos[];
    in vec4 vColor[];

    out vec4 gColor;

    void main()
    {
        gColor = vColor[0];

        gl_Position = gl_in[0].gl_Position;
        EmitVertex();

        gl_Position = trailPos[0];
        EmitVertex();

        EndPrimitive();
    }
);

const char* FRAGMENT_SHADER = CI_GLSL(150,
    in  vec4 gColor;
    out vec4 oColor;

    void main(void) {
        oColor = gColor;
    }
);

/// ============================================================================================
/// App class definition
/// ============================================================================================

class HoushouApp : public App
{
    /// Simulation configuration
    /// ========================

    // Screen width and height.
    int windowWidth, windowHeight;

    // Grid scale.
    float gridScale;

    // Number of threads used.
    int threadCount;

    // Target FPS (-1 = disable initial frame rate, use 30 FPS if enabled while running).
    float targetFPS;

    // true = benchmark mode, false = not benchmark mode.
    bool benchmark;

    // Benchmark name.
    string benchmarkLabel;

    // Benchmark start time.
    tm benchmarkStartTime;

    // (Benchmark mode) Number of records stored before exiting the benchmark.
    int benchmarkRecordLimit;

    /// Simulation
    /// ==========

    // Simulator pointer to instance.
    Simulator* simulator = nullptr;

    /// Benchmarking
    /// ============

    // List of benchmark records.
    vector<BenchmarkRecord> benchmarkRecords;

    // Timer to count time before adding benchmark record.
    Timer benchmarkRecordTimer;

    /// Fluid rendering via shader
    /// ==========================

    geom::BufferLayout particleLayout;
    gl::VboRef particleVbo;
    gl::BatchRef particleBatch;
    gl::GlslProgRef shaderProgram;

    // Terrain texture for drawing terrain.
    gl::TextureRef terrainTexture;

    /// Timer for update and draw time
    /// ==============================

    Timer avgUpdateTimer;
    double avgUpdateTime = 0;
    double sumUpdateTime = 0;
    int updateCount = 0;

    Timer avgDrawTimer;
    double avgDrawTime = 0;
    double sumDrawTime = 0;
    int drawCount = 0;

    /// Private functions
    /// =================

    // Initializes the simulation from the beginning.
    // After init, update() and draw() should work without problems.
    void init();

    // Write benchmark result and exit.
    void writeBenchmarkResult() const;

public:

    /// Override Cinder lifecycle functions
    /// ===================================

    void setup() override;
    void cleanup() override;
    void update() override;
    void draw() override;

    /// Override mouse event handling
    /// =============================

    void mouseUp(MouseEvent event) override;
    void mouseDown(MouseEvent event) override;
    void mouseMove(MouseEvent event) override;
    void mouseDrag(MouseEvent event) override;
    void mouseWheel(MouseEvent event) override;

    /// Override keyboard event handling
    /// ================================

    void keyDown(KeyEvent event) override;
    void keyUp(KeyEvent event) override;
};

/// ================================================================================================
/// Initialization
/// ================================================================================================

void HoushouApp::init()
{
    // =============================================================================================
    // Open config file stream.
    // =============================================================================================

    // If args > 1 (args[0] is program name), use args[1] as file name; otherwise, use config.txt.
    const auto& args = getCommandLineArgs();
    const auto& configPath = getAssetPath((args.size() > 1) ? (format{ "%s.txt" } % args[1]).str() : "config.txt");

    // Open config file stream.
    fstream fconfig(configPath, ios::in);

    console() << format{ "Reading configuration file from '%s'." } % configPath << endl;

    // =============================================================================================
    // Read window width and height (for grid width and height) and grid scale.
    // =============================================================================================

    fconfig
        >> windowWidth >> windowHeight
        >> gridScale;

    // Set window size to match grid width and height.
    setWindowSize(windowWidth, windowHeight);

    // Create a new simulator.
    // If simulator is not nullptr, delete the old simulator.
    if (simulator != nullptr)
    {
        delete simulator;
    }
    simulator = new Simulator(windowWidth, windowHeight, gridScale);

    console() << format{ "Simulation size: %dx%d (grid scale=%g)" } % windowWidth % windowHeight % gridScale << endl;

    // =============================================================================================
    // Read thread count.
    // =============================================================================================

    fconfig >> threadCount;

    // If thread count <= 0, use number from thread::hardware_concurrency().
    if (threadCount <= 0)
    {
        threadCount = thread::hardware_concurrency();
        console() << "Using number of threads based on thread::hardware_concurrency()." << endl;
    }

    // Set simulator thread count.
    simulator->setThreadCount(threadCount);

    console() << format{ "Thread count: %d" } % threadCount << endl;

    // =============================================================================================
    // Read target FPS.
    // =============================================================================================

    fconfig >> targetFPS;

    // If target FPS is 0 or negative, disable frame limit.
    // Otherwise, set target FPS to the provided number.
    if (targetFPS <= 0)
    {
        disableFrameRate();
        console() << format{ "Frame rate limit is disabled" } << endl;
    }
    else // (targetFPS > 0)
    {
        setFrameRate(targetFPS);
        console() << format{ "Target frame rate: %g" } % targetFPS << endl;
    }

    // =============================================================================================
    // Read benchmark mode.
    // =============================================================================================

    string benchmarkMode;
    fconfig >> benchmarkMode;

    // If benchmark mode is "on", enable benchmark mode, then read benchmark label and benchmark record limit.
    // Otherwise, disable benchmark mode.
    if (benchmarkMode == "on")
    {
        fconfig >> benchmarkLabel >> benchmarkRecordLimit;

        benchmark = true;

        // Set benchmark start time.
        const auto t = time(nullptr);
        localtime_s(&benchmarkStartTime, &t);

        console() << format{ "Benchmark mode is on (timeout: %g s)" } % benchmarkRecordLimit << endl;
    }
    else // (benchmarkMode != "on")
    {
        benchmark = false;
        console() << "Benchmark mode is off" << endl;
    }

    // =============================================================================================
    // Read material parameters.
    // =============================================================================================

    for (int n = 0; n < MATERIALS_COUNT; n++)
    {
        Material& m = simulator->materials[n];

        // Copy default material values.
        m = Material();

        // Read material values in order.
        // Ignore if value == 0 (leave material default values).
        float* values[] = {
            &m.mass,
            &m.restDensity,
            &m.stiffness,
            &m.bulkViscosity,
            &m.surfaceTension,
            &m.kElastic,
            &m.maxDeformation,
            &m.meltRate,
            &m.viscosity,
            &m.damping,
            &m.friction,
            &m.stickiness,
            &m.smoothing,
            &m.gravity
        };
        for (int i = 0; i < _countof(values); i++) // note: _countof(values) might not be portable on Linux
        {
            float val;
            fconfig >> val;

            if (val >= 0)
            {
                *(values[i]) = val;
            }
        }

        // Read color.
        float r, g, b;
        fconfig >> r >> g >> b;
        m.color = Color(clamp(r, 0.0f, 1.0f), clamp(g, 0.0f, 1.0f), clamp(b, 0.0f, 1.0f));

        console()
            << format{ "Material #%d:" } % m.index << endl
            << format{ "    Mass:                %.3g" } % m.mass << endl
            << format{ "    Rest density:        %.3g" } % m.restDensity << endl
            << format{ "    Stiffness:           %.3g" } % m.stiffness << endl
            << format{ "    Bulk velocity:       %.3g" } % m.bulkViscosity << endl
            << format{ "    Surface tension:     %.3g" } % m.surfaceTension << endl
            << format{ "    Elastic coefficient: %.3g" } % m.kElastic << endl
            << format{ "    Maximum deformation: %.3g" } % m.maxDeformation << endl
            << format{ "    Melt rate:           %.3g" } % m.meltRate << endl
            << format{ "    Viscosity:           %.3g" } % m.viscosity << endl
            << format{ "    Damping:             %.3g" } % m.damping << endl
            << format{ "    Friction:            %.3g" } % m.friction << endl
            << format{ "    Stickiness:          %.3g" } % m.stickiness << endl
            << format{ "    Smoothing:           %.3g" } % m.smoothing << endl
            << format{ "    Gravity:             %.3g" } % m.gravity << endl
            << format{ "    Color:               (%.2f, %.2f, %.2f)" } % m.color.r % m.color.g % m.color.b << endl;
    }

    // =============================================================================================
    // Read terrain texture.
    // =============================================================================================

    // Read terrain name.
    string terrainName;
    fconfig >> terrainName;

    // Load terrain texture.
    auto terrainPath = getAssetPath((format{ "terrain/%s.bmp" } % terrainName).str());
    auto terrainImage = Surface::create(loadImage(terrainPath));

    // Create terrain surface.
    auto terrainView = Surface::create(windowWidth, windowHeight, false);

    int
        imageW = terrainImage->getWidth(),
        imageH = terrainImage->getHeight();

    // Iterate through all grid cells to set static polygon matrix.
    for (int gridX = 0; gridX < simulator->gridWidth; gridX++)
    {
        for (int gridY = 0; gridY < simulator->gridHeight; gridY++)
        {
            vec2 pos(
                (float(gridX) / simulator->gridWidth) * imageW,
                (float(gridY) / simulator->gridHeight) * imageH
            );

            // If pixel color on position is black, set static polygon matrix to true.
            bool black = terrainImage->getPixel(pos) == ColorA8u::black();
            int i = GRID_INDEX(gridX, gridY, simulator->gridHeight);
            simulator->terrainMatrix[i] = black;
        }
    }

    // Iterate through all (x, y) on the terrain view.
    for (int viewX = 0; viewX < terrainView->getWidth(); viewX++)
    {
        for (int viewY = 0; viewY < terrainView->getHeight(); viewY++)
        {
            vec2 pos(
                (float(viewX) / terrainView->getWidth()) * imageW,
                (float(viewY) / terrainView->getHeight()) * imageH
            );

            // If pixel color on position is black, set terrain view pixel to black.
            bool black = terrainImage->getPixel(pos) == ColorA8u::black();
            terrainView->setPixel(vec2(viewX, viewY), black ? ColorA(TERRAIN_COLOR) : ColorA(0, 0, 0, 0));
        }
    }

    // Create a texture for drawing from surface.
    terrainTexture = gl::Texture::create(*terrainView);

    console() << format{ "Terrain texture loaded from '%s'." } % terrainPath << endl;

    // =============================================================================================
    // Read fluids.
    // =============================================================================================

    // Read fluid count.
    int fluidCount;
    fconfig >> fluidCount;

    // Loop through each fluids.
    for (int i = 0; i < fluidCount; i++)
    {
        float
            xMin, yMin,
            xMax, yMax,
            spacingX, spacingY;
        int materialIndex;

        // Read fluid values. 
        fconfig
            >> xMin >> yMin
            >> xMax >> yMax
            >> spacingX >> spacingY
            >> materialIndex;

        // Add new particles at the described fluid position.
        const Material& m = simulator->materials[materialIndex];
        for (float x = xMin; x <= xMax; x += spacingX)
        {
            for (float y = yMin; y <= yMax; y += spacingY)
            {
                simulator->particles.push_back(Particle(*simulator, m, x, y));
            }
        }

        console()
            << format{ "Fluid with material #%d created from [%.1f, %.1f] to [%.1f, %.1f] (spacing: [%.1f, %.1f])." }
            % materialIndex
            % xMin % yMin
            % xMax % yMax
            % spacingX % spacingY
            << endl;
    }

    // =============================================================================================
    // Read polygons.
    // =============================================================================================

    // Read polygon count.
    int polygonCount;
    fconfig >> polygonCount;

    // Loop through each polygons.
    for (int i = 0; i < polygonCount; i++)
    {
        Polygon p;

        // Read vertex count.
        int vertexCount;
        fconfig >> vertexCount;

        // Loop through each vertices.
        for (int j = 0; j < vertexCount; j++)
        {
            float x, y;
            fconfig >> x >> y;
            p.points.push_back(vec2(x, y));
        }

        // Push new polygon into simulator.
        simulator->polygons.push_back(p);
    }

    // =============================================================================================
    // Additional things to do after init.
    // =============================================================================================

    // Close config stream.
    fconfig.close();

    // Setup average update and draw timers.
    avgUpdateTimer.start();
    avgDrawTimer.start();

    // Start benchmark timer if in benchmark mode.
    if (benchmark)
    {
        benchmarkRecordTimer.start();
    }
}

/// ================================================================================================
/// Setup and cleanup
/// ================================================================================================

void HoushouApp::setup()
{
    // Precompile config shader.
    try {
        shaderProgram = gl::GlslProg::create(
            gl::GlslProg::Format()
            .vertex(VERTEX_SHADER)
            .geometry(GEOMETRY_SHADER)
            .fragment(FRAGMENT_SHADER)
        );
    }
    catch (gl::GlslProgCompileExc e)
    {
        console() << e.what() << endl;
        quit();
    }

    // Prepare particle memory layout for GPU.
    particleLayout.append(geom::Attrib::POSITION, 3, sizeof(Particle), offsetof(Particle, position));
    particleLayout.append(geom::Attrib::CUSTOM_9, 3, sizeof(Particle), offsetof(Particle, trail));
    particleLayout.append(geom::Attrib::COLOR, 4, sizeof(Particle), offsetof(Particle, color));

    // Perform initialization.
    init();
}

void HoushouApp::cleanup()
{
    // Destroy simulator.
    delete simulator;
}

/// ================================================================================================
/// Application loop
/// ================================================================================================

void HoushouApp::update()
{
    // =============================================================================================
    // Perform update.
    // =============================================================================================

    Timer updateTimer(true);

    simulator->update();

    auto timeToUpdate = updateTimer.getSeconds();

    // =============================================================================================
    // Calculate average update time.
    // =============================================================================================

    // If update timer exceeds 1 second, calculate average update time, then reset accumulators and timer.
    // Otherwise, increase accumulators.
    if (avgUpdateTimer.getSeconds() > 1)
    {
        avgUpdateTime = sumUpdateTime / updateCount;

        sumUpdateTime = 0;
        updateCount = 0;

        avgUpdateTimer.start();
    }
    else // (avgUpdateTimer.getSeconds() < 1)
    {
        sumUpdateTime += timeToUpdate;
        updateCount++;
    }

    // =============================================================================================
    // Benchmark checking.
    // =============================================================================================

    if (benchmark)
    {
        // If time exceeds 1 second, push a new record.
        if (benchmarkRecordTimer.getSeconds() > 1)
        {
            BenchmarkRecord br;
            br.t = getElapsedSeconds();
            br.u = avgUpdateTime;
            br.d = avgDrawTime;
            br.fps = getAverageFps();
            benchmarkRecords.push_back(br);

            // If record size exceeds record limit, write benchmark result and exit.
            // Otherwise, restart timer.
            if (int(benchmarkRecords.size()) >= benchmarkRecordLimit)
            {
                writeBenchmarkResult();
                exit(0);
            }
            else // (benchmarkRecords.size() < benchmarkRecordLimit)
            {
                benchmarkRecordTimer.start(0);
            }
        }
    }
}

void HoushouApp::draw()
{
    Timer drawTimer(true);
    // =============================================================================================
    // BEGIN DRAW
    // =============================================================================================

    gl::clear();

    // =============================================================================================
    // Draw terrain.
    // =============================================================================================

    // Set perspective matrix and color.
    gl::setMatricesWindowPersp(getWindowSize());
    gl::color(Color::white());

    // Draw terrain texture.
    gl::draw(terrainTexture);

    // =============================================================================================
    // Draw contour (TODO enable this only on interactive mode)
    // =============================================================================================

    /*
    gl::setMatricesWindowPersp(getWindowSize());
    gl::color(Color(1.0f, 0.0f, 0.0f));

    for (int x = 0; x < sim->gridW - 1; x++)
    {
        for (int y = 0; y < sim->gridH - 1; y++)
        {
            vec2
                // offset
                offset(0.5f, 0.5f),
                // top
                t((vec2(x, y - 0.5f) + offset) * sim->scale),
                // bottom
                b((vec2(x, y + 0.5f) + offset) * sim->scale),
                // left
                l((vec2(x - 0.5f, y) + offset) * sim->scale),
                // right
                r((vec2(x + 0.5f, y) + offset) * sim->scale);

            int index = sim->normalMatrix[GRID_INDEX(x, y, sim->gridH)];
            switch (index)
            {
            case 0b0000:
            case 0b1111:
                // [   ]
                // [   ]
                break;
            case 0b0001:
            case 0b1110:
                gl::drawLine(l, b);
                // [   ]
                // [\  ]
                break;
            case 0b0010:
            case 0b1101:
                gl::drawLine(r, b);
                // [   ]
                // [  /]
                break;
            case 0b0100:
            case 0b1011:
                gl::drawLine(r, t);
                // [  \]
                // [   ]
                break;
            case 0b0111:
            case 0b1000:
                // [/  ]
                // [   ]
                gl::drawLine(l, t);
                break;
            case 0b0101:
            case 0b1010:
                gl::drawLine(r, t);
                gl::drawLine(l, b);
                // [  \]
                // [\  ]
                break;
            case 0b0011:
            case 0b1100:
                gl::drawLine(l, r);
                // [___]
                // [   ]
                break;
            case 0b0110:
            case 0b1001:
                gl::drawLine(t, b);
                // [ | ]
                // [ | ]
                break;
            default:
                throw Exception((format{ "Unknown index: %d" } % index).str());
            }
        }
    }
    */

    // =============================================================================================
    // Draw polygons.
    // =============================================================================================

    // Set perspective matrix and color.
    gl::setMatricesWindowPersp(getWindowSize());
    gl::color(POLYGON_COLOR);

    // Loop through polygons and draw.
    for (auto it = simulator->polygons.cbegin(); it != simulator->polygons.cend(); it++)
    {
        const auto polygon = *it;

        // TODO Avoid recreating PolyLine2 objects: create PolyLine2 on Polygon object.
        PolyLine2 line;
        for (auto jt = polygon.points.cbegin(); jt != polygon.points.cend(); jt++)
        {
            line.push_back(*jt);
        }

        gl::drawSolid(line);
    }

    // =============================================================================================
    // Draw fluid.
    // =============================================================================================

    // Set perspective matrix to match window size.
    gl::setMatricesWindowPersp(getWindowSize());

    // Setup particle VBO.
    particleVbo = gl::Vbo::create(GL_ARRAY_BUFFER, simulator->particles, GL_STREAM_DRAW);

    // Copy particle data to GPU.
    void *gpuMem = particleVbo->mapReplace();
    memcpy(gpuMem, simulator->particles.data(), simulator->particles.size() * sizeof(Particle));
    particleVbo->unmap();

    // Create mesh by pairing our particle layout with our particle VBO.
    auto mesh = gl::VboMesh::create((int)simulator->particles.size(), GL_POINTS, { { particleLayout, particleVbo } });
    gl::Batch::AttributeMapping mapping({ { geom::Attrib::CUSTOM_9, "trailPosition" } });

    // Draw mesh with particle batch.
    particleBatch = gl::Batch::create(mesh, shaderProgram, mapping);
    particleBatch->draw();

    // =============================================================================================
    // Draw debug string.
    // =============================================================================================

    stringstream debug;

    // Basic debug data.
    debug
        << format{ "Number of particles: %d" } % simulator->particles.size() << endl
        << format{ "Average time to update: %.3f ms" } % (avgUpdateTime * 1000) << endl
        << format{ "Average time to draw: %.3f ms" } % (avgDrawTime * 1000) << endl
        << format{ "Average FPS: %d frames/second" } % int(getAverageFps()) << endl;

    // Additional data on benchmark mode.
    if (benchmark)
    {
        debug
            << endl
            << format{ "Benchmark is currently running (timeout: %d)" } % benchmarkRecordLimit << endl
            << format{ "Record size: %d" } % benchmarkRecords.size() << endl;
    }

    // Draw string at top-left.
    gl::drawString(debug.str(), vec2(10, 10), Color::white(), Font("Arial", 16));

    // =============================================================================================
    // END DRAW
    // =============================================================================================
    auto timeToDraw = drawTimer.getSeconds();

    // =============================================================================================
    // Calculate average draw time.
    // =============================================================================================

    // If draw timer exceeds 1 second, calculate average update time, then reset accumulators and timer.
    // Otherwise, increase accumulators.
    if (avgDrawTimer.getSeconds() > 1)
    {
        avgDrawTime = sumDrawTime / drawCount;

        sumDrawTime = 0;
        drawCount = 0;

        avgDrawTimer.start();
    }
    else // (avgDrawTimer.getSeconds() < 1)
    {
        sumDrawTime += timeToDraw;
        drawCount++;
    }
}

// =============================================================================================
// Handle mouse event
// =============================================================================================

void HoushouApp::mouseUp(MouseEvent event)
{ }

void HoushouApp::mouseDown(MouseEvent event)
{ }

void HoushouApp::mouseMove(MouseEvent event)
{ }

void HoushouApp::mouseDrag(MouseEvent event)
{ }

void HoushouApp::mouseWheel(MouseEvent event)
{ }

// =============================================================================================
// Handle keyboard event
// =============================================================================================

void HoushouApp::keyDown(KeyEvent event)
{ }

void HoushouApp::keyUp(KeyEvent event)
{ }

// =============================================================================================
// Main app entry point
// =============================================================================================

CINDER_APP(HoushouApp, RendererGl, [](App::Settings *settings)
{
    // Enable console (for logging).
    settings->setConsoleWindowEnabled(true);
})

/// ================================================================================================
/// Write benchmark result
/// ================================================================================================

void HoushouApp::writeBenchmarkResult() const
{
    // =============================================================================================
    // Create benchmark output file name
    // =============================================================================================

    // Benchmark file name string.
    string filename = "benchmark";

    // Add label.
    filename += "-" + benchmarkLabel;

#ifdef DISABLE_COLLISION
    filename += "-DISABLE_COLLISION";
#endif

    // Add "DEBUG"/"RELEASE" depending on the compilation flag.
#ifdef _DEBUG
    filename += "-DEBUG";
#else
    filename += "-RELEASE";
#endif

    // Add "x86"/"x64" depending on the compilation flag.
#ifdef _M_X64
    filename += "-x64";
#else
    filename += "-x86";
#endif

    // Append time and .txt extension to file name.
    filename += (format{ "-%s.txt" } % put_time(&benchmarkStartTime, "%Y%m%d-%H.%M.%S")).str();

    // Open file stream for writing.
    fstream fbenchmark(filename, ios::out);

    // Write benchmark header.
    fbenchmark
        << benchmarkLabel
        << endl
        << endl;

    // Write benchmark result.
    for (auto it = benchmarkRecords.cbegin(); it != benchmarkRecords.cend(); it++)
    {
        auto& log = *it;

        fbenchmark
            << format{ "t = %.3f s, u = %.3f ms, d = %.3f ms, fps = %d fps" }
            % log.t % (log.u * 1000) % (log.d * 1000) % log.fps
            << endl;
    }

    fbenchmark.close();
}
