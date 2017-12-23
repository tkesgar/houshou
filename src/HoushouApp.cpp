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

extern const char
    *VERTEX_SHADER,
    *GEOMETRY_SHADER,
    *FRAGMENT_SHADER;

const Color
    TERRAIN_COLOR = Color(0.64f, 0.64f, 0.16f),
    SKY_COLOR = Color::black();

struct PerformanceLog
{
    double elapsedTime;
    double avgUpdateTime;
    double avgDrawTime;
    int avgFps;
};

class HoushouApp : public App
{
    int windowW;
    int windowH;
    float gridScale;
    int threadCount;
    float targetFPS;
    bool benchmark;
    string benchmarkLabel;
    int benchmarkTimeout;

    Timer benchmarkTimer;
    vector<PerformanceLog> benchmarkRecord;

    Simulator* sim = nullptr;

    gl::GlslProgRef shaderProgram;
    geom::BufferLayout particleLayout;
    gl::VboRef particleVbo;
    gl::BatchRef particleBatch;

    gl::TextureRef terrainTexture;

    Timer lastUpdateTimer;
    double lastAvgUpdateTime = 0;
    double avgUpdateTime = 0;
    double totalUpdateTime = 0;
    int totalUpdateCount = 0;

    Timer lastDrawTimer;
    double lastAvgDrawTime = 0;
    double avgDrawTime = 0;
    double totalDrawTime = 0;
    int totalDrawCount = 0;

    void init();
    
public:
    void setup() override;
    void cleanup() override;
    void update() override;
    void draw() override;

    void mouseUp(MouseEvent event) override;
    void mouseDown(MouseEvent event) override;
    void mouseMove(MouseEvent event) override;
    void mouseDrag(MouseEvent event) override;
    void mouseWheel(MouseEvent event) override;

    void keyDown(KeyEvent event) override;
    void keyUp(KeyEvent event) override;
};

// =============================================================================================
// Initialization
// =============================================================================================

void HoushouApp::init()
{
    // ---------------------------------------------------------------------------------------------
    // Open config file stream.
    // ---------------------------------------------------------------------------------------------

    const auto& args = getCommandLineArgs();
    auto configFile = args.size() > 1 ? (format{ "%s.txt" } % args[1]).str() : "config.txt";
    auto configPath = getAssetPath(configFile);
    fstream fconfig(configPath, ios::in);
    
    console() << format{ "Reading configuration file from '%s'." } % configPath << endl;

    // ---------------------------------------------------------------------------------------------
    // Read window width and height (for grid width and height) and grid scale.
    // ---------------------------------------------------------------------------------------------

    // If simulator is not nullptr, delete the old simulator.
    if (sim != nullptr)
    {
        delete sim;
        sim = nullptr;
    }

    fconfig
        >> windowW >> windowH
        >> gridScale;
    
    setWindowSize(windowW, windowH);
    sim = new Simulator(windowW, windowH, gridScale);

    console() << format{ "MPM grid size: %dx%d (grid scale=%g)" } % windowW % windowH % gridScale << endl;

    // ---------------------------------------------------------------------------------------------
    // Read thread count.
    // ---------------------------------------------------------------------------------------------

    fconfig >> threadCount;

    if (threadCount <= 0)
    {
        threadCount = thread::hardware_concurrency();
        console() << "Using number of threads based on thread::hardware_concurrency()." << endl;
    }

    sim->setThreadCount(threadCount);

    console() << format{ "Thread count: %d" } % threadCount << endl;

    // ---------------------------------------------------------------------------------------------
    // Read target FPS.
    // ---------------------------------------------------------------------------------------------

    fconfig >> targetFPS;
    if (targetFPS > 0)
    {
        setFrameRate(targetFPS);
        console() << format{ "Target frame rate: %g" } % targetFPS << endl;
    }
    else
    {
        disableFrameRate();
        console() << format{ "Frame rate limit is disabled" } << endl;
    }

    // ---------------------------------------------------------------------------------------------
    // Set benchmark mode.
    // ---------------------------------------------------------------------------------------------

    string benchmarkMode;
    fconfig >> benchmarkMode;

    if (benchmarkMode == "on")
    {
        fconfig >> benchmarkLabel >> benchmarkTimeout;
        benchmark = true;
        console() << format{ "Benchmark mode is on (timeout: %g s)" } % benchmarkTimeout << endl;
    }
    else // (benchmarkMode != "on")
    {
        benchmark = false;
        console() << "Benchmark mode is off" << endl;
    }

    // ---------------------------------------------------------------------------------------------
    // Read material parameters.
    // ---------------------------------------------------------------------------------------------

    for (int n = 0; n < MATERIALS_COUNT; n++)
    {
        Material& m = sim->materials[n];

        // Copy default material values.
        m = Material();

        // Read values in order. Ignore if value == 0 (leave material default values).
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
        for (int i = 0; i < _countof(values); i++)
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

        // Print material data to console.
        console()
            << format{ "Values for material #%d:" } % m.index << endl
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

    // ---------------------------------------------------------------------------------------------
    // Read terrain texture.
    // ---------------------------------------------------------------------------------------------

    // Read terrain name.
    string terrainName;
    fconfig >> terrainName;

    // Load terrain texture.
    auto terrainPath = getAssetPath((format{ "terrain/%s.bmp" } % terrainName).str());
    auto terrainImage = Surface::create(loadImage(terrainPath));
    auto terrainView = Surface::create(windowW, windowH, false);

    int
        imageW = terrainImage->getWidth(),
        imageH = terrainImage->getHeight();

    // Iterate through all grid cells.
    for (int gridX = 0; gridX < sim->gridW; gridX++)
    {
        for (int gridY = 0; gridY < sim->gridH; gridY++)
        {
            vec2 pos(
                (float(gridX) / sim->gridW) * imageW,
                (float(gridY) / sim->gridH) * imageH
            );

            // Set static polygon matrix to TRUE if the pixel color is black.
            bool black = terrainImage->getPixel(pos) == ColorA8u::black();
            sim->terrainMatrix[GRID_INDEX(gridX, gridY, sim->gridH)] = black;
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

            // TRUE if image color is black, FALSE otherwise.
            bool black = terrainImage->getPixel(pos) == ColorA8u::black();
            terrainView->setPixel(vec2(viewX, viewY), ColorA(black ? TERRAIN_COLOR : SKY_COLOR));
        }
    }

    // Create a texture for drawing from view.
    terrainTexture = gl::Texture::create(*terrainView);

    console() << format{ "Terrain texture loaded from '%s'." } % terrainPath << endl;

    // ---------------------------------------------------------------------------------------------
    // Read fluid list
    // ---------------------------------------------------------------------------------------------

    int fluidCount;
    fconfig >> fluidCount;
    for (int i = 0; i < fluidCount; i++)
    {
        float
            xMin, yMin,
            xMax, yMax,
            spX, spY;
        int materialIndex;

        fconfig
            >> xMin >> yMin
            >> xMax >> yMax
            >> spX >> spY
            >> materialIndex;

        const Material& m = sim->materials[materialIndex];
        for (float x = xMin; x <= xMax; x += spX)
        {
            for (float y = yMin; y <= yMax; y += spY)
            {
                sim->particles.push_back(Particle(*sim, m, x, y));
            }
        }

        console() << format{ "Fluid with material #%d created from [%.1f, %.1f] to [%.1f, %.1f] (spacing: [%.1f, %.1f])." }
            % materialIndex
            % xMin % yMin
            % xMax % yMax
            % spX % spY
            << endl;
    }

    // ---------------------------------------------------------------------------------------------
    // Read dynamic polygon data
    // ---------------------------------------------------------------------------------------------

    int polygonCount;
    fconfig >> polygonCount;
    for (int i = 0; i < polygonCount; i++)
    {
        Polygon p;

        int vertexCount;
        fconfig >> vertexCount;
        for (int j = 0; j < vertexCount; j++)
        {
            float x, y;
            fconfig >> x >> y;
            p.points.push_back(vec2(x, y));
        }

        sim->polygons.push_back(p);
    }

    // ---------------------------------------------------------------------------------------------
    // After init
    // ---------------------------------------------------------------------------------------------

    // Close config stream.
    fconfig.close();

    // Setup last average times.
    lastAvgUpdateTime = lastAvgDrawTime = getElapsedSeconds();

    // Start benchmark timer.
    if (benchmark)
    {
        benchmarkTimer.start();
    }
}

// =============================================================================================
// Setup and cleanup
// =============================================================================================

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

    init();
}

void HoushouApp::cleanup()
{
    // Destroy simulator.
    delete sim;
    sim = nullptr;
}

// =============================================================================================
// Application loop
// =============================================================================================

void HoushouApp::update()
{
    // ---------------------------------------------------------------------------------------------
    // Calculate delta time from last update timer
    // ---------------------------------------------------------------------------------------------

    double deltaTime = lastUpdateTimer.getSeconds();
    lastUpdateTimer.start();

    // ---------------------------------------------------------------------------------------------
    // Perform update (surrounded by timer)
    // ---------------------------------------------------------------------------------------------

    Timer updateTimer(true);

    sim->update(deltaTime);

    auto timeToUpdate = updateTimer.getSeconds();

    // ---------------------------------------------------------------------------------------------
    // Calculate average update time.
    // ---------------------------------------------------------------------------------------------

    double elapsedSeconds = getElapsedSeconds();
    if (elapsedSeconds - lastAvgUpdateTime > 1)
    {
        avgUpdateTime = totalUpdateTime / totalUpdateCount;
        lastAvgUpdateTime = elapsedSeconds;
        totalUpdateTime = 0;
        totalUpdateCount = 0;
    }
    else
    {
        totalUpdateTime += timeToUpdate;
        totalUpdateCount++;
    }

    // ---------------------------------------------------------------------------------------------
    // Check benchmark.
    // ---------------------------------------------------------------------------------------------

    if (benchmark)
    {
        if (benchmarkTimer.getSeconds() > 1)
        {
            PerformanceLog log;
            log.elapsedTime = getElapsedSeconds();
            log.avgUpdateTime = avgUpdateTime;
            log.avgDrawTime = avgDrawTime;
            log.avgFps = int(getAverageFps());
            benchmarkRecord.push_back(log);

            if (benchmarkRecord.size() >= benchmarkTimeout)
            {
                // ---------------------------------------------------------------------------------
                // Write benchmark output
                // ---------------------------------------------------------------------------------

                auto t = time(nullptr);
                tm time;
                localtime_s(&time, &t);

                string fname = "benchmark";
#ifdef _DEBUG
                fname += "-DEBUG";
#else
                fname += "-RELEASE";
#endif
#ifdef _M_X64
                fname += "-x64";
#else
                fname += "-x86";
#endif
                fname += (format{ "-%s.txt" } % put_time(&time, "%Y%m%d-%H.%M.%S")).str();

                fstream fbenchmark(fname, ios::out);
                fbenchmark
                    << fname << endl
                    << benchmarkLabel << endl
                    << endl;

                for (auto it = benchmarkRecord.cbegin(); it != benchmarkRecord.cend(); it++)
                {
                    auto& log = *it;

                    fbenchmark
                        << format{ "t = %.3f s, u = %.3f ms, d = %.3f ms, fps = %d fps" }
                        % log.elapsedTime
                        % (log.avgUpdateTime * 1000)
                        % (log.avgDrawTime * 1000)
                        % log.avgFps
                        << endl;
                }

                fbenchmark.close();

                exit(0);
            }
            else
            {
                benchmarkTimer.start(0);
            }
        }
    }
}

void HoushouApp::draw()
{
    // Start draw timer.
    Timer drawTimer(true);

    // ---------------------------------------------------------------------------------------------
    // BEGIN DRAW
    // ---------------------------------------------------------------------------------------------

    gl::clear();

    // ---------------------------------------------------------------------------------------------
    // Draw terrain
    // ---------------------------------------------------------------------------------------------

    gl::setMatricesWindowPersp(getWindowSize());
    gl::color(Color::white());
    gl::draw(terrainTexture);

    // ---------------------------------------------------------------------------------------------
    // Draw contour
    // ---------------------------------------------------------------------------------------------
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
    // ---------------------------------------------------------------------------------------------
    // Draw polygon
    // ---------------------------------------------------------------------------------------------

    gl::setMatricesWindowPersp(getWindowSize());
    gl::color(Color(0.0f, 0.0f, 0.75f));
    for (auto it = sim->polygons.cbegin(); it != sim->polygons.cend(); it++)
    {
        const auto polygon = *it;
        
        PolyLine2 line;
        for (auto jt = polygon.points.cbegin(); jt != polygon.points.cend(); jt++)
        {
            line.push_back(*jt);
        }

        gl::drawSolid(line);
    }

    // ---------------------------------------------------------------------------------------------
    // Draw fluid
    // ---------------------------------------------------------------------------------------------

    // Set perspective matrix to match window size.
    gl::setMatricesWindowPersp(getWindowSize());

    // Setup particle VBO.
    particleVbo = gl::Vbo::create(GL_ARRAY_BUFFER, sim->particles, GL_STREAM_DRAW);

    // Copy particle data to GPU.
    void *gpuMem = particleVbo->mapReplace();
    memcpy(gpuMem, sim->particles.data(), sim->particles.size() * sizeof(Particle));
    particleVbo->unmap();

    // Create mesh by pairing our particle layout with our particle VBO.
    auto mesh = gl::VboMesh::create((int)sim->particles.size(), GL_POINTS, { { particleLayout, particleVbo } });
    gl::Batch::AttributeMapping mapping({ { geom::Attrib::CUSTOM_9, "trailPosition" } });

    // Draw mesh with particle batch.
    particleBatch = gl::Batch::create(mesh, shaderProgram, mapping);
    particleBatch->draw();

    // ---------------------------------------------------------------------------------------------
    // Draw debug string
    // ---------------------------------------------------------------------------------------------

    stringstream debug;
    debug
        << format{ "Number of particles: %d" } % sim->particles.size() << endl
        << format{ "Average time to update: %.3f ms" } % (avgUpdateTime * 1000) << endl
        << format{ "Average time to draw: %.3f ms" } % (avgDrawTime * 1000) << endl
        << format{ "Average FPS: %d frames/second" } % int(getAverageFps()) << endl;

    if (benchmark)
    {
        debug
            << endl
            << format{ "Benchmark is currently running (timeout: %d)" } % benchmarkTimeout << endl
            << format{ "Record size: %d" } % benchmarkRecord.size() << endl;
    }

    gl::drawString(debug.str(), vec2(10, 10), Color::white(), Font("Arial", 16));

    // ---------------------------------------------------------------------------------------------
    // END DRAW
    // ---------------------------------------------------------------------------------------------

    // End draw timer.
    auto timeToDraw = drawTimer.getSeconds();

    // ---------------------------------------------------------------------------------------------
    // Calculate average draw time.
    // ---------------------------------------------------------------------------------------------

    double elapsedSeconds = getElapsedSeconds();
    if (elapsedSeconds - lastAvgDrawTime > 1)
    {
        avgDrawTime = totalDrawTime / totalDrawCount;
        lastAvgDrawTime = elapsedSeconds;
        totalDrawTime = 0;
        totalDrawCount = 0;
    }
    else
    {
        totalDrawTime += timeToDraw;
        totalDrawCount++;
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
{
    /*
    switch (event.getCode())
    {
    // F1: disable/enable frame rate locking.
    case KeyEvent::KEY_F1:
        if (isFrameRateEnabled())
        {
            disableFrameRate();
            console() << "Frame rate is now disabled." << endl;
        }
        else // (!isFrameRateEnabled())
        {
            setFrameRate(float(targetFPS));
            console() << format{ "Target frame rate is now %d FPS." } % targetFPS << endl;
        }
        break;
    // F2: disable/enable draw contour.
    case KeyEvent::KEY_F2:
        showContour = !showContour;
        if (showContour)
        {
            console() << "Contour is now shown (note: drawing might be slowed down)." << endl;
        }
        else // (!showContour)
        {
            console() << "Contour is now not shown." << endl;
        }
        break;
    // F3: pause/resume simulation.
    case KeyEvent::KEY_F3:
        run = !run;
        if (run)
        {
            console() << "Simulation resumed." << endl;
        }
        else // (!runUpdate)
        {
            console() << "Simulation paused." << endl;
        }
        break;
    // Right/Left: increase/decrease addFluidU
    // Up/Down:    increase/decrease addFluidV
    // Backspace:  reset addFluidU-V
    case KeyEvent::KEY_BACKSPACE:
        addFluidU = addFluidV = 0;
        break;
    case KeyEvent::KEY_UP:
        addFluidV += ADD_FLUID_UV_SPEED * FPS_TIME_MULTIPLIER;
        break;
    case KeyEvent::KEY_DOWN:
        addFluidV -= ADD_FLUID_UV_SPEED * FPS_TIME_MULTIPLIER;
        break;
    case KeyEvent::KEY_LEFT:
        addFluidU -= ADD_FLUID_UV_SPEED * FPS_TIME_MULTIPLIER;
        break;
    case KeyEvent::KEY_RIGHT:
        addFluidU += ADD_FLUID_UV_SPEED * FPS_TIME_MULTIPLIER;
        break;
    default:
        // Ignore unknown keystrokes.
        break;
    }
    */
}

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

// =============================================================================================
// Shaders
// =============================================================================================

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
