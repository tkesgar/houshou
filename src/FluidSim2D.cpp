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

#define CONFIG_TXT "config.txt"
#define ADD_FLUID_UV_SPEED 100
#define FPS_TIME_MULTIPLIER (1.0f / getAverageFps())

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

class FluidSim2DApp : public App
{
    int frameRate;
    int windowW;
    int windowH;
    float scale;

    Simulator* sim;

    gl::GlslProgRef    shaderProgram;
    geom::BufferLayout particleLayout;
    gl::VboRef         particleVbo;
    gl::BatchRef       particleBatch;

    gl::TextureRef terrainTexture;

    Timer  lastUpdateTimer;
    double lastAvgUpdateTime = 0;
    double avgUpdateTime = 0;
    double totalUpdateTime = 0;
    int    totalUpdateCount = 0;

    bool run = true;

    bool showContour = false;

    int mouseX = 0;
    int mouseY = 0;

    float addFluidRadius;
    float addFluidSpacing;
    int   addFluidMaterialIndex;
    bool  addFluid  = false;
    float addFluidU = 0;
    float addFluidV = 0;

public:
    void setup()   override;
    void cleanup() override;
    void update()  override;
    void draw()    override;

    void mouseUp   (MouseEvent event) override;
    void mouseDown (MouseEvent event) override;
    void mouseMove (MouseEvent event) override;
    void mouseDrag (MouseEvent event) override;
    void mouseWheel(MouseEvent event) override;

    void keyDown(KeyEvent event) override;
    void keyUp  (KeyEvent event) override;
};

// =============================================================================================
// Setup and cleanup
// =============================================================================================

void FluidSim2DApp::setup()
{
    // ---------------------------------------------------------------------------------------------
    // Before setup
    // ---------------------------------------------------------------------------------------------

    // Log thread hardware concurrency.
    console() << format{ "Number of threads: %d" } % thread::hardware_concurrency() << endl;

    // ---------------------------------------------------------------------------------------------
    // Drawing preparation
    // ---------------------------------------------------------------------------------------------

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

    // ---------------------------------------------------------------------------------------------
    // Open config file stream.
    // ---------------------------------------------------------------------------------------------

    auto configPath = getAssetPath("config.txt");
    fstream fconfig(configPath, ios::in);

    console() << format{ "Reading configuration file from '%s'." } % configPath << endl;

    // ---------------------------------------------------------------------------------------------
    // Read initial config
    // ---------------------------------------------------------------------------------------------
    
    fconfig
        // window width and height
        >> windowW >> windowH
        // grid scale
        >> scale
        // frame rate
        >> frameRate
        // add fluid radius and spacing
        >> addFluidRadius >> addFluidSpacing >> addFluidMaterialIndex;

    setWindowSize(windowW, windowH);
    setFrameRate(float(frameRate));

    // ---------------------------------------------------------------------------------------------
    // Instantiate simulator
    // ---------------------------------------------------------------------------------------------
    
    sim = new Simulator((int)(windowW / scale), (int)(windowH / scale), scale);

    console() << format{ "Simulator created (w = %d, h = %d, scale = %.1f)." }
        % sim->gridWidth
        % sim->gridHeight
        % sim->gridScale
        << endl;

    // ---------------------------------------------------------------------------------------------
    // Read materials from config
    // ---------------------------------------------------------------------------------------------

    for (int n = 0; n < MATERIALS_COUNT; n++)
    {
        // List material values (config input must be written in this order).
        Material& m = sim->materials[n];
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

        // Read values. Ignore if value == 0 (leave material default values).
        // (Note that color is not read yet.)
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
    // Read static terrain image
    // ---------------------------------------------------------------------------------------------

    // Read terrain name.
    string terrainName;
    fconfig >> terrainName;

    // Read terrain and sky color.
    Color terrainColor, skyColor;
    fconfig
        >> terrainColor.r >> terrainColor.g >> terrainColor.b
        >> skyColor.r >> skyColor.g >> skyColor.b;

    // Load terrain texture.
    auto terrainPath = getAssetPath((format{ "terrain/%s.bmp" } % terrainName).str());
    auto terrainImage = Surface::create(loadImage(terrainPath));
    auto terrainView = Surface::create(windowW, windowH, false);

    int
        imageW = terrainImage->getWidth(),
        imageH = terrainImage->getHeight(),
        gridW = sim->gridWidth,
        gridH = sim->gridHeight;

    // Iterate through all grid cells.
    for (int x = 0; x < sim->gridWidth; x++)
    {
        for (int y = 0; y < sim->gridHeight; y++)
        {
            vec2 pos(
                (float(x) / sim->gridWidth) * imageW,
                (float(y) / sim->gridHeight) * imageH
            );

            // Set static polygon matrix to TRUE if the pixel color is black.
            bool black = terrainImage->getPixel(pos) == ColorA8u::black();
            sim->terrainMatrix[GRID_INDEX(x, y, sim->gridHeight)] = black;
        }
    }

    // Iterate through all (x, y) on the terrain view.
    for (int x = 0; x < terrainView->getWidth(); x++)
    {
        for (int y = 0; y < terrainView->getHeight(); y++)
        {
            vec2 pos(
                (float(x) / terrainView->getWidth()) * imageW,
                (float(y) / terrainView->getHeight()) * imageH
            );

            // TRUE if image color is black, FALSE otherwise.
            bool black = terrainImage->getPixel(pos) == ColorA8u::black();
            terrainView->setPixel(vec2(x, y), ColorA(black ? terrainColor : skyColor));
        }
    }

    // Create a texture for drawing from view.
    terrainTexture = gl::Texture::create(*terrainView);

    console() << format{ "Successfully loaded texture from '%s'." } % terrainPath << endl;

    // ---------------------------------------------------------------------------------------------
    // TODO Read dynamic polygon data
    // ---------------------------------------------------------------------------------------------

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
    // After setup
    // ---------------------------------------------------------------------------------------------

    // Setup update timer.
    lastAvgUpdateTime = getElapsedSeconds();
}

void FluidSim2DApp::cleanup()
{
    // Destroy simulator.
    delete sim;
}

// =============================================================================================
// Application loop
// =============================================================================================

void FluidSim2DApp::update()
{
    // ---------------------------------------------------------------------------------------------
    // Calculate delta time from last update timer
    // ---------------------------------------------------------------------------------------------

    double deltaTime = lastUpdateTimer.getSeconds();
    lastUpdateTimer.start();

    // ---------------------------------------------------------------------------------------------
    // Add fluid (if not paused)
    // ---------------------------------------------------------------------------------------------

    if (run)
    {
        if (addFluid)
        {
            for (float x = mouseX - addFluidRadius; x <= mouseX + addFluidRadius; x += addFluidSpacing)
            {
                for (float y = mouseY - addFluidRadius; y <= mouseY + addFluidRadius; y += addFluidSpacing)
                {
                    float
                        dx = x - mouseX,
                        dy = y - mouseY;

                    if (dx * dx + dy * dy < addFluidRadius * addFluidRadius)
                    {
                        const Material& m = sim->materials[addFluidMaterialIndex];
                        const float
                            px = x / sim->gridScale + 0.01f * (rand() / RAND_MAX),
                            py = y / sim->gridScale + 0.01f * (rand() / RAND_MAX),
                            u = addFluidU / sim->gridScale * FPS_TIME_MULTIPLIER + 0.01f * (rand() / RAND_MAX),
                            v = addFluidV / sim->gridScale * FPS_TIME_MULTIPLIER + 0.01f * (rand() / RAND_MAX);

                        sim->particles.push_back(Particle(*sim, m, px, py, u, v));
                    }
                }
            }
        }
    }

    // ---------------------------------------------------------------------------------------------
    // Perform update (surrounded by timer)
    // ---------------------------------------------------------------------------------------------

    Timer updateTimer(true);

    // Update fluid simulator.
    if (run)
    {
        sim->update();
    }

    auto timeToUpdate = updateTimer.getSeconds();

    // ---------------------------------------------------------------------------------------------
    // Calculate average update time
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
}

void FluidSim2DApp::draw()
{
    gl::clear();

    // ---------------------------------------------------------------------------------------------
    // Draw terrain
    // ---------------------------------------------------------------------------------------------

    gl::setMatricesWindowPersp(getWindowSize());
    gl::draw(terrainTexture);

    // ---------------------------------------------------------------------------------------------
    // Draw contour
    // ---------------------------------------------------------------------------------------------

    if (showContour)
    {
        gl::setMatricesWindowPersp(getWindowSize());
        gl::color(Color(1.0f, 1.0f, 0.5f));

        for (int x = 0; x < sim->gridWidth - 1; x++)
        {
            for (int y = 0; y < sim->gridHeight - 1; y++)
            {
                vec2
                    // offset
                    offset(0.5f, 0.5f),
                    // top
                    t((vec2(x, y - 0.5f) + offset) * sim->gridScale),
                    // bottom
                    b((vec2(x, y + 0.5f) + offset) * sim->gridScale),
                    // left
                    l((vec2(x - 0.5f, y) + offset) * sim->gridScale),
                    // right
                    r((vec2(x + 0.5f, y) + offset) * sim->gridScale);

                int index = sim->normalMatrix[GRID_INDEX(x, y, sim->gridHeight)];
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
    }

    // ---------------------------------------------------------------------------------------------
    // TODO draw polygon
    // ---------------------------------------------------------------------------------------------

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
    // Draw add fluid speed pointer
    // ---------------------------------------------------------------------------------------------

    gl::drawLine(
        vec2(mouseX, mouseY),
        vec2(mouseX + addFluidU, mouseY + addFluidV)
    );

    // ---------------------------------------------------------------------------------------------
    // Draw debug string
    // ---------------------------------------------------------------------------------------------

    stringstream debug;
    debug
        << format{ "Number of particles: %d" } % sim->particles.size() << endl
        << format{ "Average time to update: %.3f ms" } % (avgUpdateTime * 1000) << endl
        << format{ "Average FPS: %d frames/second" } % int(getAverageFps()) << endl;
    gl::drawString(debug.str(), vec2(10, 10), Color::white(), Font("Arial", 16));
}

// =============================================================================================
// Handle mouse event
// =============================================================================================

void FluidSim2DApp::mouseUp(MouseEvent event)
{
    // Disable add fluid if left mouse button is up.
    if (event.isLeft())
    {
        addFluid = false;
    }
}

void FluidSim2DApp::mouseDown(MouseEvent event)
{
    // Start add fluid if left mouse button is down.
    if (event.isLeft())
    {
        addFluid = true;
        mouseX = event.getX();
        mouseY = event.getY();
    }
}

void FluidSim2DApp::mouseMove(MouseEvent event)
{
    // Set mouse position.
    mouseX = event.getX();
    mouseY = event.getY();
}

void FluidSim2DApp::mouseDrag(MouseEvent event)
{
    // Set mouse position.
    mouseX = event.getX();
    mouseY = event.getY();
}

void FluidSim2DApp::mouseWheel(MouseEvent event)
{ }

// =============================================================================================
// Handle keyboard event
// =============================================================================================

void FluidSim2DApp::keyDown(KeyEvent event)
{
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
            setFrameRate(float(frameRate));
            console() << format{ "Target frame rate is now %d FPS." } % frameRate << endl;
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
}

void FluidSim2DApp::keyUp(KeyEvent event)
{ }

// =============================================================================================
// Main app entry point
// =============================================================================================

CINDER_APP(FluidSim2DApp, RendererGl, [](App::Settings *settings)
{
    // Enable console (for logging).
    settings->setConsoleWindowEnabled(true);
})
