# houshou

This is a repository for implementation of my college Final Project on fluid simulation with rigid body collisions.

## Usage

Run the program with:

```
houshou.exe [<configuration name>]
```

The program will look for the specified configuration text file by the name of `<configuration name>.txt` in the folder `assets`. If configuration name is not specified, `config` is used.

### Configuration

Below is the configuration file format specification. Since I use C++ `std::cin` to read the configuration, the whitespace does not actually matters; I just write them to make it more readable.

```
<window width> <window height>
<grid scale>
<thread count>
<target FPS>
<benchmark mode> [<benchmark mode timeout>]

<material parameter lines>

<terrain name>

<number of initial fluid>

<fluid #n top-left x>     <fluid #n top-left y>
<fluid #n bottom-right x> <fluid #n bottom-right y>
<fluid #n spacing x>      <fluid #n spacing y>
<fluid #n material index [0..3]>
```

- **<window width>**, **<window height>**: the window size of the simulation, in pixels.
- **<grid scale>**: the MPM grid scale. 1.0 means every grid represents 1 pixel.
- **<thread count>**: specifies the number of threads that will be used for concurrency. Use 0 or negative integers to use value from `std::thread::hardware_concurrency()`.
- **<target FPS>**: target FPS to obtain, in frames per second. Use negative values (-1) to disable target FPS (run as fast as possible).
- **<benchmark mode>**: `on`/`off`. If benchmark mode is enabled, user input will be disabled and program will automatically exit after a specified time. A benchmark result will be written at `benchmark-<time>.txt`. **<benchmark time>** must also be specified, in seconds.
- **<material parameters>**: consists of exactly 4 lines of material parameter line (see below).
- **<terrain name>>**: terrain name to be used. See the available terrain list below.
- **<number of fluid>**: number of fluids at the start of simulation
  - For each fluids: top-left coordinate, bottom-right coordinate, spacings, and material index ([0..3]) of each fluid. All coordinates and spacings are in the screen coordinate, top-left is `(0, 0)` and bottom-right is `(<window width>, <window height>)`.

Here is a sample configuration file:

```
800 480
2.0
4
30
on 30

1.0     -1      -1      -1      -1      -1      -1      -1      0.04    -1      -1      -1      -1      -1      1.0  0.5  0.5
1.0     10.0    1.0     3.0     -1      1.0     -1      1.0     1.0     -1      -1      -1      -1      -1      0.5  1.0  0.5
0.7     -1      -1      -1      -1      -1      -1      -1      0.03    -1      -1      -1      -1      -1      0.5  0.5  1.0
-1      -1      -1      -1      -1      -1      -1      -1      -1      -1      -1      -1      -1      -1      0.5  0.5  0.5

box

2

50   10
75   60
0.5  0.25
1

100  50
150  100
0.75 0.75
2
```

#### Material parameter line

```
mass    rest    stiff   bulk    surface k       max     melt    visco-  damping fric-   stick-  smooth- gravity color(rgb)
        density ness    visco-  tension elastic defor-  rate    sity            tion    iness   ing
                        sity                    mation
1st     2nd     3rd     4th     5th     6th     7th     8th     9th     10th    11th    12th    13th    14th    15th 16th 17th
1.0     2.0     1.0     1.0     0.0     0.0     0.0     0.0     0.02    0.001   0.0     0.0     0.02    0.03    -1   -1   -1
-1      -1      -1      -1      -1      -1      -1      -1      -1      -1      -1      -1      -1      -1      0.1  0.75 0.1
-1      -1      -1      -1      -1      -1      -1      -1      -1      -1      -1      -1      -1      -1      0.1  0.1  0.75
-1      -1      -1      -1      -1      -1      -1      -1      -1      -1      -1      -1      -1      -1      0.5  0.5  0.5
```

If a negative value is specified, a default value will be used; the default values for each material value matches the first line of the provided sample above. Colors are interpreted as red, green, and blue values respectively; any values will be clamped between 0 and 1.

### Controls

Note that controls are disabled in benchmark mode.

- **F1**: enable/disable frame rate limit to the specified value on the configuration. If negative value is used in the configuration, 30 FPS will be used.
- **F2**: enable/disable contour drawing. The current implementation is highly inefficient (using Cinder's `drawLine()` function), so the simulation perhaps will lag.
- **P**: pauses/resumes the simulation.
- **R**: reset the simulation, rereading the configuration file.
- **Left click**: creates a fluid in the specified position. On terrain brush mode, draws terrain in the specified position.
- **Right click**: On terrain brush mode, clears terrain in the specified position.
- **1**: remove mouse pointer polygon.
- **2**: use 'rod' mouse pointer polygon.
- **3**: use 'ladle' mouse pointer polygon.
- **4**: use terrain brush mouse pointer.
- **Q**: rotate mouse pointer counter-clockwise.
- **E**: rotate mouse pointer clockwise.
- **W**: increase mouse pointer fluid speed. On terrain brush mode, increase brush size.
- **S**: decrease mouse pinter fluid speed. On terrain brush mode, decrease brush size.

## Terrains

TBD

## Development

TBD

## Licence

Licenced under [MIT License](LICENSE).

The implementation is heavily modified from [dodydharma/Multithreaded-MPM-Fluid-Simulation][dodydharma-mpm], which itself is forked from [kotsoft/FluidCinder][kotsoft-fluidsim].

[cinder-github]:  https://github.com/cinder/cinder
[cinder-a064d9d]: https://github.com/cinder/cinder/commit/a064d9d
[dodydharma-mpm]: https://github.com/dodydharma/Multithreaded-MPM-Fluid-Simulation
[kotsoft-fluidsim]: https://github.com/kotsoft/FluidCinder
