# PathTracer
A PathTracing Engine for rendering 3D scenes with realistic light-ray behaviour. This was a project for the Computer Graphics course CSCD18 at the University of Toronto. The project was completed by myself (Kate Nickoriuk) and Isaiah Reimer. 

## Running the Program
- Compile the code using the makefile by running `make`. Explicit light sampling, importance sampling, and russian roulette ray termination can be activated while running the `make` command. Use `make IS=1` to enable importance sampling, `make ES=1` to enable explicit sampling, `make RR=1` to enable russian roulette, or any combination of the three. Run `make clean` before changing which of these flags are enabled.
- Run the program with the following parameters: `./PathTracer image_size depth num_samples output_name aperture_size`
- The `aperture_size` parameter is optional; not adding anything sets aperture size to 0 and disables DOF effects. 
- The scene to render, as well as camera setup is initialized in `buildScene.c`. Some samples for this file are in the `Scenes` folder.
- At runtime, the program will print out which features are being used in the render.

## Sample Images
```
./PathTracer 2048 12 2000 bubbletea.ppm .15
```
<img src="/Renders/BubbleTea.png" width="400"/>

```
./PathTracer 2048 12 400 scificube.ppm .05
```
<img src="/Renders/SciFiCube.png" width="400"/>

```
./PathTracer 1024 6 1000 cornell.ppm
```
<img src="/Renders/CornellBox.png" width="400"/>
