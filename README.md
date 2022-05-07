# PathTracer
A PathTracing Engine for rendering 3D scenes with realistic light-ray behaviour. This was a project for the Computer Graphics course CSCD18 at the University of Toronto. The project was completed by myself (Kate Nickoriuk) and Isaiah Reimer.  

## Running the Program
- Explicit light sampling, importance sampling, and russian roulette can be activated while running the `make` command. Use `make IS=1` to enable importance sampling, `make ES=1` to enable explicit sampling, `make RR=1` to enable russian roulette, or any combination of the three. Our renders used all three: `make IS=1 ES=1 RR=1`. Run `make clean` before changing which of these flags are enabled.
- There is an additional parameter at runtime to enable depth of field, but it is optional (`./PathTracer image_size depth num_samples output_name aperture_size`). Not adding anything sets aperture size to 0 and disables DOF effects. Additionally setting up the view requires you to specify how far away from the camera the "focus plane" will be, at which the rays through a thin-lens aperture should converge.
- Camera/view setup is initialized in `buildScene.c`
- At runtime, the program will print out which features are being used in this render.

## Sample Images

<img src="/Renders/BubbleTea_PostProcessed.jpg" width="400"/>
<img src="/Renders/SciFiCube_PostProcessed.png" width="400"/>