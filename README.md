# 2IMV20 - Visualization

## Asignment 1: Volume rendering

This repository contains the source code of the VolVis application for Volume Rendering.

The most important implementations has been done as next:

- Ray casting methods are implemented in class `RaycastRenderer.java` of package `volvis`
- Interpolation methods are coded in `InterpolationMethods.java` of the package `implementation`
- Gradient computation is done in `GradientVolume.java` of package `volume`
- Levoy's opactity and Phong's shading is computed in `LevoysIllumination.java` of package `implementation`
- Triange widget extension is done in `TransferFunction2DEditor.java` and used in method `getAlpha()` of the `LevoysIllumination.java` class.
