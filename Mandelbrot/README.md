# Mandelbrot
A basic mandelbrot implimentation with no optimizations

# Rendering in Paraview
* Open `mandelbrot.raw`
* toggle open advanced options(gear icon at top of properties tab)
* set `Data Scalar Type` to `unsigned char`
* set `File Dimensionality` to `2`
* set `Data Extent` upper bounds in X,Y to the number of pixels in the fractal image minus 1
* Click `Apply` 
* unselect `Map Scalars` under `Scalar Coloring` to use RGB specified in file
* Click `Edit` under the heading `Lights`, unselect `Light Kit` and select `Additional Headlight`
