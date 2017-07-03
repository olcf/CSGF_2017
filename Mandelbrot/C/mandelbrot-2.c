// Unoptimized mandelbrot set using escape iterations(dwell) and distance estimator
// Gray scale visulization

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Use continuous color based upon distance
void CalculateColor(int dwell, double distance, int max_iterations, double pixel_size, unsigned char *pixel) {
  if(distance <= 0.5*pixel_size) {
    *pixel = pow(distance/(0.5*pixel_size), 1.0/3.0)  * 255.0;
  } else {
    *pixel = 255;
  }
}

// Distance estimator algorithm
// https://www.mrob.com/pub/muency/distanceestimator.html
void CalculatePixel(double xO, double yO, double pixel_size, unsigned char *pixel) {
  const int max_iterations = 1000;
  const double escape_radius = 4.0;

  double x=0.0,y=0.0;
  double dx=0.0,dy=0.0;
  int dwell = 0;
  double distance = 0.0;
  while( (x*x + y*y) < (escape_radius*escape_radius) && dwell < max_iterations) {
    // Iterate orbit
    const double x_new = x*x - y*y + xO;
    const double y_new = 2.0*x*y + yO;
 
    // Calculate derivative
    const double tmp_dx = 2.0*(x*dx - y*dy) + 1.0;
    dy = 2.0*(x*dy + y*dx);
    dx = tmp_dx;

    // Update orbit
    x = x_new;
    y = y_new;
    dwell++;
  }

  // Calculate the distance if the orbit escaped
  if((x*x + y*y) >= (escape_radius*escape_radius)) {
    // Calculate distance
    const double mag_z = sqrt(x*x + y*y);
    const double mag_dz = sqrt(dx*dx + dy*dy);
    distance = log(mag_z*mag_z) * mag_z / mag_dz;
  }

  // Calculate color based on dwell and distance
  CalculateColor(dwell, distance, max_iterations, pixel_size, pixel);
}

int main(int argc, char **argv) {
  // Image bounds
  const double center_x = -0.75;
  const double center_y =  0.0;
  const double length_x =  2.75;
  const double length_y =  2.0;

  // Convenience variables based on image bounds
  const double x_min = center_x - length_x/2.0;
  const double x_max = center_x + length_x/2.0;
  const double y_min = center_y - length_y/2.0;
  const double y_max = center_y - length_y/2.0;
  const double pixel_size = 0.001;
  const int pixels_x = length_x / pixel_size;
  const int pixels_y = length_y / pixel_size; 

  // Linearized 2D image data packed in RGB format in range [0-255]
  size_t pixel_bytes = sizeof(unsigned char)*pixels_x*pixels_y;
  unsigned char *pixels = malloc(pixel_bytes);

  // Iterate over each pixel and calculate RGB color
  for(int n_y=0; n_y<pixels_y; n_y++) {
      double y = y_min + n_y * pixel_size;
    for(int n_x=0; n_x<pixels_x; n_x++) {
      double x = x_min + n_x * pixel_size;

      unsigned char *pixel = pixels + (pixels_x * n_y + n_x);
      CalculatePixel(x, y, pixel_size, pixel);
    }
  }  

  // Write pixels to PPM P6 formatted file
  FILE *file = fopen("mandelbrot.ppm", "wb");
  fprintf(file, "P5\n%d %d\n%d\n", pixels_x, pixels_y, 255);
  fwrite(pixels, sizeof(unsigned char), pixels_x*pixels_y, file);

  // Cleanup
  fclose(file);
  free(pixels);

  return 0;
}
