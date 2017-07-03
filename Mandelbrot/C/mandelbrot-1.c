// Unoptimized mandelbrot set using escape iterations(dwell)
// Binary(black and white) visulization

#include <stdio.h>
#include <stdlib.h>

// If point has escaped color it black, otherwise white
void CalculateColor(int dwell, int max_iterations, unsigned char *pixel) {
  if(dwell >= max_iterations) {
    *pixel = 0;
  } else {
    *pixel[0] = 255;
  }
}

// Escape time algorithm
// https://en.wikipedia.org/wiki/Mandelbrot_set
void CalculatePixel(double xO, double yO, double pixel_size, unsigned char *pixel) {
  const int max_iterations = 1000;
  const double escape_radius = 2.0;

  double x = 0.0;
  double y = 0.0;
  int dwell = 0;
  while( (x*x + y*y) < escape_radius*escape_radius && dwell < max_iterations) {
    const double tmp_x = x*x - y*y + xO;
    y = 2.0*x*y + yO;
    x = tmp_x;
    dwell++;
  }

  // Calculate color based upon escape iterations
  CalculateColor(dwell, max_iterations, pixel);
}

int main(int argc, char **argv) {
  // Image bounds
  const double center_x = -0.75;
  const double center_y =  0.0;
  const double length_x =  2.75;
  const double length_y =  2.0;

  // Convenience variables based on image bounds
  const double x_min = center_x - length_x/2.0;
  const double y_min = center_y - length_y/2.0;
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
