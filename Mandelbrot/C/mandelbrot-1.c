// A basic unoptimized mandelbrot set using escape time
// Binary(black and white) visulization

#include <stdio.h>
#include <stdlib.h>

// If point has escaped color it black, otherwise white
void CalculateColor(int iteration, int max_iterations, unsigned char *pixel) {
  if(iteration >= max_iterations) {
    pixel[0] = 0;
    pixel[1] = 0;
    pixel[2] = 0;
  } else {
    pixel[0] = 255;
    pixel[1] = 255;
    pixel[2] = 255;
  }
}

// Escape time algorithm
// https://en.wikipedia.org/wiki/Mandelbrot_set
void CalculatePixel(double xO, double yO, double pixel_size, unsigned char *pixel) {
  const int max_iterations = 1000;

  double x = 0.0;
  double y = 0.0;
  int iteration = 0;
  while( (x*x + y*y) < 4.0 && iteration < max_iterations) {
    const double tmp_x = x*x - y*y + xO;
    y = 2.0*x*y + yO;
    x = tmp_x;
    iteration++;
  }

  // Calculate color based upon escape iterations
  CalculateColor(iteration, max_iterations, pixel);
}

int main(int argc, char **argv) {
  // Image bounds
  const double center_x = -1.0;
  const double center_y =  0.0;
  const double length_x =  2.5;
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
  size_t pixel_bytes = sizeof(unsigned char)*3*pixels_x*pixels_y;
  unsigned char *pixels = malloc(pixel_bytes);

  // Iterate over each pixel and calculate RGB color
  for(int n_y=0; n_y<pixels_y; n_y++) {
      double y = y_min + n_y * pixel_size;
    for(int n_x=0; n_x<pixels_x; n_x++) {
      double x = x_min + n_x * pixel_size;

      unsigned char *pixel = pixels + (3 * (pixels_x * n_y + n_x));
      CalculatePixel(x, y, pixel_size, pixel);
    }
  }  

  // Dump pixels to binary file
  FILE *file = fopen("mandelbrot.raw", "wb");
  fwrite(pixels, sizeof(char), 3*pixels_x*pixels_y, file);
  printf("Saved %d by %d image\n", pixels_x, pixels_y);

  // Cleanup
  fclose(file);
  free(pixels);

  return 0;
}
