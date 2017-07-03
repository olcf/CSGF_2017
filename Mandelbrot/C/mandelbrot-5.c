// Unoptimized Mandelbrot set using continuous dwell and distance estimator
// Full MuEncy color Mandelbrot

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct pixel {
  unsigned char R;
  unsigned char G;
  unsigned char B;
} pixel_t;

// Convert HSV [0,1] to RGB [0,1]
void HSVtoRGB(float H, float S, float V, pixel_t *pixel) {
  float R = 0.0;
  float G = 0.0;
  float B = 0.0;
  float C = (V*S);
  float H_p = 360.0*H/60.0;
  float X = C * (1.0 - fabs((fmod(H_p,2.0) - 1.0)));

  if(H_p >= 0.0 && H_p <= 1.0) {
    R = C;
    G = X;
    B = 0.0;
  } else if(H_p > 1.0 && H_p <= 2.0) {
    R = X;
    G = C;
    B = 0.0;
  } else if(H_p > 2.0 && H_p <= 3.0) {
    R = 0.0;
    G = C;
    B = X;
  } else if(H_p > 3.0 && H_p <= 4.0) {
    R = 0.0;
    G = X;
    B = C;
  } else if(H_p > 4.0 && H_p <= 5.0) {
    R = X;
    G = 0.0;
    B = C;
  } else if(H_p > 5.0 && H_p <= 6.0) {
    R = C;
    G = 0.0;
    B = X;
  }

  float m = V - C;
  R += m;
  G += m;
  B += m;

  (*pixel).R = R * 255;
  (*pixel).G = G * 255;
  (*pixel).B = B * 255;
}

// Calculate HSV color in range [0,1]
// https://www.mrob.com/pub/muency/color.html
void CalculateColor(int dwell, double fractional_dwell, double distance, double final_y, int max_iterations, double pixel_spacing, pixel_t *pixel) {
  float H,S,V;

  // Point is within Mandelbrot set, color white
  if(dwell >= max_iterations) {
    H = 0.0;
    S = 0.0;
    V = 1.0;
    HSVtoRGB(H,S,V, pixel);
    return;
  }

  // log scale distance
  const double dist_scaled = log2(distance / pixel_spacing / 2.0);

  // Convert scaled distance to Value in 8 intervals
  if (dist_scaled > 0.0) {
    V = 1.0;
  } else if (dist_scaled > -8.0) {
    V = (8.0 + dist_scaled) / 8.0;
  } else {
    V = 0.0;
  }

  // log scale dwell
  const double dwell_scale = 1.0;
  double dwell_scaled = log(dwell)/log(dwell_scale * max_iterations);

  // Remap the scaled dwell onto Hue and Saturation
  if(dwell_scaled < 0.5) {
    dwell_scaled = 1.0 - 1.5*dwell_scaled;
    H = 1.0 - dwell_scaled;
    S = sqrt(dwell_scaled);
  } else {
    dwell_scaled = 1.5*dwell_scaled - 0.5;
    H = dwell_scaled;
    S = sqrt(dwell_scaled);
  }

  // Lighten every other stripe
  if(dwell%2) {
    H += 0.2;
    // Can also modify Saturation if desired S *= 0.8;
  }

  // Break the stripes up depending on angle of orbit at escape
  if(final_y < 0.0) {
    H += 0.02;
    S += 0.1;
  }

  // Break square into full color gradient
  H += 0.05 * fractional_dwell;

  // Go around color wheel 10x to cover range 1->max_iterations
  H *= 10.0;
  H -= floor(H);
  S -= floor(S);

  // Convert to RGB and set pixel value
  HSVtoRGB(H, S, V, pixel);
}

// Distance estimator and continuous dwell algorithm
// https://www.mrob.com/pub/muency/distanceestimator.html
// https://www.mrob.com/pub/muency/continuousdwell.html
void CalculatePixel(double xO, double yO, double pixel_size, pixel_t *pixel) {
  const int max_iterations = 10000;

  double x=0.0,y=0.0;
  double dx=0.0,dy=0.0;
  int dwell = 0;
  const double escape_radius = 100.0; // Increased from 2.0 to find better distance estimate
  double distance = 0.0;
  double continuous_dwell = 0.0;

  while( (x*x + y*y) < escape_radius*escape_radius && dwell < max_iterations) {
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

  double fractional_dwell = 0.0;

  // Calculate the distance if the orbit escaped
  if((x*x + y*y) >= escape_radius*escape_radius) {
    const double mag_z = sqrt(x*x + y*y);
    const double mag_dz = sqrt(dx*dx + dy*dy);
    distance = log(mag_z*mag_z) * mag_z / mag_dz;
    fractional_dwell = log2(log2(mag_z)) - log2(log2(escape_radius));
  }

  // Calculate color based on dwell and distance
  CalculateColor(dwell, fractional_dwell, distance, y, max_iterations, pixel_size, pixel);
}

int main(int argc, char **argv) {
  // Image bounds
  const double center_x = -0.745429;//-0.75;
  const double center_y =  0.113008;//0.0;
  const double length_x =  0.00001;//2.75;
  const double length_y =  0.00001;//2.0;

  // Convenience variables based on image bounds
  const double x_min = center_x - length_x/2.0;
  const double y_min = center_y - length_y/2.0;
  const double pixel_size = length_x/2000.0;
  const int pixels_x = length_x / pixel_size;
  const int pixels_y = length_y / pixel_size; 

  // Linearized 2D image data packed in RGB format in range [0-255]
  size_t pixel_bytes = sizeof(pixel_t)*pixels_x*pixels_y;
  pixel_t *pixels = malloc(pixel_bytes);

  // Iterate over each pixel and calculate RGB color
  for(int n_y=0; n_y<pixels_y; n_y++) {
      double y = y_min + n_y * pixel_size;
    for(int n_x=0; n_x<pixels_x; n_x++) {
      double x = x_min + n_x * pixel_size;

      pixel_t *pixel = pixels + (pixels_x * n_y + n_x);
      CalculatePixel(x, y, pixel_size, pixel);
    }
  }  

  // Write pixels to PPM P6 formatted file
  FILE *file = fopen("mandelbrot.ppm", "wb");
  fprintf(file, "P6\n%d %d\n%d\n", pixels_x, pixels_y, 255);
  fwrite(pixels, sizeof(pixel_t), pixels_x*pixels_y, file);

  // Cleanup
  fclose(file);
  free(pixels);

  return 0;
}
