// Unoptimized Mandelbrot set using Dwell and Distance Estimator methods
// Continuous color

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Convert HSV [0,1] to RGB [0,1]
// https://en.wikipedia.org/wiki/HSL_and_HSV#Converting_to_RGB
void HSVtoRGB(float H, float S, float V, unsigned char *pixel) {
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

  pixel[0] = R * 255;
  pixel[1] = G * 255;
  pixel[2] = B * 255;
}

// Calculate HSV color in range [0,1]
// Based off basic principles outlines in https://www.mrob.com/pub/muency/color.html
void CalculateColor(int dwell,  double distance, double mag_z, int max_iterations, double pixel_spacing, unsigned char *pixel) {
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
  double dwell_scaled = log(dwell)/log(max_iterations);

  // Remap the scaled dwell onto Hue and Saturation
  if(dwell_scaled < 0.5) {
    dwell_scaled = 1.0 - 2.0*dwell_scaled;
    H = 1.0 - dwell_scaled;
    S = sqrt(dwell_scaled);
  } else {
    dwell_scaled = 1.5*dwell_scaled - 0.5;
    H = dwell_scaled;
    S = sqrt(dwell_scaled);
  }

  // Go around color wheel 10x to cover range 1->max_iterations
  H *= 10.0;
  H -= floor(H);
  S -= floor(S);

  // log scale dwell for next band
  dwell_scaled = log(dwell+1)/log(max_iterations);
  double Hp1, Sp1;
  // Remap the scaled dwell onto Hue and Saturation
  if(dwell_scaled < 0.5) {
    dwell_scaled = 1.0 - 2.0*dwell_scaled;
    Hp1 = 1.0 - dwell_scaled;
    Sp1 = sqrt(dwell_scaled);
  } else {
    dwell_scaled = 1.5*dwell_scaled - 0.5;
    Hp1 = dwell_scaled;
    Sp1 = sqrt(dwell_scaled);
  }

  // Go around color wheel 10x to cover range 1->max_iterations
  Hp1 *= 10.0;
  Hp1 -= floor(Hp1);
  Sp1 -= floor(Sp1);

  // Continuously vary between dwell bands in range 0,1
  double frac = 1.0 - log2(log(mag_z)/log(1<<18));

  // Take the shortest path between angles on the color wheel
  // Without this crossing over 0 will give eronously large deltas
  double delta_H = (Hp1 - H);
  double mag_delta_H = fabs(delta_H);               
  delta_H = fmin(1.0 - mag_delta_H, mag_delta_H); // Find shortest path around color wheel

  double linear_interp_H = delta_H * frac;
  double linear_interp_S = (Sp1 - S) * frac;

  // Add interpolated delta value
  H += linear_interp_H;
//  S += linear_interp_S;

  // Wrap values around 0
  if(H < 0.0) {
    H = 1.0 + H;
  } else if(H > 1.0) {
    H = H - 1.0;
  }

  // Convert to RGB and set pixel value
  HSVtoRGB(H, S, V, pixel);
}

// Distance estimator and continuous dwell algorithm
// https://www.mrob.com/pub/muency/distanceestimator.html
void CalculatePixel(double xO, double yO, double pixel_size, unsigned char *pixel) {
  const int max_iterations = 10000;

  double x=0.0,y=0.0;
  double dx=0.0,dy=0.0;
  int dwell = 0;
  const double escape_radius = 1<<18; // Increased from 2.0 to find better distance estimate
  double distance = 0.0;
  double continuous_dwell = 0.0;

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
 
   const double mag_z = sqrt(x*x + y*y);
 
  // Calculate the distance if the orbit escaped
  if((x*x + y*y) >= (escape_radius*escape_radius)) {
    const double mag_dz = sqrt(dx*dx + dy*dy);
    distance = log(mag_z*mag_z) * mag_z / mag_dz;
  }

  // Calculate color based on dwell and distance
  CalculateColor(dwell, distance, mag_z,  max_iterations, pixel_size, pixel);
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
  const double pixel_size = length_x/2000.0;
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