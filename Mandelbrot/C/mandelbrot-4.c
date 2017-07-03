// Unoptimized Mandelbrot set using Dwell and Distance Estimator methods
// Continuous color

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct pixel {
  unsigned char R;
  unsigned char G;
  unsigned char B;
} pixel_t;

// Convert HSV [0,1] to RGB [0,1]
// https://en.wikipedia.org/wiki/HSL_and_HSV#Converting_to_RGB
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

double CalculateHue(double dwell_scaled) {
  double H;

  // Remap the scaled dwell onto Hue and Saturation
  if(dwell_scaled < 0.5) {
    dwell_scaled = 1.0 - 2.0*dwell_scaled;
    H = 1.0 - dwell_scaled;
  } else {
    dwell_scaled = 1.5*dwell_scaled - 0.5;
    H = dwell_scaled;
  }

  // Go around color wheel 10x to cover range 1->max_iterations
  H *= 10.0;
  H -= floor(H);

  return H;
}

double CalculateSaturation(double dwell_scaled) {
  double S;

  // Remap the scaled dwell onto Saturation
  S = sqrt(dwell_scaled);
  S -= floor(S);
  return S;
}

// Calculate HSV color in range [0,1]
// Based off basic principles outlines in https://www.mrob.com/pub/muency/color.html
void CalculateColor(int dwell,  double distance, double mag_z, double escape_radius, int max_iterations, double pixel_spacing, pixel_t *pixel) {
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
  double dwell_scaled_plus = log(dwell+1)/log(max_iterations);

  // Calculate current and next dwell band values
  H = CalculateHue(dwell_scaled);
  S = CalculateSaturation(dwell_scaled);
  double H_plus = CalculateHue(dwell_scaled_plus);
  double S_plus = CalculateSaturation(dwell_scaled_plus);

  // Continuously vary between dwell bands in range 0,1
  double escape_interpolation = 1.0 - log2(log(mag_z)/log(escape_radius));

  // Take the shortest path between angles on the color wheel
  // Without this crossing over 0 will give eronously large deltas
  double delta_H = (H_plus - H);
  double mag_delta_H = fabs(delta_H);               
  delta_H = fmin(1.0 - mag_delta_H, mag_delta_H); // Find shortest path around color wheel

  double linear_interp_H = delta_H * escape_interpolation;
  double linear_interp_S = (S_plus - S) * escape_interpolation;

  // Add interpolated delta value
  H += linear_interp_H;
  S += linear_interp_S;

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
void CalculatePixel(double xO, double yO, double pixel_size, pixel_t *pixel) {
  const int max_iterations = 10000;
  const double escape_radius = 1<<18;

  double x=0.0,y=0.0;
  double dx=0.0,dy=0.0;
  int dwell = 0;
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
  CalculateColor(dwell, distance, mag_z, escape_radius, max_iterations, pixel_size, pixel);
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
