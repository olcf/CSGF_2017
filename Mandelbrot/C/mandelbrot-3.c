// A basic unoptimized mandelbrot set using continuous dwell  and distance estimator
// Gray scale visulization

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Convert HSV [0,1] to RGB [0,1]
void HSVtoRGB(float H, float S, float V, unsigned char *pixel) {
  // convert hue to degrees
  H *= 360.0;

  if(S <= 0.0) {
    pixel[0] = (unsigned char)(V * 255.0);
    pixel[1] = (unsigned char)(V * 255.0);
    pixel[2] = (unsigned char)(V * 255.0);
    return;
  }

  float R = 0.0;
  float G = 0.0;
  float B = 0.0;
  float C = (V*S);
  float H_p = H/60.0;
  float X = C * (1-abs(((int)H_p%2)-1));

  if(H_p >= 0.0 && H_p <= 1.0) {
    R = C;
    G = X;
    B = 0.0;
  }
  else if(H_p > 1.0 && H_p <= 2.0) {
    R = X;
    G = C;
    B = 0.0;
  } 
  else if(H_p > 2.0 && H_p <= 3.0) {
    R = 0.0;
    G = C;
    B = X;
  }
  else if(H_p > 3.0 && H_p <= 4.0) {
    R = 0.0;
    G = X;
    B = C;
  }
  else if(H_p > 4.0 && H_p <= 5.0) {
    R = X;
    G = 0.0;
    B = C;
  }
  else if(H_p > 5.0 && H_p <= 6.0) {
    R = C;
    G = 0.0;
    B = X;
  }

  float m = V - C;
  R += m;
  G += m;
  B += m;

  unsigned char rc,gc,bc;

  rc = (unsigned char)(R * 255.0);
  gc = (unsigned char)(G * 255.0);
  bc = (unsigned char)(B * 255.0);

  pixel[0] = rc;
  pixel[1] = gc;
  pixel[2] = bc;
}

// Calculate HSV color in range [0,1]
// https://www.mrob.com/pub/muency/color.html
void CalculateColor(double continuous_dwell, double distance, int max_iterations, double pixel_spacing, unsigned char *pixel) {
  float H,S,V;

  // Split dwell into scalar and fractional pieces
  const int dwell = floor(continuous_dwell);
  const double final_rad = continuous_dwell - (double)dwell;

  // Point is within Mandelbrot set, color white
  if(dwell >= max_iterations) {
    H = 0.0;
    S = 0.0;
    V = 1.0;
    HSVtoRGB(H,S,V, pixel);
    return;
  }

  // log scale distance
  const double dist_scaled = log2(distance / pixel_spacing);

  // Convert scaled distance to Value in 8 intervals
  if (dist_scaled > 0.0) {
    V = 1.0;
  } else if (dist_scaled > -8.0) {
    V = (8.0 + dist_scaled) / 8.0;
  } else {
    V = 0.0;
  }

  // log scale dwell
  double dwell_scaled = log(dwell)/log(100*max_iterations);

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
    V *= 0.90;
    S *= 0.8;
  }

  // Break the stripes up depending on angle of orbit at escape
//  if(final_angle > M_PI) {
//    H += 0.02;
//   }

  // Break square into full color gradient
  H += 0.0001 * final_rad;

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
void CalculatePixel(double xO, double yO, double pixel_size, unsigned char *pixel) {
  const int max_iterations = 100;

  double x=0.0,y=0.0;
  double dx=0.0,dy=0.0;
  int iteration = 0;
  const double escape_radius = 10000; // Increased from 2.0 to find better distance estimate
  double distance = 0.0;
  double continuous_dwell = 0.0;

  while( (x*x + y*y) < escape_radius*escape_radius && iteration < max_iterations) {
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
    iteration++;
  }

  // Calculate the continuous dwell
  continuous_dwell = iteration;

  // Calculate the distance if the orbit escaped
  if((x*x + y*y) >= escape_radius*escape_radius) {
    const double mag_z = sqrt(x*x + y*y);
    const double mag_dz = sqrt(dx*dx + dy*dy);
    distance = log(mag_z*mag_z) * mag_z / mag_dz;
    continuous_dwell += log2(log2(mag_z)) - log2(log2(escape_radius));
  }

  // Calculate color based on dwell and distance
  CalculateColor(continuous_dwell, distance, max_iterations, pixel_size, pixel);
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
