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

double CalculateHue(double dwell_scaled) {
  double H;

  // Remap the scaled dwell onto Hue
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
void CalculateColor(int dwell,  double distance, double mag_z, double escape_radius, int max_iterations, double pixel_spacing, unsigned char *pixel) {
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

  // Lighten every other stripe
  if(dwell%2) {
    V *= 0.95;
  }

  // log scale dwell
  double dwell_scaled = log(dwell)/log(max_iterations);

  // Calculate current and next dwell band values
  H = CalculateHue(dwell_scaled);
  S = CalculateSaturation(dwell_scaled);

  // Convert to RGB and set pixel value
  HSVtoRGB(H, S, V, pixel);
}

// Distance estimator and continuous dwell algorithm
// https://www.mrob.com/pub/muency/distanceestimator.html
void CalculatePixel(double xO, double yO, double pixel_size, unsigned char *pixel) {
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
  // Animation constraints
  int frame_count = 200;
  const double center_x = -0.745429;
  const double center_y =  0.113008;
  const double length_x_begin = 2.0;
  const double length_y_begin = 2.0;
  const double length_x_end =  0.00001;
  const double length_y_end =  0.00001;
  const double delta_length = (length_x_end - length_x_begin)/(double)frame_count;

  double length_x = length_x_begin;
  double length_y = length_y_begin;
  double pixel_size = length_x/2000.0;
  const int pixels_x = length_x / pixel_size;
  const int pixels_y = length_y / pixel_size;

  // Linearized 2D image data packed in RGB format in range [0-255]
  size_t pixel_bytes = sizeof(unsigned char)*3*pixels_x*pixels_y;
  unsigned char *pixels = malloc(pixel_bytes);

  char file_name[100];

  for(int frame=0; frame<frame_count; frame++) {
    // Convenience variables based on image bounds
    const double x_min = center_x - length_x/2.0;
    const double x_max = center_x + length_x/2.0;
    const double y_min = center_y - length_y/2.0;
    const double y_max = center_y - length_y/2.0;

    printf("%f, %f, %f, %f\n", x_min, x_max, y_min, y_max);

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
    sprintf(file_name, "mandelbrot%03d.ppm", frame);
    FILE *file = fopen(file_name, "wb");
    fprintf(file, "P6\n%d %d\n%d\n", pixels_x, pixels_y, 255);
    fwrite(pixels, sizeof(unsigned char), 3*pixels_x*pixels_y, file);
    printf("Saved %d by %d image %d of %d\n", pixels_x, pixels_y, frame+1, frame_count);
    fclose(file);

    // "Zoom" in image
    length_x += delta_length;
    length_y += delta_length;
    pixel_size = length_x/2000.0;
  }

  free(pixels);

  return 0;
}
