#include <iostream>
#include <Kokkos_Core.hpp>
#include <Kokkos_Parallel.hpp>
#include <Kokkos_View.hpp>
#include <assert.h>
#include <limits> 

int main(int argc, char **argv) {
  // Initialize Kokkos
  Kokkos::initialize(argc, argv);

  const int dim = 1000;

  // Allocate our arrays
  Kokkos::View<double*> a("a", dim);    
  Kokkos::View<double*> b("b", dim);
  Kokkos::View<double*> c("c", dim);

  // Create host mirror of C
  auto  c_mirror = Kokkos::create_mirror_view(c);

  // Initialize a and b
  Kokkos::parallel_for(dim, KOKKOS_LAMBDA(int i) {
      a(i) = sin((double)i)*sin((double)i);
      b(i) = cos((double)i)*cos((double)i);    
  });
 
  // Preform the vector addition
  Kokkos::parallel_for(dim, KOKKOS_LAMBDA(int i) {
    c(i) = a(i) + b(i);
  });

  // Update the mirror
  deep_copy(c_mirror, c);

  // Verify the results of the host mirror
  for(int i=0; i<dim; i++) {
    double eps = abs(1.0 - c_mirror(i));
    assert(eps <= std::numeric_limits<double>::epsilon());    
  };

  std::cout<<"Sum Correct\n";

  Kokkos::finalize();
  return 0;
}
