#include "navier-stokes/centered.h"
#define FILTERED
#include "three-phase.h"
#include "tension.h"
#include "distance.h"
#include "adapt_wavelet_limited_v2.h"

#define MAXlevel 12                                             // maximum level
#define MINlevel 3                                               // maximum level
#define Ldomain 40
#define tmax 100
#define tsnap (0.05)

// Error tolerances
#define fErr (1e-3) // error tolerance in VOF
#define KErr (1e-4) // error tolerance in KAPPA
#define VelErr (1e-2) // error tolerances in velocity
#define OmegaErr (1e-3) // error tolerances in vorticity

double Oh1, Oh2, Oh3;
double RHO21, RHO31;
double BoX, BoY; 

// Boundary conditions
u.t[left] = dirichlet(0.0);
f1[left] = dirichlet(1.0);
f2[left] = dirichlet(0.0);

int main(int argc, char const *argv[]) {
  
  Oh1 = 0.1;  // Film Ohnesorge
  Oh2 = 1.0;  // Drop Ohnesorge
  Oh3 = 0.0001;  // air has a smaller density

  BoX = 1.0;  // Bond_x
  BoY = 0.5;  // Bond_y

  RHO21 = 1.0;  // rho2/rho1
  RHO31 = 0.001;  // rho3/rho1

  L0=Ldomain;  // Length of domain
  X0=0.; Y0=0.;  // Origin
  init_grid (1 << (12));

  char comm[80];
  sprintf (comm, "mkdir -p intermediate");
  system(comm);

  rho1 = 1.0; mu1 = Oh1;
  rho2 = RHO21; mu2 = Oh2;
  rho3 = RHO31; mu3 = Oh3;

  f1.sigma = 0.8;  // gamma1/gamma2
  f2.sigma = 1.0;  // gamma2/gamma2 = 1.0

  run();
}


event acceleration(i++) {
  face vector av = a;
  foreach_face(x){
    av.x[] -= BoX;
  }
  
  foreach_face(y){
    av.y[] -= BoY;
  }
}

// Refined region
int refRegion(double x, double y, double z){

  return (x < 0.2 && y < 40.0 && y > 0.0 ? MAXlevel:
          x < 1.3 && y < 40.0 && y > 0.0 ? MAXlevel-1:
          x < 2.0 && y < 42.0 && y > -2.0 ? MAXlevel-2:
          x < 3.0 && y < 43.0 && y > -3.0 ? MAXlevel-3:
          x < 4.0 && y < 45.0 && y > -5.0 ? MAXlevel-4:
          MAXlevel-5
        );
}

// Initialize drop and film profile
event init(t = 0){

  if(!restore (file = "dump")){
    char filename[60];
    sprintf(filename,"f1.dat");
    FILE * fp = fopen(filename,"rb");
    if (fp == NULL){
      fprintf(ferr, "There is no file named %s\n", filename);
      return 1;
    }
    coord* InitialShape;
    InitialShape = input_xy(fp);
    fclose (fp);
    scalar d[];
    distance (d, InitialShape);
    while (adapt_wavelet_limited ((scalar *){f1, d}, (double[]){1e-8, 1e-8}, refRegion).nf);
    vertex scalar phi[];
    foreach_vertex(){
      phi[] = (d[] + d[-1] + d[0,-1] + d[-1,-1])/4.;
    }
    fractions (phi, f1);
    sprintf(filename,"f2.dat");
    fp = fopen(filename,"rb");
    if (fp == NULL){
      fprintf(ferr, "There is no file named %s\n", filename);
      return 1;
    }
    InitialShape = input_xy(fp);
    fclose (fp);
    distance (d, InitialShape);
    
    while (adapt_wavelet_limited ((scalar *){f1, d}, (double[]){1e-8, 1e-8}, refRegion).nf);
    foreach_vertex(){
      phi[] = -(d[] + d[-1] + d[0,-1] + d[-1,-1])/4.;
    }
    fractions (phi, f2);
    
    foreach () {
      u.x[] = 0.0;
      u.y[] = 0.0;
    }
    boundary((scalar *){f1, f2, u.x, u.y});
  }
}

// Adaptive Mesh Refinement
event adapt(i++) {

  scalar KAPPA1[], KAPPA2[];
  curvature(f1, KAPPA1);
  curvature(f2, KAPPA2);
  boundary ((scalar *){KAPPA1, KAPPA2});
  adapt_wavelet_limited ((scalar *){f1, f2, u.x, u.y, KAPPA1, KAPPA2},
    (double[]){fErr, fErr, VelErr, VelErr, KErr, KErr},
    refRegion, MINlevel);
}

// Dumping files
event writingFiles (t = 0; t += tsnap; t <= tmax + tsnap) {
  dump (file = "dump");
  char nameOut[80];
  sprintf (nameOut, "intermediate/snapshot-%5.4f", t);
  dump (file = nameOut);
}

// Log files
event logWriting (i+=10) {
  
  double ke = 0.;
  foreach (reduction(+:ke)){
    ke += f2[]*(sq(u.x[]) + sq(u.y[]))*(2*pi*y)*sq(Delta);
  }
  static FILE * fp;
  if (i == 0) {
    fprintf (ferr, "i dt t ke\n");
    fp = fopen ("log", "w");
    fprintf (fp, "i dt t ke\n");
    fprintf (fp, "%d %g %g %g\n", i, dt, t, ke);
    fclose(fp);
  } else {
    fp = fopen ("log", "a");
    fprintf (fp, "%d %g %g %g\n", i, dt, t, ke);
    fclose(fp);
  }
  fprintf (ferr, "%d %g %g %g\n", i, dt, t, ke);

}
