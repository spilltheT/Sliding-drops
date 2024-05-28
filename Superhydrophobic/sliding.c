#include "navier-stokes/centered.h"
#define FILTERED
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "reduced.h"
#include "distance.h"
#include "adapt_wavelet_limited_v2.h"

int MAXlevel = 12; 
#define MINlevel 3 
#define Ldomain 40
#define tmax 200
#define tsnap (0.05)

// Error tolerances
#define fErr (1e-3) // error tolerance in VOF
#define KErr (1e-4) // error tolerance in KAPPA
#define OmegaErr (1e-3) // error tolerances in vorticity
#define velErr (1e-4) // error tolerances in velocity

#define Rho21 (1e-3) // density ratio between the gas and drop

// Boundary conditions
f[left] = dirichlet(0.0);
u.t[left] = dirichlet(0.0);

double BondX, BondY;

int main(int argc, char const *argv[]){
  periodic (bottom);
  L0=Ldomain;  // Domain length
  origin(0.,0.);  // Origin
  init_grid (1 << (5));
  BondX = 1.0;  // Bond_x
  BondY = 0.5;  // Bond_y
  G.x = -BondX;
  G.y = -BondY;
  f.sigma = 1.;
  mu1 = 1.0; mu2 = 0.0001;  // Viscosity of drop, air
  rho1 = 1.; rho2 = Rho21;
  run();
}

// Refined region
int refRegion(double x, double y, double z){

  return (x < 2.0 && y < 38.5 && y > 0.0 ? MAXlevel:
          x < 3.0 && y < 40.0 && y > 0.0 ? MAXlevel-1:
          x < 4.0 && y < 40.0 && y > 0.0 ? MAXlevel-2:
          x < 6.0 && y < 40.0 && y > 0.0 ? MAXlevel-3:
          MAXlevel-4
        );
}

// Initializing drop profile
event init(t = 0)
{
  if(!restore (file = "dump")){

    char comm[80];
    sprintf (comm, "mkdir -p intermediate");
    system(comm);

    char filename[60];
    sprintf(filename,"Bo0.1.dat");
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

    while (adapt_wavelet_limited ((scalar *){f, d}, (double[]){1e-8, 1e-8}, refRegion).nf);

    vertex scalar phi[];
    foreach_vertex(){
      phi[] = -(d[] + d[-1] + d[0,-1] + d[-1,-1])/4.;
    }
    fractions (phi, f);

    foreach () {
      u.x[] = 0.0;
      u.y[] = 0.0;
    }
    boundary((scalar *){f, u.x, u.y});
  }
}

// Log files
event logWriting (i=i+25) {
  double ke = 0.;
  foreach (reduction(+:ke)){
    ke += (0.5*rho(f[])*(sq(u.x[]) + sq(u.y[])));
  }
  static FILE * fp;
  if (i == 0) {
    fprintf (ferr, "i dt t ke MaxLevel\n");
    fp = fopen ("log", "w");
    fprintf (fp, "i dt t ke\n");
    fprintf (fp, "%d %g %g %g\n", i, dt, t, ke);
    fclose(fp);
  } else {
    fp = fopen ("log", "a");
    fprintf (fp, "%d %g %g %g\n", i, dt, t, ke);
    fclose(fp);
  }
  fprintf (ferr, "%d %g %g %g %d\n", i, dt, t, ke, MAXlevel);

}

// Adaptive Mesh Refinement
event adapt(i++){
  scalar KAPPA[], omega[], uG[];
  curvature(f, KAPPA);
  vorticity(u,omega);
  foreach(){
    omega[] *= (f[]);
  } 
  foreach(){
    uG[] = u.y[]*f[];
  }
  
  boundary ((scalar *){KAPPA, omega, uG});
  adapt_wavelet ((scalar *){f, KAPPA, omega, uG},
    (double[]){fErr, KErr, OmegaErr, velErr},
    MAXlevel, MINlevel);
  adapt_wavelet_limited ((scalar *){f, KAPPA, omega, uG},
    (double[]){fErr, KErr, OmegaErr, velErr},
    refRegion, MINlevel);
}

// Dumping files
event writingFiles (t = 0; t += tsnap; t <= tmax) {
  dump (file = "dump");
  char nameOut[80];
  sprintf (nameOut, "intermediate/snapshot-%5.4f", t);
  dump (file = nameOut);
}
