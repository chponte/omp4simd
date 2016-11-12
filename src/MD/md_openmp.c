# include <stdlib.h>
# include <stdio.h>
# include <time.h>
# include <math.h>
# include <omp.h>
# include <string.h>

// Include appropiate library depending on compiler
#if defined(__INTEL_COMPILER)
#include <malloc.h>
#else
#include <mm_malloc.h>
#endif

# define NP 1024
# define ALIGNMENT 64

int main ( int argc, char *argv[] );
void compute ( int np, int nd, const int nppadded, double pos[], double vel[],
  double mass, double f[], double *pot, double *kin );
double pt_dist(int nd, int nppadded, double pos[], int p1, int p2, double dr[]);
double min(double val[], int pos, double value);
double f_calc(int nppadded, double dr[], double d[], double d2[], int i, int j, int k);
void initialize ( int np, int nd, int nppadded, double box[], int *seed, double pos[],
  double vel[], double acc[] );
double r8_uniform_01 ( int *seed );
void timestamp ( void );
void update ( int np, int nd, int nppadded, double pos[], double vel[], double f[],
  double acc[], double mass, double dt );
void write_to_file(const int np, const int nd, int nppadded, double *acc,
  double *force, double *pos, double *vel, char *outfile);

/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for MD_OPENMP.

  Discussion:

    MD implements a simple molecular dynamics simulation.

    The program uses Open MP directives to allow parallel computation.

    The velocity Verlet time integration scheme is used. 

    The particles interact with a central pair potential.

  Licensing:

    This code is distributed under the GNU GPLv3 license. 

  Modified:

    10 November 2016

  Author:

    Original FORTRAN77 version by Bill Magro.
    C version by John Burkardt.
    OpenMP SIMD version by Christian Ponte.

  Parameters:

    None
*/
{
  const int divisor = ALIGNMENT>sizeof(double) ? ALIGNMENT / sizeof(double) : 1;
  const int nppadded = (NP+divisor-1)/divisor * divisor;
  double *acc;
  double *box;
  double dt = 0.0001;
  double e0;
  double *force;
  int i;
  double kinetic;
  double mass = 1.0;
  int nd = 3;
  int np = NP;
  double *pos;
  double potential;
  int proc_num;
  int seed = 123456789;
  int step;
  int step_num = 400;
  int step_print;
  int step_print_index;
  int step_print_num;
  double *vel;
  double wtime;
  char * outfile = NULL;

  timestamp ( );

  proc_num = omp_get_num_procs ( );

  acc = (double *) _mm_malloc(nd*nppadded * sizeof(double), ALIGNMENT);
  box = (double *) malloc(nd * sizeof(double));
  force = (double *) _mm_malloc(nd*nppadded * sizeof(double), ALIGNMENT);
  pos = (double *) _mm_malloc(nd*nppadded * sizeof(double), ALIGNMENT);
  vel = (double *) _mm_malloc(nd*nppadded * sizeof(double), ALIGNMENT);

  printf ( "\n" );
  printf ( "MD_OPENMP\n" );
  printf ( "  C/OpenMP version\n" );
  printf ( "\n" );
  printf ( "  A molecular dynamics program.\n" );

  printf ( "\n" );
  printf ( "  NP, the number of particles in the simulation is %d\n", np );
  printf ( "  STEP_NUM, the number of time steps, is %d\n", step_num );
  printf ( "  DT, the size of each time step, is %f\n", dt );

  printf ( "\n" );
  printf ( "  Number of processors available = %d\n", omp_get_num_procs ( ) );
  printf ( "  Number of threads =              %d\n", omp_get_max_threads ( ) );

  if (argc > 1){
    outfile = (char *) malloc(sizeof(char) * strlen(argv[1]));
    strcpy(outfile, argv[1]);
  }

/*
  Set the dimensions of the box.
*/
  for ( i = 0; i < nd; i++ )
  {
    box[i] = 10.0;
  }

  printf ( "\n" );
  printf ( "  Initializing positions, velocities, and accelerations.\n" );
/*
  Set initial positions, velocities, and accelerations.
*/
  initialize ( np, nd, nppadded, box, &seed, pos, vel, acc );
/*
  Compute the forces and energies.
*/
  printf ( "\n" );
  printf ( "  Computing initial forces and energies.\n" );

  compute ( np, nd, nppadded, pos, vel, mass, force, &potential, &kinetic );

  e0 = potential + kinetic;
/*
  This is the main time stepping loop:
    Compute forces and energies,
    Update positions, velocities, accelerations.
*/
  printf ( "\n" );
  printf ( "  At each step, we report the potential and kinetic energies.\n" );
  printf ( "  The sum of these energies should be a constant.\n" );
  printf ( "  As an accuracy check, we also print the relative error\n" );
  printf ( "  in the total energy.\n" );
  printf ( "\n" );
  printf ( "      Step      Potential       Kinetic        (P+K-E0)/E0\n" );
  printf ( "                Energy P        Energy K       Relative Energy Error\n" );
  printf ( "\n" );

  step_print = 0;
  step_print_index = 0;
  step_print_num = 10;

  step = 0;
  printf ( "  %8d  %14f  %14f  %14e\n",
    step, potential, kinetic, ( potential + kinetic - e0 ) / e0 );
  step_print_index = step_print_index + 1;
  step_print = ( step_print_index * step_num ) / step_print_num;

  wtime = omp_get_wtime ( );

  for ( step = 1; step <= step_num; step++ )
  {
    compute ( np, nd, nppadded, pos, vel, mass, force, &potential, &kinetic );

    if ( step == step_print )
    {
      printf ( "  %8d  %14f  %14f  %14e\n", step, potential, kinetic,
       ( potential + kinetic - e0 ) / e0 );
      step_print_index = step_print_index + 1;
      step_print = ( step_print_index * step_num ) / step_print_num;
    }
    update ( np, nd, nppadded, pos, vel, force, acc, mass, dt );
  }
  wtime = omp_get_wtime ( ) - wtime;

  printf ( "\n" );
  printf ( "  Elapsed seconds = %f\n", wtime );

  if (outfile != NULL){
    printf( "  Writting result to %s\n", outfile);
    write_to_file(np, nd, nppadded, acc, force, pos, vel, outfile);
    free(outfile);
  }

  _mm_free(acc);
  free(box);
  _mm_free(force);
  _mm_free(pos);
  _mm_free(vel);
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "MD_OPENMP\n" );
  printf ( "  Normal end of execution.\n" );

  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void compute ( int np, int nd, const int nppadded, double pos[], double vel[],
  double mass, double f[], double *pot, double *kin )

/******************************************************************************/
/*
  Purpose:

    COMPUTE computes the forces and energies.

  Discussion:

    The computation of forces and energies is fully parallel.

    The potential function V(X) is a harmonic well which smoothly
    saturates to a maximum value at PI/2:

      v(x) = ( sin ( min ( x, PI2 ) ) )**2

    The derivative of the potential is:

      dv(x) = 2.0 * sin ( min ( x, PI2 ) ) * cos ( min ( x, PI2 ) )
            = sin ( 2.0 * min ( x, PI2 ) )

  Licensing:

    This code is distributed under the GNU GPLv3 license. 

  Modified:

    10 November 2016

  Author:

    Original FORTRAN77 version by Bill Magro.
    C version by John Burkardt.
    OpenMP SIMD version by Christian Ponte.

  Parameters:

    Input, int NP, the number of particles.

    Input, int ND, the number of spatial dimensions.

    Input, const int NPPADDED, NP padded to a multiple of ALIGNMENT.

    Input, double POS[ND*NP], the position of each particle.

    Input, double VEL[ND*NP], the velocity of each particle.

    Input, double MASS, the mass of each particle.

    Output, double F[ND*NP], the forces.

    Output, double *POT, the total potential energy.

    Output, double *KIN, the total kinetic energy.
*/
{
  double *d;
  double *d2;
  double *rij;
  int i;
  int j;
  int k;
  double ke;
  double pe;
  double PI2 = 3.141592653589793 / 2.0;

  pe = 0.0;
  ke = 0.0;

  # pragma omp parallel \
  shared ( f, nd, np, pos, vel ) \
  private ( i, j, k, rij, d, d2 )
  {
    d = (double *) _mm_malloc(nppadded * sizeof(double), ALIGNMENT);
    d2 = (double *) _mm_malloc(nppadded * sizeof(double), ALIGNMENT);
    rij = (double *) _mm_malloc(3*nppadded * sizeof(double), ALIGNMENT);

    # pragma omp for reduction ( + : pe, ke )
    for ( k = 0; k < np; k++ )
    {
  /*
    Compute the potential energy and forces.
  */
      # pragma omp simd reduction(+:pe) aligned(d:ALIGNMENT, pos:ALIGNMENT, rij:ALIGNMENT, d2:ALIGNMENT)
      for (j=0; j<np; j++){
        d[j] = pt_dist(nd, nppadded, pos, k, j, rij);
        d2[j] = min(d, j, PI2);
        /*
          Attribute half of the potential energy to particle J.
        */
        pe += 0.5 * pow(sin(d2[j]), 2);
      }

      # pragma omp simd reduction(+:ke) aligned(f:ALIGNMENT, vel:ALIGNMENT)
      for ( i = 0; i < nd; i++ ){
        f[i*nppadded+k] = 0.0;
        /*
          Compute the kinetic energy.
        */
        ke += vel[i*nppadded+k] * vel[i*nppadded+k];
      }

      for ( i = 0; i < nd; i++ ){
        # pragma omp simd aligned(f:ALIGNMENT, rij:ALIGNMENT, d:ALIGNMENT, d2:ALIGNMENT)
        for ( j = 0; j < np; j++ ){
          f[i*nppadded+k] -= f_calc(nppadded, rij, d, d2, i, j, k);
        }
      }
    }

    _mm_free(d);
    _mm_free(d2);
    _mm_free(rij);
  }

  ke = ke * 0.5 * mass;

  *pot = pe;
  *kin = ke;

  return;
}
/******************************************************************************/

# pragma omp declare simd notinbranch uniform(nd, nppadded, pos, dr, p1) linear(p2:1) \
  aligned(pos:ALIGNMENT, dr:ALIGNMENT)
double pt_dist(int nd, int nppadded, double pos[], int p1, int p2, double dr[])

/******************************************************************************/
/*
  Purpose:

    PT_DIST computes the displacement (and its norm) between two particles.

  Licensing:

    This code is distributed under the GNU GPLv3 license. 

  Modified:

    10 November 2016

  Author:

    Original FORTRAN77 version by Bill Magro.
    C version by John Burkardt.
    OpenMP SIMD version by Christian Ponte.

  Parameters:

    Input, int ND, the number of spatial dimensions.

    Input, const int NPPADDED, NP padded to a multiple of ALIGNMENT. 

    Input, double POS[], address of the position array.

    Input, int P1, P2, offset of the particles inside the position array.

    Output, double DR[ND], the displacement vector.

    Output, double D, the Euclidean norm of the displacement.
*/
{
  double d;
  int i;

  d = 0.0;
  for ( i = 0; i < nd; i++ )
  {
    dr[i*nppadded+p2] = pos[i*nppadded+p1] - pos[i*nppadded+p2];
    d = d + dr[i*nppadded+p2] * dr[i*nppadded+p2];
  }
  return sqrt(d);
}
/******************************************************************************/

# pragma omp declare simd notinbranch uniform(val, value) linear(pos:1) \
  aligned(val:ALIGNMENT)
double min(double val[], int pos, double value)

/******************************************************************************/
/*
  Purpose:

    MIN calculates the minimum between two values.

  Licensing:

    This code is distributed under the GNU GPLv3 license. 

  Modified:

    10 November 2016

  Author: Christian Ponte.

  Parameters:

    Input, double VAL[], the address of the value array.

    Input, int POS, the offset of the value inside the array. 

    Input, double VALUE, the second value to compare.

    Output, double, the minimum between the two values.
*/
{
  return val[pos] < value ? val[pos] : value;
}
/******************************************************************************/

# pragma omp declare simd notinbranch uniform(nppadded, dr, d, d2, i, k) linear(j:1) \
  aligned(dr:ALIGNMENT, d:ALIGNMENT, d2:ALIGNMENT)
double f_calc(int nppadded, double dr[], double d[], double d2[], int i, int j, int k)
/******************************************************************************/
/*
  Purpose:

    F_CALC calculates the force induced by one particle onto another based on their
    distance.

  Licensing:

    This code is distributed under the GNU GPLv3 license. 

  Modified:

    10 November 2016

  Author: Christian Ponte.

  Parameters:

    Input, int NPPADDED, NP padded to a multiple of ALIGNMENT.

    Input, double DR[ND], the displacement vector. 

    Input, double d[NPPADDED], the distances vector.

    Input, double d2[NPPADDED], the distances vector considering a distance threshold.

    Input, int I, int J, int K, the dimension number and the two particles number.

    Output, double, the resulting force.
*/
{
  return j==k ? 0 : (dr[i*nppadded+j] * sin (2.0 * d2[j]) / d[j]);
}

/******************************************************************************/

void initialize ( int np, int nd, int nppaddded, double box[], int *seed,
  double pos[], double vel[], double acc[] )

/******************************************************************************/
/*
  Purpose:

    INITIALIZE initializes the positions, velocities, and accelerations.

  Licensing:

    This code is distributed under the GNU GPLv3 license. 

  Modified:

    10 November 2016

  Author:

    Original FORTRAN77 version by Bill Magro.
    C version by John Burkardt.
    OpenMP SIMD version by Christian Ponte.

  Parameters:

    Input, int NP, the number of particles.

    Input, int ND, the number of spatial dimensions.

    Input, const int NPPADDED, NP padded to a multiple of ALIGNMENT. 

    Input, double BOX[ND], specifies the maximum position
    of particles in each dimension.

    Input, int *SEED, a seed for the random number generator.

    Output, double POS[ND*NP], the position of each particle.

    Output, double VEL[ND*NP], the velocity of each particle.

    Output, double ACC[ND*NP], the acceleration of each particle.
*/
{
  int i;
  int j;
/*
  Give the particles random positions within the box.
*/
  for ( i = 0; i < nd; i++ )
  {
    for ( j = 0; j < np; j++ )
    {
      pos[i*nppaddded+j] = box[i] * r8_uniform_01 ( seed );
    }
  }

  for ( i = 0; i < nd; i++ )
  {
    for ( j = 0; j < np; j++ )
    {
      vel[i*nppaddded+j] = 0.0;
    }
  }
  for ( i = 0; i < nd; i++ )
  {
    for ( j = 0; j < np; j++ )
    {
      acc[i*nppaddded+j] = 0.0;
    }
  }
  return;
}
/******************************************************************************/

double r8_uniform_01 ( int *seed )

/******************************************************************************/
/*
  Purpose:

    R8_UNIFORM_01 is a unit pseudorandom R8.

  Discussion:

    This routine implements the recursion

      seed = 16807 * seed mod ( 2**31 - 1 )
      unif = seed / ( 2**31 - 1 )

    The integer arithmetic never requires more than 32 bits,
    including a sign bit.

  Licensing:

    This code is distributed under the GNU GPLv3 license. 

  Modified:

    11 August 2004

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Springer Verlag, pages 201-202, 1983.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, pages 362-376, 1986.

  Parameters:

    Input/output, int *SEED, a seed for the random number generator.

    Output, double R8_UNIFORM_01, a new pseudorandom variate, strictly between
    0 and 1.
*/
{
  int k;
  double r;

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + 2147483647;
  }

  r = ( double ) ( *seed ) * 4.656612875E-10;

  return r;
}
/******************************************************************************/

void timestamp ( void )

/******************************************************************************/
/*
  Purpose:

    TIMESTAMP prints the current YMDHMS date as a time stamp.

  Example:

    31 May 2001 09:45:54 AM

  Licensing:

    This code is distributed under the GNU GPLv3 license. 

  Modified:

    24 September 2003

  Author:

    John Burkardt

  Parameters:

    None
*/
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  printf ( "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}
/******************************************************************************/

void update ( int np, int nd, int nppadded, double pos[], double vel[],
  double f[], double acc[], double mass, double dt )

/******************************************************************************/
/*
  Purpose:

    UPDATE updates positions, velocities and accelerations.

  Discussion:

    The time integration is fully parallel.

    A velocity Verlet algorithm is used for the updating.

    x(t+dt) = x(t) + v(t) * dt + 0.5 * a(t) * dt * dt
    v(t+dt) = v(t) + 0.5 * ( a(t) + a(t+dt) ) * dt
    a(t+dt) = f(t) / m

  Licensing:

    This code is distributed under the GNU GPLv3 license. 

  Modified:

    10 November 2016

  Author:

    Original FORTRAN77 version by Bill Magro.
    C version by John Burkardt.
    OpenMP SIMD version by Christian Ponte.

  Parameters:

    Input, int NP, the number of particles.

    Input, int ND, the number of spatial dimensions.

    Input, const int NPPADDED, NP padded to a multiple of ALIGNMENT. 

    Input/output, double POS[ND*NP], the position of each particle.

    Input/output, double VEL[ND*NP], the velocity of each particle.

    Input, double F[ND*NP], the force on each particle.

    Input/output, double ACC[ND*NP], the acceleration of each particle.

    Input, double MASS, the mass of each particle.

    Input, double DT, the time step.
*/
{
  int i;
  int j;
  double rmass;

  rmass = 1.0 / mass;

# pragma omp parallel \
  shared ( acc, dt, f, nd, np, pos, rmass, vel ) \
  private ( i, j )

# pragma omp for
  for ( i = 0; i < nd; i++ )
  {
    # pragma omp simd aligned(pos:ALIGNMENT, vel:ALIGNMENT, acc:ALIGNMENT, f:ALIGNMENT)
    for ( j = 0; j < np; j++ )
    {
      pos[i*nppadded+j] = pos[i*nppadded+j] + vel[i*nppadded+j] * dt + 0.5 * acc[i*nppadded+j] * dt * dt;
      vel[i*nppadded+j] = vel[i*nppadded+j] + 0.5 * dt * ( f[i*nppadded+j] * rmass + acc[i*nppadded+j] );
      acc[i*nppadded+j] = f[i*nppadded+j] * rmass;
    }
  }

  return;
}
/******************************************************************************/

void write_to_file(const int np, const int nd, int nppadded, double *acc,
  double *force, double *pos, double *vel, char *outfile)

/******************************************************************************/
/*
  Purpose:

    WRITE_TO_FILE writes the arrays ACC, FORCE, POS and VEL into the specified file in OUTFILE.

  Licensing:

    This code is distributed under the GNU GPLv3 license. 

  Modified:

    4 October 2016

  Author:

    Christian Ponte

  Parameters:

    Input, int np, the number of particles represented on the array.

    Input, int nd, the number of dimensions for each particle represented in the array.

    Input, int npadded, the padded horizontal dimension of the array.

    Input, double *acc, the accelerations array.

    Input, double *force, the forces array.

    Input, double *pos, the positions array.

    Input, double *vel, the velocities array.

    Input, char *outfile, the output file name.
*/
{
  int i,j;

  FILE *f = fopen(outfile, "w");
  if (f==NULL) return;

  for (i=0; i<np; i++){
    for (j=0; j<nd; j++){
      fprintf(f, "%f %f %f %f\n",
        acc[j*nppadded+i], force[j*nppadded+i], pos[j*nppadded+i], vel[j*nppadded+i]);
    }
  }

  fclose(f);
}
# undef NP