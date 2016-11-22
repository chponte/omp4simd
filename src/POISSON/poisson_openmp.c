# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <omp.h>
# include <string.h>

// Include appropiate library depending on compiler
#if defined(__INTEL_COMPILER)
# include <malloc.h>
#else
# include <mm_malloc.h>
#endif

#define NX 1024
#define NY 1024
#define ALIGNMENT 64

int main ( int argc, char *argv[] );
double r8mat_rms ( int nx, int ny, int npadded, double *a );
void rhs ( int nx, int ny, int npadded, double *f );
void sweep ( int nx, int ny, int npadded, double dx, double dy, double *f,
  int itold, int itnew, double *u, double *unew );
void timestamp ( void );
double u_exact ( double x, double y );
double uxxyy_exact ( double x, double y );
void write_to_file(const int nx, const int ny, int npadded, double *u, char *outfile);

/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for POISSON_OPENMP.

  Discussion:

    POISSON_OPENMP is a program for solving the Poisson problem.

    This program uses OpenMP for parallel execution.

    The Poisson equation

      - DEL^2 U(X,Y) = F(X,Y)

    is solved on the unit square [0,1] x [0,1] using a grid of NX by
    NX evenly spaced points.  The first and last points in each direction
    are boundary points.

    The boundary conditions and F are set so that the exact solution is

      U(x,y) = sin ( pi * x * y )

    so that

      - DEL^2 U(x,y) = pi^2 * ( x^2 + y^2 ) * sin ( pi * x * y )

    The Jacobi iteration is repeatedly applied until convergence is detected.

    For convenience in writing the discretized equations, we assume that NX = NY.

  Licensing:

    This code is distributed under the GNU GPLv3 license. 

  Modified:

    10 November 2016

  Author:

    John Burkardt
    Christian Ponte
*/
{
  const int divisor = ALIGNMENT<sizeof(double) ? 1 : ALIGNMENT/sizeof(double);
  const int npadded = (NY + divisor -1)/divisor * divisor;
  int converged;
  double diff;
  double dx;
  double dy;
  double error;
  double *f = (double *) _mm_malloc(NX*npadded*sizeof(double), ALIGNMENT);
  int i;
  int id;
  int itnew;
  int itold;
  int j;
  int nx = NX;
  int ny = NY;
  double tolerance = 0.000001;
  double *u = (double *) _mm_malloc(NX*npadded*sizeof(double), ALIGNMENT);
  double u_norm;
  double *udiff = (double *) _mm_malloc(NX*npadded*sizeof(double), ALIGNMENT);
  double *uexact = (double *) _mm_malloc(NX*npadded*sizeof(double), ALIGNMENT);
  double *unew = (double *) _mm_malloc(NX*npadded*sizeof(double), ALIGNMENT);
  double unew_norm;
  double wtime;
  double x;
  double y;
  char * outfile = NULL;

  dx = 1.0 / ( double ) ( nx - 1 );
  dy = 1.0 / ( double ) ( ny - 1 );
/*
  Print a message.
*/
  timestamp ( );
  printf ( "\n" );
  printf ( "POISSON_OPENMP:\n" );
  printf ( "  C version\n" );
  printf ( "  A program for solving the Poisson equation.\n" );
  printf ( "\n" );
  printf ( "  Use OpenMP for parallel execution.\n" );
  printf ( "  The number of processors is %d\n", omp_get_num_procs ( ) );
# pragma omp parallel
{
  id = omp_get_thread_num ( );
  if ( id == 0 )
  {
    printf ( "  The maximum number of threads is %d\n", omp_get_num_threads ( ) );
  }
}
  printf ( "\n" );
  printf ( "  -DEL^2 U = F(X,Y)\n" );
  printf ( "\n" );
  printf ( "  on the rectangle 0 <= X <= 1, 0 <= Y <= 1.\n" );
  printf ( "\n" );
  printf ( "  F(X,Y) = pi^2 * ( x^2 + y^2 ) * sin ( pi * x * y )\n" );
  printf ( "\n" );
  printf ( "  The number of interior X grid points is %d\n", nx );
  printf ( "  The number of interior Y grid points is %d\n", ny );
  printf ( "  The X grid spacing is %f\n", dx );
  printf ( "  The Y grid spacing is %f\n", dy );

  if (argc > 1){
    outfile = (char *) malloc(sizeof(char) * strlen(argv[1]));
    strcpy(outfile, argv[1]);
  }

/*
  Set the right hand side array F.
*/
  rhs ( nx, ny,npadded, f );
/*
  Set the initial solution estimate UNEW.
  We are "allowed" to pick up the boundary conditions exactly.
*/
  for ( i = 0; i < nx; i++ )
  {
    for ( j = 0; j < ny; j++ )
    {
      if ( i == 0 || i == nx - 1 || j == 0 || j == ny - 1 )
      {
        unew[i*npadded+j] = f[i*npadded+j];
      }
      else
      {
        unew[i*npadded+j] = 0.0;
      }
    }
  }
  unew_norm = r8mat_rms ( nx, ny, npadded, unew );
/*
  Set up the exact solution UEXACT.
*/
  for ( i = 0; i < nx; i++ )
  {
    x = ( double ) ( i ) / ( double ) ( nx - 1 );
    for ( j = 0; j < ny; j++ )
    {
      y = ( double ) ( j ) / ( double ) ( ny - 1 );
      uexact[i*npadded+j] = u_exact ( x, y );
    }
  }
  u_norm = r8mat_rms ( nx, ny, npadded, uexact );
  printf ( "  RMS of exact solution = %g\n", u_norm );
/*
  Do the iteration.
*/
  converged = 0;

  printf ( "\n" );
  printf ( "  Step    ||Unew||     ||Unew-U||     ||Unew-Exact||\n" );
  printf ( "\n" );

  for ( i = 0; i < nx; i++ )
  {
    for ( j = 0; j < ny; j++ )
    {
      udiff[i*npadded+j] = unew[i*npadded+j] - uexact[i*npadded+j];
    }
  }
  error = r8mat_rms ( nx, ny, npadded, udiff );
  printf ( "  %4d  %14g                  %14g\n", 0, unew_norm, error );

  wtime = omp_get_wtime ( );

  itnew = 0;

  for ( ; ; )
  {
    itold = itnew;
    itnew = itold + 500;
/*
  SWEEP carries out 500 Jacobi steps in parallel before we come
  back to check for convergence.
*/
    sweep ( nx, ny, npadded, dx, dy, f, itold, itnew, u, unew );
/*
  Check for convergence.
*/
    u_norm = unew_norm;
    unew_norm = r8mat_rms ( nx, ny, npadded, unew );

    # pragma omp parallel for shared(udiff, unew, u, nx, ny) private(i, j)
    for ( i = 0; i < nx; i++ )
    {
      # pragma omp simd aligned(udiff:ALIGNMENT, unew:ALIGNMENT, u:ALIGNMENT)
      for ( j = 0; j < ny; j++ )
      {
        udiff[i*npadded+j] = unew[i*npadded+j] - u[i*npadded+j];
      }
    }
    diff = r8mat_rms ( nx, ny, npadded, udiff );

    # pragma omp parallel for shared(udiff, unew, uexact, nx, ny) private(i, j)
    for ( i = 0; i < nx; i++ )
    {
      # pragma omp simd aligned(udiff:ALIGNMENT, unew:ALIGNMENT, uexact:ALIGNMENT)
      for ( j = 0; j < ny; j++ )
      {
        udiff[i*npadded+j] = unew[i*npadded+j] - uexact[i*npadded+j];
      }
    }
    error = r8mat_rms ( nx, ny, npadded, udiff );

    printf ( "  %4d  %14g  %14g  %14g\n", itnew, unew_norm, diff, error );

    if ( diff <= tolerance )
    {
      converged = 1;
      break;
    }

  }

  if ( converged )
  {
    printf ( "  The iteration has converged.\n" );
  }
  else
  {
    printf ( "  The iteration has NOT converged.\n" );
  }

  wtime = omp_get_wtime ( ) - wtime;
  printf ( "\n" );
  printf ( "  Elapsed seconds = %g\n", wtime );
/*
  Terminate.
*/

  if (outfile != NULL){
    printf( "  Writting result to %s\n", outfile);
    write_to_file(nx, ny, npadded, unew, outfile);
    free(outfile);
  }

  _mm_free(f);
  _mm_free(u);
  _mm_free(udiff);
  _mm_free(uexact);
  _mm_free(unew);

  printf ( "\n" );
  printf ( "POISSON_OPENMP:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

double r8mat_rms ( int nx, int ny, int npadded, double *a )

/******************************************************************************/
/*
  Purpose:

    R8MAT_RMS returns the RMS norm of a vector stored as a matrix.

  Licensing:

    This code is distributed under the GNU GPLv3 license.

  Modified:

    10 November 2016

  Author:

    John Burkardt
    Christian Ponte

  Parameters:

    Input, int NX, NY, the number of rows and columns in A.

    Input, int NPADDED, the width of the array padded to a multiple of ALIGNMENT.

    Input, double *A, the vector.

    Output, double R8MAT_RMS, the root mean square of the entries of A.
*/
{
  int i;
  int j;
  double v;

  v = 0.0;

  # pragma omp parallel for shared(a) private(j,i) reduction(+:v)
  for ( i = 0; i < nx; i++ )
  {
    # pragma omp simd aligned(a:ALIGNMENT) reduction(+:v)
    for ( j = 0; j < ny; j++ )
    {
      v = v + a[i*npadded+j] * a[i*npadded+j];
    }
  }
  v = sqrt ( v / ( double ) ( nx * ny )  );

  return v;
}
/******************************************************************************/

void rhs ( int nx, int ny, int npadded, double *f )

/******************************************************************************/
/*
  Purpose:

    RHS initializes the right hand side "vector".

  Discussion:

    It is convenient for us to set up RHS as a 2D array.  However, each
    entry of RHS is really the right hand side of a linear system of the
    form

      A * U = F

    In cases where U(I,J) is a boundary value, then the equation is simply

      U(I,J) = F(i,j)

    and F(I,J) holds the boundary data.

    Otherwise, the equation has the form

      (1/DX^2) * ( U(I+1,J)+U(I-1,J)+U(I,J-1)+U(I,J+1)-4*U(I,J) ) = F(I,J)

    where DX is the spacing and F(I,J) is the value at X(I), Y(J) of

      pi^2 * ( x^2 + y^2 ) * sin ( pi * x * y )

  Licensing:

    This code is distributed under the GNU GPLv3 license.

  Modified:

    10 November 2016

  Author:

    John Burkardt
    Christian Ponte

  Parameters:

    Input, int NX, NY, the X and Y grid dimensions.

    Input, int NPADDED, the width of the array padded to a multiple of ALIGNMENT.

    Output, double *F, the initialized right hand side data.
*/
{
  double fnorm;
  int i;
  int j;
  double x;
  double y;
/*
  The "boundary" entries of F store the boundary values of the solution.
  The "interior" entries of F store the right hand sides of the Poisson equation.
*/
  for ( i = 0; i < nx; i++ )
  {
    x = ( double ) ( i ) / ( double ) ( nx - 1 );
    for ( j = 0; j < ny; j++ )
    {
      y = ( double ) ( j ) / ( double ) ( ny - 1 );
      if ( i == 0 || i == nx - 1 || j == 0 || j == ny - 1 )
      {
        f[i*npadded+j] = u_exact ( x, y );
      }
      else
      {
        f[i*npadded+j] = - uxxyy_exact ( x, y );
      }
    }
  }

  fnorm = r8mat_rms ( nx, ny, npadded, f );

  printf ( "  RMS of F = %g\n", fnorm );

  return;
}
/******************************************************************************/

void sweep ( int nx, int ny, int npadded, double dx, double dy, double *f,
  int itold, int itnew, double *u, double *unew )

/******************************************************************************/
/*
  Purpose:

   SWEEP carries out one step of the Jacobi iteration.

  Discussion:

    Assuming DX = DY, we can approximate

      - ( d/dx d/dx + d/dy d/dy ) U(X,Y)

    by

      ( U(i-1,j) + U(i+1,j) + U(i,j-1) + U(i,j+1) - 4*U(i,j) ) / dx / dy

    The discretization employed below will not be correct in the general
    case where DX and DY are not equal.  It's only a little more complicated
    to allow DX and DY to be different, but we're not going to worry about
    that right now.

  Licensing:

    This code is distributed under the GNU GPLv3 license.

  Modified:

    10 November 2016

  Author:

    John Burkardt
    Christian Ponte

  Parameters:

    Input, int NX, NY, the X and Y grid dimensions.

    Input, int NPADDED, the width of the array padded to a multiple of ALIGNMENT.

    Input, double DX, DY, the spacing between grid points.

    Input, double *F, the right hand side data.

    Input, int ITOLD, the iteration index on input.

    Input, int ITNEW, the desired iteration index
    on output.

    Input, double *U, the solution estimate on
    iteration ITNEW-1.

    Input/output, double *UNEW, on input, the solution
    estimate on iteration ITOLD.  On output, the solution estimate on
    iteration ITNEW.
*/
{
  int i;
  int it;
  int j;

# pragma omp parallel shared ( dx, dy, f, itnew, itold, nx, ny, u, unew ) private ( i, it, j )

  for ( it = itold + 1; it <= itnew; it++ )
  {
/*
  Save the current estimate.
*/
# pragma omp for
    for ( i = 0; i < nx; i++ )
    {
      # pragma omp simd aligned(u:ALIGNMENT, unew:ALIGNMENT)
      for ( j = 0; j < ny; j++ )
      {
        u[i*npadded+j] = unew[i*npadded+j];
      }
    }
/*
  Compute a new estimate.
*/

    # pragma omp for simd aligned(unew:ALIGNMENT, f:ALIGNMENT)
    for (i=0; i<ny; i++){
      unew[i] = f[i];
    }

    # pragma omp for simd aligned(unew:ALIGNMENT, f:ALIGNMENT)
    for (i=0; i<ny; i++){
      unew[(nx-1)*npadded + i] = f[(nx-1)*npadded + i];
    }

    # pragma omp for
    for ( i = 1; i < nx-1; i++ )
    {
      unew[i*npadded] = f[i*npadded];
      # pragma omp simd aligned(unew:ALIGNMENT, u:ALIGNMENT, f:ALIGNMENT)
      for ( j = 1; j < ny-1; j++ )
      {
        unew[i*npadded+j] = 0.25 * ( u[(i-1)*npadded+j] + u[i*npadded+j+1] + u[i*npadded+j-1] + u[(i+1)*npadded+j] + f[i*npadded+j] * dx * dy );
      }
      unew[i*npadded + ny-1] = f[i*npadded + ny-1];
    }
  }
  return;
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
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  printf ( "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}
/******************************************************************************/

double u_exact ( double x, double y )

/******************************************************************************/
/*
  Purpose:

    U_EXACT evaluates the exact solution.

  Licensing:

    This code is distributed under the GNU GPLv3 license.

  Modified:

    25 October 2011

  Author:

    John Burkardt

  Parameters:

    Input, double X, Y, the coordinates of a point.

    Output, double U_EXACT, the value of the exact solution
    at (X,Y).
*/
{
  double pi = 3.141592653589793;
  double value;

  value = sin ( pi * x * y );

  return value;
}
/******************************************************************************/

double uxxyy_exact ( double x, double y )

/******************************************************************************/
/*
  Purpose:

    UXXYY_EXACT evaluates ( d/dx d/dx + d/dy d/dy ) of the exact solution.

  Licensing:

    This code is distributed under the GNU GPLv3 license.

  Modified:

    25 October 2011

  Author:

    John Burkardt

  Parameters:

    Input, double X, Y, the coordinates of a point.

    Output, double UXXYY_EXACT, the value of
    ( d/dx d/dx + d/dy d/dy ) of the exact solution at (X,Y).
*/
{
  double pi = 3.141592653589793;
  double value;

  value = - pi * pi * ( x * x + y * y ) * sin ( pi * x * y );

  return value;
}
/******************************************************************************/

void write_to_file(const int nx, const int ny, int npadded, double *u, char *outfile)

/******************************************************************************/
/*
  Purpose:

    WRITE_TO_FILE writes the array U into the specified file in OUTFILE.

  Licensing:

    This code is distributed under the GNU GPLv3 license. 

  Modified:

    10 November 2016

  Author:

    Christian Ponte

  Parameters:

    Input, int nx, the horizontal size of the array.

    Input, int ny, the vertical size of the array.

    Input, int npadded, the padded horizontal dimension of the array.

    Input, double *u, the input array.

    Input, char *outfile, the output file name.
*/
{
  int i,j;

  FILE *f = fopen(outfile, "w");
  if (f==NULL) return;

  for (i=0; i<NX; i++){
    for (j=0; j<NY; j++){
      fprintf(f, "%f ", u[i*npadded+j]);
    }
    fprintf(f, "\n");
  }

  fclose(f);
}


# undef NX
# undef NY