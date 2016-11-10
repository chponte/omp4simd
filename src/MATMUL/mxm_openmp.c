# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <omp.h>
# include <string.h>

// Include appropiate library depending on compiler
#if defined(__INTEL_COMPILER)
#include <malloc.h>
#else
#include <mm_malloc.h>
#endif

# define N 1024
# define ALIGNMENT 64

int main ( int argc, char *argv[] );
void timestamp ( void );
void write_to_file(const int n, const int npadded, double *array, char *outfile);

/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for MXM_OPENMP.

  Licensing:

    This code is distributed under the GNU GPLv3 license. 

  Modified:

    10 November 2016

  Author:

    John Burkardt
    Christian Ponte
*/
{
  const int divider = ALIGNMENT<sizeof(double) ? 1 : ALIGNMENT/sizeof(double);
  const int npadded = (N+divider-1)/divider * divider;
  double *a = (double *) _mm_malloc(N*N*sizeof(double), ALIGNMENT);
  double *trans = (double *) _mm_malloc(N*N*sizeof(double), ALIGNMENT);
  double *b = (double *) malloc(N*N*sizeof(double));
  double *c = (double *) malloc(N*N*sizeof(double));
  double angle;
  double temp;
  int i;
  int j;
  int k;
  int n = N;
  double pi = 3.141592653589793;
  double s;
  int thread_num;
  double wtime;
  char * outfile = NULL;

  timestamp ( );

  printf ( "\n" );
  printf ( "MXM_OPENMP:\n" );
  printf ( "  C/OpenMP version\n" );
  printf ( "  Compute matrix product C = A * B.\n" );

  thread_num = omp_get_max_threads ( );

  printf ( "\n" );
  printf ( "  The number of processors available = %d\n", omp_get_num_procs ( ) );
  printf ( "  The number of threads available    = %d\n", thread_num );

  printf ( "  The matrix order N                 = %d\n", n );
  printf ( "  Padded matrix order                = %d\n", npadded );

  if (argc > 1){
    outfile = (char *) malloc(sizeof(char) * strlen(argv[1]));
    strcpy(outfile, argv[1]);
  }

/*
  Loop 1: Evaluate A.
*/
  s = 1.0 / sqrt ( ( double ) ( n ) );

# pragma omp parallel shared ( wtime, a, b, c, npadded, pi, s ) private ( angle, i, j, k )
{
  # pragma omp for
  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      angle = 2.0 * pi * i * j / ( double ) n;
      a[i*npadded+j] = s * ( sin ( angle ) + cos ( angle ) );
    }
  }
/*
  Loop 2: Copy A into B.
*/
  # pragma omp for
  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      b[i*npadded+j] = a[i*npadded+j];
    }
  }
/*
  Loop 3: Compute C = A * B.
*/

  # pragma omp single
  {
    wtime = omp_get_wtime ( );
  }

  // First transpose B
  # pragma omp for
  for (i=0; i<n; i++){
    for (j=0; j<n; j++){
      trans[i*npadded+j] = b[j*npadded+i];
    }
  }

  # pragma omp for
  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      temp = 0;
      # pragma omp simd aligned(a,trans:ALIGNMENT) reduction(+:temp)
      for ( k = 0; k < n; k++ )
      {
        temp += a[i*npadded+k] * trans[j*npadded+k];
      }
      c[i*npadded+j] = temp;
    }
  }

}
  wtime = omp_get_wtime ( ) - wtime;
  printf ( "  Elapsed seconds = %g\n", wtime );

  if (outfile != NULL){
    printf( "  Writting result to %s\n", outfile);
    write_to_file(n, npadded, c, outfile);
    free(outfile);
  }

/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "MXM_OPENMP:\n" );
  printf ( "  Normal end of execution.\n" );

  printf ( "\n" );
  timestamp ( );

  return 0;
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

void write_to_file(const int n, const int npadded, double *array, char *outfile)

/******************************************************************************/
/*
  Purpose:

    WRITE_TO_FILE writes the array ARRAY into the specified file in OUTFILE.

  Licensing:

    This code is distributed under the GNU GPLv3 license. 

  Modified:

    10 November 2016

  Author:

    Christian Ponte

  Parameters:

    Input, int N, the size of the array.

    Input, int NPADDED, the width of the array padded to a multiple of ALIGNMENT.

    Input, double *array, the input array.

    Input, char *outfile, the output file name.
*/
{
  int i,j;

  FILE *f = fopen(outfile, "w");
  if (f==NULL) return;

  for (i=0; i<n; i++){
    for (j=0; j<n; j++){
      fprintf(f, "%f ", array[i*npadded+j]);
    }
    fprintf(f, "\n");
  }

  fclose(f);
}
# undef N