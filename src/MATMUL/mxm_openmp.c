# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <omp.h>

#define N 500

int main ( void );
void timestamp ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for MXM_OPENMP.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 October 2011

  Author:

    John Burkardt
*/
{
  double *a = (double *) malloc(N*N*sizeof(double));
  double angle;
  double *b = (double *) malloc(N*N*sizeof(double));
  double *c = (double *) malloc(N*N*sizeof(double));
  int i;
  int j;
  int k;
  int n = N;
  double pi = 3.141592653589793;
  double s;
  int thread_num;
  double wtime;

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
/*
  Loop 1: Evaluate A.
*/
  s = 1.0 / sqrt ( ( double ) ( n ) );

  wtime = omp_get_wtime ( );

# pragma omp parallel shared ( a, b, c, n, pi, s ) private ( angle, i, j, k )
{
  # pragma omp for
  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      angle = 2.0 * pi * i * j / ( double ) n;
      a[i*n+j] = s * ( sin ( angle ) + cos ( angle ) );
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
      b[i*n+j] = a[i*n+j];
    }
  }
/*
  Loop 3: Compute C = A * B.
*/
  # pragma omp for
  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      c[i*n+j] = 0.0;
      for ( k = 0; k < n; k++ )
      {
        c[i*n+j] = c[i*n+j] + a[i*n+k] * b[k*n+j];
      }
    }
  }

}
  wtime = omp_get_wtime ( ) - wtime;
  printf ( "  Elapsed seconds = %g\n", wtime );
/*
  Terminate.
*/
  free(a);
  free(b);
  free(c);

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

    This code is distributed under the GNU LGPL license. 

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
