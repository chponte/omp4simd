# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <omp.h>
# include <string.h>

# define N 1024

int main ( int argc, char *argv[] );
void timestamp ( void );
void write_to_file(const int nwidth, double *array, char *outfile);

/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for MXM_OPENMP.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    4 October 2016

  Author:

    John Burkardt
    Christian Ponte
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

  if (argc > 1){
    outfile = (char *) malloc(sizeof(char) * strlen(argv[1]));
    strcpy(outfile, argv[1]);
  }

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

  if (outfile != NULL){
    printf( "  Writting result to %s\n", outfile);
    write_to_file(n, c, outfile);
    free(outfile);
  }

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
/******************************************************************************/

void write_to_file(const int nwidth, double *array, char *outfile)

/******************************************************************************/
/*
  Purpose:

    WRITE_TO_FILE writes the array ARRAY into the specified file in OUTFILE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    4 October 2016

  Author:

    Christian Ponte

  Parameters:

    Input, int nwidth, the size of the array.

    Input, double *array, the input array.

    Input, char *outfile, the output file name.
*/
{
  int i,j;

  FILE *f = fopen(outfile, "w");
  if (f==NULL) return;

  for (i=0; i<N; i++){
    for (j=0; j<N; j++){
      fprintf(f, "%f ", array[i*nwidth+j]);
    }
    fprintf(f, "\n");
  }

  fclose(f);
}
# undef N