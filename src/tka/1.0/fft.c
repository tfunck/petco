#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include "nrutil.h"
#include "fft.h"

#ifndef TRUE
#  define TRUE 1
#  define FALSE 0
#endif

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr


/* Globals */
static double PI = 3.1415926535897932384626433832795;

/**
 * Returns 1 if n is a power of 2, 0 otherwise
 */
int isPowerOf2(unsigned long n)
{
  int val=1;
  float temp_n = (float)n;

  if ((unsigned long)(temp_n) != n)
  {
    fprintf(stderr,"isPowerOf2(%d) failed, possibly out of range\n",n);
    exit(1);
  }

  while (temp_n != 1.0)
  {
    printf("temp_n is %f\n",temp_n);
    if (temp_n != (int)temp_n || temp_n == 0.0)
    {
      /* If we have a decimal or zero it is not a power of two */
      val=0;
      break;
    }
    else
    {
      /* Divide by two for the next iteration */
      temp_n=temp_n/2.0;
    }
  }
  return val;
}

/**
 * Returns the next power of 2, such that it is => n
 */
int nextPowerOf2(unsigned long n)
{
  unsigned int temp_n = 2; 

  while (temp_n < n)
  {
    temp_n *= 2;
  }
  return temp_n;
}

/**
 * Replaces data by its n-dimensional discrete Fourier transform, if isign is
 * input as 1.  nn[1..ndim] is an integer array containing the lengths of each
 * dimension (number of complex values), which MUST all be powers of 2.  data
 * is a real array of length twice the product of these lengths, in which the
 * data are stored as in a multidimensional complex array: real and imaginary
 * parts of each element are in consecutive locations, and the right-most index
 * of the array increases most rapidly as one proceeds along data.  For a
 * two-dimensional array, this is equivalent to storing the array by rows.
 * If isign is input as -1, data is replaced by its inverse transform times the
 * product of the lengths of all dimensions.
 *
 * Note that after a FFT the output spectrum is returned as a logically packaged
 * complex 3 dimensional array, SPEC, of dimensions [1..nn1][1..nn2][1..nn3/2+1].
 * The third of the three dimensions returns only the positive half of the frequency
 * spectrum.
 *
 * The physical storage is somewhat different to the logical storage, with the bulk
 * of the data replacing the original data and the remainder in the speq array:
 *
 *   Re(SPEC[i1][i2][i3])=data[i1][i2][2*i3-1]
 *   Im(SPEC[i1][i2][i3])=data[i1][i2][2*i3]
 *
 *   Re(SPEC[i1][i2][nn3/2+1])=speq[i1][2*i2-1]
 *   Im(SPEC[i1][i2][nn3/2+1])=speq[i1][2*i2]
 *
 * See Numerical Recipies in C book for more information if required.
 *
 * Thanks to authors of Numerical Recipies in C for the original implementation.
 * This code is unmodified from the original implementation.
 */
void fourn(float data[], unsigned long nn[], int ndim, int isign)
{

  int idim;
  unsigned long i1,i2,i3,i2rev,i3rev,ip1,ip2,ip3,ifp1,ifp2;
  unsigned long ibit,k1,k2,n,nprev,nrem,ntot;
  float tempi,tempr;
  double theta,wi,wpi,wpr,wr,wtemp; /* Double precision for trigonometric recurrences */

  /* Compute total number of complex values */
  for (ntot=1,idim=1;idim<=ndim;idim++) ntot *= nn[idim];
  nprev=1;

  /* Main loop over the dimensions. */
  for (idim=ndim;idim>=1;idim--)
  {
    n=nn[idim];
    nrem=ntot/(n*nprev);
    ip1=nprev << 1;
    ip2=ip1*n;
    ip3=ip2*nrem;
    i2rev=1;
    for (i2=1;i2<=ip2;i2+=ip1)
    {
      /* This is the bit-reversal section of the routine. */
      if (i2 < i2rev)
      {
        for (i1=i2;i1<=i2+ip1-2;i1+=2)
        {
          for (i3=i1;i3<=ip3;i3+=ip2)
          {
            i3rev=i2rev+i3-i2;
            SWAP(data[i3],data[i3rev]);
            SWAP(data[i3+1],data[i3rev+1]);
          }
        }
      }
      ibit=ip2 >> 1;
      while (ibit >= ip1 && i2rev > ibit)
      {
        i2rev -= ibit;
        ibit >>= 1;
      }
      i2rev += ibit;
    }
    ifp1=ip1;

    /* Here begins the Danielson-Lanczos section of the routine */
    while (ifp1 < ip2)
    {
      ifp2=ifp1 << 1;
      theta=isign*6.28318530717959/(ifp2/ip1);

      /* Initialize for the trig. recurrence */
      wtemp=sin(0.5*theta);
      wpr = -2.0*wtemp*wtemp;
      wpi=sin(theta);
      wr=1.0;
      wi=0.0;
      for (i3=1;i3<=ifp1;i3+=ip1)
      {
        for (i1=i3;i1<=i3+ip1-2;i1+=2)
        {
          for (i2=i1;i2<=ip3;i2+=ifp2)
          {
            /* Danielson-Lanczos formula: */
            k1=i2;
            k2=k1+ifp1;
            tempr=(float)wr*data[k2]-(float)wi*data[k2+1];
            tempi=(float)wr*data[k2+1]+(float)wi*data[k2];
            data[k2]=data[k1]-tempr;
            data[k2+1]=data[k1+1]-tempi;
            data[k1] += tempr;
            data[k1+1] += tempi;
          }
        }
        /* Trigonometric recurrence. */
        wr=(wtemp=wr)*wpr-wi*wpi+wr;
        wi=wi*wpr+wtemp*wpi+wi;
      }
      ifp1=ifp2;
    }
    nprev *= n;
  }
}


/**
 * Summary: performs a FFT in 2d or 3d. For a 2d transform, the first dimension (nn1)
 * must be of length 1 - the mechanics of the routine are identical for both 2 and 3d.
 *
 * Given a three-dimensional real array iFFTdata[1..nn1][1..nn2][1..nn3] (where nn1=1 for the
 * case of a logically two-dimensional array), this routine returns (for isign=1) the complex
 * fast Fourier transform as two complex arrays: On output, iFFTdata contains the zero and
 * positive frequency values of the third frequency component, while speq[1..nn1][1..2*nn2]
 * contains the Nyquist critical frequency values of the third frequency component. First
 * (and second) frequency components are stored for zero, positive, and negative frequencies,
 * in standard wraparound order.
 *
 * Complex values are arranged as a pair of sequential floats (real, imaginary) within the
 * iFFTdata array. As the array is indexed from 1, odd indexes are for the real components whilst
 * the even indexes are for the imaginary components.
 *
 * For isign=-1,the inverse transform (times nn1*nn2*nn3/2 as a constant multiplicative
 * factor) is performed, with output data (viewed as a real array) deriving from input data
 * (viewed as complex) and speq. For inverse transforms on data not generated first by a
 * forward transform, make sure the complex input data array satisfies property (12.5.2)
 * (see Chapter 12.5 in Numerical Receipes in C for more information).
 *
 * The dimensions nn1, nn2, nn3 must always be integer powers of 2. Conceptually:
 *    nn1 represents the depth
 *    nn2 represents the rows
 *    nn3 represents the columns
 *
 * The number of rows need not be equal to the number of columns. The space for the imaginary
 * component of the data must be included in the number of columns specified here, so for
 * example, for a 2d matrix of 256x256, you would use nn1=1, nn2=256, nn3=512.
 *
 * Thanks to authors of Numerical Recipies in C for the original implementation.
 * This code is unmodified from the original implementation.
 *
 **/
void rlft3(float ***FFTdata, float **speq, unsigned long nn1, unsigned long nn2,
                  unsigned long nn3, int isign)
{
  void nrerror(char error_text[]);
  unsigned long i1,i2,i3,j1,j2,j3,nn[4],ii3;
  double theta,wi,wpi,wpr,wr,wtemp;
  float c1,c2,h1r,h1i,h2r,h2i;

  if (1+&FFTdata[nn1][nn2][nn3]-&FFTdata[1][1][1] != nn1*nn2*nn3)
  {
    nrerror("rlft3: problem with dimensions or contiguity of data array\n");
  }

  c1=0.5;
  c2 = -0.5*isign;
  theta=isign*(6.28318530717959/nn3);
  wtemp=sin(0.5*theta);
  wpr = -2.0*wtemp*wtemp;
  wpi=sin(theta);
  nn[1]=nn1;
  nn[2]=nn2;
  nn[3]=nn3 >> 1;

  if (isign == 1)
  {
    /* Case of forward transform in 2 dimensions */
    fourn(&FFTdata[1][1][1]-1,nn,3,isign);

    /* Computationally expensive part */
    for (i1=1;i1<=nn1;i1++)
    {
      for (i2=1,j2=0;i2<=nn2;i2++)
      {
        /* Extend FFTdata periodically into speq */
        speq[i1][++j2]=FFTdata[i1][i2][1];
        speq[i1][++j2]=FFTdata[i1][i2][2];
      }
    }
  }
  for (i1=1;i1<=nn1;i1++)
  {
    j1=(i1 != 1 ? nn1-i1+2 : 1);
    wr=1.0;

    /* Initialize trigonometric recurrence. */
    wi=0.0;

    /* Zero frequency is its own reflection , otherwise local corresponding negative *
     * frequency in wrap around order.                                               */
    for (ii3=1,i3=1;i3<=(nn3>>2)+1;i3++,ii3+=2)
    {
      for (i2=1;i2<=nn2;i2++)
      {
        if (i3 == 1)
        {
          j2=(i2 != 1 ? ((nn2-i2)<<1)+3 : 1);
          h1r=c1*(FFTdata[i1][i2][1]+speq[j1][j2]);
          h1i=c1*(FFTdata[i1][i2][2]-speq[j1][j2+1]);
          h2i=c2*(FFTdata[i1][i2][1]-speq[j1][j2]);
          h2r= -c2*(FFTdata[i1][i2][2]+speq[j1][j2+1]);
          FFTdata[i1][i2][1]=h1r+h2r;
          FFTdata[i1][i2][2]=h1i+h2i;
          speq[j1][j2]=h1r-h2r;
          speq[j1][j2+1]=h2i-h1i;
        }
        else
        {
          j2=(i2 != 1 ? nn2-i2+2 : 1);
          j3=nn3+3-(i3<<1);
          h1r=c1*(FFTdata[i1][i2][ii3]+FFTdata[j1][j2][j3]);
          h1i=c1*(FFTdata[i1][i2][ii3+1]-FFTdata[j1][j2][j3+1]);
          h2i=c2*(FFTdata[i1][i2][ii3]-FFTdata[j1][j2][j3]);
          h2r= -c2*(FFTdata[i1][i2][ii3+1]+FFTdata[j1][j2][j3+1]);
          FFTdata[i1][i2][ii3]=h1r+wr*h2r-wi*h2i;
          FFTdata[i1][i2][ii3+1]=h1i+wr*h2i+wi*h2r;
          FFTdata[j1][j2][j3]=h1r-wr*h2r+wi*h2i;
          FFTdata[j1][j2][j3+1]= -h1i+wr*h2i+wi*h2r;
        }
      }
      /* Do the recurrence */
      wr=(wtemp=wr)*wpr-wi*wpi+wr;
      wi=wi*wpr+wtemp*wpi+wi;
    }
  }
  if (isign == -1)
  {
    /* Case of inverse transform in 2 dimensions */
    fourn(&FFTdata[1][1][1]-1,nn,3,isign);
  }
}


/**
 * After a forward FFT and its inverse, the original matrix is restored, however, with a
 * mulitplication factor of nn1 * nn2 * nn3. This routine simply divides each value by
 * this factor to restore the real values.
 *
 * Note there are no imaginary values at this point as the matrix is back into its original
 * form of real values.
 *
 **/
void normaliseInverse(float ***FFTdata, unsigned long nn1, unsigned long nn2,
                             unsigned long nn3)
{
  int i, j, k;
  float factor = (nn1 * nn2 * nn3) / 2.0;
  for (i=1; i<=nn1; i++)
    for (j=1; j<=nn2; j++)
      for (k=1; k<=nn3; k++)
        FFTdata[i][j][k]=FFTdata[i][j][k]/factor;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : loadImageData3d
@INPUT      : imagedata - 3d array of shorts of dims [i_dim][j_dim][k_dim]
              i_dim - the number of slices of image data
              j_dim - the number of rows of image data
              k_dim - the number of cols of image data
              i_pad - the number of slices of image data after padding
              j_pad - the number of rows of image data after padding
              k_pad - the number of cols of image data after padding
@OUTPUT     : iFFTdata - the manipulated image data
@RETURNS    : (nothing)
@DESCRIPTION: Routine to load the 3 dimensional image data into the data
              structure required for compatibility with the FFT algorithm
              This means the data is converted into floats and any dimensions
              which are not powers of 2 are zero padded to make them so
@METHOD     :
@GLOBALS    :
@CALLS      :
@CREATED    : March 2008, Paul Lyon
@MODIFIED   :
---------------------------------------------------------------------------- */
void loadImageData3d(float ***imagedata,
                     long i_dim, long j_dim, long k_dim,
                     long i_pad, long j_pad, long k_pad,
                     float ***iFFTdata)
{
  int i, j, k;
  unsigned int i_start_offset, j_start_offset, k_start_offset;
  unsigned int i_stop, j_stop, k_stop;

  /* Calculate the start offsets and stops for non-padded image data */
  i_start_offset = (int)((float)(i_pad-i_dim)/2.0);
  j_start_offset = (int)((float)(j_pad-j_dim)/2.0);
  k_start_offset = (int)((float)(k_pad-k_dim)/2.0);
  i_stop = i_start_offset+i_dim;
  j_stop = j_start_offset+j_dim;
  k_stop = k_start_offset+k_dim;



  /* Fill in the padded image data */
  for (i=1; i<=i_pad; i++)
  {
    for (j=1; j<=j_pad; j++)
    {
      for (k=1; k<=k_pad; k++)
      {
        /* Fill the iFFTdata array with the real or padded image values */
        /* Note we are not using imaginary values at this time          */
        if ( i > i_start_offset && i <= i_stop && j > j_start_offset && j <= j_stop && k > k_start_offset && k <= k_stop )
        {
		    //fprintf(stderr, "%f\t%f\t%f\n", (double) i-i_start_offset-1,(double) j-j_start_offset-1,(double) k-k_start_offset-1 );
          iFFTdata[i][j][k] = (float)imagedata[i-i_start_offset-1][j-j_start_offset-1][k-k_start_offset-1];
        }
        else
        {
          iFFTdata[i][j][k] = (float) 0.0;
        }
      }
    }
  }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : loadFilterData3d
@INPUT      : filterdata - 3d array of floats of dims [i_dim][j_dim][k_dim]
              i_dim - the number of slices of filter data
              j_dim - the number of rows of filter data
              k_dim - the number of cols of filter data
              i_pad - the number of slices of filter data after padding
              j_pad - the number of rows of filter data after padding
              k_pad - the number of cols of filter data after padding
@OUTPUT     : fFFTdata - the manipulated filter data
@RETURNS    : (nothing)
@DESCRIPTION: Routine to load the 3 dimensional filter data into the data
              structure required for compatibility with the FFT algorithm
              This means the data is converted into floats and any dimensions
              which are not powers of 2 are zero padded to make them so
@METHOD     :
@GLOBALS    :
@CALLS      :
@CREATED    : March 2008, Paul Lyon
@MODIFIED   :
---------------------------------------------------------------------------- */
void loadFilterData3d(float ***filterdata,
                      long i_dim, long j_dim, long k_dim,
                      long i_pad, long j_pad, long k_pad,
                      float ***fFFTdata)
{
  int i, j, k;
  unsigned int i_start_offset, j_start_offset, k_start_offset;
  unsigned int i_stop, j_stop, k_stop;

  /* Calculate the start offsets and stops for non-padded gaussian data */
  i_start_offset = (int)((float)(i_pad-i_dim)/2.0);
  j_start_offset = (int)((float)(j_pad-j_dim)/2.0);
  k_start_offset = (int)((float)(k_pad-k_dim)/2.0);
  i_stop = i_start_offset+i_dim;
  j_stop = j_start_offset+j_dim;
  k_stop = k_start_offset+k_dim;

  /* Fill in the padded filter data */
  for (i=1; i<=i_pad; i++)
  {
    for (j=1; j<=j_pad; j++)
    {
      for (k=1; k<=k_pad; k++)
      {
        /* Fill the fFFTdata array with the real or padded guass values */
        /* Note we are not using imaginary values at this time          */
        if ( i > i_start_offset && i <= i_stop &&
             j > j_start_offset && j <= j_stop &&
             k > k_start_offset && k <= k_stop )
        {
          fFFTdata[i][j][k] = (float) filterdata[i-i_start_offset-1]
                                               [j-j_start_offset-1]
                                               [k-k_start_offset-1];
        }
        else
        {
          fFFTdata[i][j][k] = (float) 0.0;
        }
      }
    }
  }
}

/**
 * Displays the 3/2D matrix of unsigned shorts
 */
void displayImage(float ***iFFTdata, int nrows, int ncols)
{
  int i,j;

  /* List format for big images arrays *
  for (i=1; i<=nrows; i++)
  {
    for (j=1; j<=ncols; j++)
    {
      printf("iFFTdata(%d,%d)=%f\n",i,j,iFFTdata[1][i][j]/(nrows*ncols/2));
    }
    printf("\n");
  }
  */

  /* Matrix format for small image arrays *
  for (i=1; i<=nrows; i++)
  {
    for (j=1; j<=ncols; j++)
    {
      printf("%3.3f ",iFFTdata[1][i][j]/(nrows*ncols/2));
    }
    printf("\n");
  }
  */

  /* Matrix format for small image arrays */
  for (i=1; i<=nrows; i++)
  {
    for (j=1; j<=ncols; j++)
    {
      printf("%2.3f ",iFFTdata[1][i][j]);
    }
    printf("\n");
  }
}

/**
 * Invokes the FFT method in 3 dimensional mode on our 3D image data.
 *
 * Filter data is a 3d array of appropriate dimensions, if NULL no
 * filtering will be performed, just a forward and a inverse FFT.
 *
 */
void applyFFT(float ***imagedata, float*** filterdata, int nslices, int nrows, int ncols)
{
  float ***iFFTdata;             /* Holds the padded image FFT data    */
  float ***fFFTdata;             /* Holds the padded gaussian FFT data */
  float **ispeq;                 /* Holds the image speq data          */
  float **fspeq;                 /* Holds the gaussian speq data       */
  unsigned int pad_nslices,      /* Variable for the padded dimensions */
               pad_nrows,
               pad_ncols;
//printf("1\n");
  /* Obtain the new padded dims as all dimensions must be a power of 2 */
  pad_nslices = nextPowerOf2(nslices);
  pad_nrows   = nextPowerOf2(nrows);
  pad_ncols   = nextPowerOf2(ncols);


//printf("2\n");
  /* f3tensor() allocates a 3rd rank (3d) tensor (array) of size:        *
   *                                                                     *
   *      pad_nslices x pad_nrows x pad_ncols                            *
   *                                                                     *
   * and with indexing starting from 1.                                  *
   *                                                                     *
   * Note that as the ncols definition does not require any factor of 2  *
   * in order to store the imaginary data. This is because the output    *
   * array which is generated only includes from [1..nn3/2] in the third *
   * dimension.                                                          */
  iFFTdata=f3tensor(1,pad_nslices,1,pad_nrows,1,pad_ncols);
  fFFTdata=f3tensor(1,pad_nslices,1,pad_nrows,1,pad_ncols);
//printf("Imagedata test: %f\n", imagedata[5][10][10] );
//printf("3\n");
  /* Load the image data into the iFFTdata variable                      */
  loadImageData3d(imagedata,nslices,nrows,ncols, pad_nslices,pad_nrows,pad_ncols,iFFTdata);

//printf("4\n");
  /* Load the filter data into the fFFTdata variable                     */
  if (filterdata != NULL)
    loadFilterData3d(filterdata,nslices,nrows,ncols,
                    pad_nslices,pad_nrows,pad_ncols,fFFTdata);
//printf("5\n");
  /* Creates a 2 dimensional matrix of size pad_nslices x 2*pad_nrows    *
   * This the matrix used to hold the spectrum of frequencies for the    *
   * indexes SPEC[1..nn1][1..nn2][nn3/2+1] of the image data.            *
   * See the fourn comments for a more detailed description on the data  *
   * structures used.                                                    */
  ispeq=matrix(1,pad_nslices,1,2*pad_nrows);
//printf("6\n");
  /* Creates a 2 dimensional matrix of size pad_nslices x 2*pad_nrows    *
   * This the matrix used to hold the spectrum of frequencies for the    *
   * indexes SPEC[1..nn1][1..nn2][nn3/2+1] of the Filter data.           *
   * See the fourn comments for a more detailed description on the data  *
   * structures used.                                                    */
  fspeq=matrix(1,pad_nslices,1,2*pad_nrows);
//printf("7\n");
  /* Perform a forward FFT to get the image data in the frequency domain */
  rlft3(iFFTdata,ispeq,pad_nslices,pad_nrows,pad_ncols,1);
//printf("8\n");
  if (filterdata != NULL){
    /* Perform a forward FFT to get the filter data in the freq. domain  */
    rlft3(fFFTdata,fspeq,pad_nslices,pad_nrows,pad_ncols,1);
 
    /* Perform matrix convultion in order to filter the image            */
    applyFilter3d(iFFTdata,fFFTdata,ispeq,fspeq,
                  pad_nslices,pad_nrows,pad_ncols);
  }
//printf("9\n");
  /* Perform an inverse FFT to get the data back into the spacial domain */
  rlft3(iFFTdata,ispeq,pad_nslices,pad_nrows,pad_ncols,-1);

  /* Normalise the inverse to eliminate any scaling (intensity change)   */
  normaliseInverse(iFFTdata,pad_nslices,pad_nrows,pad_ncols);

  /* Unload the filtered image from iFFTdata back into imagedata array   */
  unloadImageData3d(imagedata,nslices,nrows,ncols,
                    pad_nslices,pad_nrows,pad_ncols,iFFTdata);

  /* Clean up after ourselves */
  free_matrix(ispeq,1,pad_nslices,1,2*pad_nrows);
  free_matrix(fspeq,1,pad_nslices,1,2*pad_nrows);
  free_f3tensor(iFFTdata,1,pad_nslices,1,pad_nrows,1,pad_ncols);
  if (filterdata != NULL)
    free_f3tensor(fFFTdata,1,pad_nslices,1,pad_nrows,1,pad_ncols);
}

/**
 * Applies the transformed filter to the 3d image frequency data
 * as a complete 3d matrix multiplication
 */
void applyFilter3d(float*** iFFTdata, float*** fFFTdata,
                   float** ispeq,  float** fspeq,
                   int nn1, int nn2, int nn3)
{
  int j;
  float fac, r, i;
  float *sp1, *sp2;

  /* Factor used to normalise the filter, which has a muliplicative       */
  /* factor of nn1 * nn2 * nn3/2  (in the 2 dimensional case nn1 = 1)     */
  fac=2.0/(nn1*nn2*nn3);

  /* Assign temporary pointers for dealing with FFT arrays */
  sp1 = &iFFTdata[1][1][1];
  sp2 = &fFFTdata[1][1][1];
  for (j=1; j<= nn1*nn2*nn3/2; j++)
  {
    r = sp1[0]*sp2[0] - sp1[1]*sp2[1];
    i = sp1[0]*sp2[1] + sp1[1]*sp2[0];

    /* We don't need to multiply by fac as per original code */
    /* as we used the normaliseInverse() method to do this   */
    /* Best to keep separate so the NO_FILTER can be tested  */
    /* Have tested both ways and give very similar values    */
    /*
    sp1[0] = fac*r;
    sp1[1] = fac*i;
    */

    sp1[0] = r;
    sp1[1] = i;
    sp1 += 2;
    sp2 += 2;
  }

  /* Assign temporary pointers for dealing with speq arrays */
  sp1 = &ispeq[1][1];
  sp2 = &fspeq[1][1];
  for (j=1; j<=nn1*nn2; j++)
  {
    r = sp1[0]*sp2[0] - sp1[1]*sp2[1];
    i = sp1[0]*sp2[1] + sp1[1]*sp2[0];

    /* We don't need to multiply by fac as per original code */
    /* as we used the normaliseInverse() method to do this   */
    /* Best to keep separate so the NO_FILTER can be tested  */
    /* Have tested both ways and give very similar values    */
    /*
    sp1[0]=fac*r;
    sp1[1]=fac*i;
    */

    sp1[0]=r;
    sp1 += 2;
    sp2 += 2; 
  }
}

/**
 * Given the dimensions, allocates memory for e.g. gaussian kernel, indexed
 * from zero. Note that loadFilterData3d is required to copy the array
 * into a new memory space indexed from 1 as required by the FFT algorithm.
 */
float*** allocateKernelMatrix(unsigned int n1, unsigned int n2, unsigned int n3)
{
  int i, j;
  float*** matrix;
  /* Dynamically allocate the memory for the Gaussian matrix */
  matrix = malloc(n1 * sizeof(float**));

  if (matrix == NULL)
  {
    fprintf(stderr,"Unable to allocate enough memory for Gaussian PSF!\n");
    exit(1);
  }
  for (i=0; i<n1; i++)
  {
    matrix[i] = malloc(n2 * sizeof(float*));
    if (matrix[i] == NULL)
    {
      fprintf(stderr,"Unable to allocate enough memory for Gaussian PSF!\n");
      exit(1);
    }
    for (j=0; j<n2; j++)
    {
      matrix[i][j] = malloc(n3 * sizeof(float));
      if (matrix[i][j] == NULL)
      {
        fprintf(stderr,"Unable to allocate enough memory for Gaussian PSF!\n");
        exit(1);
      }
    }
  }


  return matrix;
}

/**
 * Given the dimensions generates a 3d gaussian kernel indexed from zero.
 * Note that loadFilterData3d is required to copy the array
 * into a new memory space indexed from 1 as required by the FFT algorithm.
 *
 * step[] is an array of doubles with the 3 elements representing the
 * step size in the x, y and z dimensions respectively.
 */
float*** genGaussianKernel(float* step,
                           unsigned int n1, unsigned int n2, unsigned int n3)
{
  int i, j, k;
  float*** matrix = allocateKernelMatrix(n1,n2,n3);
  /* Generate normalised Gaussian matrix with stdev=GAUSS_STD_DEV with */
  /* the default centroid coordinates (note FWHM is overridden if the  */
  /* stdev is zero (0.0))                                              */
  psfGaussian(matrix, n1, n2, n3, step[2], step[1], step[0],
              GAUSS_FWHM_MM, GAUSS_FWHM_MM, GAUSS_FWHM_MM,
              GAUSS_STD_DEV, GAUSS_STD_DEV, GAUSS_STD_DEV,
              0.0, 0.0, 0.0, TRUE);

  /* FIXME - display the guts of the matrix */
  /*
  for (i=0; i<n1; i++)
  {
    for (j=0; j<n2; j++)
    {
      for (k=0; k<n3; k++)
      {
        if (matrix[i][j][k] > 0.0001) printf("%f %d %d %d\n",matrix[i][j][k],i,j,k);
      }
    }
  }
  */
  return matrix;
}

/**
 * Given the dimensions generates a 3d identity kernel indexed from zero.
 * Note that loadFilterData3d is required to copy the array
 * into a new memory space indexed from 1 as required by the FFT algorithm.
 */
float*** genIdentityKernel(unsigned int n1, unsigned int n2, unsigned int n3)
{
  int i, j, k;
  float*** matrix = allocateKernelMatrix(n1,n2,n3);

  /* Generate normalised Identity matrix */
  for (i=0; i<n1; i++)
  {
    for (j=0; j<n2; j++)
    {
      for (k=0; k<n3; k++)
      {
        matrix[i][j][k] = 1.0;
      }
    }
  }
  return matrix;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : unloadImageData3d
@INPUT      : iFFTdata - 3 dimensional array of image FFT data
              i_dim - the number of slices of image FFT data
              j_dim - the number of rows of image FFT data
              k_dim - the number of cols of image FFT data
              i_pad - the number of slices of image FFT data after padding
              j_pad - the number of rows of image FFT data after padding
              k_pad - the number of cols of image FFT data after padding
@OUTPUT     : imagedata - 3d array of image data of dims [i_dim][j_dim][k_dim]
@RETURNS    : (nothing)
@DESCRIPTION: Routine to unload the 3 dimensional image FFT data back into the
              3 dimensional array as required by loadImageData3d()
              This is required for compatibility with the FFT algorithm
@METHOD     :
@GLOBALS    :
@CALLS      :
@CREATED    : March 2008, Paul Lyon
@MODIFIED   :
---------------------------------------------------------------------------- */
void unloadImageData3d(float ***imagedata,
                  long i_dim, long j_dim, long k_dim,
                  long i_pad, long j_pad, long k_pad,
                  float ***iFFTdata)
{
  int i, j, k;
  unsigned int i_start_offset, j_start_offset, k_start_offset;
  unsigned int i_stop, j_stop, k_stop;

  /* Calculate the start offsets and stops for non-padded image data */
  i_start_offset = (int)((float)(i_pad-i_dim)/2.0);
  j_start_offset = (int)((float)(j_pad-j_dim)/2.0);
  k_start_offset = (int)((float)(k_pad-k_dim)/2.0);
  i_stop = i_start_offset+i_dim;
  j_stop = j_start_offset+j_dim;
  k_stop = k_start_offset+k_dim;

  for (i=1; i<=i_pad; i++)
  {
    for (j=1; j<=j_pad; j++)
    {
      for (k=1; k<=k_pad; k++)
      {
        /* Fill the imagedata array with the real values ignoring padding */
        /* Note we are not using imaginary values at this time            */
        if ( i > i_start_offset && i <= i_stop &&
             j > j_start_offset && j <= j_stop &&
             k > k_start_offset && k <= k_stop )
        { 
          imagedata[i-i_start_offset-1]
                   [j-j_start_offset-1]
                   [k-k_start_offset-1] = (signed short) iFFTdata[i][j][k];
        }
      }
    }
  }
}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : freeKernelMatrix
@INPUT      : gaussMatrix - dynamically allocated 3d array to be freed up
              n1 - the x dimension
              n2 - the y dimension
@OUTPUT     :
@RETURNS    : (nothing)
@DESCRIPTION: Frees the gaussian matrix 3d array
@METHOD     :
@GLOBALS    :
@CALLS      :
@CREATED    : March 2008, Paul Lyon
@MODIFIED   :
---------------------------------------------------------------------------- */
void freeKernelMatrix(float*** matrix, unsigned int n1, unsigned int n2)
{
  int i, j;

  for (i=0; i<n1; i++)
  {
    for (j=0; j<n2; j++)
    {
      free(matrix[i][j]);
    }
    free(matrix[i]);
  }
  free(matrix);
}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : psfGaussian
@INPUT      : psf    - pointer to an empty 3d matrix to hold the gaussian data
                       (memory must already be allocated)
              xdim   - number pixels for x dimension
              ydim   - number pixels for y dimension
              zdim   - number pixels for z dimension
              xstep  - the size of the step in the x dimension (mm)
              ystep  - the size of the step in the y dimension (mm)
              zstep  - the size of the step in the z dimension (mm)
              xfwhm  - the Full-Width Half-Max (pixels) for x dimension (mm)
              yfwhm  - the Full-Width Half-Max (pixels) for y dimension (mm)
              zfwhm  - the Full-Width Half-Max (pixels) for z dimension (mm)
              xstdev - the standard dev. for x dim (overrides xfwhm if not 0)
              ystdev - the standard dev. for y dim (overrides yfwhm if not 0)
              zstdev - the standard dev. for z dim (overrides zfwhm if not 0)
              xcntrd - x coordinate of the centroid (pixel no of PSF maximum)
                       0.5 is centre of a pixel, defaults to centre if zero
              ycntrd - y coordinate of the centroid (pixel no of PSF maximum)
                       0.5 is centre of a pixel, defaults to centre if zero
              zcntrd - z coordinate of the centroid (pixel no of PSF maximum)
                       0.5 is centre of a pixel, defaults to centre if zero
              normalise - causes resulting PSF to be normalised so psf sum = 1
@OUTPUT     :
@RETURNS    : (nothing)
@DESCRIPTION: Create a 3d Gaussian with specified fwhm, center
              Return a point spread function having Gaussian profiles,
              as 3D volumetric-data.
              For example, to create a 31x31x31 array containing a normalised
              centered gaussian with an X, Y and Z fwhm's of 4.3, 3.6 and 4.0
              respectively, with std dev = 1.0 and default centering, call as:
                 psfGaussian(31, 31, 31, 0, 0, 0, 1, 1, 1, 0, 0, 0, TRUE)

              If specifing the filter size using FWHM (Full Width Half Max),
              the voxel size is taken into account when constructing the
              filter-matrix. Usual filter sizes for PET-studies are 8mm FWHM.
              If specifing the filter size using SD (standard deviation),
              values of 0.5-1.5 are sensible.
@METHOD     :
@GLOBALS    :
@CALLS      :
@CREATED    : March 2008, Paul Lyon
@MODIFIED   : Written, Frank Varosi NASA/GSFC 1991.
              Changed GAUSSIAN to ASTRO_GAUSSIAN to avoid confusion with our
              routine of same name, Alex Schuster, Aug. 98.
              Ported to C by Paul Lyon 2008 for 3D (volume) case only
---------------------------------------------------------------------------- */

void psfGaussian(float*** psf, int xdim, int ydim, int zdim,
                 double xstep, double ystep, double zstep,
                 float xfwhm, float yfwhm, float zfwhm,
                 float xstdev, float ystdev, float zstdev,
                 float xcntrd, float ycntrd, float zcntrd, int normalise)
{
  /* Initialisation... */
  int i, j, k;
  float total;
  float* x = malloc(xdim * sizeof(float));
  float* y = malloc(ydim * sizeof(float));
  float* z = malloc(zdim * sizeof(float));
  float* psfx = malloc(xdim * sizeof(float));
  float* psfy = malloc(ydim * sizeof(float));
  float* psfz = malloc(zdim * sizeof(float));

  if (psf == NULL)
  {
    fprintf(stderr,"psfGaussian was called with NULL matrix variable!");
    exit(1);
  }

  /* If the centre has not been specified, use the default centre */
  if (xcntrd == 0 && ycntrd == 0 && zcntrd == 0)
  {
    /* FIXME - image is scued if I use this code? 
    xcntrd = (xdim-1.0)/2.0;
    ycntrd = (ydim-1.0)/2.0;
    zcntrd = (zdim-1.0)/2.0;
    */
    xcntrd = 0.0;
    ycntrd = 0.0;
    zcntrd = 0.0;
  }
  else
  {
    /* Remove 0.5 to centre on the middle of a pixel */
    xcntrd-=0.5;
    ycntrd-=0.5;
    zcntrd-=0.5;
  }

  /* Set up the standard deviations */
  if (xstdev == 0.0)
  {
    if (xfwhm == 0.0)
    {
      fprintf(stderr,"Both xfwhm and xstdev are unset!");
      exit(1);
    }
    else
    {
      xstdev=fabs(xfwhm / xstep) / (2.0 * sqrt (2.0 * log(2.0) ));
    }
  }
  if (ystdev == 0.0)
  {
    if (yfwhm == 0.0)
    {
      fprintf(stderr,"Both yfwhm and ystdev are unset!");
      exit(1);
    }
    else
    {
      ystdev=fabs(yfwhm / ystep) / (2.0 * sqrt (2.0 * log(2.0) ));
    }
  }
  if (zstdev == 0.0)
  {
    if (zfwhm == 0.0)
    {
      fprintf(stderr,"Both zfwhm and zstdev are unset!");
      exit(1);
    }
    else
    {
      zstdev=fabs(zfwhm / zstep) / (2.0 * sqrt (2.0 * log(2.0) ));
    }
  }
  printf("Gaussian filter parameters are:\n");
  printf("  step (x,y,t)=(%2.2f,%2.2f,%2.2f)\n",xstep,  ystep,  zstep);
  printf("  fwhm (x,y,t)=(%2.2f,%2.2f,%2.2f)\n",xfwhm,  yfwhm,  zfwhm);
  printf("  stdev(x,y,t)=(%2.2f,%2.2f,%2.2f)\n",xstdev, ystdev, zstdev);
 printf("  cntrd(x,y,t)=(%2.2f,%2.2f,%2.2f)\n",xcntrd, ycntrd, zcntrd);

  /* Set up the x, y, z temporary arrays */
  for (i=0; i<xdim; i++) x[i] = (float)i-xcntrd;
  for (j=0; j<ydim; j++) y[j] = (float)j-ycntrd;
  for (k=0; k<zdim; k++) z[k] = (float)k-zcntrd;

  /* Generate the 3D Gaussian matrix */
  astroGaussian(x, xdim, psfx, GAUSS_MAX, GAUSS_MEAN, xstdev);
  astroGaussian(y, ydim, psfy, GAUSS_MAX, GAUSS_MEAN, ystdev);
  astroGaussian(z, zdim, psfz, GAUSS_MAX, GAUSS_MEAN, zstdev);

  for (k=0; k<zdim; k++)
  {
    for (j=0; j<ydim; j++)
    {
      for (i=0; i<xdim; i++)
      {
        /* FIXME - this is an assumption - look at original code, IDL shorthand?? */
        psf[i][j][k] = psfx[i] * psfy[j] * psfz[k];
      }
    }
  }

  /* Normalise such that the sum of matrix values = 1 */
  if (normalise == TRUE)
  {
    total=0.0;
    for (i=0; i<xdim; i++)
      for (j=0; j<ydim; j++)
        for (k=0; k<zdim; k++)
          total+=psf[i][j][k];
    for (i=0; i<xdim; i++)
      for (j=0; j<ydim; j++)
        for (k=0; k<zdim; k++)
          psf[i][j][k]/=total;
  }

  /* Clean up... */
  free(psfx);
  free(psfy);
  free(psfz);
  free(x);
  free(y);
  free(z);
}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : astroGaussian
@INPUT      : xi    - array, independent variable of Gaussian function.
              xidim - the dimension of xi
              gauss - array to hold the result
              max   - maximum value (factor) of Gaussian,
              mean  - mean value (center) of Gaussian,
              sigma - standard deviation (sigma) of Gaussian.
@OUTPUT     : gauss - 1d array
@RETURNS    : (nothing)
@DESCRIPTION: Compute the 1-d Gaussian function calculated at indices of xi
              For example, to evaulate a Gaussian centered at x=0, with sigma=1,
              and a peak value of 10 at the points 0.5 and 1.5 call with the
              following parameters:
                astroGaussian([0.5,1.5], 2, gauss, 10.0, 0.0, 1.0)
                Result ==> [8.825,3.25]
@METHOD     :
@GLOBALS    :
@CALLS      :
@CREATED    : March 2008, Paul Lyon
@MODIFIED   : Written, Frank Varosi NASA/GSFC 1991.
              Changed GAUSSIAN to ASTRO_GAUSSIAN to avoid confusion with our
              routine of same name, Alex Schuster, Aug. 98.
              Ported to C by Paul Lyon 2008 for 3D (volume) case only
---------------------------------------------------------------------------- */

void astroGaussian(float* xi, int xidim, float* gauss,
                   float max, float mean, float sigma)
{
  int i;

  for (i=0; i<xidim; i++)
  {
    gauss[i]=(xi[i]-mean)/sigma;
    gauss[i]=gauss[i]*gauss[i];

    /* We do not have to worry about overflow here ...               */
    /* exp(x) computes e to the power of x (base e exponential of x) */
    /* As x increases, exp(-x) tends to zero                         */
    gauss[i]=exp(-(double)gauss[i]/2.0) * max;
  }
}

