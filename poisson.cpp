#include "slab.hpp"

//Globals so only generated once. Can be cleaned up later.
float* formphi = nullptr;
Array3D_Old<float> formphi3D;
//FFT logic can probably move into a class.
bool fftwPlan = false;
fftwf_plan phiFFT = nullptr; //Note, use extra f in fftwf (everywhere!) for floating point. (And in makefile!)
fftwf_plan phiIFFT = nullptr;
fftwf_complex* phik = nullptr; //Complex output array (single precision)
Array3D_Old<std::complex<float>> phik3D;

void poisson(Array3D<float>& den, Array3D<float>& phi, const int im, const int jm, const int km,
             const float lx, const float ly, const float lz) {

   //Calculation of 1/k^2 only needs to happen once I believe?
   if (!formphi) {
       //Note this logic can all be handled in the 3D array class constructor.
      float* formphi = new float[im*jm*km];
      std::fill(formphi, formphi + im*jm*km, 0);
      formphi3D.CreateArray3D(formphi,im,jm,km);

      std::cout << "Calculate ksq." << std::endl;

      //Are these the correct ranges? Assuming in ORB indices were 0:imx,0:jmx,0:kmx inclusive in the Fortran allocations. This means
      //C++ arrays should be (imx+1)*(jmx+1)*(kmx+1) not imx*jmx*kmx to use this logic? Same for all other loop/indexing logic below...
      //From the fourier space loop below it seems phi and formphi are indexed differently though so its confusing.
      for (int l = 0; l <= im-1; ++l) {
         for (int m = 0; m <= jm-1; ++m) {
            for (int n = 0; n <= km-1; ++n) {
               int l1 = l;
               int m1 = m;
               int n1 = n;

               if (l1 > (im/2)) l1 = im - l1;
               if (m1 > (jm/2)) m1 = jm - m1;
               if (n1 > (km/2)) n1 = km - n1;

               float kx = 2*M_PI*static_cast<float>(l1)/lx;
               float ky = 2*M_PI*static_cast<float>(m1)/ly;
               float kz = 2*M_PI*static_cast<float>(n1)/lz;

               //Note these are not different than the above?
               float skx = 2*M_PI*static_cast<float>(l1)/lx;
               float sky = 2*M_PI*static_cast<float>(m1)/ly;
               float skz = 2*M_PI*static_cast<float>(n1)/lz;

               //Not sure what shapes should be. Assume 1.0 for slab?
               float xshape = 1.0;
               float yshape = 1.0;
               float zshape = 1.0;

               //Kind of annoying there is no ** operator in C++. But there is longer pow(kx,2),pow(xshape*skx,2),etc.
               float b = kx*kx + ky*ky + kz*kz; //Note, there was no kz here in ORB code but b2 had kz.
               float b2 = xshape*xshape*skx*skx + yshape*yshape*sky*sky + zshape*zshape*skz*skz;
               b2 = exp(-b2);

               float gam0 = 0.0;
               float gam1 = 0.0;
               srcbes(b, gam0, gam1);

               if ((b == 0.0) || (kz == 0.0)) {
                  formphi3D(l,m,n) = 0.0;
               } else {
                  //formphi3D(l,m,n) = b2/(1 + (1-gam0)*static_cast<float>(im*jm*km)); Old ORB logic.
                  formphi3D(l,m,n) = 1.0/b; //Use to divide rho/k^2 to calculate phi further on.
               }
            }
         }
      }
   }

   std::cout << "Load density" << std::endl;

   //Simply load den (since q=1) into phi.
   for (int i = 0; i <= im-1; ++i) {
      for (int j = 0; j <= jm-1; ++j) {
         for (int k = 0; k <= km-1; ++k) {
            //phi(i+1,j+1,k+1) = den(i,j,k); ORB version
            phi(i,j,k) = den(i,j,k);
         }
      }
   }

   //This fft init logic can be moved into a class or something.
   //Use FFTW_MEASURE rather than FFTW_ESTIMATE to rigorously calculate a fast option.
   //I found FFTW_ESTIMATE was slower than python FFTs with the Hasegawa-Mima code. But lots of FFTW docs use ESTIMATE.
   //Also have FFTW_PATIENT and FFTW_EXHAUSTIVE for more intense analysis.
   if (!fftwPlan) {
      std::cout << "Creating FFT plans." << std::endl;
      fftwPlan = true;
      phik = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * im * jm * (km / 2 + 1));
      phiFFT  = fftwf_plan_dft_r2c_3d(im, jm, km, phi.start(), phik, FFTW_MEASURE);
      phiIFFT = fftwf_plan_dft_c2r_3d(im, jm, km, phik, phi.start(), FFTW_MEASURE);
      phik3D.CreateArray3D(reinterpret_cast<std::complex<float>*>(phik),im,jm,km/2+1);
      std::cout << "Created FFT plans." << std::endl;
   }
   
   //std::cout << phik3D.start() << std::endl;

   //Transform to Fourier space.
   std::cout << "Transform to fourier space." << std::endl;
   fftwf_execute(phiFFT);

   //Apply form factor and ifft scaling. 
   for (int i = 1; i <= im; ++i) {
      for (int j = 1; j <= jm; ++j) {
         for (int k = 1; k <= (km/2 + 1); ++k) {
            //Scale the data by array size (im*jm*km) for ifft to give back correct values in real space. Assuming this can be done early.
            //phik3D(i,j,k) = formphi3D(i-1,j-1,k-1)*phik3D(i,j,k)*static_cast<float>(im*jm*km); ORB version. Different scaling in fft?
            phik3D(i-1,j-1,k-1) = formphi3D(i-1,j-1,k-1)*phik3D(i-1,j-1,k-1)/static_cast<float>(im*jm*km);
         }
      }
   }

   //Transform back to real space.
   std::cout << "Transform to real space." << std::endl;
   fftwf_execute(phiIFFT);

   std::cout << "Poisson solve complete." << std::endl;

   //Temporary reverse copy to use Cray fft. I imagine this was used for further code in ORB here and can be deleted?
   //for (int i = 1; i <= im; ++i) {
   //   for (int j = 1; j <= jm; ++j) {
   //      for (int k = 1; k <= km; ++k) {
   //         phi(i,j,k) = ez1(i,j,k); 
   //      }
   //   }
   //}
}

void CleanupFFT() {
   // Clean up FFT plan and allocated memory
   if (phiFFT) fftwf_destroy_plan(phiFFT);
   if (phiIFFT) fftwf_destroy_plan(phiIFFT);
   if (phik) fftwf_free(phik);
   if (formphi) delete[] formphi;
}

//Calculates gamma nought and gamma 1. (Abramowitz and Stegun).
void srcbes(const float biz, float& gam0, float& gam1) {
   double t1, t2;

   if (biz > 3.75) {
      t2 = 1/sqrt(biz);
      t1 = 3.75/biz;
      gam0 = t2*((((((((.00392377*t1 - 0.01647633)*t1 + 0.02635537)
                *t1 - 0.02057706)*t1 + 0.00916281)*t1 - 0.00157565)
                *t1 + 0.00225319)*t1 + 0.01328592)*t1 + 0.39894228);
      gam1 = t2*((((((((-0.00420059*t1 + 0.01787654)*t1 - 0.02895312)
                  *t1 + 0.02282967)*t1 - 0.01031555)*t1 + 0.00163801)
                  *t1 - 0.00362018)*t1 - 0.03988024)*t1 + 0.39894228);
   } else {
      t1 = (biz/3.75);
      t1 = t1*t1;
      t2 = exp(-biz);
      gam0 = t2*((((((0.0045813*t1 + 0.0360768)*t1 + 0.2659732)
               *t1 + 1.2067492)*t1 + 3.0899424)
               *t1 + 3.5156229)*t1 + 1.0);
      gam1 = t2*biz*((((((0.00032411*t1 + 0.00301532)*t1 + 0.02658733)
                   *t1 + 0.15084934)*t1 + 0.51498869)
                   *t1 + 0.87890594)*t1 + 0.5);
   }
}