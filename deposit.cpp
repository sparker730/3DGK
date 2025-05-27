
#include "slab.hpp"

// deposits density from particles on the grid
void deposit(float x[],float y[],float z[],float mu[],float w[],int mm,Array3D<float> &den,int im,int jm,int km,float lx,float ly,float lz,float e,float Mi,float B) {      
// assumes particles are in bounds. no bounds checking. x>=0, x<lx, etc.

  int m,i,j,k,ip,jp,kp;
  float wghtnorm,wght; // particle weight factor
  float dx,dy,dz,wx,wy,wz,wxp,wyp,wzp;
  float rwx[4]={-1.,1.,0.,0.};
  float rwy[4]={0.,0.,-1.,1.};
  float xt,yt,rho;
  int l;
  
  dx=lx/(float)im; dy=ly/(float)jm; dz=lz/(float)km;
  wghtnorm=0.25/(dx*dy*dz);  
  for  (i=0;i<im;++i){
    for  (j=0;j<jm;++j){
      for  (k=0;k<km;++k){
        den(i, j, k)=0;
      }
    }
  }
  for (m=0; m<mm; m++){
	wght=wghtnorm*w[m];
	k=(int)(z[m]/dz);
	kp=(k+1) % km;
	wz=(float)(k+1)-z[m]/dz;
	wzp=1.-wz;
	rho=sqrt(2.*Mi*mu[m]/(e*e*B)); // assuming B=1, mu0=1, q=1
	for(l=0;l<4;l++){  			// 4-point gyro-average
		xt=x[m]+rwx[l]*rho;
		yt=y[m]+rwy[l]*rho;
		if(xt>=lx){xt=xt-lx;};
		if(xt<0.){xt=xt+lx;};
		if(yt>=lx){yt=yt-ly;};
		if(yt<0.){yt=yt+ly;};
		i=(int)(xt/dx); j=(int)(yt/dy);
		ip=(i+1) % im; jp=(j+1) % jm;
		wx=(float)(i+1)-xt/dx; wy=(float)(j+1)-yt/dy;
		wxp=1.-wx; wyp=1.-wy;
    		ip=(i+1) % im; jp=(j+1) % jm; kp=(k+1) % km;
    		den(i, j, k)=den(i, j, k)+wght*wx*wy*wz;
    		den(ip, j, k)=den(ip, j, k)+wght*wxp*wy*wz;
    		den(i, jp, k)=den(i, jp, k)+wght*wx*wyp*wz;
    		den(ip, jp, k)=den(ip, jp, k)+wght*wxp*wyp*wz;
    		den(i, j, kp)=den(i, j, kp)+wght*wx*wy*wzp;
    		den(ip, j, kp)=den(ip, j, kp)+wght*wxp*wy*wzp;
    		den(i, jp, kp)=den(i, jp, kp)+wght*wx*wyp*wzp;
    		den(ip, jp, kp)=den(ip, jp, kp)+wght*wxp*wyp*wzp;
	}
  }
}
