// MAKING CHANGE TO SEE IF GITHUB WORKS!!
#include "slab.hpp"

// pushers for 3d slab gk simulation. 
// predictor push - push(flag='p'), corrector push - push(flag='c')

void push(char flag,
float xm[],float ym[],float zm[],float mu[],float vzm[],float wm[],
float xn[],float yn[],float zn[],float vzn[],float wn[],
float xp[],float yp[],float zp[],float vzp[],float wp[],
int mm,Array3D<float> Ex, Array3D<float> Ey,Array3D<float> Ez,
int im,int jm,int km,float lx,float ly,float lz,float dt,
float e, float Mi,float Ti,float B) { 
     
// Argument list is a mess.  We need to fix this without a lot of class structure.  
// More discussion needed.  I am going to use a simple if/than for predictor corrector push.
// To combine predictor-push and corrector-push into one function would make things tricky
// to understand. So, I just have it broken into two using if/than.

// xm is x_n-1, xp is x_n+1, xn is x_n, etc...

  int m,i,j,k,ip,jp,kp,l;
  float dx,dy,dz,wx,wy,wz,wxp,wyp,wzp;
  float xt,yt,zt,Exp,Eyp,Ezp;
  float kapn,kap,kapt,vfac,rho;
  float rwx[4]={-1.,1.,0.,0.};
  float rwy[4]={0.,0.,-1.,1.};
  
  dx=lx/(float)im; dy=ly/(float)jm; dz=lz/(float)km;

// predictor step begins
if(flag=='p'){
  for (m=0; m<mm; m++){
	k=(int)(zn[m]/dz);
	 kp=(k+1) % km;
	wz=(float)(k+1)-zn[m]/dz;
	wzp=1.-wz;
	Exp=0.;Eyp=0.;Ezp=0.;
	rho=sqrt(2.*Mi*mu[m]/(e*e*B)); // assuming B=1, mu0=1, q=1
	for(l=0;l<4;l++){  			// 4-point gyro-average
		xt=xn[m]+rwx[l]*rho;
		yt=yn[m]+rwy[l]*rho;
		if(xt>=lx){xt=xt-lx;};
		if(xt<0.){xt=xt+lx;};
		if(yt>=lx){yt=yt-ly;};
		if(yt<0.){yt=yt+ly;};
		i=(int)(xt/dx); j=(int)(yt/dy);
		ip=(i+1) % im; jp=(j+1) % jm;
		wx=(float)(i+1)-xt/dx; wy=(float)(j+1)-yt/dy;
		wxp=1.-wx; wyp=1.-wy;
		Exp=Exp+wx*wy*wz*Ex(i,j,k)+wxp*wy*wz*Ex(ip,j,k)+wx*wyp*wz*Ex(i,jp,k)
			+wxp*wyp*wz*Ex(ip,jp,k)+wx*wy*wzp*Ex(i,j,kp)+wxp*wy*wzp*Ex(ip,j,kp)
			+wx*wyp*wzp*Ex(i,jp,kp)+wxp*wyp*wzp*Ex(ip, jp,kp);
		Eyp=Eyp+wx*wy*wz*Ey(i,j,k)+wxp*wy*wz*Ey(ip,j,k)+wx*wyp*wz*Ey(i,jp,k)
			+wxp*wyp*wz*Ey(ip,jp,k)+wx*wy*wzp*Ey(i,j,kp)+wxp*wy*wzp*Ey(ip,j,kp)
			+wx*wyp*wzp*Ey(i,jp,kp)+wxp*wyp*wzp*Ey(ip, jp,kp);
		Ezp=Ezp+wx*wy*wz*Ez(i,j,k)+wxp*wy*wz*Ez(ip,j,k)+wx*wyp*wz*Ez(i,jp,k)
			+wxp*wyp*wz*Ez(ip,jp,k)+wx*wy*wzp*Ez(i,j,kp)+wxp*wy*wzp*Ez(ip,j,kp)
			+wx*wyp*wzp*Ez(i,jp,kp)+wxp*wyp*wzp*Ez(ip, jp,kp);
//do not need to gyro-average Ez, but it makes simper code
	}
//	printf("after 4-pt avg\n");
	Exp=0.25*Exp; Eyp=0.25*Eyp; Ezp=0.25*Ezp;
	vzp[m]=vzm[m]+2.*dt*e*Ezp;
	vfac=0.5*vzn[m]*vzn[m]+ (e*B/Mi)*mu[m];
// normalization in slab code is hard to figure out in dimensional units.  I will come 
// back to this later... I will need to define Ti at least.
	kap=kapn-(1.5-vfac)*kapt;
	xp[m]=xm[m]+2.*dt*Eyp; // B is positive and in the z-direction
	yp[m]=ym[m]+2.*dt*(-Exp);
	zp[m]=zm[m]+2.*dt*vzn[m];

	wp[m] = wm[m]+2.*dt*(kap*Eyp + e*Exp*vzn[m]);

// Some trickiness here.  We no longer need x_n-1, so we store 
// x_n + 0.5*dt*Eyp_n, that we can use going into the corrector-push to speed it up
// otherwise we would have to calculate the 4 pt avg field at time level n again...
	xm[m]=xn[m]+0.25*(xp[m]-xm[m]);
	ym[m]=yn[m]+0.25*(yp[m]-ym[m]);
	zm[m]=zn[m]+0.25*(zp[m]-zm[m]);
	vzm[m]=vzn[m]+0.25*(vzp[m]-vzm[m]);
	wm[m]=wn[m]+0.25*(wp[m]-wm[m]);

// Now enforce periodicity, so that all particles are in the domain going into deposit()
// May only need to do this in the z direction if 4-pt avg is done...
// using fmod this way, we are assuming a particle does not move more than -Lz in one timestep
	xp[m]=fmod(xp[m]+lx,lx);
	yp[m]=fmod(yp[m]+ly,ly);
	zp[m]=fmod(zp[m]+lz,lz);


  } // end of particle loop m
} // end of predictor step if/than

//corrector step begins
if(flag=='c'){
	printf("corrector-push not written");
} // end of corrector step if/than

} 
