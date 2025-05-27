#include "slab.hpp"

using namespace std; 
// declare global variables, these are fixed physical constants
   // float mu0,eps0,e, mi, me, dt; //added some of the slab.in vars, makes my life simpler and sets it up in case we want to use it

// All variables in next 5 lines read in by init.
    // float kappaT, kappan; 
    // int im, jm, km;
    // int nm;
    // int gyropts;
    // int nsnap;

//--------------------------------------------------------------------------------------------------------------------------------------------------------
int main() { 
  auto start = std::chrono::high_resolution_clock::now();

  int m,mm;//,im,jm,km;
  int n;
  int im, jm, km;
  int nm;
  int gyropts;
  int nsnap;
// all variables in next 5 lines read in by init.  
  float lx,ly,lz;
  float vt,B=1.,Mi=1.,Me=0.01;
  float Ti=1.,Te=1.;
  float mu0,eps0,e, mi, me, dt; 
  float kappaT, kappan; 
    
  // set parameters
  FieldConst f;
  NumericalConst numc; //use these in init instead?
  ParticleConst p;

  init(B, nm, dt, nsnap, mm, lx, ly, lz, mu0, eps0, e, mi, me, im, jm, km, kappan, kappaT,gyropts);
  f = {B, mu0, eps0, lx, ly, lz, nsnap};
  numc = {im, jm, km, nm, dt, mm};
  p = {e, gyropts, kappaT, kappan, Mi, Me};
  vt=sqrt(Ti/Mi);
  
//  float x[mm],y[mm],z[mm];
//  float vx[mm],vpar[mm],mu[mm];
  float xn[mm],yn[mm],zn[mm],vzn[mm],wn[mm],mu[m];
  float xm[mm],ym[mm],zm[mm],vzm[mm],wm[mm];
  float xp[mm],yp[mm],zp[mm],vzp[mm],wp[mm];

  Array3D<float> den(im,jm,km);
  Array3D<float> phi(im,jm,km);
  Array3D<float> Ex(im,jm,km); 
  Array3D<float> Ey(im,jm,km); 
  Array3D<float> Ez(im,jm,km);

// load particles
  printf("before load\n");
  load(xn, yn, zn, vzn, mu, xp, yp, zp, vzp, xm, ym, zm, vzm, vt, f, numc, p);
  printf("before deposit\n");
//  deposit the particles on the grid
  deposit(xn,yn,zn,den,f,numc);
  printf("after deposit\n");
  
  print3DArray("den.csv", den);
  printParticleValues(mm, xn, yn, zn, vzn, mu);
  
//    exit(0);
// Calculate E from phi 
  grad(phi, Ex, Ey, Ez, f, numc, p); //E,phi,im,jm,km,dx,dy,dz
  
// main time loop
  for (n = 0; n <= nm; n++) {
    printf("timestep n=%d\n", n);
  // pre-push
	  printf("before ppush\n");
    push('p',xm,ym,zm,mu,vzm,wm,xn,yn,zn,vzn,wn,
    xp,yp,zp,vzp,wp, mm, Ex, Ey, Ez, im, jm, km, lx, ly, lz, dt, e, Mi, Ti, B);
    printf("after ppush\n");
    deposit(xp, yp, zp, den, f, numc);
    //poisson(den, phi, im, jm, km, lx, ly, lz);
  // grid
  // corrector push
    push('c',xm,ym,zm,mu,vzm,wm,xn,yn,zn,vzn,wn,
    xp,yp,zp,vzp,wp,mm,Ex, Ey, Ez, im,jm,km,lx,ly,lz,dt,e,Mi,Ti,B);	
    deposit(xn, yn, zn, den, f, numc);
    //poisson(den, phi, im, jm, km, lz, ly, lz);
    grad(phi, Ex, Ey, Ez, f, numc, p);
  }

   CleanupFFT();

  auto end = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

  std::cout << "Execution time: " << duration.count() << " milliseconds" << std::endl;

  return 0;
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------

void load(float xn[],float yn[],float zn[],float vzn[],float mu[],
	float xp[],float yp[],float zp[],float vzp[],
	float xm[],float ym[],float zm[],float vzm[],
	float vt, FieldConst &f, NumericalConst &numc, ParticleConst &p) {
  int m; 
  float r1,r2,r3,vperp2;
        
  for (m = 0; m < numc.mm; ++m){
//     xn[m]=lx*(float)rand()/(float)RAND_MAX;
//     yn[m]=ly*(float)rand()/(float)RAND_MAX;
//     zn[m]=lz*(float)rand()/(float)RAND_MAX;
     xn[m]=f.lx*revers(m+1,2);
     yn[m]=f.ly*revers(m+1,3);
     zn[m]=f.lz*revers(m+1,5);
  
    // for the velocity distribution, use Box-Muller transform (BMT). The BMT gives two gaussian #'s from two uniform #'s.
    // for vperp we can use one random #.  for v_parallel, we will use two and throw away one.
//    r1=(float)rand()/(float)RAND_MAX;
//    r2=(float)rand()/(float)RAND_MAX;
//    r3=(float)rand()/(float)RAND_MAX;
    r1=revers(m+1,7);
    r2=revers(m+1,11);
    r3=revers(m+1,13);
    vzn[m] = vt*sqrt(-2.0 * log(r1)) * cos(2.0 * M_PI * r2);
    vperp2=vt*vt*(-2.*log(r3));
    mu[m]=0.5*p.Mi*vperp2/f.B;

    xp[m]=xn[m];yp[m]=yn[m];zp[m]=zn[m];vzp[m]=vzn[m];
    xm[m]=xn[m];ym[m]=yn[m];zm[m]=zn[m];vzm[m]=vzn[m];
  }
}

void grad(Array3D<float> &Ex, Array3D<float> &Ey, Array3D<float> &Ez, Array3D<float> &phi, 
		  FieldConst &f, NumericalConst &numc, ParticleConst &p){
    float dx=f.lx/static_cast<float>(numc.im);
    float dy=f.ly/static_cast<float>(numc.jm); 
    float dz=f.lz/static_cast<float>(numc.km);

  for(int i = 0; i < numc.im; ++i){
    for(int j = 0; j < numc.jm; ++j){
      for(int k = 0; k < numc.km; ++k){
        if(i == 0){
          Ex(0, j, k) = (phi(numc.im-1, j, k) - phi(1, j, k))/(2.*dx);
        }
        else if(i == numc.im-1){
          Ex(numc.im, j, k) = (phi(numc.im-2, j, k)-phi(0, j, k))/(2.*dx);
        }
        else{
          Ex(i,j,k) = (phi(i-1, j, k) - phi(i+1, j, k))/(2.*dx);
        }
       
        if(j == 0){
          Ey(i,0,k) = (phi(i, numc.jm-1, k) - phi(i, 1, k))/(2.*dy);
        }
        else if(j == numc.jm-1){
          Ey(i, numc.jm-1, k) = (phi(i, numc.jm-2, k) - phi(i, 0, k))/(2.*dy);
        } else{
          Ey(i,j,k) = (phi(i, j-1, k) - phi(i, j+1, k))/(2.*dy);  
        }

        if(k == 0){
          Ez(i, j, 0) = (phi(i, j, numc.km-1) - phi(i, j, 1))/(2.*dz);
        }
        else if(k == numc.km-1){
          Ez(i,j,numc.km-1) = (phi(i,j,numc.km-1) - phi(i,j,0))/(2.*dz);
        }
        else{
          Ez(i,j,k) =  (phi(i, j, k-1) - phi(i, j, k+1))/(2.*dz);
        }
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------- 

// deposits density from particles on the grid
void deposit(float x[],float y[],float z[], Array3D<float> &den, FieldConst &f, NumericalConst &numc) {      
// assumes particles are in bounds. no bounds checking. x>=0, x<lx, etc.

  int m,i,j,k,ip,jp,kp;
  float wght; // particle weight factor
  float dx,dy,dz,wx,wy,wz,wxp,wyp,wzp;
  
  dx=f.lx/(float)numc.im; dy=f.ly/(float)numc.jm; dz=f.lz/(float)numc.km;
  wght=1./(dx*dy*dz);  
  for  (i=0;i<numc.im;++i){
    for  (j=0;j<numc.jm;++j){
      for  (k=0;k<numc.km;++k){
        den(i, j, k)=0;
      }
    }
  }
  for (m=0; m<numc.mm; m++){
    i=(int)(x[m]/dx); 
    j=(int)(y[m]/dy); 
    k=(int)(z[m]/dz);
    wx=(float)(i+1)-x[m]/dx; wy=(float)(j+1)-y[m]/dy; wz=(float)(k+1)-z[m]/dz;
    wxp=1.-wx; wyp=1.-wy; wzp=1.-wz;
    ip=(i+1) % numc.im; jp=(j+1) % numc.jm; kp=(k+1) % numc.km;
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

//--------------------------------------------------------------------------------------------------------------------------------------------------------

//Iniailizes Slab for run (All this does at the moment is call the function readParams)
void init(float &B, int &nm, float &dt, int &nsnap, int &mm, float &lx, float &ly, float &lz,
                float &mu0, float &eps0, float &e, float &Mi, float &Me, int &im, int &jm, int &km, float &kappaT, float &kappan, int &gyropts){
  readParams(B,nm, dt, nsnap, mm, lx, ly, lz, mu0, eps0, e, Mi, Me, im, jm, km, kappaT, kappan, gyropts);
}

//--------------------------------------------------------------------------------------------------------------------------------------------------------

//reads and initializes parameter values from the file slab.in
void readParams(float &B, int &nm, float &dt, int &nsnap, int &mm, float &lx, float &ly, float &lz,
                float &mu0, float &eps0, float &e, float &Mi, float &Me, int &im, int &jm, int &km, float &kappaT, float &kappan, int &gyropts){
  FILE *in_file = fopen("slab.in", "r");
  if (in_file == NULL)
  {  
    printf("Error! Could not open file: slab\n");
    exit(-1); 
  }else
  {
    //Doing these by line: %*x will ignore the value, %f will map to variable in order below.
    //As example, first scan below ignores the 6 var names and reads the 6 floats. Each deals with these kinds of lines
    fscanf(in_file, "%*[^\n]\n");
    fscanf(in_file, "%f %f %f %f %f %f", &mu0, &eps0, &e, &B, &Mi, &Me);
    fscanf(in_file, " %*[^\n]\n"); //space at front of comment is intentional - ignores empty line
    fscanf(in_file, "%d %f %f %d", &nm, &dt, &B, &nsnap);
    fscanf(in_file, " %*[^\n]\n");
    fscanf(in_file, "%d %d %d %f %f %f", &im, &jm, &km, &lx, &ly, &lz);
    fscanf(in_file, " %*[^\n]\n");
    fscanf(in_file, "%d %f %f %d", &mm, &kappaT, &kappan, &gyropts);
  }
  fclose(in_file);
}

//--------------------------------------------------------------------------------------------------------------------------------------------------------

//prints values to arg1 (file name) using arg2 (array)
void print3DArray(string fileName, Array3D<float> &array){
  // Open file and set premission to write only
  ofstream outFile(fileName); 
  if (outFile.fail())
    {  
      printf("Error! Could not open file: den.csv\n");
      exit(-1);
    }
  int imax = array.getX();
  int jmax = array.getY();
  int kmax = array.getZ();
  for(int i = 0; i < imax; ++i){
    for(int j = 0; j < jmax; ++j){
      for(int k = 0; k < kmax; ++k){
          outFile << setprecision(5) << array(i, j, k) << "\n";
      }
    }
  }
  outFile.close();
}

//--------------------------------------------------------------------------------------------------------------------------------------------------------

//this function will print the values: x,y,z,vpar, and mu (in that order) to the file particle.csv
void printParticleValues(int mm, float x[], float y[], float z[], float vpar[], float mu[]){
  ofstream outFile("particle.dat"); 
  if (outFile.fail())
    {  
      printf("Error! Could not open file: particle.dat\n");
      exit(-1);
    }
  outFile << "x,   y,    z,     vpar,     mu\n";
  for(int i = 0; i < mm; ++i){
    outFile << setprecision(5) << x[i] << "  " <<  y[i] << "  " <<  z[i] << "  " <<  vpar[i] << "  " <<  mu[i];
    if(i+1 < mm){
      outFile << "\n";
    }
  }
  outFile.close();
}

//------------------------------------------------------------------------------------------------------------

double revers(const int& num, const int& n){
   double rev = 0.0;
   double power = 1.0;
   int inum = num;
   int iquot = 0;
   int irem = 0;

   while(inum > 0){
      iquot = int(inum/n);
      irem = inum - n*iquot;
      power = power/n;
      rev = rev + irem*power;
      inum = iquot;
   }
   return rev;
}
