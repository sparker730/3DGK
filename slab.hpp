#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <chrono>
#include <string>
#include <complex>
#include <vector>
#include <fftw3.h>
#include "MultiArrays.hpp"
#pragma once
struct FieldConst{
  float B;
  float mu0; 
  float eps0;
  float lx;
  float ly;
  float lz;
  int nsnap;
};

struct NumericalConst{
  int im;
  int jm;
  int km;
  int nm; 
  float dt;
  int mm;
};

struct ParticleConst{      
  float e;
  int gyropts;
  float kappaT;
  float kappan;
  float Mi;
  float Me;
};
void load(float xn[],float yn[],float zn[],float vzn[],float mu[],
	float xp[],float yp[],float zp[],float vzp[],
	float xm[],float ym[],float zm[],float vzm[],
	float vt, FieldConst &f, NumericalConst &numc, ParticleConst &p);
void deposit(float x[],float y[],float z[], Array3D<float> &den, FieldConst &f, NumericalConst &numc);
void init(float &B, int &nm, float &dt, int &nsnap, int &mm, float &lx, float &ly, float &lz,
                float &mu0, float &eps0, float &e, float &Mi, float &Me, int &im, int &jm, int &km, float &kappaT, float &kappan, int &gyropts);
void readParams(float &B, int &nm, float &dt, int &nsnap, int &mm, float &lx, float &ly, float &lz,
                float &mu0, float &eps0, float &e, float &Mi, float &Me, int &im, int &jm, int &km, float &kappaT, float &kappan, int &gyropts);
void print3DArray(std::string fileName, Array3D<float> &den);
void printParticleValues(const int mm,float x[],float y[], float z[], float vpar[], float mu[]);
double revers(const int& num, const int& n);
void push(char flag,
float xm[],float ym[],float zm[],float mu[],float vzm[],float wm[],
float xn[],float yn[],float zn[],float vzn[],float wn[],
float xp[],float yp[],float zp[],float vzp[],float wp[],
int mm,Array3D<float> &Ex, Array3D<float> &Ey,Array3D<float> &Ez,
int im,int jm,int km,float lx,float ly,float lz,float dt,
float e, float Mi,float Ti,float B);
void grad(Array3D<float> &Ex, Array3D<float> &Ey, Array3D<float> &Ez, Array3D<float> &phi, 
		  FieldConst &f, NumericalConst &numc, ParticleConst &p);
// void poisson(Array3D<float> &den, Array3D<float> &phi, const int im, const int jm, const int km,
//              const float lx, const float ly, const float lz);
// void CleanupFFT();
// void srcbes(const float biz, float& gam0, float& gam1);