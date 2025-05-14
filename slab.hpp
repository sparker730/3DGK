#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <chrono>
#include <string>
#include "MultiArrays.hpp"
#pragma once

void load(float xn[],float yn[],float zn[],float vzn[],float mu[],
	float xp[],float yp[],float zp[],float vzp[],
	float xm[],float ym[],float zm[],float vzm[],
	int mm,float lx,float ly,float lz,float vt, float B, float Mi);
void deposit(float x[],float y[],float z[],int mm, Array3D<float> den,int im,int jm,int km,float lx,float ly,float lz);
void init(float &B, int &nm, float &dt, int &nsnap, int &mm, float &lx, float &ly, float &lz);
void readParams(float &B, int &nm, float &dt, int &nsnap, int &mm, float &lx, float &ly, float &lz);
void print3DArray(std::string fileName, Array3D<float> den);
void printParticleValues(const int mm,float x[],float y[], float z[], float vpar[], float mu[]);
double revers(const int& num, const int& n);
void push(char flag,
float xm[],float ym[],float zm[],float mu[],float vzm[],float wm[],
float xn[],float yn[],float zn[],float vzn[],float wn[],
float xp[],float yp[],float zp[],float vzp[],float wp[],
int mm,Array3D<float> Ex, Array3D<float> Ey,Array3D<float> Ez,
int im,int jm,int km,float lx,float ly,float lz,float dt,
float e, float Mi,float Ti,float B);
