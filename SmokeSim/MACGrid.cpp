#include "MACGrid.h"
#include "GL/glut.h"
#include "camera.h"
#include "ConjGrad.h"
#include <math.h>
#include <map>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include "boost/numeric/ublas/matrix.hpp"
#include <map>

// Globals
MACGrid target;
extern int theDim[3];
extern double theCellSize;

// NOTE: x -> cols, z -> rows, y -> stacks
MACGrid::RenderMode MACGrid::theRenderMode = SHEETS;
bool MACGrid::theDisplayVel = true;

#define FOR_EACH_CELL \
   for(int k = 0; k < theDim[MACGrid::Z]; k++)  \
      for(int j = 0; j < theDim[MACGrid::Y]; j++) \
         for(int i = 0; i < theDim[MACGrid::X]; i++) 

#define FOR_EACH_FACE \
   for(int k = 0; k < theDim[MACGrid::Z]+1; k++) \
      for(int j = 0; j < theDim[MACGrid::Y]+1; j++) \
         for(int i = 0; i < theDim[MACGrid::X]+1; i++) 

#ifdef _DEBUGA
asdg
#endif

MACGrid::MACGrid()
{
   initialize();
}

MACGrid::MACGrid(const MACGrid& orig)
{
   mU = orig.mU;
   mV = orig.mV;
   mW = orig.mW;
   mP = orig.mP;
   mD = orig.mD;
   mT = orig.mT;
   mS = orig.mS;
}

MACGrid& MACGrid::operator=(const MACGrid& orig)
{
   if (&orig == this)
   {
      return *this;
   }
   mU = orig.mU;
   mV = orig.mV;
   mW = orig.mW;
   mP = orig.mP;
   mD = orig.mD;
   mT = orig.mT;
   mS = orig.mS;

   return *this;
}

MACGrid::~MACGrid()
{
}

void MACGrid::reset()
{
   mU.initialize();
   mV.initialize();
   mW.initialize();
   mP.initialize();
   mD.initialize();
   mS.initialize();
   mT.initialize(0.0);
}

void MACGrid::initialize()
{
   reset();
}

void MACGrid::updateSources()
{
    //TODO: Set initial values for density, temperature, velocity
#ifdef _DEBUG
	mV(2,0,0) = 0.1;
	mD(2,0,0) = 1.0;
	mT(2,0,0) = 200.0;
#else
	mV(5,1,0) += 5.0;
	mD(5,0,0) += 0.5;
	
	mV(10,29,0) -= 5.0;
	mT(10,29,0) = -900;
	mD(10,29,0) += 0.5;


	mS(1,7,0) = 1.0;
	mS(2,7,0) = 1.0;
	mS(3,7,0) = 1.0;
	mS(4,7,0) = 1.0;
	mS(5,7,0) = 1.0;
	mS(6,7,0) = 1.0;
	mS(7,7,0) = 1.0;
	mS(10,7,0) = 1.0;
	mS(11,7,0) = 1.0;
	mS(12,7,0) = 1.0;
	mS(13,7,0) = 1.0;
	mS(1,7,1) = 1.0;
	mS(2,7,1) = 1.0;
	mS(3,7,1) = 1.0;
	mS(4,7,1) = 1.0;
	mS(5,7,1) = 1.0;
	mS(6,7,1) = 1.0;
	mS(7,7,1) = 1.0;
	mS(10,7,1) = 1.0;
	mS(11,7,1) = 1.0;
	mS(12,7,1) = 1.0;
	mS(13,7,1) = 1.0;
	mT(5,0,0) = 900.0;
#endif
}

void MACGrid::advectVelocity(double dt)
{
    // TODO: Calculate new velocities and store in target
	//Advect Each Velocity on Each Face
	
	//Solid Cell Boundaries
	//FOR_EACH_CELL {
	//	if(mS(i,j,k) > EPSILON){
	//		//6 Sides -> 0
	//		if(i > 0) {mU(i,j,k)=0.0;}
	//		if(i < theDim[0]) {mU(i+1,j,k)=0.0;}
	//		if(j > 0) {mV(i,j,k)=0.0;}
	//		if(j < theDim[1]) {mV(i,j+1,k)=0.0;}
	//		if(k > 0) {mW(i,j,k)=0.0;}
	//		if(k < theDim[2]) {mW(i,j,k+1)=0.0;}
	//	}
	//}

	//Advect Velocity on Each Face
	FOR_EACH_FACE
	{
		if (isFace(i,j,k,X)){
			if((i==0 || i==theDim[0]) || mS(i,j,k) > EPSILON || mS(i-1,j,k) > EPSILON) {target.mU(i,j,k)=0.0;}
			else{
				vec3 traceBack = getXFace(i,j,k) - dt*getVelocity(getXFace(i,j,k));
				traceBack = clampPoint(traceBack); // Clamp point to outside of solid cell
				target.mU(i,j,k) = getVelocityX(traceBack);
			}
		}
		if (isFace(i,j,k,Y)){
			if((j==0 || j==theDim[1])|| mS(i,j,k) > EPSILON || mS(i,j-1,k) > EPSILON) {target.mV(i,j,k)=0.0;}
			else{
				vec3 traceBack = getYFace(i,j,k) - dt*getVelocity(getYFace(i,j,k));
				traceBack = clampPoint(traceBack); // Clamp point to outside of solid cell
				target.mV(i,j,k) = getVelocityY(traceBack);
			}
		}
		if (isFace(i,j,k,Z)){
			if((k==0 || k==theDim[2])|| mS(i,j,k) > EPSILON || mS(i,j,k-1) > EPSILON) {target.mW(i,j,k)=0.0;}
			else{
				vec3 traceBack = getZFace(i,j,k) - dt*getVelocity(getZFace(i,j,k));
				traceBack = clampPoint(traceBack); // Clamp point to outside of solid cell
				target.mW(i,j,k) = getVelocityZ(traceBack);
			}
		}
	}

	////Solid Cell Boundaries
	//FOR_EACH_CELL {
	//	if(mS(i,j,k) > EPSILON){
	//		//6 Sides -> 0
	//		if(i > 0) {target.mU(i,j,k)=0.0;}
	//		if(i < theDim[0]) {target.mU(i+1,j,k)=0.0;}
	//		if(j > 0) {target.mV(i,j,k)=0.0;}
	//		if(j < theDim[1]) {target.mV(i,j+1,k)=0.0;}
	//		if(k > 0) {target.mW(i,j,k)=0.0;}
	//		if(k < theDim[2]) {target.mW(i,j,k+1)=0.0;}
	//	}
	//}

	//Then save the result to our object
    mU = target.mU;
    mV = target.mV;
    mW = target.mW;
}

vec3 MACGrid::clampPoint(vec3 p) {
	int i,j,k;
	mP.getCell(p,i,j,k);

	//If Inside Solid Cell...
	if(mS(i,j,k) > EPSILON)
	{	
		//std::cout << "Advecting From Solid Cell" << std::endl;
		//Center of Solid Cell
		vec3 center = getCenter(i,j,k);
		double half_cell = theCellSize/2.0;
		vector<double> diff_vec;
		vector<int> diff_vec_i;
		
		//Push Difference ONLY if Neighboring Cell in that Direction is Fluid Cell 
		if(mS(i+1,j,k) < 1.0){diff_vec.push_back(fabs((center[X]+half_cell)-p[X]));diff_vec_i.push_back(0);} //X Hi
		if(mS(i-1,j,k) < 1.0){diff_vec.push_back(fabs((center[X]-half_cell)-p[X]));diff_vec_i.push_back(1);} //X Lo
		if(mS(i,j+1,k) < 1.0){diff_vec.push_back(fabs((center[Y]+half_cell)-p[Y]));diff_vec_i.push_back(2);} //Y Hi
		if(mS(i,j-1,k) < 1.0){diff_vec.push_back(fabs((center[Y]-half_cell)-p[Y]));diff_vec_i.push_back(3);} //Y Lo
		if(mS(i,j,k+1) < 1.0){diff_vec.push_back(fabs((center[Z]+half_cell)-p[Z]));diff_vec_i.push_back(4);} //Z Hi
		if(mS(i,j,k-1) < 1.0){diff_vec.push_back(fabs((center[Z]-half_cell)-p[Z]));diff_vec_i.push_back(5);} //Z Lo

		//Get Minimum Difference of Possible Neighboring Cells
		double min = diff_vec[0];
		int min_i = diff_vec_i[0];
		for(int i=1;i<diff_vec.size();i++) {
			if(diff_vec[i] < min){
				min = diff_vec[i];
				min_i = diff_vec_i[i];
			}
		}
		
		//Clamp Value to Min Distance Value
		min += EPSILON; //Avoid Floating Point Boundary Errors
		switch(min_i){
			case 0: //X Hi
				return p+vec3(min,0,0);
				break;
			case 1: //X Lo
				return p-vec3(min,0,0);
				break;
			case 2://Y Hi
				return p+vec3(0,min,0);
				break;
			case 3: //Y Lo
				return p-vec3(0,min,0);
				break;
			case 4: //Z Hi
				return p+vec3(0,0,min);
				break;
			case 5: //Z Lo
				return p-vec3(0,0,min);
				break;
		}
	}
	return p;
}

bool MACGrid::isFace(int i, int j, int k, MACGrid::Direction d)
{
	switch (d)
	{
		case X: return (i >= 0 && i < theDim[X]+1 &&
			j >= 0 && j < theDim[Y] &&
			k >= 0 && k < theDim[Z]);
		case Y: return (i >= 0 && i < theDim[X] &&
			j >= 0 && j < theDim[Y]+1 &&
			k >= 0 && k < theDim[Z]);
		case Z: return (i >= 0 && i < theDim[X] &&
			j >= 0 && j < theDim[Y] &&
			k >= 0 && k < theDim[Z]+1);
	}
	printf("Error: bad direction passed to isFace\n");
	return false;
}

void MACGrid::advectTemperature(double dt)
{
    // TODO: Calculate new temp and store in target
	FOR_EACH_CELL {
		if(mS(i,j,k) > EPSILON)
			target.mT(i,j,k) = 0.0;
		else{
			vec3 traceBack = getCenter(i,j,k) - dt*getVelocity(getCenter(i,j,k));
			traceBack = clampPoint(traceBack); // Clamp point to outside of solid cell
			target.mT(i,j,k) = getTemperature(traceBack);
		}
	}
	// Then save the result to our object
    mT = target.mT;
}

void MACGrid::advectDensity(double dt)
{
    // TODO: Calculate new densitities and store in target
	FOR_EACH_CELL {
		if(mS(i,j,k) > EPSILON)
			target.mD(i,j,k) = 0.0;
		else{
			vec3 traceBack = getCenter(i,j,k) - dt*getVelocity(getCenter(i,j,k));
			traceBack = clampPoint(traceBack); // Clamp point to outside of solid cell
			target.mD(i,j,k) = getDensity(traceBack);
		}
	}
	// Then save the result to our object
    mD = target.mD;
}

void MACGrid::computeBouyancy(double dt)
{
   // TODO: Calculate bouyancy and store in target
   double a = 1.0;
   double b = 1.0;
   double T_amb = 0.0;
   
   FOR_EACH_FACE{
	   if(isFace(i,j,k,Y)){
		   if(j==0 || j==theDim[1] || mS(i,j,k) > EPSILON || mS(i,j-1,k) > EPSILON) {target.mV(i,j,k) = 0.0;}
		   else
			   target.mV(i,j,k) = mV(i,j,k)+dt*(-a*getDensity(getYFace(i,j,k))+b*(getTemperature(getYFace(i,j,k))-T_amb));
	   }
   }

   //Solid Cell Boundaries
	//FOR_EACH_CELL {
	//	if(mS(i,j,k) > EPSILON){
	//		//6 Sides -> 0
	//		if(i > 0) {target.mU(i,j,k)=0.0;}
	//		if(i < theDim[0]) {target.mU(i+1,j,k)=0.0;}
	//		if(j > 0) {target.mV(i,j,k)=0.0;}
	//		if(j < theDim[1]) {target.mV(i,j+1,k)=0.0;}
	//		if(k > 0) {target.mW(i,j,k)=0.0;}
	//		if(k < theDim[2]) {target.mW(i,j,k+1)=0.0;}
	//	}
	//}

   // and then save the result to our object
   mV = target.mV;
}

vec3 MACGrid::getOmega(int i,int j,int k)
{
	vec3 c_i1jk = getCenter(i+1,j,k);
	vec3 c_ij1k = getCenter(i,j+1,k);
	vec3 c_ijk1 = getCenter(i,j,k+1);
	vec3 c_n_i1jk = getCenter(i-1,j,k);
	vec3 c_n_ij1k = getCenter(i,j-1,k);
	vec3 c_n_ijk1 = getCenter(i,j,k-1);

	vec3 omega_ijk(
		   (getVelocityZ(c_ij1k)-getVelocityZ(c_n_ij1k))/(2.0*theCellSize)-(getVelocityY(c_ijk1)-getVelocityY(c_n_ijk1))/(2.0*theCellSize),
		   (getVelocityX(c_ijk1)-getVelocityX(c_n_ijk1))/(2.0*theCellSize)-(getVelocityZ(c_i1jk)-getVelocityZ(c_n_i1jk))/(2.0*theCellSize),
		   (getVelocityY(c_i1jk)-getVelocityY(c_n_i1jk))/(2.0*theCellSize)-(getVelocityX(c_ij1k)-getVelocityX(c_n_ij1k))/(2.0*theCellSize));
	return omega_ijk;
}

void MACGrid::computeVorticityConfinement(double dt)
{
   // TODO: Calculate vorticity confinement forces
	GridData fU;	//<--This made me chuckle and is probably inappropriate, but i left it in ;)
	GridData fV;
	GridData fW;
	fU.initialize(0.0);
	fV.initialize(0.0);
	fW.initialize(0.0);

   double epsilon = 40.0;
   FOR_EACH_CELL{
	   vec3 omega_grad_ijk(
		   (getOmega(i+1,j,k).Length()-getOmega(i-1,j,k).Length())/(2*theCellSize),
		   (getOmega(i,j+1,k).Length()-getOmega(i,j-1,k).Length())/(2*theCellSize),
		   (getOmega(i,j,k+1).Length()-getOmega(i-1,j,k-1).Length())/(2*theCellSize));
	   omega_grad_ijk.Normalize();
	   vec3 f_conf = epsilon*theCellSize*(omega_grad_ijk.Cross(getOmega(i,j,k)));
	   //std::cout << f_conf << std::endl;	   
	   //Add to Appropriate f-vector
	   fU(i,j,k) = f_conf[0];
	   fV(i,j,k) = f_conf[1];
	   fW(i,j,k) = f_conf[2];
   }

   // Apply the forces to the current velocity and store the result in target
   FOR_EACH_FACE
	{
		if (isFace(i,j,k,X)){
			if((i==0 || i==theDim[0]) || mS(i,j,k) > EPSILON || mS(i-1,j,k) > EPSILON) {target.mU(i,j,k)=0.0;}
			else
				target.mU(i,j,k) = mU(i,j,k)+dt*fU.interpolate(getXFace(i,j,k));
		}
		if (isFace(i,j,k,Y)){
			if((j==0 || j==theDim[1]) || mS(i,j,k) > EPSILON || mS(i,j-1,k) > EPSILON) {target.mV(i,j,k)=0.0;}
			else
				target.mV(i,j,k) = mV(i,j,k)+dt*fV.interpolate(getYFace(i,j,k));
		}
		if (isFace(i,j,k,Z)){
			if((k==0 || k==theDim[2]) || mS(i,j,k) > EPSILON || mS(i,j,k-1) > EPSILON) {target.mW(i,j,k)=0.0;}
			else 
				target.mW(i,j,k) = mW(i,j,k)+dt*fW.interpolate(getZFace(i,j,k));
		}
   }

   //Solid Cell Boundaries
 //  FOR_EACH_CELL {
	//	if(mS(i,j,k) > EPSILON){
	//		//6 Sides -> 0
	//		if(i > 0) {target.mU(i,j,k)=0.0;}
	//		if(i < theDim[0]) {target.mU(i+1,j,k)=0.0;}
	//		if(j > 0) {target.mV(i,j,k)=0.0;}
	//		if(j < theDim[1]) {target.mV(i,j+1,k)=0.0;}
	//		if(k > 0) {target.mW(i,j,k)=0.0;}
	//		if(k < theDim[2]) {target.mW(i,j,k+1)=0.0;}
	//	}
	//}

   // Then save the result to our object
   mU = target.mU;
   mV = target.mV;
   mW = target.mW;
}

void MACGrid::addExternalForces(double dt)
{
   computeBouyancy(dt);
   computeVorticityConfinement(dt);
}

int MACGrid::getIndex(int i, int j, int k)
{
	int col = i;
	int row = k*theDim[0];
	int stack = j*theDim[0]*theDim[2];
	return col+row+stack;
}

void MACGrid::project(double dt)
{
	//Solid Cell Boundaries -> 0 (Final Check Before Projection!);
	//Outer Boundaries
	FOR_EACH_FACE {
		if (isFace(i,j,k,X)){
			if((i==0 || i==theDim[0])) {mU(i,j,k)=0.0;}
		}
		if (isFace(i,j,k,Y)){
			if((j==0 || j==theDim[1])) {mV(i,j,k)=0.0;}
		}
		if (isFace(i,j,k,Z)){
			if((k==0 || k==theDim[2])) {mW(i,j,k)=0.0;}
		}
	}

	//Solid Cell Boundaries
	FOR_EACH_CELL {
		if(mS(i,j,k) > EPSILON){
			//6 Sides -> 0
			if(i > 0) {mU(i,j,k)=0.0;}
			if(i < theDim[0]) {mU(i+1,j,k)=0.0;}
			if(j > 0) {mV(i,j,k)=0.0;}
			if(j < theDim[1]) {mV(i,j+1,k)=0.0;}
			if(k > 0) {mW(i,j,k)=0.0;}
			if(k < theDim[2]) {mW(i,j,k+1)=0.0;}
		}
	}

	// TODO: Solve Ap = d for pressure
    // 1. Contruct d
	GridData d;
	d.initialize();
	double rho = 1.2;
	
	//Calculate Divergence
	//std::ofstream outPreDiv("pre_div.txt");
	FOR_EACH_CELL
	{
		d(i,j,k) = -1*((rho*theCellSize)/dt)*((mU(i+1,j,k)-mU(i,j,k)) + (mV(i,j+1,k)-mV(i,j,k)) + (mW(i,j,k+1)-mW(i,j,k))); //Discrete Divergence
		//outPreDiv << (1.0/theCellSize)*((mU(i+1,j,k)-mU(i,j,k)) + (mV(i,j+1,k)-mV(i,j,k)) + (mW(i,j,k+1)-mW(i,j,k))) << std::endl;
		//outPreDiv << -1*((rho*theCellSize)/dt)*((mU(i+1,j,k)-mU(i,j,k)) + (mV(i,j+1,k)-mV(i,j,k)) + (mW(i,j,k+1)-mW(i,j,k))) << std::endl;
	}
	//outPreDiv.close();

	// 2. Construct A
	ublas::matrix<double> A(theDim[0]*theDim[1]*theDim[2],theDim[0]*theDim[1]*theDim[2]);

	//Initialize Matrix to 0's
	A.clear();

	//Determine Row Coefficients
	FOR_EACH_CELL {
		//If Solid Cell -> Ignore and Go to Next Cell
		if(mS(i,j,k) > EPSILON)
			continue;
		//Pijk Index
		int p_ijk = getIndex(i,j,k);
		
		//Neighbors Index
		int p_i1jk =  getIndex(i+1,j,k);
		int p_ij1k =  getIndex(i,j+1,k);
		int p_ijk1 =  getIndex(i,j,k+1);
		int p_n_i1jk = getIndex(i-1,j,k);
		int p_n_ij1k =  getIndex(i,j-1,k);
		int p_n_ijk1 = getIndex(i,j,k-1);

		//Coefficients
		A(p_ijk,p_ijk) = 6.0;
		
		//Decrease Count for Neighboring Boundary Cells
		if(i+1 >= theDim[0] || mS(i+1,j,k) > EPSILON)
			A(p_ijk,p_ijk) -= 1.0;
		if(i-1 < 0 || mS(i-1,j,k) > EPSILON)
			A(p_ijk,p_ijk) -= 1.0;
		if(j+1 >= theDim[1] || mS(i,j+1,k) > EPSILON)
			A(p_ijk,p_ijk) -= 1.0;
		if(j-1 < 0 || mS(i,j-1,k) > EPSILON)
			A(p_ijk,p_ijk) -= 1.0;
		if(k+1 >= theDim[2] || mS(i,j,k+1) > EPSILON)
			A(p_ijk,p_ijk) -= 1.0;
		if(k-1 < 0 || mS(i,j,k-1) > EPSILON)
			A(p_ijk,p_ijk) -= 1.0;

		//Update Other Coefficients
		if(i+1 < theDim[0] && mS(i+1,j,k) < EPSILON)
			A(p_ijk,p_i1jk) = -1.0;
		if(j+1 < theDim[1] && mS(i,j+1,k) < EPSILON)
			A(p_ijk,p_ij1k) = -1.0;
		if(k+1 < theDim[2] && mS(i,j,k+1) < EPSILON) 
			A(p_ijk,p_ijk1) = -1.0;
		if(i-1 >= 0 && mS(i-1,j,k) < EPSILON)
			A(p_ijk,p_n_i1jk) = -1.0;
		if(j-1 >= 0 && mS(i,j-1,k) < EPSILON)
			A(p_ijk,p_n_ij1k) = -1.0;
		if(k-1 >= 0 && mS(i,j,k-1) < EPSILON)
			A(p_ijk,p_n_ijk1) = -1.0;
	}

	
	// 3. Solve for p
	//Output and Check Matrix A
	//std::ofstream outFile("a_matrix.txt");
	//for(int i=0;i<A.size1();i++){
	//	for(int j=0;j<A.size2();j++){
	//		outFile << A(i,j) << " ";
	//	}
	//	outFile << std::endl;
	//}
	//outFile.close();

	//Compute MIC(0) Preconditioner
	double tau = 0.97; //Tuning Constant
	ublas::vector<double> E_inv(theDim[0]*theDim[1]*theDim[2]);
	E_inv.clear();

	FOR_EACH_CELL {
		if(mS(i,j,k) < EPSILON){
			//Indices
			int p_ijk = getIndex(i,j,k);
			int p_i1jk =  getIndex(i+1,j,k);
			int p_ij1k =  getIndex(i,j+1,k);
			int p_ijk1 =  getIndex(i,j,k+1);
			int p_n_i1jk = getIndex(i-1,j,k);
			int p_n_ij1k =  getIndex(i,j-1,k);
			int p_n_ijk1 = getIndex(i,j,k-1);

			//Values (Only If Indices Work Out Though)
			double a_ijk_ijk= A(p_ijk,p_ijk);
		    double a_in1jk_ijk=0.0;
			double e_in1jk=0.0;
			if(i-1 >= 0){
				a_in1jk_ijk=A(p_n_i1jk,p_ijk);
				e_in1jk=E_inv(p_n_i1jk);
			}
			double a_ijn1k_ijk=0.0;
			double e_ijn1k=0.0;
			if(j-1 >= 0){
				a_ijn1k_ijk=A(p_n_ij1k,p_ijk);
				e_ijn1k=E_inv(p_n_ij1k);
			}
			double a_ijkn1_ijk=0.0;
			double e_ijkn1=0.0;
			if(k-1 >= 0){
				a_ijkn1_ijk=A(p_n_ijk1,p_ijk);
				e_ijkn1=E_inv(p_n_ijk1);
			}

			//Modified Incomplete Cholesky Additions
			double A_in1jk_in1j1k = 0.0;
			double A_in1jk_in1jk1 = 0.0;
			double A_ijn1k_i1jn1k = 0.0;
			double A_ijn1k_ijn1k1 = 0.0;
			double A_ijkn1_i1jkn1 = 0.0;
			double A_ijkn1_ij1kn1 = 0.0;
			if(i-1 >= 0 && j+1 < theDim[1])
				A_in1jk_in1j1k = A(p_n_i1jk,getIndex(i-1,j+1,k));
			if(i-1 >= 0 && k+1 < theDim[2])
				A_in1jk_in1jk1 = A(p_n_i1jk,getIndex(i-1,j,k+1));
			if(j-1 >= 0 && i+1 < theDim[0])
				A_ijn1k_i1jn1k = A(p_n_ij1k,getIndex(i+1,j-1,k));
			if(j-1 >= 0 && k+1 < theDim[2])
				A_ijn1k_ijn1k1 = A(p_n_ij1k,getIndex(i,j-1,k+1));
			if(k-1 >= 0 && i+1 < theDim[0])
				A_ijkn1_i1jkn1 = A(p_n_ijk1,getIndex(i+1,j,k-1));
			if(k-1 >= 0 && j+1 < theDim[1])
				A_ijkn1_ij1kn1 = A(p_n_ijk1,getIndex(i,j+1,k-1));

			double e = a_ijk_ijk - (a_in1jk_ijk*e_in1jk)*(a_in1jk_ijk*e_in1jk) - 
									(a_ijn1k_ijk*e_ijn1k)*(a_ijn1k_ijk*e_ijn1k)-
									(a_ijkn1_ijk*e_ijkn1)*(a_ijkn1_ijk*e_ijkn1);
			
			e = e - tau*(a_in1jk_ijk*(A_in1jk_in1j1k+A_in1jk_in1jk1)*e_in1jk*e_in1jk +
						 a_ijn1k_ijk*(A_ijn1k_i1jn1k+A_ijn1k_ijn1k1)*e_ijn1k*e_ijn1k +
						 a_ijkn1_ijk*(A_ijkn1_i1jkn1+A_ijkn1_ij1kn1)*e_ijkn1*e_ijkn1);

			E_inv(p_ijk) = 1.0/sqrt(e + 0.000000001); 
		}
	}


	//Un-Preconditioned
	//if(!cg_solve(A,d.data(),mP.data(),200,EPSILON))
	//	mP.data().clear();
	//Preconditioned
	if(!cg_preconditioned_solve(A,d.data(),E_inv,mP.data(),200,EPSILON))
		mP.data().clear();

	//std::ofstream outPressure("pressure_matrix.txt");
	//FOR_EACH_CELL
	//	outPressure << mP(i,j,k) << std::endl;
	//outPressure.close();

	// Subtract pressure from our velocity and save in target
	FOR_EACH_FACE
	{
		if (isFace(i,j,k,X)){
			if(i==0 || i==theDim[0] || mS(i,j,k) > EPSILON || mS(i-1,j,k) > EPSILON)
				target.mU(i,j,k) = 0.0;
			else
				target.mU(i,j,k) = mU(i,j,k)-(dt/rho)*((mP(i,j,k)-mP(i-1,j,k))/theCellSize);
		}
		if (isFace(i,j,k,Y)){
			if(j==0 || j==theDim[1] || mS(i,j,k) > EPSILON || mS(i,j-1,k) > EPSILON)
				target.mV(i,j,k) = 0.0;
			else
				target.mV(i,j,k) = mV(i,j,k)-(dt/rho)*((mP(i,j,k)-mP(i,j-1,k))/theCellSize);
		}
		if (isFace(i,j,k,Z)){
			if(k==0 || k==theDim[2] || mS(i,j,k) > EPSILON || mS(i,j,k-1) > EPSILON)
				target.mW(i,j,k) = 0.0;
			else
				target.mW(i,j,k) = mW(i,j,k)-(dt/rho)*((mP(i,j,k)-mP(i,j,k-1))/theCellSize);
		}
	}
   
   // Then save the result to our object
   mU = target.mU;
   mV = target.mV;
   mW = target.mW;

   // IMPLEMENT THIS IS A SANITY CHECK: assert (checkDivergence());
   checkDivergence();
}

bool MACGrid::checkDivergence() {
	//Check Divergence of Each Cell
	//std::ofstream outU("div_out.txt");
	FOR_EACH_CELL
	{
		double divergence = (1.0/theCellSize)*((mU(i+1,j,k)-mU(i,j,k)) + (mV(i,j+1,k)-mV(i,j,k)) + (mW(i,j,k+1)-mW(i,j,k)));
		if(fabs(divergence) > EPSILON){
			//outU << divergence << std::endl;
			return false;
		}
	}
	//outU.close();
	return true;
}

vec3 MACGrid::getVelocity(const vec3& pt)
{
   vec3 vel;
   vel[0] = getVelocityX(pt); 
   vel[1] = getVelocityY(pt); 
   vel[2] = getVelocityZ(pt); 
   return vel;
}

double MACGrid::getVelocityX(const vec3& pt)
{
   return mU.interpolate(pt);
}

double MACGrid::getVelocityY(const vec3& pt)
{
   return mV.interpolate(pt);
}

double MACGrid::getVelocityZ(const vec3& pt)
{
   return mW.interpolate(pt);
}

double MACGrid::getTemperature(const vec3& pt)
{
   return mT.interpolate(pt);
}

double MACGrid::getDensity(const vec3& pt)
{
   return mD.interpolate(pt);
}

vec3 MACGrid::getCenter(int i, int j, int k)
{
   double xstart = theCellSize/2.0;
   double ystart = theCellSize/2.0;
   double zstart = theCellSize/2.0;

   double x = xstart + i*theCellSize;
   double y = ystart + j*theCellSize;
   double z = zstart + k*theCellSize;
   return vec3(x, y, z);
}

vec3 MACGrid::getZFace(int i, int j, int k)
{
   return getCenter(i,j,k)-vec3(0,0,theCellSize/2.0);
}

vec3 MACGrid::getYFace(int i, int j, int k)
{
   return getCenter(i,j,k)-vec3(0,theCellSize/2.0,0);
}

vec3 MACGrid::getXFace(int i, int j, int k)
{
   return getCenter(i,j,k)-vec3(theCellSize/2.0,0,0);
}

void MACGrid::draw(const Camera& c)
{   
   drawWireGrid();
   if (theDisplayVel) drawVelocities();   
   if (theRenderMode == CUBES) drawSmokeCubes(c);
   else drawSmoke(c);
   drawSolids(c);
}

void MACGrid::drawVelocities()
{
   // draw line at each center
   glColor4f(0.0, 1.0, 0.0, 1.0);
   glBegin(GL_LINES);
      FOR_EACH_CELL
      {
         vec3 pos = getCenter(i,j,k);
         vec3 vel = getVelocity(pos);
         if (vel.Length() > 0.0001)
         {
           vel.Normalize();
           vel *= theCellSize/2.0;
           vel += pos;
           glVertex3dv(pos.n);
           glVertex3dv(vel.n);
         }
      }
   glEnd();
}

void MACGrid::drawSolids(const Camera& c) 
{
   std::multimap<double, MACGrid::Cube, std::greater<double>> cubes;
   FOR_EACH_CELL
   {
	   if(mS(i,j,k) > EPSILON){
		   MACGrid::Cube cube;
		   cube.color = vec4(0.5,0.5,0.5,1.0);
		   cube.pos = getCenter(i,j,k);
		   cube.dist = DistanceSqr(cube.pos, c.getPosition());
		   cubes.insert(make_pair(cube.dist, cube));
	   }
   } 

   // Draw Solid cubes from back to front
   std::multimap<double, MACGrid::Cube, std::greater<double>>::const_iterator it;
   for (it = cubes.begin(); it != cubes.end(); ++it)
   {
      drawCube(it->second);
   }
}

vec4 MACGrid::getRenderColor(int i, int j, int k)
{
    double value = mD(i, j, k); 
    return vec4(1.0, 0.9, 1.0, value);
}

vec4 MACGrid::getRenderColor(const vec3& pt)
{
    double value = getDensity(pt); 
    return vec4(1.0, 1.0, 1.0, value);
}

void MACGrid::drawZSheets(bool backToFront)
{
   // Draw K Sheets from back to front
   double back =  (theDim[2])*theCellSize;
   double top  =  (theDim[1])*theCellSize;
   double right = (theDim[0])*theCellSize;
  
   double stepsize = theCellSize*0.25;

   double startk = back - stepsize;
   double endk = 0;
   double stepk = -theCellSize;

   if (!backToFront)
   {
      startk = 0;
      endk = back;   
      stepk = theCellSize;
   }

   for (double k = startk; backToFront? k > endk : k < endk; k += stepk)
   {
     for (double j = 0.0; j < top; )
      {
         glBegin(GL_QUAD_STRIP);
         for (double i = 0.0; i <= right; i += stepsize)
         {
            vec3 pos1 = vec3(i,j,k); 
            vec3 pos2 = vec3(i, j+stepsize, k); 

            vec4 color1 = getRenderColor(pos1);
            vec4 color2 = getRenderColor(pos2);

            glColor4dv(color1.n);
            glVertex3dv(pos1.n);

            glColor4dv(color2.n);
            glVertex3dv(pos2.n);
         } 
         glEnd();
         j+=stepsize;

         glBegin(GL_QUAD_STRIP);
         for (double i = right; i >= 0.0; i -= stepsize)
         {
            vec3 pos1 = vec3(i,j,k); 
            vec3 pos2 = vec3(i, j+stepsize, k); 

            vec4 color1 = getRenderColor(pos1);
            vec4 color2 = getRenderColor(pos2);

            glColor4dv(color1.n);
            glVertex3dv(pos1.n);

            glColor4dv(color2.n);
            glVertex3dv(pos2.n);
         } 
         glEnd();
         j+=stepsize;
      }
   }
}

void MACGrid::drawXSheets(bool backToFront)
{
   // Draw K Sheets from back to front
   double back =  (theDim[2])*theCellSize;
   double top  =  (theDim[1])*theCellSize;
   double right = (theDim[0])*theCellSize;
  
   double stepsize = theCellSize*0.25;

   double starti = right - stepsize;
   double endi = 0;
   double stepi = -theCellSize;

   if (!backToFront)
   {
      starti = 0;
      endi = right;   
      stepi = theCellSize;
   }

   for (double i = starti; backToFront? i > endi : i < endi; i += stepi)
   {
     for (double j = 0.0; j < top; )
      {
         glBegin(GL_QUAD_STRIP);
         for (double k = 0.0; k <= back; k += stepsize)
         {
            vec3 pos1 = vec3(i,j,k); 
            vec3 pos2 = vec3(i, j+stepsize, k); 

            vec4 color1 = getRenderColor(pos1);
            vec4 color2 = getRenderColor(pos2);

            glColor4dv(color1.n);
            glVertex3dv(pos1.n);

            glColor4dv(color2.n);
            glVertex3dv(pos2.n);
         } 
         glEnd();
         j+=stepsize;

         glBegin(GL_QUAD_STRIP);
         for (double k = back; k >= 0.0; k -= stepsize)
         {
            vec3 pos1 = vec3(i,j,k); 
            vec3 pos2 = vec3(i, j+stepsize, k); 

            vec4 color1 = getRenderColor(pos1);
            vec4 color2 = getRenderColor(pos2);

            glColor4dv(color1.n);
            glVertex3dv(pos1.n);

            glColor4dv(color2.n);
            glVertex3dv(pos2.n);
         } 
         glEnd();
         j+=stepsize;
      }
   }
}


void MACGrid::drawSmoke(const Camera& c)
{
   vec3 eyeDir = c.getBackward();
   double zresult = fabs(Dot(eyeDir, vec3(1,0,0)));
   double xresult = fabs(Dot(eyeDir, vec3(0,0,1)));
   double yresult = fabs(Dot(eyeDir, vec3(0,1,0)));

   if (zresult < xresult)
   {      
      drawZSheets(c.getPosition()[2] < 0);
   }
   else 
   {
      drawXSheets(c.getPosition()[0] < 0);
   }
}

void MACGrid::drawSmokeCubes(const Camera& c)
{
   std::multimap<double, MACGrid::Cube, std::greater<double>> cubes;
   FOR_EACH_CELL
   {
      MACGrid::Cube cube;
      cube.color = getRenderColor(i,j,k);
      cube.pos = getCenter(i,j,k);
      cube.dist = DistanceSqr(cube.pos, c.getPosition());
      cubes.insert(make_pair(cube.dist, cube));
   } 

   // Draw cubes from back to front
   std::multimap<double, MACGrid::Cube, std::greater<double>>::const_iterator it;
   for (it = cubes.begin(); it != cubes.end(); ++it)
   {
      drawCube(it->second);
   }
}

void MACGrid::drawWireGrid()
{
   // Display grid in light grey, draw top & bottom

   double xstart = 0.0;
   double ystart = 0.0;
   double zstart = 0.0;
   double xend = theDim[0]*theCellSize;
   double yend = theDim[1]*theCellSize;
   double zend = theDim[2]*theCellSize;

   glPushAttrib(GL_LIGHTING_BIT | GL_LINE_BIT);
      glDisable(GL_LIGHTING);
      glColor3f(0.25, 0.25, 0.25);

      glBegin(GL_LINES);
      for (int i = 0; i <= theDim[0]; i++)
      {
         double x = xstart + i*theCellSize;
         glVertex3d(x, ystart, zstart);
         glVertex3d(x, ystart, zend);

         glVertex3d(x, yend, zstart);
         glVertex3d(x, yend, zend);
      }

      for (int i = 0; i <= theDim[2]; i++)
      {
         double z = zstart + i*theCellSize;
         glVertex3d(xstart, ystart, z);
         glVertex3d(xend, ystart, z);

         glVertex3d(xstart, yend, z);
         glVertex3d(xend, yend, z);
      }

      glVertex3d(xstart, ystart, zstart);
      glVertex3d(xstart, yend, zstart);

      glVertex3d(xend, ystart, zstart);
      glVertex3d(xend, yend, zstart);

      glVertex3d(xstart, ystart, zend);
      glVertex3d(xstart, yend, zend);

      glVertex3d(xend, ystart, zend);
      glVertex3d(xend, yend, zend);
      glEnd();
   glPopAttrib();

   glEnd();
}

#define LEN 0.5
void MACGrid::drawFace(const MACGrid::Cube& cube)
{
   glColor4dv(cube.color.n);
   glPushMatrix();
      glTranslated(cube.pos[0], cube.pos[1], cube.pos[2]);      
      glScaled(theCellSize, theCellSize, theCellSize);
      glBegin(GL_QUADS);
         glNormal3d( 0.0,  0.0, 1.0);
         glVertex3d(-LEN, -LEN, LEN);
         glVertex3d(-LEN,  LEN, LEN);
         glVertex3d( LEN,  LEN, LEN);
         glVertex3d( LEN, -LEN, LEN);
      glEnd();
   glPopMatrix();
}

void MACGrid::drawCube(const MACGrid::Cube& cube)
{
   glColor4dv(cube.color.n);
   glPushMatrix();
      glTranslated(cube.pos[0], cube.pos[1], cube.pos[2]);      
      glScaled(theCellSize, theCellSize, theCellSize);
      glBegin(GL_QUADS);
         glNormal3d( 0.0, -1.0,  0.0);
         glVertex3d(-LEN, -LEN, -LEN);
         glVertex3d(-LEN, -LEN,  LEN);
         glVertex3d( LEN, -LEN,  LEN);
         glVertex3d( LEN, -LEN, -LEN);         

         glNormal3d( 0.0,  0.0, -0.0);
         glVertex3d(-LEN, -LEN, -LEN);
         glVertex3d(-LEN,  LEN, -LEN);
         glVertex3d( LEN,  LEN, -LEN);
         glVertex3d( LEN, -LEN, -LEN);

         glNormal3d(-1.0,  0.0,  0.0);
         glVertex3d(-LEN, -LEN, -LEN);
         glVertex3d(-LEN, -LEN,  LEN);
         glVertex3d(-LEN,  LEN,  LEN);
         glVertex3d(-LEN,  LEN, -LEN);

         glNormal3d( 0.0, 1.0,  0.0);
         glVertex3d(-LEN, LEN, -LEN);
         glVertex3d(-LEN, LEN,  LEN);
         glVertex3d( LEN, LEN,  LEN);
         glVertex3d( LEN, LEN, -LEN);

         glNormal3d( 0.0,  0.0, 1.0);
         glVertex3d(-LEN, -LEN, LEN);
         glVertex3d(-LEN,  LEN, LEN);
         glVertex3d( LEN,  LEN, LEN);
         glVertex3d( LEN, -LEN, LEN);

         glNormal3d(1.0,  0.0,  0.0);
         glVertex3d(LEN, -LEN, -LEN);
         glVertex3d(LEN, -LEN,  LEN);
         glVertex3d(LEN,  LEN,  LEN);
         glVertex3d(LEN,  LEN, -LEN);
      glEnd();
   glPopMatrix();
}