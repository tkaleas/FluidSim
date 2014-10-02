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

#define GRAVITY -9.81 * 2.0
#define FLUID 0.0
#define AIR -1.0
#define SOLID 1.0

// Globals
MACGrid target;
extern int theDim[3];
extern double theCellSize;

// NOTE: x -> cols, z -> rows, y -> stacks
MACGrid::RenderMode MACGrid::theRenderMode = SPHERES;
bool MACGrid::theDisplayVel = false;
bool MACGrid::theDisplayFluReg = false;
bool MACGrid::theDisplayLevSet = false;
bool MACGrid::theDisplaySol = false;
bool MACGrid::theDisplayMarkMesh = false;


#define FOR_EACH_CELL \
   for(int k = 0; k < theDim[MACGrid::Z]; k++)  \
      for(int j = 0; j < theDim[MACGrid::Y]; j++) \
         for(int i = 0; i < theDim[MACGrid::X]; i++) 

#define FOR_EACH_FACE \
   for(int k = 0; k < theDim[MACGrid::Z]+1; k++) \
      for(int j = 0; j < theDim[MACGrid::Y]+1; j++) \
         for(int i = 0; i < theDim[MACGrid::X]+1; i++) 

#ifdef _DEBUGA

#endif

MarkerParticle::MarkerParticle(vec3 p){
	position = p;
	sign = -1;
}

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

   return *this;
}

MACGrid::~MACGrid()
{
}

void MACGrid::reset()
{
   //Intialize Marching Cubes
 mMarchCube = new MarchCube();
 mSurfaceLS = new IsoSurface(&mL);
 mSurfaceM = new IsoSurface(this); 
 //Initialize GridData
   mU.initialize();
   mV.initialize();
   mW.initialize();
   mP.initialize();
   mS.initialize(AIR);
   mPG_u.initialize();
   mPG_v.initialize();
   mPG_w.initialize();
   mL.initialize();

   //Clear Particles
   mParticles.clear();

   initializeSources();

   //Add Marker Particles Based on Initial Fluid Cells
   FOR_EACH_CELL 
   {
	   //If Fluid Cell mS = 0.0
	   if(mS(i,j,k) < EPSILON && mS(i,j,k) > N_EPSILON){
		   //Create 8 Particles in Random Places within the Cell
		   for(int ii=0; ii < 32;ii++){
			   double x_offset  = (theCellSize/2.0)*(double)((rand() % 180)-90)/100.0;
			   double y_offset  =(theCellSize/2.0)*(double)((rand() % 180)-90)/100.0;
			   double z_offset  = (theCellSize/2.0)*(double)((rand() % 180)-90)/100.0;
			   mParticles.push_back(MarkerParticle(getCenter(i,j,k)+vec3(x_offset,y_offset,z_offset)));
			   int x;
			   int y;
			   int z;
			   //mS.getCell(mParticles.back().position,x,y,z);
			  // std::cout << mParticles.back().position << std::endl;
			   //std::cout <<"IS:" << x << ", " << y << ", "<< z << std::endl;
			   //std::cout << "SHOULE BE: " << i << ", " << j << ", "<< k << std::endl;
		   }
	   }
   }

   initializeLevelSet();
}

void MACGrid::initialize()
{
   reset();
}

void MACGrid::updateSources()
{
    //TODO: Set initial values for velocity
	//mV(1,1,0) = 0.1;
}

void MACGrid::initializeSources()
{
	//TODO: Set Initial Fluid Cells
	//for(int i=0; i < theDim[0];i++){
	//	mS(i,1,0) = FLUID;
	//}
	#ifdef _DEBUG
	mS(1,1,0) = FLUID;
	#else
	for(int i=1; i < 8;i++){
		for(int j=4; j < 10;j++)
			mS(i,j,0) = FLUID;
	}

	//for(int i=9; i < 14;i++){
	//	for(int j=9; j < 14;j++){
	//		for(int k=3; k < 8;k++)
	//			mS(i,j,k) = FLUID;
	//	}
	//}
	//for(int i=40; i < 50;i++){
	//	for(int j=40; j < 50;j++)
	//		mS(i,j,0) = FLUID;
	//}
	//
	//for(int i=0; i < 60;i++){
	//	for(int j=0; j < 20;j++)
	//		mS(i,j,0) = FLUID;
	//}

	#endif
	mV(1,1,0) = 0.1;
}

void MACGrid::initializeLevelSet(){
	double particleRadius = 0.12f;
	//std::ofstream outLS("cell_init.txt");
	FOR_LS 
		if(mS(i-1,j-1,k-1) >= N_EPSILON){
			//Fluid Cells (Set to -cellSize/2.0 since all that is important is the surface on the outside here)
			mL(i,j,k) = -theCellSize;
		} 
		else if (mS(i-1,j-1,k-1) < N_EPSILON) {
			//Air Cells (Distance to the Cloest Particle is the Surface)
			double closest = std::numeric_limits<double>::infinity();
			
			for(int ii=0; ii < mParticles.size(); ii++)
			{
				//Distance from Grid Center -> Particle
				vec3 to_particle = getLSCenter(i,j,k) - mParticles[ii].position;
				double phi_temp = to_particle.Length()-particleRadius;
				if(phi_temp < closest){ closest = phi_temp;}
			}

			//outLS << "GRID POINT:" << i << "," << j << "," << k << std::endl;
			//outLS << closest << std::endl;
			mL(i,j,k) = closest;
		}
	END_FOR_THREE
	
	mL.ReInitialize();
	mL.SetBoundarySignedDist();
	//outLS.close();
}

void MACGrid::updateGrid(){
	//Update Grid based on Particle Locations
	mS.initialize(AIR);
	
	//std::ofstream out("Cell_Grid_Update.txt");
	//Update Fluid Cells -> Fluid
	for(int i=0; i < mParticles.size();i++){
		int x;
		int y;
		int z;
		mS.getCell(mParticles[i].position,x,y,z);
		mS(x,y,z) = FLUID;
		//out << x << ", " << y << ", "<< z << std::endl;
	}
	//out.close();

	//Update Fluid Cell Regions from Level Set
	//FOR_EACH_CELL {
	//	if(mL(i+1,j+1,k+1) < 0)
	//		mS(i,j,k) = FLUID;
	//}
}

void MACGrid::advectVelocity(double dt)
{
    // TODO: Calculate new velocities and store in target

	//Advect Velocity on Each Face (RK2 Integration)
	//Euler y = x + dt*u(x)
	//RK2 y = x + dt*u(x+(dt/2)*u(x))
	FOR_EACH_FACE
	{
		if (isFace(i,j,k,X)){
			if((i==0 || i==theDim[0]) || mS(i,j,k) > EPSILON || mS(i-1,j,k) > EPSILON) {target.mU(i,j,k)=0.0;}
			else{
				//vec3 traceBack1 = getXFace(i,j,k) - (dt/2.0)*getVelocity(getXFace(i,j,k));
				//traceBack2 = clampPoint(traceBack2); // Clamp point to outside of solid cell
				//vec3 traceBack2 = getXFace(i,j,k) - dt*getVelocity(traceback1);
				//traceBack2 = clampPoint(traceBack2); // Clamp point to outside of solid cell
				//target.mU(i,j,k) = getVelocityX(traceBack2);

				vec3 traceBack = getXFace(i,j,k) - dt*getVelocity(getXFace(i,j,k));
				traceBack = clampPoint(traceBack); // Clamp point to outside of solid cell
				target.mU(i,j,k) = getVelocityX(traceBack);
			}
		}
		if (isFace(i,j,k,Y)){
			if((j==0 || j==theDim[1])|| mS(i,j,k) > EPSILON || mS(i,j-1,k) > EPSILON) {target.mV(i,j,k)=0.0;}
			else{
				//vec3 traceBack1 = getYFace(i,j,k) - (dt/2.0)*getVelocity(getYFace(i,j,k));
				//traceBack2 = clampPoint(traceBack2); // Clamp point to outside of solid cell
				//vec3 traceBack2 = getYFace(i,j,k) - dt*getVelocity(traceback1);
				//traceBack2 = clampPoint(traceBack2); // Clamp point to outside of solid cell
				//target.mV(i,j,k) = getVelocityY(traceBack2);

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

   double epsilon = 1.0;
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

   // Then save the result to our object
   mU = target.mU;
   mV = target.mV;
   mW = target.mW;
}

void MACGrid::computeGravity(double dt)
{
   // TODO: Calculate Gravity and store in target
   FOR_EACH_FACE{
	   if(isFace(i,j,k,Y)){
		   if(j==0 || j==theDim[1] || mS(i,j,k) > EPSILON || mS(i,j-1,k) > EPSILON) {target.mV(i,j,k) = 0.0;}
		   else if(mS(i,j,k) > N_EPSILON || mS(i,j-1,k) > N_EPSILON)
			   target.mV(i,j,k) = mV(i,j,k)+dt*GRAVITY;
	   }
   }

   // and then save the result to our object
   mV = target.mV;
}

void MACGrid::addExternalForces(double dt)
{
   computeGravity(dt);
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
	double rho = 1000.0;
	
	//Calculate Divergence
	std::ofstream outPreDiv("pre_div.txt");
	FOR_EACH_CELL
	{
		//FLUID
		if(mS(i,j,k) > N_EPSILON){
			d(i,j,k) = -1*((rho*theCellSize)/dt)*((mU(i+1,j,k)-mU(i,j,k)) + (mV(i,j+1,k)-mV(i,j,k)) + (mW(i,j,k+1)-mW(i,j,k))); //Discrete Divergence
			//outPreDiv << (1.0/theCellSize)*((mU(i+1,j,k)-mU(i,j,k)) + (mV(i,j+1,k)-mV(i,j,k)) + (mW(i,j,k+1)-mW(i,j,k))) << std::endl;
			outPreDiv << -1*((rho*theCellSize)/dt)*((mU(i+1,j,k)-mU(i,j,k)) + (mV(i,j+1,k)-mV(i,j,k)) + (mW(i,j,k+1)-mW(i,j,k))) << std::endl;
		}
		else {
			outPreDiv << 0 << std::endl;
		}
	}
	outPreDiv.close();

	// 2. Construct A
	ublas::matrix<double> A(theDim[0]*theDim[1]*theDim[2],theDim[0]*theDim[1]*theDim[2]);

	//Initialize Matrix to 0's
	A.clear();

	//Determine Row Coefficients
	FOR_EACH_CELL {
		//If Solid or Air Cell -> Ignore and Go to Next Cell
		if(mS(i,j,k) > EPSILON || mS(i,j,k) < N_EPSILON)
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
		if(i+1 < theDim[0] && mS(i+1,j,k) < EPSILON && mS(i+1,j,k) > N_EPSILON)
			A(p_ijk,p_i1jk) = -1.0;
		if(j+1 < theDim[1] && mS(i,j+1,k) < EPSILON && mS(i,j+1,k) > N_EPSILON)
			A(p_ijk,p_ij1k) = -1.0;
		if(k+1 < theDim[2] && mS(i,j,k+1) < EPSILON && mS(i,j,k+1) > N_EPSILON) 
			A(p_ijk,p_ijk1) = -1.0;
		if(i-1 >= 0 && mS(i-1,j,k) < EPSILON && mS(i-1,j,k) > N_EPSILON)
			A(p_ijk,p_n_i1jk) = -1.0;
		if(j-1 >= 0 && mS(i,j-1,k) < EPSILON && mS(i,j-1,k) > N_EPSILON)
			A(p_ijk,p_n_ij1k) = -1.0;
		if(k-1 >= 0 && mS(i,j,k-1) < EPSILON && mS(i,j,k-1) > N_EPSILON)
			A(p_ijk,p_n_ijk1) = -1.0;
	}

	
	// 3. Solve for p
	//Output and Check Matrix A
	std::ofstream outFile("a_matrix.txt");
	for(int i=0;i<A.size1();i++){
		for(int j=0;j<A.size2();j++){
			outFile << A(i,j) << " ";
		}
		outFile << std::endl;
	}
	outFile.close();

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

	std::ofstream outPressure("pressure_matrix.txt");
	FOR_EACH_CELL
		outPressure << mP(i,j,k) << std::endl;
	outPressure.close();

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
	//Check Divergence of Each Fluid Cell
	std::ofstream outU("div_out.txt");
	FOR_EACH_CELL
	{
		if(mS(i,j,k) > N_EPSILON){
			double divergence = (1.0/theCellSize)*((mU(i+1,j,k)-mU(i,j,k)) + (mV(i,j+1,k)-mV(i,j,k)) + (mW(i,j,k+1)-mW(i,j,k)));
			if(fabs(divergence) > EPSILON){	
				outU << "INCORRECT DIVERGENCE" << std::endl;
			}
			outU << divergence << std::endl;
		}
	}
	outU.close();
	return true;
}

vec3 MACGrid::getPhiGradient(int i,int j,int k)
{
	double phi_x,phi_y,phi_z;
	vec3 vel = getVelocity(getLSCenter(i,j,k));
	if(vel[0] < 0.0)
		phi_x = (mL(i+1,j,k) - mL(i,j,k))/theCellSize;
	else
		phi_x = (mL(i,j,k)-mL(i-1,j,k))/theCellSize;
	if(vel[1] < 0.0)
		phi_y = (mL(i,j+1,k) - mL(i,j,k))/theCellSize;
	else
		phi_y = (mL(i,j,k)-mL(i,j-1,k))/theCellSize;
	if(vel[2] < 0.0)
		phi_z = (mL(i,j,k+1) - mL(i,j,k))/theCellSize;
	else
		phi_z = (mL(i,j,k)-mL(i,j,k-1))/theCellSize;
	vec3 phi_gradient(phi_x,phi_y,phi_z);
	phi_gradient.Normalize();
	return phi_gradient;
}

void MACGrid::extrapolateVelocities(double dt)
{
	//Update PhiGradient Temporary Grid
	FOR_EACH_CELL {
		vec3 phi_grad = getPhiGradient(i,j,k);
		mPG_u(i,j,k) = phi_grad[0];
		mPG_v(i,j,k) = phi_grad[1];
		mPG_w(i,j,k) = phi_grad[2];
	}

	//Extrapolate Velocity Based on the Gradient Level Set (Closet Point on Fluid Surface)
	FOR_EACH_FACE {
		if (isFace(i,j,k,X)){
			if((i==0 || i==theDim[0])) {mU(i,j,k)=0.0;}
			else if(mS(i,j,k) < N_EPSILON && mS(i-1,j,k) < N_EPSILON && getSignedDistance(getXFace(i,j,k)) < theCellSize * 6.0) {
				vec3 x_face = getXFace(i,j,k);
				mU(i,j,k) = getVelocityX(x_face + -getSignedDistance(x_face)*vec3(mPG_u.interpolate(x_face),mPG_v.interpolate(x_face),mPG_w.interpolate(x_face)).Normalize());
			}
		}
		if (isFace(i,j,k,Y)){
			if((j==0 || j==theDim[1])) {mV(i,j,k)=0.0;}
			else if(mS(i,j,k) < N_EPSILON && mS(i,j-1,k) < N_EPSILON && getSignedDistance(getYFace(i,j,k)) < theCellSize * 6.0) {
				vec3 y_face = getYFace(i,j,k);
				mV(i,j,k) = getVelocityY(y_face + -getSignedDistance(y_face)*vec3(mPG_u.interpolate(y_face),mPG_v.interpolate(y_face),mPG_w.interpolate(y_face)).Normalize());
			}
		}
		if (isFace(i,j,k,Z)){
			if((k==0 || k==theDim[2])) {mW(i,j,k)=0.0;}
			else if(mS(i,j,k) < N_EPSILON && mS(i,j,k-1) < N_EPSILON && getSignedDistance(getZFace(i,j,k)) < theCellSize * 6.0) {
				vec3 z_face = getZFace(i,j,k);
				mW(i,j,k) = getVelocityZ(z_face + -getSignedDistance(z_face)*vec3(mPG_u.interpolate(z_face),mPG_v.interpolate(z_face),mPG_w.interpolate(z_face)).Normalize());
			}
		}
	}
}

void MACGrid::fixLS(const MarkerParticle &p, int i, int j, int k)
{
	//Passing In Particle, LevelSet Index (offset by +1 for ijk
	static double particlePhi;
	for(int dx = -1; dx < 2; dx++) {
		for(int dy = -1; dy < 2; dy++) {
            for(int dz = -1; dz < 2; dz++) {
				vec3 cell_center = getCenter(i+dx-1,j+dy-1,k+dz-1);
				particlePhi = -1*(0.1-(cell_center-p.position).Length());
				mL(i+dx,j+dy,k+dz) = min(mL(i+dx,j+dy,k+dz),particlePhi);
			}
		}
	}

}

void MACGrid::errorCorrectionLS(double dt)
{
	for(int i=0;i < mParticles.size();i++){
		double phi = getSignedDistance(mParticles[i].position);
		double particlePhi;
		if(phi > 0.0) {
			//Particle Crossed the Level Set Boundary
			mParticles[i].sign = 1;
			vec3 pos = mParticles[i].position;
			//Pass In Level Set Index
			fixLS(mParticles[i],int(pos[0]/theCellSize)+1,int(pos[1]/theCellSize)+1,int(pos[2]/theCellSize)+1);

		} else {
			mParticles[i].sign = -1;
		}
	}
	mL.SetBoundarySignedDist();
}

void MACGrid::updateLevelSet(double dt)
{
	//Update the Level Set to Reflect the True Fluid Boundary
	advectLevelSet(dt);
	errorCorrectionLS(dt);
	//ReInitialize Level Set Values
	mL.ReInitialize();
	errorCorrectionLS(dt);
}

void MACGrid::advectLevelSet(double dt)
{
	// TODO: Calculate new level set and store in target
	for(int k = 1; k <= theDim[MACGrid::Z]; k++)  {
		for(int j = 1; j <= theDim[MACGrid::Y]; j++) {
			for(int i = 1; i <= theDim[MACGrid::X]; i++)  {

				//Upwind Differencing Method (In Level Set Primer Notes)
				double phi_x,phi_y,phi_z;
				vec3 vel = getVelocity(getLSCenter(i,j,k));

				if(vel[0] < 0.0)
					phi_x = (mL(i+1,j,k) - mL(i,j,k))/theCellSize;
				else
					phi_x = (mL(i,j,k)-mL(i-1,j,k))/theCellSize;
				
				if(vel[1] < 0.0)
					phi_y = (mL(i,j+1,k) - mL(i,j,k))/theCellSize;
				else
					phi_y = (mL(i,j,k)-mL(i,j-1,k))/theCellSize;
				
				if(vel[2] < 0.0)
					phi_z = (mL(i,j,k+1) - mL(i,j,k))/theCellSize;
				else
					phi_z = (mL(i,j,k)-mL(i,j,k-1))/theCellSize;

				vec3 phi_grad(phi_x,phi_y,phi_z);

				//phi_n1 = phi_n - dt*(V * phi_grad)
				target.mL(i,j,k) = mL(i,j,k) - dt*Dot(vel,phi_grad);

				////Semi-Langrangian v1
				//target.mL(i,j,k) = getSignedDistance(getLSCenter(i,j,k) - dt*getVelocity(getLSCenter(i,j,k)));

				////Semi-Lagrangian Phi Advection v2
				//    int r,s,t;
				//    double a,b,c;
				//	vec3 u = getVelocity(getLSCenter(i,j,k));

				//	r = i - int(ceil((u[0] * dt) /theCellSize));
				//	s = j - int(ceil((u[1] * dt) /theCellSize));
				//	t = k - int(ceil((u[2] * dt) /theCellSize));

				//	//doesn't get from boundary
				//	r = Clamp(r,1,theDim[0]+1);
				//	s = Clamp(s,1,theDim[1]+1);
				//	t = Clamp(t,1,theDim[2]+1);

				//	a = (double(i - r) * theCellSize - u[0] * dt) / theCellSize;
				//	b = (double(j - s) * theCellSize - u[1] * dt) / theCellSize;
				//	c = (double(k - t) * theCellSize - u[2] * dt) / theCellSize;

				//	target.mL(i,j,k) =    a  *    b  *    c  * mL(r+1, s+1, t+1) +
				//					 (1-a) *    b  *    c  * mL(r  , s+1, t+1) +
				//						a  * (1-b) *    c  * mL(r+1, s  , t+1) +
				//						a  *    b  * (1-c) * mL(r+1, s+1, t  ) +
				//					 (1-a) * (1-b) *    c  * mL(r  , s  , t+1) +
				//					 (1-a) *    b  * (1-c) * mL(r  , s+1, t  ) +
				//						a  * (1-b) * (1-c) * mL(r+1, s  , t  ) +
				//					 (1-a) * (1-b) * (1-c) * mL(r  , s  , t  );
			}
		}
	}

	// Then save the result to our object
    mL = target.mL;
	mL.SetBoundarySignedDist();
}



void MACGrid::moveParticles(double dt)
{

	for(int i=0;i<mParticles.size();i++)
	{
		//RK2 Integration
		vec3 half_step_p = mParticles[i].position + (dt/2.0)*getVelocity(mParticles[i].position);
		mParticles[i].position = mParticles[i].position + dt*getVelocity(half_step_p);
	}

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

double MACGrid::getSignedDistance(const vec3& pt)
{
   return mL.interpolate(pt);
}


vec3 MACGrid::getLSCenter(int i, int j, int k)
{
   double xstart = -theCellSize/2.0;
   double ystart = -theCellSize/2.0;
   double zstart = -theCellSize/2.0;

   double x = xstart + i*theCellSize;
   double y = ystart + j*theCellSize;
   double z = zstart + k*theCellSize;
   return vec3(x, y, z);
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
   if (theRenderMode == BILLBOARD) drawWaterBillboards(c);
   else drawWaterParticles(c);
   if (theDisplayFluReg) drawFluidAirRegions(c);
   if (theDisplaySol) drawSolids(c);
   if (theDisplayLevSet) drawLevelSet(c);
   if (theDisplayMarkMesh) drawMarkerMesh(c);
}

Double	MACGrid::eval(const Point3d& location)
{
	//Based on Particles Directly
	double MARCH_CUBE_RADIUS_SQR = 0.027;
	vec3 loc_vec(location[0],location[1],location[2]);
	double min_dist = INFINITY;

	for(int i=0;i<mParticles.size();i++){
		double min_dist_tmp = (loc_vec-mParticles[i].position).SqrLength();
		
		if(min_dist_tmp < MARCH_CUBE_RADIUS_SQR)
			return -(loc_vec-mParticles[i].position).Length();
		
		else if(min_dist_tmp < min_dist)
			min_dist = min_dist_tmp;
	}

	return min_dist;

	////Based On Fluid Region Interpolation
	//vec3 loc_vec(location[0],location[1],location[2]);
	//return mS.interpolate(loc_vec);
	
	////Based on Fluid Region 0-1
	//int i,j,k;
	//mS.getCell(loc_vec,i,j,k);
	//if(mS(i,j,k) < N_EPSILON){
	//	return theCellSize/2.0;
	//} else {
	//	return -theCellSize/2.0;
	//}
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

void MACGrid::drawFluidAirRegions(const Camera& c) 
{
   std::multimap<double, MACGrid::Cube, std::greater<double>> cubes;
   FOR_EACH_CELL
   {
	   //Draw Fluid Regions
	   if(mS(i,j,k) < EPSILON && mS(i,j,k) > N_EPSILON){
		   MACGrid::Cube cube;
		   cube.color = vec4(1.0,0.1,0.1,0.2);
		   cube.pos = getCenter(i,j,k);
		   cube.dist = DistanceSqr(cube.pos, c.getPosition());
		   cubes.insert(make_pair(cube.dist, cube));
	   }

	   //Draw Air Regions
	   if(mS(i,j,k) < N_EPSILON){
		   MACGrid::Cube cube;
		   cube.color = vec4(0.1,1.0,0.1,0.2);
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

void MACGrid::drawMarkerMesh(const Camera& c)
{
	mMarchCube->setThreshold(0.0);
	std::cout << theCellSize*theDim[0] << std::endl;
	mMarchCube->setSize(theCellSize*(theDim[0]+2),theCellSize*(theDim[1]+2),theCellSize*(theDim[2]+2));
	mMarchCube->setRes(200,200,5);
	mMarchCube->setCenter(theCellSize*(theDim[0])/2.0,theCellSize*(theDim[1])/2.0,theCellSize*(theDim[2])/2.0);
	mMarchCube->march(*mSurfaceM);
	
	glColor4f(1.0, 0.0, 0.0, 0.5);
	mSurfaceM->glDraw(c);
}

void MACGrid::drawLevelSet(const Camera& c) 
{
 //  std::multimap<double, MACGrid::Cube, std::greater<double>> cubes;
 //  
 //  double max_dist = 3.0* theCellSize * max(theDim[0],max(theDim[1],theDim[2]));
 //  
 //  FOR_ALL_LS	   
	//	   MACGrid::Cube cube;
	//	   double opacity = mL(i,j,k)/max_dist;
	//	   if(opacity < 0)
	//		   cube.color = vec4(0.0,1.0,0.0,1.0+5.0*opacity);
	//	   else
	//		   cube.color = vec4(0.1,0.1,1.0,opacity);
	//	   cube.pos = getLSCenter(i,j,k);
	//	   cube.dist = DistanceSqr(cube.pos, c.getPosition());
	//	   cubes.insert(make_pair(cube.dist, cube));
	//END_FOR_THREE

 //  // Draw Solid cubes from back to front
 //  std::multimap<double, MACGrid::Cube, std::greater<double>>::const_iterator it;
 //  for (it = cubes.begin(); it != cubes.end(); ++it)
 //  {
 //     drawCube(it->second);
 //  }
 //  cubes.clear();

   //Implicit Surface Rendering of Marching Cubes
	
	mMarchCube->setThreshold(0);
	std::cout << theCellSize*theDim[0] << std::endl;
	mMarchCube->setSize(theCellSize*(theDim[0]+2),theCellSize*(theDim[1]+2),theCellSize*(theDim[2]+2));
	mMarchCube->setRes(100,100,100);
	mMarchCube->setCenter(theCellSize*(theDim[0])/2.0,theCellSize*(theDim[1])/2.0,theCellSize*(theDim[2])/2.0);
	mMarchCube->march(*mSurfaceLS);
	
	glColor4f(0.0, 1.0, 0.0, 0.5);
	mSurfaceLS->glDraw(c);

}

void MACGrid::drawWaterParticles(const Camera& c)
{
   std::multimap<double, MACGrid::Sphere, std::greater<double>> spheres;
   
   //Insert Particles into Sphere Map
   for(int i=0;i<mParticles.size();i++){
	    MACGrid::Sphere sphere;
		if(mParticles[i].sign < 0)
			sphere.color = vec4(0.0,0.25,0.75,1.0);
		else
			sphere.color = vec4(0.75,0.0,0.2,1.0);
		sphere.radius = 0.025;
		sphere.pos = mParticles[i].position;
		//std::cout << "MP" << std::endl;
		//std::cout << sphere.pos << std::endl;
		sphere.dist = DistanceSqr(sphere.pos, c.getPosition());
		spheres.insert(make_pair(sphere.dist,sphere));
   }

   // Draw spheres from back to front
   std::multimap<double, MACGrid::Sphere, std::greater<double>>::const_iterator it;
   for (it = spheres.begin(); it != spheres.end(); ++it)
   {
      drawSphere(it->second);
   }

}

void MACGrid::drawWaterBillboards(const Camera& c)
{
   //TODO: Change to Draw Billboards for Better Rendering Performance
   std::multimap<double, MACGrid::Sphere, std::greater<double>> spheres;
   //Insert Spheres Into the Map

   // Draw spheres from back to front
   std::multimap<double, MACGrid::Sphere, std::greater<double>>::const_iterator it;
   for (it = spheres.begin(); it != spheres.end(); ++it)
   {
      drawSphere(it->second);
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

void MACGrid::drawSphere(const MACGrid::Sphere& sphere)
{
	glColor4dv(sphere.color.n);
	glPushMatrix();
		glTranslated(sphere.pos[0], sphere.pos[1], sphere.pos[2]);      
		glutSolidSphere(sphere.radius,5,5);
	glPopMatrix();
}