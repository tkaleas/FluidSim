#ifndef MACGrid_H_
#define MACGrid_H_

#pragma warning(disable: 4244 4267 4996)
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
using namespace boost::numeric;

#include <windows.h>
#include "GL/gl.h"
#include "vec.h"
#include "GridData.h"

class Camera;
class MACGrid
{
   friend MACGrid;
public:
   MACGrid();
   ~MACGrid();
   MACGrid(const MACGrid& orig);
   MACGrid& operator=(const MACGrid& orig);

   void reset();

   void draw(const Camera& c);
   void updateSources();
   void advectVelocity(double dt);
   void addExternalForces(double dt);
   void project(double dt);
   void advectTemperature(double dt);
   void advectDensity(double dt);
   bool checkDivergence(); //Sanity Check?

protected:

   // Setup
   void initialize();

   // Simulation
   void computeBouyancy(double dt);
   void computeVorticityConfinement(double dt);

   // Rendering
   struct Cube { vec3 pos; vec4 color; double dist; };
   void drawWireGrid();
   void drawSmokeCubes(const Camera& c);
   void drawSmoke(const Camera& c);
   void drawSolids(const Camera& c);
   void drawCube(const MACGrid::Cube& c);
   void drawFace(const MACGrid::Cube& c);
   void drawVelocities();
   vec4 getRenderColor(int i, int j, int k);
   vec4 getRenderColor(const vec3& pt);
   void drawZSheets(bool backToFront);
   void drawXSheets(bool backToFront);

   // GridData accessors
   enum Direction { X, Y, Z };
   vec3 getVelocity(const vec3& pt);
   double getVelocityX(const vec3& pt);
   double getVelocityY(const vec3& pt);
   double getVelocityZ(const vec3& pt);
   double getTemperature(const vec3& pt);
   double getDensity(const vec3& pt);
   vec3 getCenter(int i, int j, int k);
   bool isFace(int i, int j, int k,Direction d);
   vec3 getXFace(int i, int j, int k);
   vec3 getYFace(int i, int j, int k);
   vec3 getZFace(int i, int j, int k);
   int getIndex(int i, int j, int k); //Helper method for getting Matrix/Vector index in ublas::vector/matrix;
   vec3 getOmega(int i, int j, int k); //Get Omega for Vorticity Confinement

   vec3 clampPoint(vec3 p);

   GridDataX mU; // X component of velocity, stored on X faces, size is (dimX+1)*dimY*dimZ
   GridDataY mV; // Y component of velocity, stored on Y faces, size is dimX*(dimY+1)*dimZ
   GridDataZ mW; // W component of velocity, stored on Z faces, size is dimX*dimY*(dimZ+1)
   GridData mP;  // Pressure, stored at grid centers, size is dimX*dimY*dimZ
   GridData mD;  // Density, stored at grid centers, size is dimX*dimY*dimZ
   GridData mT;  // Temperature, stored at grid centers, size is dimX*dimY*dimZ
   GridData mS;  // Solid, marker stored at grid centers; 1.0 if Solid, 0.0 if Fluid

   template<typename Matrix, typename Vector>
   bool applyPreconditioner(const Matrix &A, const Vector &precon, const Vector &r, Vector &z)
   {

	   Vector q(theDim[0]*theDim[1]*theDim[2]);
	   q.clear();
	   double t;

	   //Solve Lq=r
	   for(int k = 0; k < theDim[MACGrid::Z]; k++){ 
		   for(int j = 0; j < theDim[MACGrid::Y]; j++){
			   for(int i = 0; i < theDim[MACGrid::X]; i++){
				   if(mS(i,j,k) < EPSILON){
					   int p_ijk = getIndex(i,j,k);
					   int p_in1jk = getIndex(i-1,j,k);
					   int p_ijn1k = getIndex(i,j-1,k);
					   int p_ijkn1 = getIndex(i,j,k-1);
					   double t1=0.0;
					   if(i-1 >= 0)
						   t1 = A(p_in1jk,p_ijk)*precon(p_in1jk)*q(p_in1jk);
					   double t2=0.0;
					   if(j-1 >= 0)
						   t2 = A(p_ijn1k,p_ijk)*precon(p_ijn1k)*q(p_ijn1k);
					   double t3=0.0;
					   if(k-1 >= 0)
						   t3 = A(p_ijkn1,p_ijk)*precon(p_ijkn1)*q(p_ijkn1);

					   t = r(p_ijk)-t1-t2-t3;
					   q(p_ijk) = t * precon(p_ijk);
				   }
			   }
		   }
	   }

	   //Solve transpose(L)z=q
	   for(int k=theDim[MACGrid::Z]-1; k >= 0; k--){ 
		   for(int j=theDim[MACGrid::Y]-1; j >= 0; j--){
			   for(int i=theDim[MACGrid::X]-1; i >= 0; i--){
				    if(mS(i,j,k) < EPSILON){
					   int p_ijk = getIndex(i,j,k);
					   int p_i1jk = getIndex(i+1,j,k);
					   int p_ij1k = getIndex(i,j+1,k);
					   int p_ijk1 = getIndex(i,j,k+1);
					   double t1=0.0;
					   if(i+1 < theDim[MACGrid::X])
						   t1 = A(p_ijk,p_i1jk)*precon(p_ijk)*z(p_i1jk);
					   double t2=0.0;
					   if(j+1 < theDim[MACGrid::Y])
						   t2 = A(p_ijk,p_ij1k)*precon(p_ijk)*z(p_ij1k);
					   double t3=0.0;
					   if(k+1 < theDim[MACGrid::Z])
						   t3 = A(p_ijk,p_ijk1)*precon(p_ijk)*z(p_ijk1);

					   t = q(p_ijk)-t1-t2-t3;
					   z(p_ijk) = t * precon(p_ijk);
				   }
			   }
		   }
	   }

	   return true;
   }
   
   template<typename Matrix, typename Vector>
   bool cg_preconditioned_solve(const Matrix &A, const Vector &b, const Vector &precon, Vector &x, int max_iter, double tol)
   {
	    Vector z(theDim[0]*theDim[1]*theDim[2]);
		std::fill(x.begin(), x.end(), 0);
		std::fill(z.begin(), z.end(), 0);
		Vector r = b;
		applyPreconditioner(A,precon,r,z); 
		Vector d = z;
		Vector temp;
		double resign;
		for(int niter = 0; niter < max_iter; niter++)
		{
		  temp = prod(A, d);
		  double dotzr = inner_prod(z, r);
		  double alpha = dotzr/inner_prod(d, temp);
		  x += alpha*d;

		  r -= alpha*temp;
		  resign = norm_2(r);
		  if(resign < tol) 
		  { 
			std::cout << "numiters: "<< niter <<std::endl;
			return true; 
		  }

		  applyPreconditioner(A,precon,r,z);
		  double beta = inner_prod(z, r)/dotzr;
		  d = z + beta*d;
	   }

	   std::cout << "WARNING: cg_preconditioned_solve did not converge: " << resign << std::endl;
	   return false;
	}

public:
   enum RenderMode { CUBES, SHEETS };
   static RenderMode theRenderMode;
   static bool theDisplayVel;
};

#endif