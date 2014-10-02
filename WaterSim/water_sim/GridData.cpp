#include "GridData.h"
#include "LevelSet.h"
#include "main.h"

#define LERP(a,b,t) (1-t)*a + t*b

#ifdef _DEBUG
int theDim[3] = {5,5, 1};
#else
int theDim[3] = {15, 15, 1};
#endif
double theCellSize = 0.5;

GridData::GridData() :
   mDfltValue(0.0), mMax(0.0,0.0,0.0)
{
}

GridData::GridData(const GridData& orig) :
   mDfltValue(orig.mDfltValue)
{
   mData = orig.mData;
   mMax = orig.mMax;
}

GridData::~GridData() 
{
}

ublas::vector<double>& GridData::data()
{
   return mData;
}

GridData& GridData::operator=(const GridData& orig)
{
   if (this == &orig)
   {
      return *this;
   }
   mDfltValue = orig.mDfltValue;
   mData = orig.mData;
   mMax = orig.mMax;
   return *this;
}

void GridData::initialize(double dfltValue)
{
   mDfltValue = dfltValue;
   mMax[0] = theCellSize*theDim[0];
   mMax[1] = theCellSize*theDim[1];
   mMax[2] = theCellSize*theDim[2];
   mData.resize(theDim[0]*theDim[1]*theDim[2], false);
   std::fill(mData.begin(), mData.end(), mDfltValue);
}

double& GridData::operator()(int i, int j, int k)
{
	////Default Value Method
	//if (i< 0 || j<0 || k<0 || 
	//	i > theDim[0]-1 || 
	//	j > theDim[1]-1 || 
	//	k > theDim[2]-1) return mDfltValue;
	
	//Clamp Value Method
	if(i < 0){i=0;}
	if(i > theDim[0]-1) {i=theDim[0]-1;}
	if(j < 0){j=0;}
	if(j > theDim[1]-1) {j=theDim[1]-1;} 
	if(k < 0){k=0;}
	if(k > theDim[2]-1) {k=theDim[2]-1;} 

	int col = i;
	int row = k*theDim[0];
	int stack = j*theDim[0]*theDim[2];
	return mData(col+row+stack);
}

void GridData::getCell(const vec3& pt, int& i, int& j, int& k)
{
   vec3 pos = pt;  //= worldToSelf(pt); 
   i = (int) (pos[0]/theCellSize);
   j = (int) (pos[1]/theCellSize);
   k = (int) (pos[2]/theCellSize);   
}

double GridData::interpolate(const vec3& pt)
{
   vec3 pos = worldToSelf(pt);

   int i = (int) (pos[0]/theCellSize);
   int j = (int) (pos[1]/theCellSize);
   int k = (int) (pos[2]/theCellSize);

   double scale = 1.0/theCellSize;  
   double fractx = scale*(pos[0] - i*theCellSize);
   double fracty = scale*(pos[1] - j*theCellSize);
   double fractz = scale*(pos[2] - k*theCellSize);

   assert (fractx < 1.0 && fractx >= 0);
   assert (fracty < 1.0 && fracty >= 0);
   assert (fractz < 1.0 && fractz >= 0);

   double tmp1 = (*this)(i,j,k);
   double tmp2 = (*this)(i,j+1,k);
   double tmp3 = (*this)(i+1,j,k);
   double tmp4 = (*this)(i+1,j+1,k);

   double tmp5 = (*this)(i,j,k+1);
   double tmp6 = (*this)(i,j+1,k+1);
   double tmp7 = (*this)(i+1,j,k+1);
   double tmp8 = (*this)(i+1,j+1,k+1);

   double tmp12 = LERP(tmp1, tmp2, fracty);
   double tmp34 = LERP(tmp3, tmp4, fracty);

   double tmp56 = LERP(tmp5, tmp6, fracty);
   double tmp78 = LERP(tmp7, tmp8, fracty);

   double tmp1234 = LERP (tmp12, tmp34, fractx);
   double tmp5678 = LERP (tmp56, tmp78, fractx);

   double tmp = LERP(tmp1234, tmp5678, fractz);
   return tmp;
}

vec3 GridData::worldToSelf(const vec3& pt) const
{
   vec3 out;
   out[0] = min(max(0.0, pt[0] - theCellSize*0.5), mMax[0]);
   out[1] = min(max(0.0, pt[1] - theCellSize*0.5), mMax[1]);
   out[2] = min(max(0.0, pt[2] - theCellSize*0.5), mMax[2]);
   return out;
}

GridDataX::GridDataX() : GridData()
{
}

GridDataX::~GridDataX()
{
}

void GridDataX::initialize(double dfltValue)
{
   GridData::initialize(dfltValue);
   mMax[0] = theCellSize*(theDim[0]+1);
   mMax[1] = theCellSize*theDim[1];
   mMax[2] = theCellSize*theDim[2];
   mData.resize((theDim[0]+1)*theDim[1]*theDim[2], false);
   std::fill(mData.begin(), mData.end(), mDfltValue);
}

double& GridDataX::operator()(int i, int j, int k)
{
   //if (i < 0 || i > theDim[0]) return mDfltValue;
   if(i < 0) {i = 0;} 
   if(i > theDim[0]) {i = theDim[0];}

   if (j < 0) j = 0;
   if (j > theDim[1]-1) j = theDim[1]-1;
   if (k < 0) k = 0;
   if (k > theDim[2]-1) k = theDim[2]-1;

   int col = i;
   int row = k*(theDim[0]+1);
   int stack = j*(theDim[0]+1)*theDim[2];
   return mData(stack + row + col);
}

vec3 GridDataX::worldToSelf(const vec3& pt) const
{   
   vec3 out;
   out[0] = min(max(0.0, pt[0]), mMax[0]);
   out[1] = min(max(0.0, pt[1]-theCellSize*0.5), mMax[1]);
   out[2] = min(max(0.0, pt[2]-theCellSize*0.5), mMax[2]);
   return out;
}

GridDataY::GridDataY() : GridData()
{

}

GridDataY::~GridDataY()
{
}

void GridDataY::initialize(double dfltValue)
{
   GridData::initialize(dfltValue);
   mMax[0] = theCellSize*theDim[0];
   mMax[1] = theCellSize*(theDim[1]+1);
   mMax[2] = theCellSize*theDim[2];
   mData.resize(theDim[0]*(theDim[1]+1)*theDim[2], false);
   std::fill(mData.begin(), mData.end(), mDfltValue);
}

double& GridDataY::operator()(int i, int j, int k)
{
   //if (j < 0 || j > theDim[1]) return mDfltValue;
   if(j < 0) {j = 0;} 
   if(j > theDim[1]) {j = theDim[1];}

   if (i < 0) i = 0;
   if (i > theDim[0]-1) i = theDim[0]-1;
   if (k < 0) k = 0;
   if (k > theDim[2]-1) k = theDim[2]-1;

   int col = i;
   int row = k*theDim[0];
   int stack = j*theDim[0]*theDim[2];
   return mData(stack + row + col);
}

vec3 GridDataY::worldToSelf(const vec3& pt) const
{
   vec3 out;
   out[0] = min(max(0.0, pt[0]-theCellSize*0.5), mMax[0]);
   out[1] = min(max(0.0, pt[1]), mMax[1]);
   out[2] = min(max(0.0, pt[2]-theCellSize*0.5), mMax[2]);
   return out;
}

GridDataZ::GridDataZ() : GridData()
{
}

GridDataZ::~GridDataZ()
{
}

void GridDataZ::initialize(double dfltValue)
{
   GridData::initialize(dfltValue);
   mMax[0] = theCellSize*theDim[0];
   mMax[1] = theCellSize*theDim[1];
   mMax[2] = theCellSize*(theDim[2]+1);
   mData.resize(theDim[0]*theDim[1]*(theDim[2]+1), false);
   std::fill(mData.begin(), mData.end(), mDfltValue);
}

double& GridDataZ::operator()(int i, int j, int k)
{
   //if (k < 0 || k > theDim[2]) return mDfltValue;
   if(k < 0) {k = 0;} 
   if(k > theDim[2]) {k = theDim[2];}
   
   if (i < 0) i = 0;
   if (i > theDim[0]-1) i = theDim[0]-1;
   if (j < 0) j = 0;
   if (j > theDim[1]-1) j = theDim[1]-1;

   int col = i;
   int row = k*theDim[0];
   int stack = j*theDim[0]*(theDim[2]+1);

   return mData(stack + row + col);
}

vec3 GridDataZ::worldToSelf(const vec3& pt) const
{
   vec3 out;
   out[0] = min(max(0.0, pt[0]-theCellSize*0.5), mMax[0]);
   out[1] = min(max(0.0, pt[1]-theCellSize*0.5), mMax[1]);
   out[2] = min(max(0.0, pt[2]), mMax[2]);
   return out;
}

LevelSet::LevelSet() : GridData()
{
	mFM = new FastMarch(theDim[0],theDim[1],theDim[2],theCellSize);
}

LevelSet::~LevelSet()
{
	delete mFM;
}

void LevelSet::initialize(double dfltValue)
{
   GridData::initialize(dfltValue);
   mMax[0] = theCellSize*(theDim[0]+2);
   mMax[1] = theCellSize*(theDim[1]+2);
   mMax[2] = theCellSize*(theDim[2]+2);
   mData.resize((theDim[0]+2)*(theDim[1]+2)*(theDim[2]+2), false);
   std::fill(mData.begin(), mData.end(), mDfltValue);
}

double& LevelSet::operator()(int i, int j, int k)
{
   //Clamp Value Method
   if(i < 0){i=0;}
   if(i > theDim[0]+1) {i=theDim[0]+1;}
   if(j < 0){j=0;}
   if(j > theDim[1]+1) {j=theDim[1]+1;} 
   if(k < 0){k=0;}
   if(k > theDim[2]+1) {k=theDim[2]+1;}

   int col = i;
   int row = k*(theDim[0]+2);
   int stack = j*(theDim[0]+2)*(theDim[2]+2);

   return mData(stack + row + col);
}

vec3 LevelSet::worldToSelf(const vec3& pt) const
{
   vec3 out;
   out[0] = min(max(0.0, pt[0]+theCellSize*0.5), mMax[0]);
   out[1] = min(max(0.0, pt[1]+theCellSize*0.5), mMax[1]);
   out[2] = min(max(0.0, pt[2]+theCellSize*0.5), mMax[2]);
   return out;
}

Double	LevelSet::eval(const Point3d& location)
{
		//vec3 pos((location[0] + 1) * (theDim[0]-1) * 0.5 + 1, 
		//		   (location[1] + 1) * (theDim[1]-1) * 0.5 + 1,
		//		   (location[2] + 1) * (theDim[2]-1) * 0.5 + 1);
		vec3 pos(location[0],location[1],location[2]);
		//std::cout << pos << std::endl;
		return interpolate(pos);
}

void LevelSet::ReInitialize(){
	mFM->Reinitialize(*this);
	SetBoundarySignedDist();
}

inline int GILS(int i, int j, int k) { return i + (theDim[0]+2)*k + (theDim[0]+2)*(theDim[2]+2)*j; }

void LevelSet::SetBoundarySignedDist()
{
    static double phi = 3.0 * theCellSize;
	FOR_GRIDZY {
		data()[GILS(0,j,k)] = data()[GILS(theDim[0]+1,j,k)] = phi; 
	}
	END_FOR_TWO
	
	FOR_GRIDZX {
		data()[GILS(i,0,k)] = data()[GILS(i,theDim[1]+1,k)] = phi; 
	}
	END_FOR_TWO
	FOR_GRIDYX { 
		data()[GILS(i,j,0)] = data()[GILS(i,j,theDim[2]+1)] = phi; 
	}
	END_FOR_TWO
	for(int k=1 ; k<=theDim[2] ; k++ ) {
		data()[GILS(0,0   ,k)] = data()[GILS(theDim[0]+1,0   ,k)] = phi;
        data()[GILS(0,theDim[1]+1,k)] = data()[GILS(theDim[0]+1,theDim[1]+1,k)] = phi;
	}
	for(int j=1 ; j<=theDim[1] ; j++ ) {
		data()[GILS(0,j,0   )] = data()[GILS(theDim[0]+1,j,0   )] = phi;
        data()[GILS(0,j,theDim[2]+1)] = data()[GILS(theDim[0]+1,j,theDim[2]+1)] = phi;
	}
	for(int i=1 ; i<=theDim[0]; i++ ) { 
		data()[GILS(i,0,0   )] = data()[GILS(i,theDim[1]+1,0   )] = phi;
        data()[GILS(i,0,theDim[2]+1)] = data()[GILS(i,theDim[1]+1,theDim[2]+1)] = phi;
	}
    data()[GILS(0   ,0   ,0   )] = data()[GILS(theDim[0]+1,0   ,0   )] = phi;
	data()[GILS(0   ,theDim[1]+1,0   )] = data()[GILS(0   ,0   ,theDim[2]+1)] = phi;
    data()[GILS(theDim[0]+1,theDim[1]+1,0   )] = data()[GILS(theDim[0]+1,0   ,theDim[2]+1)] = phi;
    data()[GILS(0   ,theDim[1]+1,theDim[2]+1)] = data()[GILS(theDim[0]+1,theDim[1]+1,theDim[2]+1)] = phi;
}