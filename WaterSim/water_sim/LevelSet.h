
#ifndef LEVELSET_H
#define LEVELSET_H

#include "FastMarch.h"
#include "MarchCubes.h"
#include "GridData.h"

#include <vector>

class LevelSet: public GridData, public ImpSurface
{
public:
   LevelSet();
   virtual ~LevelSet();
   virtual void initialize(double dfltValue = 0.0);
   virtual double& operator()(int i, int j, int k);
   virtual vec3 worldToSelf(const vec3& pt) const;

   virtual Double	eval	(const Point3d& location);
   void SetBoundarySignedDist();
   void ReInitialize();
   //Level Set Variables
   FastMarch *mFM;
};

#endif