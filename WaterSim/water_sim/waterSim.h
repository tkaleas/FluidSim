// ==========================================================================
// Copyright (C) 2009 Aline Normoyle
// ==========================================================================
#ifndef waterSim_H_
#define waterSim_H_

#include "MACGrid.h"

class Camera;
class WaterSim
{
public:
   WaterSim();
   virtual ~WaterSim();

   virtual void reset();
   virtual void step();
   virtual void draw(const Camera& c);
   virtual void setGridDimensions(int x, int y, int z);
   virtual void setRecording(bool on);
   virtual bool isRecording();

protected:
   virtual void drawAxes();
   virtual void grabScreen();
   virtual void outputOBJ();

protected:
   MACGrid mGrid;
   bool mRecordEnabled;
   int mFrameNum;
};

#endif