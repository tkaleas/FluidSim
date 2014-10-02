#include "waterSim.h"
#include <GL/glut.h>
#include <IL/il.h>
#include <IL/ilu.h>
#include <IL/ilut.h>

WaterSim::WaterSim() : mFrameNum(0), mRecordEnabled(false)
{
   ilInit();
   iluInit();
   ilEnable(IL_FILE_OVERWRITE);
   ilutRenderer(ILUT_OPENGL);

   reset();
}

WaterSim::~WaterSim()
{
}

void WaterSim::reset()
{
   mGrid.reset();
}

void WaterSim::setGridDimensions(int x, int y, int z)
{
   extern int theDim[3]; // Naughty globals...
   theDim[0] = x;
   theDim[1] = y;
   theDim[2] = z;
   reset();
}

void WaterSim::step()
{
   double dt = 0.01;

   //Step00: Extrapolate Divergence Free Velocity Field into Air Region
   mGrid.extrapolateVelocities(dt);
   
   //Step0: Gather User Forces
   mGrid.updateSources();

   //Step1a: Update the Grid Based on Marker Particle Locations (Fluid, Air Solid)
   mGrid.updateGrid();
   
   // Step2: Calculate new velocities
   mGrid.addExternalForces(dt);
   mGrid.advectVelocity(dt);
   
   //Step1b: Update Level Sets To Reflect True Fluid Boundary
   mGrid.updateLevelSet(dt);
   //mGrid.initializeLevelSet();
   
   //Step 2b: Project Velocities to Ensure Incompressibility
   mGrid.project(dt);
   
   //Step3: Move Marker Particles Through Time
   mGrid.moveParticles(dt);
}

void WaterSim::setRecording(bool on)
{
   if (on && ! mRecordEnabled)  // reset counter
   {
      mFrameNum = 0;
   }
   mRecordEnabled = on;
}

bool WaterSim::isRecording()
{
   return mRecordEnabled;
}

void WaterSim::draw(const Camera& c)
{
   drawAxes(); 
   mGrid.draw(c);
   if (mRecordEnabled){ 
	   grabScreen();
	   outputOBJ();
   }
}

void WaterSim::drawAxes()
{
  glPushAttrib(GL_LIGHTING_BIT | GL_LINE_BIT);
      glDisable(GL_LIGHTING);

      glLineWidth(2.0); 
      glBegin(GL_LINES);
         glColor3f(1.0, 0.0, 0.0);
         glVertex3f(0.0, 0.0, 0.0);
         glVertex3f(1.0, 0.0, 0.0);

         glColor3f(0.0, 1.0, 0.0);
         glVertex3f(0.0, 0.0, 0.0);
         glVertex3f(0.0, 1.0, 0.0);

         glColor3f(0.0, 0.0, 1.0);
         glVertex3f(0.0, 0.0, 0.0);
         glVertex3f(0.0, 0.0, 1.0);
      glEnd();
  glPopAttrib();
}

void WaterSim::grabScreen()  // Code adapted from asst#1
{
	unsigned int image;
   ilGenImages(1, &image);
	ilBindImage(image);

	ILenum error = ilGetError();
	assert(error == IL_NO_ERROR);

	ilTexImage(640, 480, 1, 3, IL_RGB, IL_UNSIGNED_BYTE, NULL);

	error = ilGetError();
	assert(error == IL_NO_ERROR);

	unsigned char* data = ilGetData();

	error = ilGetError();
	assert(error == IL_NO_ERROR);

	for (int i=479; i>=0; i--) 
	{
		glReadPixels(0,i,640,1,GL_RGB, GL_UNSIGNED_BYTE, 
			data + (640 * 3 * i));
	}

	char anim_filename[2048];
	sprintf_s(anim_filename, 2048, "output/water_%04d.png", mFrameNum++); 

	ilSave(IL_PNG, anim_filename);

	error = ilGetError();
	assert(error == IL_NO_ERROR);

	ilDeleteImages(1, &image);

	error = ilGetError();
	assert(error == IL_NO_ERROR);
}

void WaterSim::outputOBJ() //For Outputting OBJ
{
	//Marker
	char anim_filename[2048];
	sprintf_s(anim_filename, 2048, "output_obj_m/water_%d.obj", mFrameNum);

	std::ofstream out;

	if(MACGrid::theDisplayMarkMesh) {
		out.open(anim_filename);
		out << *(mGrid.mSurfaceM) << endl;
		out.close();
	}
	

	//Level Set
	sprintf_s(anim_filename, 2048, "output_obj_ls/water_%d.obj", mFrameNum);
	
	if(MACGrid::theDisplayLevSet){
		out.open(anim_filename);
		out << *(mGrid.mSurfaceLS) << endl;
		out.close();
	}
	
}