/**************************************************************************
	ORIGINAL AUTHOR: 
		Emud Mokhberi (emud@ucla.edu)
	MODIFIED BY:
	
	CONTRIBUTORS:
		

-----------------------------------------------
	
 ***************************************************************
 ******General License Agreement and Lack of Warranty ***********
 ****************************************************************

 This software is distributed for noncommercial use in the hope that it will 
 be useful but WITHOUT ANY WARRANTY. The author(s) do not accept responsibility
 to anyone for the consequences of using it or for whether it serves any 
 particular purpose or works at all. No guarantee is made about the software 
 or its performance.

 You are allowed to modify the source code, add your name to the
 appropriate list above and distribute the code as long as 
 this license agreement is distributed with the code and is included at 
 the top of all header (.h) files.

 Commercial use is strictly prohibited.
***************************************************************************/

/*
	FastMarch2D : A class for resetting the grid values to be the signed distance from the interface
	Inputs: grid and cell size
			
	The class maintains its own special grid of type FMContainer and in addition to two lists, one
	of ClosePoints which is used to determine the initial set of close points and a second list for
	maintaining the min heap which is used to determine the signed distance values
	
	The class works by first resetting all the negative signed distance values and then setting all
	of the positive signed ditance values 

	Public Functions:
	Set				- initialized the value of the FMContainer grid for the specified index
	Reinitialize	- Performs the fastmarching method on the FMContainer grid. It is assumed
					  that the values to be reset are in the FMContainer grid and that is where
					  the updated grid values will be when the function is done.
	
	Private Functions:
	PopHeap			- Called by FastMarch, it pops the current closest grid value from the heap 
					  and cosidered that cell value to be done. It then restores the heap to be a min heap
	FastMarch		- Called by ReinitHalf. While there are close points in the heap, it pops the heap, 
					  determines if the threshold for the fastmarching has been met (we only need to restore
					  the signed distance function to within a finite number of cell away from the
					  interface) and if it hasn't, it updates the neighboring cells if they are not
					  considered done points
	AddToHeap		- Called by FindPhi. Adds a former far point which has become a close point to the 
				      heap and updates the heap to remain a min heap
	UpdateHeap		- Called by FindPhi. This function is called if a close point's (which is in the heap) 
					  value is changed as a result of a neighboring cell being 'done'. It updates the heap 
					  to remain a min heap
	CheckMax2		- Called by FindPhi. This function is used if a close point has two neighboring done
					  points. It is used to make sure that the distance between the two neighboring done
					  points is not greater than the cell size. If it is, it only uses the smaller of the
					  two values in updating the signed distance of the close point
	CheckMax3		- Called by FindPhi. This function is used if a close point has three neighboring done
					  points. It is used to make sure that the distance between the three neighboring done
					  points is not greater than the cell size. If it is, it only uses the two smaller
					  values in updating the signed distance of the close point
	CheckFront		- Called by FindPhi. This function checks to see if the cell in front of the current
					  close point is done
	CheckBehind		- Called by FindPhi. This function checks to see if the cell behind of the current
					  close point is done. This function also handles the case where a close points has
					  done points on both sides along one axis
	FindPhi			- Called by FastMarch. This function updates the value os the specified cell using 
					  the values of neighboring 'done' cells. If the cell is already a close point, it 
					  calls UpdateHeap. Otherwise, it changes the cell to be a close point and calls
					  AddToHeap.
	SetBoundary		- Called by ReinitHalf. Sets the boundary grid done flags to -1
	ReinitHalf		- Called by Reinitialize. Resets either the exterior or interior of the interface
					  to be a signed distance function.
	Initialize		- Called by ReinitHalf. Steps through the grid along each axis and determines which 
					  cells are 'done', and which are 'close'. It uses linear interpolation to set the 
					  values of those cells adjacent to the interface.
	InitHeap		- Called by ReinitHalf. After Initialize has run and determined the first set of
					  'close' points, this function will create a min heap of those points. This is
					  done by stepping through the ClosePoints list and calling FindPhi for each cell
					  that isn't already in the heap
	AddClose		- Adds cell to ClosePoints list.
					  

	Created by Emud Mokhberi: UCLA : 09/04/04
*/

#ifndef FastMarch_H
#define FastMarch_H

#include <vector>
#include "main.h"
class LevelSet;


class FastMarch
{
public:
	FastMarch(int nx, int ny, int nz, double hi);
	~FastMarch() { delete [] grid; }

    inline double& operator[] (int index) { return grid[index].value; }
	inline const double& operator[] (int index) const { return grid[index].value; }
	inline double& operator() (int i, int j, int k) { return grid[GI(i,j,k)].value; }
	inline const double& operator() (int i, int j, int k) const { return grid[GI(i,j,k)].value; }

    inline void Set(int index, const double &value);
	void Reinitialize(LevelSet &lset);

private:
	int PopHeap();
    void March();
    void AddToHeap(int index);
	void UpdateHeap(int index);
    inline void CheckMax2(int& a, double& phi1, const double &phi2);
    inline void CheckMax3(int& a, bool& flag, double& phi1, 
                          const double &phi2, const double &phi3);
    inline void CheckFront(double& phi, int& a, bool& flag, int index);
	inline void CheckBehind(double& phi, int& a, bool& flag, int index);
    void FindPhi(int index, int x, int y, int z);
    void SetBoundary();
    inline void ReinitHalf();
	void Initialize();		
	void InitHeap();
    inline void AddClose(int index);

    inline int GI(int i, int j, int k) const { return i + dj*j + dk*k; }
	inline void GIJK(int index, int& i, int& j, int& k)
        { k = int(index * NyNxInv); j = Mod(index,dk,NyNxInv);
          i = Mod(j,dj,NxInv);      j = int(j * NxInv); }

	struct FMContainer {
		int DoneFlag;
		int HeapPosition;
		double value;
	};

	int Nx, Ny, Nz, size, dj, dk;
    double h, hInv, NxInv, NyNxInv;

	int heapSize;
	int closeSize;
	FMContainer* grid;
	std::vector<int> FMHeap;
	std::vector<int> ClosePoints;
};

#endif