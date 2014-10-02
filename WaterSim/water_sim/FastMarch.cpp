#include "FastMarch.h"
#include "LevelSet.h"
#include "Grid.h"

#define FOR_GRID   for(int i=0; i < size; i++)
#define FOR_GRIDZY for(int k=0; k <= Nz; k++) { for(int j=0; j <= Ny; j++) { 
#define FOR_GRIDZX for(int k=0; k <= Nz; k++) { for(int i=0; i <= Nx; i++) { 
#define FOR_GRIDYX for(int j=0; j <= Ny; j++) { for(int i=0; i <= Nx; i++) {
#define END_FOR_THREE }}}
#define END_FOR_TWO }}

FastMarch::FastMarch(int nx,int ny, int nz, double hi) 
    : Nx(nx), Ny(ny), Nz(nz), h(hi), hInv(1./hi), size((nx+2)*(ny+2)*(nz+2)), 
      dj(nx+2), dk((nx+2)*(ny+2)), NxInv(1./(nx+2)), NyNxInv(1./((nx+2)*(ny+2)))
{
	FMHeap.reserve(size);
	ClosePoints.reserve(size);
	grid = new FMContainer[size];
}

void FastMarch::Reinitialize(LevelSet &lset) {
	Grid temp_grid(Nx,Ny,Nz);
	//To Temp Grid
	for(int i=0; i <= Nx+1; i++){
		for(int j=0; j <= Ny+1; j++){
			for(int k=0; k <= Nz+1; k++){
				temp_grid(i,j,k) = lset(i,j,k);
			}
		}
	}
    
	//Negative Phi first
	FOR_GRID Set(i, temp_grid[i]);
	ReinitHalf();
    //Then Positive Phi
    FOR_GRID Set(i, grid[i].value);
    ReinitHalf();
    FOR_GRID temp_grid[i] = grid[i].value;
	
	//Transfer Back to Level Set
	for(int i=0; i <= Nx+1; i++){
		for(int j=0; j <= Ny+1; j++){
			for(int k=0; k <= Nz+1; k++){
				lset(i,j,k) = temp_grid(i,j,k); 
			}
		}
	}

}

inline void FastMarch::ReinitHalf() {
    heapSize = 0;
	closeSize = 0;
    SetBoundary();
	Initialize();
    //PrintFlags();
	InitHeap();
	March();
}

inline void FastMarch::Set(int index, const double &value) {
	grid[index].value = -value; // - is to swap values for initializing negative phi
	grid[index].HeapPosition = -1;
	if(grid[index].value < 0.) grid[index].DoneFlag = -1;
	else                       grid[index].DoneFlag = 0;
}

void FastMarch::SetBoundary()
{
    FOR_GRIDZY grid[GI(0,j,k)].DoneFlag = grid[GI(Nx+1,j,k)].DoneFlag = -1; END_FOR_TWO
    FOR_GRIDZX grid[GI(i,0,k)].DoneFlag = grid[GI(i,Ny+1,k)].DoneFlag = -1; END_FOR_TWO
    FOR_GRIDYX grid[GI(i,j,0)].DoneFlag = grid[GI(i,j,Nz+1)].DoneFlag = -1; END_FOR_TWO
    for(int k=1 ; k<=Nz ; k++ ) {
		grid[GI(0,0   ,k)].DoneFlag = grid[GI(Nx+1,0   ,k)].DoneFlag = -1;
        grid[GI(0,Ny+1,k)].DoneFlag = grid[GI(Nx+1,Ny+1,k)].DoneFlag = -1;
	}
	for(int j=1 ; j<=Ny ; j++ ) {
		grid[GI(0,j,0   )].DoneFlag = grid[GI(Nx+1,j,0   )].DoneFlag = -1;
        grid[GI(0,j,Nz+1)].DoneFlag = grid[GI(Nx+1,j,Nz+1)].DoneFlag = -1;
	}
	for(int i=1 ; i<=Nx; i++ ) { 
        grid[GI(i,0,0   )].DoneFlag = grid[GI(i,Ny+1,0   )].DoneFlag = -1;
        grid[GI(i,0,Nz+1)].DoneFlag = grid[GI(i,Ny+1,Nz+1)].DoneFlag = -1;
	}
    grid[GI(0   ,0   ,0   )].DoneFlag = grid[GI(Nx+1,0   ,0   )].DoneFlag = -1;
	grid[GI(0   ,Ny+1,0   )].DoneFlag = grid[GI(0   ,0   ,Nz+1)].DoneFlag = -1;
    grid[GI(Nx+1,Ny+1,0   )].DoneFlag = grid[GI(Nx+1,0   ,Nz+1)].DoneFlag = -1;
    grid[GI(0   ,Ny+1,Nz+1)].DoneFlag = grid[GI(Nx+1,Ny+1,Nz+1)].DoneFlag = -1;
}

void FastMarch::Initialize() {
    static int ci, pi, flag;
	for(int i = 1; i <= Nx; i++)
	{
        for(int j = 1; j <= Ny; j++)
		{
            //if(Done == 0)  f= 0: change --> -1^f = -1;           NO -->  0^f = 0; 1^f = 1;
            //if(Done == -1) f=-1: change -->  0^f = -1; 1^f = -2; NO --> -1^f = 0;
            //if(Done == 1)  f= 1: change --> -1^f = -2;           NO -->  0^f = 1; 1^f = 0;
            flag = grid[GI(i,j,1)].DoneFlag;
            for(int k = 2; k <= Nz; k++)
            {
                ci = GI(i,j,k);
                if((flag ^ grid[ci].DoneFlag) < 0)
                {
                    pi = GI(i,j,k-1);
                    flag = grid[ci].DoneFlag;
                    if(flag >= 0) {
                        grid[ci].DoneFlag = 1;
                        AddClose(GI(i,j,k+1));
                        grid[ci].value = min(grid[ci].value, abs(h + grid[pi].value));
                    }
                    else {
                        grid[pi].DoneFlag = 1;
                        AddClose(GI(i,j,k-2));
                        grid[pi].value = min(grid[pi].value, abs(h + grid[ci].value));
                    }
                }
            }
        }
    }
    for(int i = 1; i <= Nx; i++)
	{
        for(int k = 1; k <= Nz; k++)
		{            
            flag = grid[GI(i,1,k)].DoneFlag;
            for(int j = 2; j <= Ny; j++)
            {
                ci = GI(i,j,k);
                if((flag ^ grid[ci].DoneFlag) < 0)
                {
                    pi = GI(i,j-1,k);
                    flag = grid[ci].DoneFlag;
                    if(flag >= 0) {
                        grid[ci].DoneFlag = 1;
                        AddClose(GI(i,j+1,k));
                        grid[ci].value = min(grid[ci].value, abs(h + grid[pi].value));
                    }
                    else {
                        grid[pi].DoneFlag = 1;
                        AddClose(GI(i,j-2,k));
                        grid[pi].value = min(grid[pi].value, abs(h + grid[ci].value));
                    }
                }
            }
        }
    }
    
	for(int j = 1; j <= Ny; j++)
	{
        for(int k = 1; k <= Nz; k++)
		{            
            flag = grid[GI(1,j,k)].DoneFlag;
            for(int i = 2; i <= Nx; i++)
            {
                ci = GI(i,j,k);
                if((flag ^ grid[ci].DoneFlag) < 0)
                {
                    pi = GI(i-1,j,k);
                    flag = grid[ci].DoneFlag;
                    if(flag >= 0) {
                        grid[ci].DoneFlag = 1;
                        AddClose(GI(i+1,j,k));
                        grid[ci].value = min(grid[ci].value, abs(h + grid[pi].value));
                    }
                    else {
                        grid[pi].DoneFlag = 1;
                        AddClose(GI(i-2,j,k));
                        grid[pi].value = min(grid[pi].value, abs(h + grid[ci].value));
                    }
                }
            }
        }
    }
}

inline void FastMarch::AddClose(int index)
{
	if(grid[index].DoneFlag == 0) {
		//Add index to list for initialization of close band
		ClosePoints[closeSize] = index;
		closeSize++;
	}
}

void FastMarch::InitHeap() {
    static int x,y,z;
	for(int i = 0; i < closeSize; i++) {
		if(grid[ClosePoints[i]].HeapPosition == -1 && grid[ClosePoints[i]].DoneFlag == 0) {
			GIJK(ClosePoints[i], x, y, z);
			FindPhi(ClosePoints[i], x, y,z);
		}
	}
}

void FastMarch::FindPhi(int index, int x, int y, int z) {
	static double phiX, phiY, phiZ, b, quotient, phi;
    static int a;
    static bool flagX, flagY, flagZ;

    phiX = phiY = phiZ = 0.;
    a = 0;
    flagX = flagY = flagZ = 0;

	//Find The phiS
	CheckFront (phiX, a, flagX, GI(x+1,  y,  z));
	CheckBehind(phiX, a, flagX, GI(x-1,  y,  z));
	CheckFront (phiY, a, flagY, GI(x  ,y+1,  z));
	CheckBehind(phiY, a, flagY, GI(x  ,y-1,  z));
    CheckFront (phiZ, a, flagZ, GI(x  ,y  ,z+1));
	CheckBehind(phiZ, a, flagZ, GI(x  ,y  ,z-1));

    //Max Tests
	if(a == 3) {
		if     ((phiX >= phiY) && (phiX >= phiZ)) CheckMax3(a, flagX, phiX, phiY, phiZ);
		else if((phiY >= phiX) && (phiY >= phiZ)) CheckMax3(a, flagY, phiY, phiX, phiZ);
		else									  CheckMax3(a, flagZ, phiZ, phiX, phiY);
	}
	if(a == 2) {
		if(!flagX) {
			if(phiY >= phiZ) CheckMax2(a, phiY, phiZ);
			else			 CheckMax2(a, phiZ, phiY);
		}
		else if(!flagY){
			if(phiX >= phiZ) CheckMax2(a, phiX, phiZ);
			else			 CheckMax2(a, phiZ, phiX);
		}
		else {
			if(phiX >= phiY) CheckMax2(a, phiX, phiY);
			else			 CheckMax2(a, phiY, phiX);
		}
	}

	b = phiX + phiY + phiZ;
	quotient = square(b) - 
		double(a) * (square(phiX) + square(phiY) + square(phiZ) - square(h));
	if(quotient < 0.) cout << "0 ";
	else {
		phi = b + sqrt(quotient);
		phi /= double(a);
		grid[index].value = phi;
		if(grid[index].HeapPosition == -1) AddToHeap(index);
	    else                               UpdateHeap(index); 
	}
}

inline void FastMarch::CheckFront(double& phi, int& a, bool& flag, int index) 
{
	if(grid[index].DoneFlag == 1) {
		phi = grid[index].value;
		flag = 1;
		a++;
	}
}

inline void FastMarch::CheckBehind(double& phi, int& a, bool& flag, int index)
{
	if(grid[index].DoneFlag == 1) {
		if(!flag) { phi = grid[index].value; a++; }
		else phi = min(grid[index].value, phi);
        flag = 1;
	}
}

inline void FastMarch::CheckMax2(int& a, double& phi1, const double &phi2) {
	if(square((phi1 - phi2)*hInv) > 1.) { phi1 = 0; a = 1; }
}

inline void FastMarch::CheckMax3(int& a, bool& flag, double& phi1, 
                                 const double &phi2, const double &phi3) {
	if((square((phi1-phi2)*hInv) + square((phi1-phi3)*hInv)) > 1.) 
    {   phi1 = 0; a = 2; flag = 0;  }
}

void FastMarch::AddToHeap(int index) {
	FMHeap[heapSize] = index;
	grid[index].HeapPosition = heapSize;
	int j, i = heapSize;
	for(i; i > 0; i = j) {
		j = (i-1)/2;
		if (grid[ (FMHeap[i]) ].value < grid[ (FMHeap[j]) ].value ){
			FMHeap[i] = FMHeap[j];
			grid[ (FMHeap[j]) ].HeapPosition = i;
			FMHeap[j] = index;
			grid[index].HeapPosition = j;
		}
		else
			break;
	}
	heapSize++;
}

void FastMarch::UpdateHeap(int index) {
	int j, i = grid[index].HeapPosition;
	for(i; i > 0; i = j) {
		j = (i-1)/2;
		if (grid[ (FMHeap[i]) ].value < grid[ (FMHeap[j]) ].value ) {
			FMHeap[i] = FMHeap[j];
			grid[ (FMHeap[j]) ].HeapPosition = i;
			FMHeap[j] = index;
			grid[index].HeapPosition = j;
		}
		else
			break;
	}
}

void FastMarch::March() {
	static int x, y, z; 
	for(int index = PopHeap(); index != -1; index = PopHeap()) {
		if(grid[index].value > FASTMARCH_LIMIT) return;
		GIJK(index, x, y, z);
        if(grid[GI(x-1,y,z)].DoneFlag == 0) FindPhi(GI(x-1,y,z),x-1,y,z);
        if(grid[GI(x+1,y,z)].DoneFlag == 0) FindPhi(GI(x+1,y,z),x+1,y,z);
        if(grid[GI(x,y-1,z)].DoneFlag == 0) FindPhi(GI(x,y-1,z),x,y-1,z);
        if(grid[GI(x,y+1,z)].DoneFlag == 0) FindPhi(GI(x,y+1,z),x,y+1,z);
        if(grid[GI(x,y,z-1)].DoneFlag == 0) FindPhi(GI(x,y,z-1),x,y,z-1);
        if(grid[GI(x,y,z+1)].DoneFlag == 0) FindPhi(GI(x,y,z+1),x,y,z+1);
	}
}

int FastMarch::PopHeap() {
	if(heapSize == 0)
		return -1;
	int j, index = FMHeap[0];
	grid[index].DoneFlag = 1;
	heapSize--;
	FMHeap[0] = FMHeap[heapSize];
	grid[FMHeap[heapSize]].HeapPosition = 0;
	for(int i = 0; i < (heapSize-1); i = j) {
		int lc = 2*i+1;
		int rc = 2*i+2;
		double current = grid[ (FMHeap[i]) ].value;
		double lv, rv;
		if(lc < heapSize) {
			lv = grid[ (FMHeap[lc]) ].value;
			if(rc < heapSize) {
				rv = grid[ (FMHeap[rc]) ].value;
				if(lv > rv) {
					lc = rc;
					lv = rv;
				}
			}
			if(current > lv) {
				FMHeap[i] = FMHeap[lc];
				grid[ FMHeap[i] ].HeapPosition = i;
				FMHeap[lc] = FMHeap[heapSize];
				grid[ FMHeap[heapSize] ].HeapPosition = lc;
				j = lc;
			}
			else
				break;
		}
		else
			break;
	}
	return index;
}