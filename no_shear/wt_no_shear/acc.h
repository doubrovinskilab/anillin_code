	const int NumXPoints = 128*2;
	const int NumYPoints = NumXPoints;

	const double Lx = 128.0*2;
	const double Ly = Lx;
	const double dx = Lx/(NumXPoints-1);
	const double dy = Ly/(NumYPoints-1);

	int i=0, j=0, k=0, m=0, n=0;

	int idx1=0, idx2=0;

	double r = 0;
	double rx=0, ry=0, r2=0;
	double x1=0, x2=0, y1=0, y2=0;

	int temp1=0, temp2=0;

	#include "inpt.h"

	double* Bndry_fx = new double[NumBndryPoints*20];
	double* Bndry_fy = new double[NumBndryPoints*20];

	//--------------------------------------------------------------------------------------------------------------------------------------------------------------
	//
	double** xs = new double*[NumXPoints];
	double** ys = new double*[NumXPoints];

	for(i=0; i<NumXPoints; i++)
	{
		xs[i] = new double[NumYPoints];
		ys[i] = new double[NumYPoints];
	}

	for(i=0; i<NumXPoints; i++)
	{
		for(j=0; j<NumYPoints; j++)
		{
			xs[i][j] = i*dx;
			ys[i][j] = j*dy;
		}
	}

	double* ApicalBndryX = new double[NumApicalPoints];
	double* ApicalBndryY = new double[NumApicalPoints];

	for(i=0; i<NumApicalPoints; i++)
	{
		ApicalBndryX[i] = BndryX[ApicalPoints[i]];
		ApicalBndryY[i] = BndryY[ApicalPoints[i]];
	}
	//--------------------------------------------------------------------------------------------------------------------------------------------------------------