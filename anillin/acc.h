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

	double** rho = new double*[NumXPoints];
	double** ux = new double*[NumXPoints];
	double** uy = new double*[NumXPoints];
	double** p = new double*[NumXPoints];
	double** Etas = new double*[NumXPoints];

	double** drho = new double*[NumXPoints];
	double** dux = new double*[NumXPoints];
	double** duy = new double*[NumXPoints];

	double** PhiX = new double*[NumXPoints];
	double** PhiY = new double*[NumXPoints];

	for(i=0; i<NumXPoints; i++)
	{
		rho[i] = new double[NumYPoints];
		ux[i] = new double[NumYPoints];
		uy[i] = new double[NumYPoints];

		p[i] = new double[NumYPoints];

		drho[i] = new double[NumYPoints];
		dux[i] = new double[NumYPoints];
		duy[i] = new double[NumYPoints];

		PhiX[i] = new double[NumYPoints];
		PhiY[i] = new double[NumYPoints];

		Etas[i] = new double[NumYPoints];
	}

	#include "inpt.h"

	double* Bndry_vx = new double[NumBndryPoints*20];
	double* Bndry_vy = new double[NumBndryPoints*20];

	double* Bndry_fx = new double[NumBndryPoints*20];
	double* Bndry_fy = new double[NumBndryPoints*20];

	//--------------------------------------------------------------------------------------------------------------------------------------------------------------
	double** q = new double* [NumXPoints];
	double** p1 = new double* [NumXPoints];
	double** px = new double* [NumXPoints];
	double** py = new double* [NumXPoints];

	double** uxProj = new double* [NumXPoints];
	double** uyProj = new double* [NumXPoints];

	for(i=0; i<NumXPoints; i++)
	{
		q[i] = new double[NumYPoints];
		p1[i] = new double[NumYPoints];
		px[i] = new double[NumYPoints];
		py[i] = new double[NumYPoints];

		uxProj[i] = new double[NumYPoints];
		uyProj[i] = new double[NumYPoints];
	}

	double* qHat = new double[2*(NumXPoints)*(NumYPoints)+1];

	double** aHat = new double*[NumXPoints];
	for(i=0; i<NumXPoints; i++)
	{
		aHat[i] = new double[NumYPoints];
	}

	double* gammaHat = new double[NumXPoints];
	double* omegaHat = new double[NumXPoints];

	const double Trshld = 1.0e-40;



	for(i=0; i<NumXPoints; i++)
	{
		gammaHat[i] = sin(2.0*Pi*i/NumXPoints)*( sqrt2/2+(1.0-sqrt2/2)*cos(2.0*Pi*i/NumXPoints) ) /h;
		omegaHat[i] = pow(cos(Pi*i/NumXPoints),2)*( 1.0-2.0*(1.0-2.0*sqrt2/Pi)*pow(sin(Pi*i/NumXPoints),2) );
	}

	for(i=0; i<NumXPoints; i++)
	{
		for(j=0; j<NumYPoints; j++)
		{
			aHat[i][j] = -(gammaHat[i]*gammaHat[i])*(omegaHat[j]*omegaHat[j])-(omegaHat[i]*omegaHat[i])*(gammaHat[j]*gammaHat[j]) ;

			if( i==0 && j==0 )
			{
				aHat[i][j] = 0;
				continue;
			}
			else if( i==NumXPoints/2  || j==NumXPoints/2 )
			{
				aHat[i][j] = 0;
				continue;
			}

			r = aHat[i][j] ;

			if( r <0 ) r = - r ;

			if( r < Trshld  )
			{
				aHat[i][j] = 0;
				continue;
			}

			aHat[i][j] = 1.0 / aHat[i][j] ;
		}
	}
	//--------------------------------------------------------------------------------------------------------------------------------------------------------------

	//--------------------------------------------------------------------------------------------------------------------------------------------------------------
	//
	double** xs = new double*[NumXPoints];
	double** ys = new double*[NumXPoints];

	double** Dist = new double*[NumXPoints];
	double** d_x = new double*[NumXPoints];
	double** d_y = new double*[NumXPoints];
	int** ClsdPnt = new int*[NumXPoints];
	int** ClsdPntId = new int*[NumXPoints];

	double* Vxs = new double[NumApicalPoints];
	double* Vys = new double[NumApicalPoints];
	double* Vs = new double[NumApicalPoints];
	
	double* Nxs = new double[NumApicalPoints];
	double* Nys = new double[NumApicalPoints];

	for(i=0; i<NumXPoints; i++)
	{
		xs[i] = new double[NumYPoints];
		ys[i] = new double[NumYPoints];
		
		Dist[i] = new double[NumYPoints];
		
		d_x[i] = new double[NumYPoints];
		d_y[i] = new double[NumYPoints];
		
		ClsdPnt[i] = new int[NumYPoints];
		ClsdPntId[i] = new int[NumYPoints];
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

	//==============================================================================================================================================================
	//
	for(i=0; i<NumXPoints; i++)
	{
		for(j=0; j<NumYPoints; j++)
		{
			ux[i][j] = 0;
			uy[i][j] = 0;

			rho[i][j] = 1.0;
		}
	}
	//==============================================================================================================================================================