/*
2D Immersed Boundary simulation of Drosophila ventral furrow formation
Initial conditions are determined from input files in ./meshfiles containing initial configuration of the cellular boundaries
Output is written to ./outpt
*/

#include <iostream>
#include <fstream>
using namespace std;
#include <math.h>
#include <string.h>

char* itoa(int val, int base);

#include "accss.h"

// void chorinProj(double** ux, double** uy, double** p, double** px, double** py, double** aHat, int NumXPoints, double h, double** q, 
//				double* qHat)

int main()
{
	//==============================================================================================================================================================
	// Parameters

	// Young's modulus of cellular boundaries
	const double E = .00010 *12*2/3;

	// Contractile stress magnitude (stressess carrid by cell boundaries scale with phi)
	const double phi = 0.0010 *2;

	// Cytoplasm density
	const double rho0 = 1.0;
	 // Cytoplasm compressibility
	const double gmma = 0.10;

	// Cytoplasm viscosity
	const double eta = 10* 1.0/6;

	// ("Small") parameter determining the rate of convergence to Stokes flow dynamics
	const double tau = 1.0*30;

	// Lattice spacing
	const double h = 1.0;

	// Length of edges discretizing solid boundaries
	const double L0 = 1.70/2 ;

	// Euler time-step
	const double dt = 3.0/3;
	//==============================================================================================================================================================

	//==============================================================================================================================================================
	// Accessories

	#include "acc.h"

	double** dx_ux = new double*[NumXPoints];
	double** dy_ux = new double*[NumXPoints];

	double** dx_uy = new double*[NumXPoints];
	double** dy_uy = new double*[NumXPoints];

	for(i=0; i<NumXPoints; i++)
	{
		dx_ux[i] = new double[NumYPoints];
		dy_ux[i] = new double[NumYPoints];

		dx_uy[i] = new double[NumYPoints];
		dy_uy[i] = new double[NumYPoints];
	}
	//==============================================================================================================================================================

	//==============================================================================================================================================================
	// Time loop
	for(i=0; ;i++)
	{

		if( !(i/100.0-i/100) ) { cout<<"Iteration # "<<i<<endl;
		}

		//----------------------------------------------------------------------------------------------------------------------------------------------------------
		//
		/*
		Determine the distribution of viscosity based on the profile of the embryonic surface
		Viscosity set to eta in the interior of the embryo, and to 1/10 that value outside
		Distribution is a function of closest distance between a point and polygon discretizing embryonic surface
		*/
		for(j=0; j<NumApicalPoints; j++)
		{
			ApicalBndryX[j] = BndryX[ApicalPoints[j]];
			ApicalBndryY[j] = BndryY[ApicalPoints[j]];
		}

		updateDistances(ApicalBndryX,ApicalBndryY,NumApicalPoints,xs,ys,NumXPoints,NumYPoints,Dist,d_x,d_y,ClsdPnt,ClsdPntId,Vxs,Vys,Vs,Nxs,Nys);

		for(j=0; j<NumXPoints; j++)
		{
			for(k=0; k<NumYPoints; k++)
			{
				dux[j][k] = 0;
				duy[j][k] = 0;

				drho[j][k] = 0;

				p[j][k] = gmma*(rho[j][k]-rho0);

				Etas[j][k] = eta*( (1.0/2+tanh(-Dist[j][k]/5.0)/2.0)*(1.0-0.10) + 0.10 );
			}
		}
		//----------------------------------------------------------------------------------------------------------------------------------------------------------

		//----------------------------------------------------------------------------------------------------------------------------------------------------------
		// Initialize force density on the fluid 
		for(j=0; j<NumXPoints; j++)
		{
			for(k=0; k<NumYPoints; k++)
			{
				PhiX[j][k] = 0;
				PhiY[j][k] = 0;
			}
		}
		//----------------------------------------------------------------------------------------------------------------------------------------------------------

		//----------------------------------------------------------------------------------------------------------------------------------------------------------
		// Initialize force density on the (solid) boundary
		for(j=0; j<NumBndryPoints; j++)
		{
			Bndry_fx[j] = 0;
			Bndry_fy[j] = 0;
		}
		//----------------------------------------------------------------------------------------------------------------------------------------------------------

		//----------------------------------------------------------------------------------------------------------------------------------------------------------
		// Compute force density on cell boundaries
		for(j=0; j<NumEdges; j++)
		{
			idx1 = Edge1[j];
			idx2 = Edge2[j];

			rx = BndryX[Edge2[j]] - BndryX[Edge1[j]];
			ry = BndryY[Edge2[j]] - BndryY[Edge1[j]];
			r2 = rx*rx + ry*ry;
			r = sqrt(r2);
			if(r2>0) { rx/=r; ry/=r; }

			//------------------------------------------------------------------------------------------------------------------------------------------------------
			// Force contribution due to elastic stresses

			// if( EdgeId[j]!=3 )
			{
				Bndry_fx[Edge1[j]] +=  Es[j]*(r-L0s[j])*rx;
				Bndry_fy[Edge1[j]] +=  Es[j]*(r-L0s[j])*ry;

				Bndry_fx[Edge2[j]] += -Es[j]*(r-L0s[j])*rx;
				Bndry_fy[Edge2[j]] += -Es[j]*(r-L0s[j])*ry;
			}
			//------------------------------------------------------------------------------------------------------------------------------------------------------

			//------------------------------------------------------------------------------------------------------------------------------------------------------
			// Force contribution due to active (myosin-generated) stresses

			//
			// Constitutive prestress (all edges)
			/*
			Bndry_fx[Edge1[j]] +=  0.25*phi*rx;
			Bndry_fy[Edge1[j]] +=  0.25*phi*ry;

			Bndry_fx[Edge2[j]] += -0.25*phi*rx;
			Bndry_fy[Edge2[j]] += -0.25*phi*ry;
			*/
			//

			if( EdgeId[j]==1 ) // apical edges
			{
				Bndry_fx[idx1] +=  12.0*phi*rx;
				Bndry_fy[idx1] +=  12.0*phi*ry;

				Bndry_fx[idx2] += -12.0*phi*rx;
				Bndry_fy[idx2] += -12.0*phi*ry;
			}
			else if( EdgeId[j]==2 ) // lateral edges
			{
				Bndry_fx[idx1] +=  1.0*phi*rx;
				Bndry_fy[idx1] +=  1.0*phi*ry;

				Bndry_fx[idx2] += -1.0*phi*rx;
				Bndry_fy[idx2] += -1.0*phi*ry;
			}
			else if( EdgeId[j]!=3 ) // Constitutive prestress (all but basal edges)
			{
				Bndry_fx[idx1] +=  0.50*phi*rx;
				Bndry_fy[idx1] +=  0.50*phi*ry;

				Bndry_fx[idx2] += -0.50*phi*rx;
				Bndry_fy[idx2] += -0.50*phi*ry;
			}

			else if( EdgeId[j]==3 ) // basal edges
			{
				Bndry_fx[idx1] +=  2.0*phi*rx;
				Bndry_fy[idx1] +=  2.0*phi*ry;

				Bndry_fx[idx2] += -2.0*phi*rx;
				Bndry_fy[idx2] += -2.0*phi*ry;

			}
			//------------------------------------------------------------------------------------------------------------------------------------------------------
		}
		//----------------------------------------------------------------------------------------------------------------------------------------------------------

		//----------------------------------------------------------------------------------------------------------------------------------------------------------
		// (Soft) constraint from the vitelline membrane
		for(j=0; j<NumBndryPoints; j++)
		{
			rx = BndryX[j] - 128.0;
			ry = BndryY[j] - 128.0;

			r2 = rx*rx + ry*ry;

			r = sqrt(r2);

			if( r>1.5*75.0+2.0 )
			{
				Bndry_fx[j] += -(0.0040)*rx/r ;
				Bndry_fy[j] += -(0.0040)*ry/r ;
			}
		}
		//------------------------------------------------------------------------------------------------------------------------------------------------------

		//------------------------------------------------------------------------------------------------------------------------------------------------------
		// Correct incompressibility using Chorin projection
		chorinProj(ux,uy,p1,px,py,aHat,NumXPoints,h,q,qHat);

		for(j=0; j<NumXPoints; j++)
		{
			for(k=0; k<NumYPoints; k++)
			{
				uxProj[j][k] = ux[j][k] - px[j][k] ;
				uyProj[j][k] = uy[j][k] - py[j][k] ;
			}
		}
		//------------------------------------------------------------------------------------------------------------------------------------------------------

		//------------------------------------------------------------------------------------------------------------------------------------------------------
		// Update force densities on cellular boundaries and on the fluid
		for(j=0; j<NumBndryPoints; j++)
		{
			int IdxX = floor(BndryX[j]);
			int IdxY = floor(BndryY[j]);

			Bndry_vx[j] = 0;
			Bndry_vy[j] = 0;

			for(m=-3; m<=3; m++)
			{
				for(n=-3; n<=3; n++)
				{
					int m_ = IdxX + m;
					int n_ = IdxY + n;

					if(m_<0)
						m_ += NumXPoints;
					if( m_>(NumXPoints-1) ) 
						m_ += NumXPoints;

					if(n_<0)
						n_ += NumYPoints;
					if( n_>(NumYPoints-1) ) 
						n_ += NumYPoints;

					rx = xs[m_][n_] - BndryX[j];
					ry = ys[m_][n_] - BndryY[j];

					Bndry_vx[j] += uxProj[m_][n_]*Delta(rx,ry,h);
					Bndry_vy[j] += uyProj[m_][n_]*Delta(rx,ry,h);

					PhiX[m_][n_] += Bndry_fx[j] * Delta(rx,ry,h);
					PhiY[m_][n_] += Bndry_fy[j] * Delta(rx,ry,h);
				}
			}
		}
		//------------------------------------------------------------------------------------------------------------------------------------------------------

		//----------------------------------------------------------------------------------------------------------------------------------------------------------
		// Pre-compute spatial derivatives of velocity
		for(j=0; j<NumXPoints; j++)
		{
			for(k=0; k<NumYPoints; k++)
			{
				if( j==0 )
				{
					dx_ux[j][k] = (ux[j+1][k]-ux[j][k])/dx;
					dx_uy[j][k] = (uy[j+1][k]-uy[j][k])/dx;
				}
				else if( j==NumXPoints-1 )
				{
					dx_ux[j][k] = (ux[j][k]-ux[j-1][k])/dx;
					dx_uy[j][k] = (uy[j][k]-uy[j-1][k])/dx;
				}
				else
				{
					dx_ux[j][k] = (ux[j+1][k]-ux[j-1][k])/(2*dx);
					dx_uy[j][k] = (uy[j+1][k]-uy[j-1][k])/(2*dx);
				}

				if( k==0 )
				{
					dy_ux[j][k] = (ux[j][k+1]-ux[j][k])/dy;
					dy_uy[j][k] = (uy[j][k+1]-uy[j][k])/dy;
				}
				else if( k==NumYPoints-1 )
				{
					dy_ux[j][k] = (ux[j][k]-ux[j][k-1])/dy;
					dy_uy[j][k] = (uy[j][k]-uy[j][k-1])/dy;
				}
				else
				{
					dy_ux[j][k] = (ux[j][k+1]-ux[j][k-1])/(2*dy);
					dy_uy[j][k] = (uy[j][k+1]-uy[j][k-1])/(2*dy);
				}
			}
		}
		//----------------------------------------------------------------------------------------------------------------------------------------------------------

		//----------------------------------------------------------------------------------------------------------------------------------------------------------
		// Stokes flow dynamics
		for(j=0+1; j<NumXPoints-1; j++)
		{
			for(k=0+1; k<NumYPoints-1; k++)
			{
				dux[j][k] += -(p[j+1][k]-p[j-1][k])/(2*dx);
				duy[j][k] += -(p[j][k+1]-p[j][k-1])/(2*dy);

				dux[j][k] += (Etas[j+1][k]*dx_ux[j+1][k]-Etas[j-1][k]*dx_ux[j-1][k])/(2*dx);
				dux[j][k] += (Etas[j][k+1]*dy_ux[j][k+1]-Etas[j][k-1]*dy_ux[j][k-1])/(2*dy);

				duy[j][k] += (Etas[j+1][k]*dx_uy[j+1][k]-Etas[j-1][k]*dx_uy[j-1][k])/(2*dx);
				duy[j][k] += (Etas[j][k+1]*dy_uy[j][k+1]-Etas[j][k-1]*dy_uy[j][k-1])/(2*dy);


				dux[j][k] += (Etas[j+1][k]*dx_ux[j+1][k]-Etas[j-1][k]*dx_ux[j-1][k])/(2*dx);
				dux[j][k] += (Etas[j][k+1]*dx_uy[j][k+1]-Etas[j][k-1]*dx_uy[j][k-1])/(2*dy);

				duy[j][k] += (Etas[j+1][k]*dy_ux[j+1][k]-Etas[j-1][k]*dy_ux[j-1][k])/(2*dx);
				duy[j][k] += (Etas[j][k+1]*dy_uy[j][k+1]-Etas[j][k-1]*dy_uy[j][k-1])/(2*dy);

				dux[j][k] += PhiX[j][k];
				duy[j][k] += PhiY[j][k];

				drho[j][k] += -(ux[j+1][k]*rho[j+1][k]-ux[j-1][k]*rho[j-1][k])/(2*dx);
				drho[j][k] += -(uy[j][k+1]*rho[j][k+1]-uy[j][k-1]*rho[j][k-1])/(2*dy);
			}
		}
		//----------------------------------------------------------------------------------------------------------------------------------------------------------

		//----------------------------------------------------------------------------------------------------------------------------------------------------------
		// Zero flux of density imposed at the domain boundaries
		for(j=0; j<NumXPoints; j++)
		{
			drho[j][0] += -(rho[j][0]*uy[j][0]+rho[j][1]*uy[j][1])/(2*dy);
			drho[j][NumYPoints-1] += (rho[j][NumYPoints-1]*uy[j][NumYPoints-1]+rho[j][NumYPoints-1-1]*uy[j][NumYPoints-1-1])/(2*dy);

			if( j>0 && j<NumXPoints-1 )
			{
				drho[j][0] += -(ux[j+1][0]*rho[j+1][0]-ux[j-1][0]*rho[j-1][0])/(2*dx);
				drho[j][NumYPoints-1] += -(ux[j+1][NumYPoints-1]*rho[j+1][NumYPoints-1]-ux[j-1][NumYPoints-1]*rho[j-1][NumYPoints-1])/(2*dx);
			}
		}

		for(k=0; k<NumYPoints; k++)
		{
			drho[0][k] += -(rho[0][k]*ux[0][k]+rho[1][k]*ux[1][k])/(2*dx);
			drho[NumXPoints-1][k] += (rho[NumXPoints-1][k]*ux[NumXPoints-1][k]+rho[NumXPoints-1-1][k]*ux[NumXPoints-1-1][k])/(2*dx);

			if( k>0 && k<NumYPoints-1 )
			{
				drho[0][k] += -(uy[0][k+1]*rho[0][k+1]-uy[0][k-1]*rho[0][k-1])/(2*dy);
				drho[NumXPoints-1][k] += -(uy[NumXPoints-1][k+1]*rho[NumXPoints-1][k+1]-uy[NumXPoints-1][k-1]*rho[NumXPoints-1][k-1])/(2*dy);
			}
		}
		//----------------------------------------------------------------------------------------------------------------------------------------------------------

		//----------------------------------------------------------------------------------------------------------------------------------------------------------
		// Euler time-step
		for(j=0; j<NumXPoints; j++)
		{
			for(k=0; k<NumYPoints; k++)
			{
				rho[j][k] += dt*drho[j][k];

				ux[j][k] += dt*dux[j][k]/tau;
				uy[j][k] += dt*duy[j][k]/tau;
			}
		}

		for(j=0; j<NumBndryPoints; j++)
		{
			BndryX[j] += dt*Bndry_vx[j];
			BndryY[j] += dt*Bndry_vy[j];
		}
		//----------------------------------------------------------------------------------------------------------------------------------------------------------

		//----------------------------------------------------------------------------------------------------------------------------------------------------------
		// If length of an edge exceeds 2*L0, split it in two halves
		temp2 = NumEdges;
		for(j=0; j<temp2; j++)
		{
			idx1 = Edge1[j];
			idx2 = Edge2[j];

			rx = BndryX[idx2] - BndryX[idx1];
			ry = BndryY[idx2] - BndryY[idx1];
			r2 = rx*rx + ry*ry;
			r = sqrt(r2);

			if( r>L0*2 )
			{
				BndryX[NumBndryPoints] = (BndryX[idx2] + BndryX[idx1])/2;
				BndryY[NumBndryPoints] = (BndryY[idx2] + BndryY[idx1])/2;

				NumBndryPoints++;


				temp1 = Edge2[j];
				Edge2[j] = NumBndryPoints-1;
				Edge1[NumEdges] = NumBndryPoints-1;
				Edge2[NumEdges] = temp1;
				EdgeId[NumEdges] = EdgeId[j];
				L0s[j] /= 2.0;
				Es[j] *= 2.0;
				L0s[NumEdges] = L0s[j];
				Es[NumEdges] = Es[j];

				NumEdges++;
			}
		}
		//----------------------------------------------------------------------------------------------------------------------------------------------------------

		//----------------------------------------------------------------------------------------------------------------------------------------------------------
		// Write output
		int qq = 2000;
		if( !(i/double(qq)-i/qq) )
		{
			char str1[280]="./outpt/out";

			char str3[280]=".m";

			strncat(str1,itoa(i/qq,10),10);
			strncat(str1,str3,10);

			ofstream * file = new ofstream(str1);

			file->precision(20);

			(*file)<<"xs=[];"<<endl;
			for(j=0; j<NumXPoints; j+=1)
			{
				(*file)<<"xs=[xs; ";
				for(k=0; k<NumYPoints; k+=1)
				{
					(*file)<<xs[j][k]<<" ";
				}
				(*file)<<"];"<<endl;
			}

			(*file)<<"ys=[];"<<endl;
			for(j=0; j<NumXPoints; j+=1)
			{
				(*file)<<"ys=[ys; ";
				for(k=0; k<NumYPoints; k+=1)
				{
					(*file)<<ys[j][k]<<" ";
				}
				(*file)<<"];"<<endl;
			}

			(*file)<<"Dist=[];"<<endl;
			for(j=0; j<NumXPoints; j+=1)
			{
				(*file)<<"Dist=[Dist; ";
				for(k=0; k<NumYPoints; k+=1)
				{
					(*file)<<Dist[j][k]<<" ";
				}
				(*file)<<"];"<<endl;
			}

			(*file)<<"rho=[];"<<endl;
			for(j=0; j<NumXPoints; j+=1)
			{
				(*file)<<"rho=[rho; ";
				for(k=0; k<NumYPoints; k+=1)
				{
					(*file)<<rho[j][k]<<" ";
				}
				(*file)<<"];"<<endl;
			}

			(*file)<<"ux=[];"<<endl;
			for(j=0; j<NumXPoints; j+=1)
			{
				(*file)<<"ux=[ux; ";
				for(k=0; k<NumYPoints; k+=1)
				{
					(*file)<<ux[j][k]<<" ";
				}
				(*file)<<"];"<<endl;
			}

			(*file)<<"uy=[];"<<endl;
			for(j=0; j<NumXPoints; j+=1)
			{
				(*file)<<"uy=[uy; ";
				for(k=0; k<NumYPoints; k+=1)
				{
					(*file)<<uy[j][k]<<" ";
				}
				(*file)<<"];"<<endl;
			}

			(*file)<<"BndryX=[";
			for(j=0; j<NumBndryPoints; j++)
			{
				(*file)<<BndryX[j]<<" ";
			}
			(*file)<<"];"<<endl;

			(*file)<<"BndryY=[";
			for(j=0; j<NumBndryPoints; j++)
			{
				(*file)<<BndryY[j]<<" ";
			}
			(*file)<<"];"<<endl;

			(*file)<<"Edge1=[";
			for(j=0; j<NumEdges; j++)
			{
				(*file)<<Edge1[j]+1<<" ";
			}
			(*file)<<"];"<<endl;

			(*file)<<"Edge2=[";
			for(j=0; j<NumEdges; j++)
			{
				(*file)<<Edge2[j]+1<<" ";
			}
			(*file)<<"];"<<endl;

                        (*file)<<"EdgeId=[";
                        for(j=0; j<NumEdges; j++)   
                        {
                                (*file)<<EdgeId[j]<<" ";
                        }
                        (*file)<<"];"<<endl;

			file->close();
			delete file;
		}
		//----------------------------------------------------------------------------------------------------------------------------------------------------------
	}
	//==============================================================================================================================================================

}

char* itoa(int val, int base) {
     static char buf[32] = {0};
     int i = 30;
     	if(val==0)
     	{
     	buf[i]='0';
	return &buf[i];
	}
     for (; val && i; --i, val /= base)
         buf[i] = "0123456789abcdef"[val % base];
     return &buf[i+1];
}
