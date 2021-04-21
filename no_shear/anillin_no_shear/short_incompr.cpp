#include <iostream>
#include <fstream>
using namespace std;
#include <math.h>
#include <string.h>
#include <stdlib.h>

char* itoa(int val, int base);

int main()
{
	//==============================================================================================================================================================
	// Parameters

	// Young's modulus of cellular boundaries
	const double E = .00010 *12*2/3;

	// Contractile stress magnitude (stresses carried by cell boundaries scale with phi0)
	/*const*/ double phi0 = 0.0010 *2;

	// Time step
	const double dt = 1.0;
	//==============================================================================================================================================================

	//==============================================================================================================================================================
	// Auxiliaries
	#include "acc.h"

	double phi = 0.0 ;
	double fx=0, fy=0, f=0;

	//--------------------------------------------------------------------------------------------------------------------------------------------------------------
	// Initial total volume of the embryo (enclosed by the apical surface)
	double V0 = 0;
	for(i=0; i<NumApicalPoints; i++)               
	{
		int Idx = i;
		int Nxt = i+1;
		if( Nxt==NumApicalPoints ) Nxt = 0;

		V0 += (BndryX[Nxt]+BndryX[Idx])*(BndryY[Nxt]-BndryY[Idx])/4 - (BndryY[Nxt]+BndryY[Idx])*(BndryX[Nxt]-BndryX[Idx])/4 ;
	}
	//--------------------------------------------------------------------------------------------------------------------------------------------------------------

	double V = 0.0;
	//==============================================================================================================================================================

	//==============================================================================================================================================================
	// Time loop
	for(i=0; ;i++)
	{

		if( !(i/100.0-i/100) ) { cout<<"Iteration # "<<i<<endl; }

		//----------------------------------------------------------------------------------------------------------------------------------------------------------
		// Time-evolve active stress: zero initially, increases linearly until time T, remains constant thereafter
		double T = 2000*50;
		if( i*dt/T > 1.0 )
		{
			phi = phi0* 1.0 ;
		}
		else
		{
			phi = phi0* i*dt/T;
		}
		//----------------------------------------------------------------------------------------------------------------------------------------------------------

		//----------------------------------------------------------------------------------------------------------------------------------------------------------
		// "Internal" loop: relax the system to equilibrium for a given imposed distribution of active stress phi
		for(int ii=0; ii<1000; ii++)
		{
			//------------------------------------------------------------------------------------------------------------------------------------------------------
			// Initialize force density at Lagrangian boundaries to zero
			for(j=0; j<NumBndryPoints; j++)
			{
				Bndry_fx[j] = 0;
				Bndry_fy[j] = 0;
			}
			//------------------------------------------------------------------------------------------------------------------------------------------------------

			//------------------------------------------------------------------------------------------------------------------------------------------------------
			// Compute force density at the boundary: contribution from active ("myosin-generated") and elastic stress
			for(j=0; j<NumEdges; j++)
			{
				idx1 = Edge1[j];
				idx2 = Edge2[j];

				rx = BndryX[Edge2[j]] - BndryX[Edge1[j]];
				ry = BndryY[Edge2[j]] - BndryY[Edge1[j]];
				r2 = rx*rx + ry*ry;
				r = sqrt(r2);
				if(r2>0) { rx/=r; ry/=r; }

				//--------------------------------------------------------------------------------------------------------------------------------------------------
				// Elastic stress (on all but basal membranes)
				if( EdgeId[j]!=3 )
				{
					Bndry_fx[Edge1[j]] +=  Es[j]*(r-L0s[j])*rx;
					Bndry_fy[Edge1[j]] +=  Es[j]*(r-L0s[j])*ry;

					Bndry_fx[Edge2[j]] += -Es[j]*(r-L0s[j])*rx;
					Bndry_fy[Edge2[j]] += -Es[j]*(r-L0s[j])*ry;
				}
				//--------------------------------------------------------------------------------------------------------------------------------------------------

				//--------------------------------------------------------------------------------------------------------------------------------------------------
				// Active stress

				//
				// Constitutive stress at all boundaries
				/*
				Bndry_fx[Edge1[j]] +=  0.25*phi*rx;
				Bndry_fy[Edge1[j]] +=  0.25*phi*ry;

				Bndry_fx[Edge2[j]] += -0.25*phi*rx;
				Bndry_fy[Edge2[j]] += -0.25*phi*ry;
				*/
				//

				if( EdgeId[j]==1 ) // stress at the apical mesodermal membranes
				{
					Bndry_fx[idx1] +=  12.0*phi*rx;
					Bndry_fy[idx1] +=  12.0*phi*ry;

					Bndry_fx[idx2] += -12.0*phi*rx;
					Bndry_fy[idx2] += -12.0*phi*ry;
				}
				else if( EdgeId[j]==2 ) // stress at the lateral mesodermal membranes
				{
					Bndry_fx[idx1] +=  2*1.0*phi*rx;
					Bndry_fy[idx1] +=  2*1.0*phi*ry;

					Bndry_fx[idx2] += -2*1.0*phi*rx;
					Bndry_fy[idx2] += -2*1.0*phi*ry;
				}
				else if( EdgeId[j]!=3 ) // Constitutive pre-stress at all (except basal) membranes
				{
					Bndry_fx[idx1] +=  0.50*phi*rx;
					Bndry_fy[idx1] +=  0.50*phi*ry;

					Bndry_fx[idx2] += -0.50*phi*rx;
					Bndry_fy[idx2] += -0.50*phi*ry;
				}

				/*
				else if( EdgeId[j]==3 ) // basal pulled
				{
				Bndry_fx[idx1] +=  2.0*phi*rx;
				Bndry_fy[idx1] +=  2.0*phi*ry;

				Bndry_fx[idx2] += -2.0*phi*rx;
				Bndry_fy[idx2] += -2.0*phi*ry;

				}
				*/
				//--------------------------------------------------------------------------------------------------------------------------------------------------
			}
			//------------------------------------------------------------------------------------------------------------------------------------------------------

			//------------------------------------------------------------------------------------------------------------------------------------------------------
			// Soft constraint due to the vitelline membrane
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
			// Force of constraint to maintain total volume of the embryo (enclosed by the apical surface)
			V = 0;
			for(j=0; j<NumApicalPoints; j++)
			{
				int Idx = j;
				int Nxt = j+1;
					if( Nxt==NumApicalPoints ) Nxt = 0;

				V += (BndryX[Nxt]+BndryX[Idx])*(BndryY[Nxt]-BndryY[Idx])/4 - (BndryY[Nxt]+BndryY[Idx])*(BndryX[Nxt]-BndryX[Idx])/4 ;
			}

			for(j=0; j<NumApicalPoints; j++)
			{
				int Prv = j-1;
				if( Prv<0 ) Prv = NumApicalPoints-1;
				int Nxt = j+1;
				if( Nxt==NumApicalPoints ) Nxt = 0;

				Bndry_fx[j] += -(2.0e-4)*(BndryY[Nxt]-BndryY[Prv])*(V-V0);
				Bndry_fy[j] += -(2.0e-4)*(BndryX[Prv]-BndryX[Nxt])*(V-V0);
			}
			//------------------------------------------------------------------------------------------------------------------------------------------------------

			//------------------------------------------------------------------------------------------------------------------------------------------------------
			// Relax towards steady state
			for(j=0; j<NumBndryPoints; j++)
			{
				BndryX[j] += dt*Bndry_fx[j];
				BndryY[j] += dt*Bndry_fy[j];
			}
			//------------------------------------------------------------------------------------------------------------------------------------------------------
		}
		//----------------------------------------------------------------------------------------------------------------------------------------------------------

		//----------------------------------------------------------------------------------------------------------------------------------------------------------
		// Write output (to ./outpt folder)
		int qq = 200;
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
