#include "rlft3.h"

const double Pi = 3.1415926535897932385;
const double sqrt2 = sqrt(2.0);

inline double gamma(int i, double h)
{
	if(i==0) return 0;

	if( (i==1) ) return -1.0/(2*h)*sqrt2/2 ;
	if( (i==2) ) return -1.0/(4*h)*(1-sqrt2/2) ;

	if( (i==-1) ) return 1.0/(2*h)*sqrt2/2 ;
	if( (i==-2) ) return 1.0/(4*h)*(1-sqrt2/2) ;

	return 0;
}

inline double omega(int i, double h)
{
	if(i==0) return 1.0/4 + sqrt2/2/Pi;

	if( (i==1) || (i==-1) ) return 1.0/4 ;

	if( (i==2) || (i==-2) ) return 1.0/8-sqrt2/4/Pi ;

	return 0;
}

void chorinProj(double** ux, double** uy, double** p, double** px, double** py, double** aHat, int NumXPoints, double h, double** q, 
				double* qHat)
{
	int i=0, j=0, k=0, m=0;
	int k_=0, m_=0;

	unsigned long* nn = new unsigned long[2+1];

	int NumYPoints = NumXPoints ;

	nn[1] = NumXPoints;
	nn[2] = NumYPoints;

	for(i=0; i<NumXPoints; i++)
	{
		for(j=0; j<NumYPoints; j++)
		{
			q[i][j] = 0;

			for(k=-2; k<3; k++)
			{
				for(m=-2; m<3; m++)
				{
					k_ = i+k;
					m_ = j+m;

					if(k_<0) k_ += (NumXPoints-1*0);
					if(k_>(NumXPoints-1)) k_ += -(NumXPoints-1*0);

					if(m_<0) m_ += (NumYPoints-1*0);
					if(m_>(NumYPoints-1)) m_ += -(NumYPoints-1*0);

					q[i][j] += ux[k_][m_]*gamma(-k,h)*omega(-m,h) + uy[k_][m_]*omega(-k,h)*gamma(-m,h) ;
				}
			}
		}
	}

	for(i=0; i<NumXPoints; i++)
	{
		for(j=0; j<NumYPoints; j++)
		{
			qHat[ i*(2*NumYPoints) + (2*j+1)] = q[i][j];
			qHat[ i*(2*NumYPoints) + (2*j+1) +1] = 0;
		}
	}

	fourn(qHat,nn,2,-1); // FFT

	for(i=0; i<NumXPoints; i++)
	{
		for(j=0; j<NumYPoints; j++)
		{
			qHat[ i*(2*NumYPoints) + (2*j+1)] *= aHat[i][j];
			qHat[ i*(2*NumYPoints) + (2*j+1) +1] *= aHat[i][j];
		}
	}

	fourn(qHat,nn,2,1);  // inverse FFT

	for(i=0; i<NumXPoints; i++)
	{
		for(j=0; j<NumYPoints; j++)
		{
			p[i][j] = qHat[ i*(2*NumYPoints) + (2*j+1)] /(NumXPoints*NumYPoints);
		}
	}

	for(i=0; i<NumXPoints; i++)
	{
		for(j=0; j<NumYPoints; j++)
		{
			px[i][j] = 0;
			py[i][j] = 0;

			for(k=-2; k<3; k++)
			{
				for(m=-2; m<3; m++)
				{
					k_ = i+k;
					m_ = j+m;

					if(k_<0) k_ += (NumXPoints-1*0);
					if(k_>(NumXPoints-1)) k_ += -(NumXPoints-1*0);

					if(m_<0) m_ += (NumYPoints-1*0);
					if(m_>(NumYPoints-1)) m_ += -(NumYPoints-1*0);

					// q[i][j] += ux[k_][m_]*gamma(-k,h)*omega(-m,h) + uy[k_][m_]*omega(-k,h)*gamma(-m,h) ;

					px[i][j] += p[k_][m_]*gamma(-k,h)*omega(-m,h) ;
					py[i][j] += p[k_][m_]*omega(-k,h)*gamma(-m,h) ;
				}
			}
		}
	}

	delete [] nn;
}

/*
Peskin discrete delta function
*/
double Delta(double x, double y, double h)
{
	if( (x*x<=4) && (y*y<4) )
		return ( 1.0+cos(Pi*x/2/h) )*( 1.0+cos(Pi*y/2/h) )/(4*h)/(4*h);

	return 0;
}

inline double sign(double d)
{
	if(d<0) return -1;
	else return 1;
}


double round(double d)
{
	if(d>0)
	{
		double cl = ceil(d);
		if( cl-d>0.5 ) return cl-1;
		else
			return cl;
	}
	else if(d<0)
	{
		double ps = -d;
		
		double cl = ceil(ps);
		if( cl-ps>0.5 ) return -(cl-1);
		else
			return -cl;
	}
	else
	{
		return 0;
	}
}

/*
Closest distance transform.

Inputs:
BndryX, BndryY: coordinates of an oriented (counterclockwise orientation) (discretized, closed) curve, i.e. a polygon, with NumBndryPoints nodes
xs, ys: coordinates of a regular grid of dimensions NumXPoints x NumYPoints

Output:
Dist: array of dimensions NumXPoints x NumYPoints. Element with indices (i,j) contains closest distance from (xs[i][j],y[i][j]) to the discretized boundary


d_x, d_y, ClsdPnt, ClsdPntId: auxiliary arrays of dimensions NumXPoints x NumYPoints
Vxs, Vys, Vs, Nxs, Nys: auxiliary arrays of dimension NumBndryPoints
*/
void updateDistances(double * BndryX, double * BndryY, int NumBndryPoints, double * * xs, double * * ys, int NumXPoints, int NumYPoints,
 double * * Dist, double * * d_x, double * * d_y, int * * ClsdPnt, int * * ClsdPntId,  double * Vxs, double * Vys, double * Vs, double * Nxs, double * Nys)
{
	int i=0, j=0, k=0;
	
	double Lx = xs[NumXPoints-1][0]-xs[0][0];
	double Ly = ys[0][NumYPoints-1]-ys[0][0];
	
	double dx = xs[1][0]-xs[0][0];
	double dy = ys[0][1]-ys[0][0];
	
	for(i=0; i<NumXPoints; i++)
	{
		for(j=0; j<NumYPoints; j++)
		{
			Dist[i][j] = 4.0*Lx;
			
			ClsdPntId[i][j] = 0;
		}
	}

	for(i=0; i<NumBndryPoints; i++)
	{
		int PrvPnt = i-1;
		if(PrvPnt<0) PrvPnt += NumBndryPoints;
		int NxtPnt = i+1;
		if(NxtPnt>(NumBndryPoints-1)) NxtPnt += -NumBndryPoints;
		
		double px = BndryX[i];
		double py = BndryY[i];

		double px1 = BndryX[i]-xs[0][0];
		double py1 = BndryY[i]-ys[0][0];

		double vx1 = BndryX[i]-BndryX[PrvPnt];
		double vy1 = BndryY[i]-BndryY[PrvPnt];
		double v1 = sqrt(vx1*vx1+vy1*vy1);
		
		double vx2 = BndryX[NxtPnt]-BndryX[i];
		double vy2 = BndryY[NxtPnt]-BndryY[i];
		double v2 = sqrt(vx2*vx2+vy2*vy2);
		
		double nx1 = vy1/v1;
		double ny1 = -vx1/v1;
		
		double nx2 = vy2/v2;
		double ny2 = -vx2/v2;
		
		if( nx1*ny2-nx2*ny1>0 ) // convex
		{
			if( nx1>1.0/sqrt2-1e-8 ) 	// I
			{
				double eps = ( ceil(px/dx)*dx-px );					
				
				for(j=0; ;j++)
				{
					int CrntX = int(round( (px1 + sign(nx1)*(eps/nx1)*nx1 + sign(nx1)*j*dx)/dx ));
					int CrntY = int(ceil( (py1+ sign(nx1)*(eps/nx1)*ny1 + sign(nx1)*j*dx*ny1/nx1)/dy ));
					
					if(CrntX>NumXPoints-1) break;
					
					for(k=0; ;k++)
					{
						if(CrntY+k<0) continue;
						if(CrntY+k>NumYPoints-1) break;
						if( (xs[CrntX][CrntY+k]-px)*ny2-(ys[CrntX][CrntY+k]-py)*nx2 < 0 ) break;
						
						double dst = sqrt( (xs[CrntX][CrntY+k]-px)*(xs[CrntX][CrntY+k]-px)+(ys[CrntX][CrntY+k]-py)*(ys[CrntX][CrntY+k]-py) );
						if( dst*dst<Dist[CrntX][CrntY+k]*Dist[CrntX][CrntY+k] ) {
							Dist[CrntX][CrntY+k] = dst;
							d_x[CrntX][CrntY+k] = (xs[CrntX][CrntY+k]-px)/dst;
							d_y[CrntX][CrntY+k] = (ys[CrntX][CrntY+k]-py)/dst;
							ClsdPnt[CrntX][CrntY+k] = i;
							ClsdPntId[CrntX][CrntY+k] = -1;
						}
					}
				}
				if(nx2<0.0)
				{
					eps = -( floor(px/dx)*dx-px );
					
					for(j=0; ;j++)
					{
						int CrntX = int(round( (px1 + sign(nx2)*(eps/nx2)*nx2 + sign(nx2)*j*dx)/dx ));
						int CrntY = int(ceil( (py1 + sign(nx2)*(eps/nx2)*ny2 + sign(nx2)*j*dx*ny2/nx2)/dy ));
						
						if(CrntX<0 || CrntX>NumXPoints-1) break;
						
						for(k=0; ;k++)
						{
							if(CrntY+k<0) continue;
							if( CrntY+k>NumYPoints-1 ) break;
							if( (xs[CrntX][CrntY+k]-px)*ny2-(ys[CrntX][CrntY+k]-py)*nx2 < 0 ) break;
							
							double dst = sqrt( (xs[CrntX][CrntY+k]-px)*(xs[CrntX][CrntY+k]-px)+(ys[CrntX][CrntY+k]-py)*(ys[CrntX][CrntY+k]-py) );
							if( dst*dst<Dist[CrntX][CrntY+k]*Dist[CrntX][CrntY+k] ) {
								Dist[CrntX][CrntY+k] = dst;
								d_x[CrntX][CrntY+k] = (xs[CrntX][CrntY+k]-px)/dst;
								d_y[CrntX][CrntY+k] = (ys[CrntX][CrntY+k]-py)/dst;
								ClsdPnt[CrntX][CrntY+k] = i;
								ClsdPntId[CrntX][CrntY+k] = -1;
							}
						}
					}
				}
			}
			else if(nx1<-1.0/sqrt2+1e-8)	// II
			{
				double eps = -( floor(px/dx)*dx-px );
				
				for(j=0; ;j++)
				{
					int CrntX = int(round( (px1 + sign(nx1)*(eps/nx1)*nx1 + sign(nx1)*j*dx)/dx ));
					int CrntY = int(floor( (py1 + sign(nx1)*(eps/nx1)*ny1 + sign(nx1)*j*dx*ny1/nx1)/dy ));
					
					if(CrntX<0) break;
					
					for(k=0; ;k--)
					{
						if(CrntY+k>(NumYPoints-1)) continue;
						if(CrntY+k<0) break;
						if( (xs[CrntX][CrntY+k]-px)*ny2-(ys[CrntX][CrntY+k]-py)*nx2 < 0 ) break;
						
						double dst = sqrt( (xs[CrntX][CrntY+k]-px)*(xs[CrntX][CrntY+k]-px)+(ys[CrntX][CrntY+k]-py)*(ys[CrntX][CrntY+k]-py) );
						if( dst*dst<Dist[CrntX][CrntY+k]*Dist[CrntX][CrntY+k] ) { 
							Dist[CrntX][CrntY+k] = dst;
							d_x[CrntX][CrntY+k] = (xs[CrntX][CrntY+k]-px)/dst;
							d_y[CrntX][CrntY+k] = (ys[CrntX][CrntY+k]-py)/dst;
							ClsdPnt[CrntX][CrntY+k] = i;
							ClsdPntId[CrntX][CrntY+k] = -1;
						}
					}
				}
				if(nx2>0.0)
				{
					eps = ( ceil(px/dx)*dx-px );
					
					for(j=0; ;j++)
					{
						int CrntX = int(round( (px1 + sign(nx2)*(eps/nx2)*nx2 + sign(nx2)*j*dx)/dx ));
						int CrntY = int(floor( (py1 + sign(nx2)*(eps/nx2)*ny2 + sign(nx2)*j*dx*ny2/nx2)/dy ));
						
						if(CrntX<0 || CrntX>NumXPoints-1) break;
						
						for(k=0; ;k--)
						{
							if( CrntY+k>NumYPoints-1 ) continue;
							if( CrntY+k<0 ) break;
							if( (xs[CrntX][CrntY+k]-px)*ny2-(ys[CrntX][CrntY+k]-py)*nx2 < 0 ) break;
							
							double dst = sqrt( (xs[CrntX][CrntY+k]-px)*(xs[CrntX][CrntY+k]-px)+(ys[CrntX][CrntY+k]-py)*(ys[CrntX][CrntY+k]-py) );
							if( dst*dst<Dist[CrntX][CrntY+k]*Dist[CrntX][CrntY+k] ) { 
								Dist[CrntX][CrntY+k] = dst;
								d_x[CrntX][CrntY+k] = (xs[CrntX][CrntY+k]-px)/dst;
								d_y[CrntX][CrntY+k] = (ys[CrntX][CrntY+k]-py)/dst;
								ClsdPnt[CrntX][CrntY+k] = i;
								ClsdPntId[CrntX][CrntY+k] = -1;
							}
						}
					}
				}
			}
			else if(ny1>=1.0/sqrt2-1e-8)	// III
			{
				double eps = ( ceil(py/dy)*dy-py );
				
				for(j=0; ;j++)
				{
					int CrntY = int(round( (py1+sign(ny1)*(eps/ny1)*ny1 + sign(ny1)*j*dy)/dy ));
					int CrntX = int(floor( (px1+sign(ny1)*(eps/ny1)*nx1 + sign(ny1)*j*dy*nx1/ny1)/dx ));
					
					if(CrntY>NumYPoints-1) break;
					
					for(k=0; ;k--)
					{
						if( CrntX+k>NumXPoints-1 ) continue;
						if( CrntX+k<0 ) break;
						if( (xs[CrntX+k][CrntY]-px)*ny2-(ys[CrntX+k][CrntY]-py)*nx2 < 0 ) break;
						
						double dst = sqrt( (xs[CrntX+k][CrntY]-px)*(xs[CrntX+k][CrntY]-px)+(ys[CrntX+k][CrntY]-py)*(ys[CrntX+k][CrntY]-py) );
						if( dst*dst<Dist[CrntX+k][CrntY]*Dist[CrntX+k][CrntY] ) { 
							Dist[CrntX+k][CrntY] = dst;
							d_x[CrntX+k][CrntY] = (xs[CrntX+k][CrntY]-px)/dst;
							d_y[CrntX+k][CrntY] = (ys[CrntX+k][CrntY]-py)/dst;
							ClsdPnt[CrntX+k][CrntY] = i;
							ClsdPntId[CrntX+k][CrntY] = -1;
						}
					}
				}
				if(ny2<0.0)
				{
					eps = -( floor(py/dy)*dy-py );
					
					for(j=0; ;j++)
					{
						int CrntY = int(round( (py1+sign(ny2)*(eps/ny2)*ny2 + sign(ny2)*j*dy)/dy ));
						int CrntX = int(floor( (px1+sign(ny2)*(eps/ny2)*nx2 + sign(ny2)*j*dy*nx2/ny2)/dx ));
						
						if(CrntY<0 || CrntY>NumYPoints-1) break;
						
						for(k=0; ;k--)
						{
							if( CrntX+k>NumXPoints-1 ) continue;
							if( CrntX+k<0 ) break;
							if( (xs[CrntX+k][CrntY]-px)*ny2-(ys[CrntX+k][CrntY]-py)*nx2 < 0 ) break; //continue;
							
							double dst = sqrt( (xs[CrntX+k][CrntY]-px)*(xs[CrntX+k][CrntY]-px)+(ys[CrntX+k][CrntY]-py)*(ys[CrntX+k][CrntY]-py) );
							if( dst*dst<Dist[CrntX+k][CrntY]*Dist[CrntX+k][CrntY] ) {
								Dist[CrntX+k][CrntY] = dst;
								d_x[CrntX+k][CrntY] = (xs[CrntX+k][CrntY]-px)/dst;
								d_y[CrntX+k][CrntY] = (ys[CrntX+k][CrntY]-py)/dst;
								ClsdPnt[CrntX+k][CrntY] = i;
								ClsdPntId[CrntX+k][CrntY] = -1;
							}
						}
					}
				}
			}
			else	// IV
			{
				double eps = -( floor(py/dy)*dy-py );
				
				for(j=0; ;j++)
				{
					int CrntY = int(round( (py1+sign(ny1)*(eps/ny1)*ny1 + sign(ny1)*j*dy)/dy ));
					int CrntX = int(ceil( (px1+sign(ny1)*(eps/ny1)*nx1 + sign(ny1)*j*dy*nx1/ny1)/dx ));
					
					if(CrntY<0) break;
					
					for(k=0; ;k++)
					{
						if( CrntX+k<0 ) continue;
						if(  CrntX+k>NumXPoints-1) break;
						if( (xs[CrntX+k][CrntY]-px)*ny2-(ys[CrntX+k][CrntY]-py)*nx2 < 0 )break;
						
						double dst = sqrt( (xs[CrntX+k][CrntY]-px)*(xs[CrntX+k][CrntY]-px)+(ys[CrntX+k][CrntY]-py)*(ys[CrntX+k][CrntY]-py) );
						if( dst*dst<Dist[CrntX+k][CrntY]*Dist[CrntX+k][CrntY] ) {
							Dist[CrntX+k][CrntY] = dst;
							d_x[CrntX+k][CrntY] = (xs[CrntX+k][CrntY]-px)/dst;
							d_y[CrntX+k][CrntY] = (ys[CrntX+k][CrntY]-py)/dst;
							ClsdPnt[CrntX+k][CrntY] = i;
							ClsdPntId[CrntX+k][CrntY] = -1;
						}
					}
				}
				if( ny2>0.0 )
				{
					eps = ( ceil(py/dy)*dy-py );

					for(j=0; ;j++)
					{
						int CrntY = int(round( (py1+sign(ny2)*(eps/ny2)*ny2 + sign(ny2)*j*dy)/dy ));
						int CrntX = int(ceil( (px1+sign(ny2)*(eps/ny2)*nx2 + sign(ny2)*j*dy*nx2/ny2)/dx ));
						
						if(CrntY<0 || CrntY>NumYPoints-1) break;
						
						for(k=0; ;k++)
						{
							if( CrntX+k<0 ) continue;
							if( CrntX+k>NumXPoints-1 ) break;
							if( (xs[CrntX+k][CrntY]-px)*ny2-(ys[CrntX+k][CrntY]-py)*nx2 < 0 ) break;
							
							double dst = sqrt( (xs[CrntX+k][CrntY]-px)*(xs[CrntX+k][CrntY]-px)+(ys[CrntX+k][CrntY]-py)*(ys[CrntX+k][CrntY]-py) );
							if( dst*dst<Dist[CrntX+k][CrntY]*Dist[CrntX+k][CrntY] ) {
								Dist[CrntX+k][CrntY] = dst;
								d_x[CrntX+k][CrntY] = (xs[CrntX+k][CrntY]-px)/dst;
								d_y[CrntX+k][CrntY] = (ys[CrntX+k][CrntY]-py)/dst;
								ClsdPnt[CrntX+k][CrntY] = i;
								ClsdPntId[CrntX+k][CrntY] = -1;
							}
						}
					}
				}
			}
		}
		else
		{
			double tmp_vx1 = vx1;
			double tmp_vy1 = vy1;
			
			double tmp_nx1 = nx1;
			double tmp_ny1 = ny1;
			
			vx1 = vx2;
			vy1 = vy2;
			
			nx1 = -nx2;
			ny1 = -ny2;
			
			vx2 = tmp_vx1;
			vy2 = tmp_vy1;
			
			nx2 = -tmp_nx1;
			ny2 = -tmp_ny1;
			
			if( nx1>1.0/sqrt2 ) 	// I
			{
				double eps = ( ceil(px/dx)*dx-px );
				
				for(j=0; ;j++)
				{
					int CrntX = int(round( (px1 + sign(nx1)*(eps/nx1)*nx1 + sign(nx1)*j*dx)/dx ));
					int CrntY = int(ceil( (py1 + sign(nx1)*(eps/nx1)*ny1 + sign(nx1)*j*dx*ny1/nx1)/dy ));
					
					if(CrntX>NumXPoints-1) break;
					
					for(k=0; ;k++)
					{
						if(CrntY+k<0) continue;
						if(CrntY+k>NumYPoints-1) break;
						if( (xs[CrntX][CrntY+k]-px)*ny2-(ys[CrntX][CrntY+k]-py)*nx2 < 0 ) break;
						
						double dst = sqrt( (xs[CrntX][CrntY+k]-px)*(xs[CrntX][CrntY+k]-px)+(ys[CrntX][CrntY+k]-py)*(ys[CrntX][CrntY+k]-py) );
						if( dst*dst<Dist[CrntX][CrntY+k]*Dist[CrntX][CrntY+k] ) {
							Dist[CrntX][CrntY+k] = -dst;
							d_x[CrntX][CrntY+k] = -(xs[CrntX][CrntY+k]-px)/dst;
							d_y[CrntX][CrntY+k] = -(ys[CrntX][CrntY+k]-py)/dst;
							ClsdPnt[CrntX][CrntY+k] = i;
							ClsdPntId[CrntX][CrntY+k] = -1;
						}
					}
				}
				if(nx2<0.0)
				{
					eps = -( floor(px/dx)*dx-px );
					
					for(j=0; ;j++)
					{
						int CrntX = int(round( (px1 + sign(nx2)*(eps/nx2)*nx2 + sign(nx2)*j*dx)/dx ));
						int CrntY = int(ceil( (py1 + sign(nx2)*(eps/nx2)*ny2 + sign(nx2)*j*dx*ny2/nx2)/dy ));
						
						if(CrntX<0 || CrntX>NumXPoints-1) break;
						
						for(k=0; ;k++)
						{
							if(CrntY+k<0) continue;
							if( CrntY+k>NumYPoints ) break;
							if( (xs[CrntX][CrntY+k]-px)*ny2-(ys[CrntX][CrntY+k]-py)*nx2 < 0 ) break;
							
							double dst = sqrt( (xs[CrntX][CrntY+k]-px)*(xs[CrntX][CrntY+k]-px)+(ys[CrntX][CrntY+k]-py)*(ys[CrntX][CrntY+k]-py) );
							if( dst*dst<Dist[CrntX][CrntY+k]*Dist[CrntX][CrntY+k] ) {
								Dist[CrntX][CrntY+k] = -dst;
								d_x[CrntX][CrntY+k] = -(xs[CrntX][CrntY+k]-px)/dst;
								d_y[CrntX][CrntY+k] = -(ys[CrntX][CrntY+k]-py)/dst;
								ClsdPnt[CrntX][CrntY+k] = i;
								ClsdPntId[CrntX][CrntY+k] = -1;
							}
						}
					}
				}
			}
			else if(nx1<-1.0/sqrt2+1e-8)	// II
			{
				double eps = -( floor(px/dx)*dx-px );
				
				for(j=0; ;j++)
				{
					int CrntX = int(round( (px1 + sign(nx1)*(eps/nx1)*nx1 + sign(nx1)*j*dx)/dx ));
					int CrntY = int(floor( (py1 + sign(nx1)*(eps/nx1)*ny1 + sign(nx1)*j*dx*ny1/nx1)/dy ));
					
					if(CrntX<0) break;
					
					for(k=0; ;k--)
					{
						if(CrntY+k>(NumYPoints-1)) continue;
						if(CrntY+k<0) break;
						if( (xs[CrntX][CrntY+k]-px)*ny2-(ys[CrntX][CrntY+k]-py)*nx2 < 0 ) break;
						
						double dst = sqrt( (xs[CrntX][CrntY+k]-px)*(xs[CrntX][CrntY+k]-px)+(ys[CrntX][CrntY+k]-py)*(ys[CrntX][CrntY+k]-py) );
						if( dst*dst<Dist[CrntX][CrntY+k]*Dist[CrntX][CrntY+k] ) { 
							Dist[CrntX][CrntY+k] = -dst;
							d_x[CrntX][CrntY+k] = -(xs[CrntX][CrntY+k]-px)/dst;
							d_y[CrntX][CrntY+k] = -(ys[CrntX][CrntY+k]-py)/dst;
							ClsdPnt[CrntX][CrntY+k] = i;
							ClsdPntId[CrntX][CrntY+k] = -1;
						}
					}
				}
				if(nx2>0.0)
				{
					eps = ( ceil(px/dx)*dx-px );
					
					for(j=0; ;j++)
					{
						int CrntX = int(round( (px1 + sign(nx2)*(eps/nx2)*nx2 + sign(nx2)*j*dx)/dx ));
						int CrntY = int(floor( (py1 + sign(nx2)*(eps/nx2)*ny2 + sign(nx2)*j*dx*ny2/nx2)/dy ));

						if(CrntX<0 || CrntX>NumXPoints-1) break;
						
						for(k=0; ;k--)
						{
							if( CrntY+k>NumYPoints-1 ) continue;
							if( CrntY+k<0 ) break;
							if( (xs[CrntX][CrntY+k]-px)*ny2-(ys[CrntX][CrntY+k]-py)*nx2 < 0 ) break;
							
							double dst = sqrt( (xs[CrntX][CrntY+k]-px)*(xs[CrntX][CrntY+k]-px)+(ys[CrntX][CrntY+k]-py)*(ys[CrntX][CrntY+k]-py) );
							if( dst*dst<Dist[CrntX][CrntY+k]*Dist[CrntX][CrntY+k] ) { 
								Dist[CrntX][CrntY+k] = -dst;
								d_x[CrntX][CrntY+k] = -(xs[CrntX][CrntY+k]-px)/dst;
								d_y[CrntX][CrntY+k] = -(ys[CrntX][CrntY+k]-py)/dst;
								ClsdPnt[CrntX][CrntY+k] = i;
								ClsdPntId[CrntX][CrntY+k] = -1;
							}
						}
					}
				}
			}
			else if(ny1>=1.0/sqrt2-1e-8)	// III
			{
				double eps = ( ceil(py/dy)*dy-py );
				
				for(j=0; ;j++)
				{
					int CrntY = int(round( (py1+sign(ny1)*(eps/ny1)*ny1 + sign(ny1)*j*dy)/dy ));
					int CrntX = int(floor( (px1+sign(ny1)*(eps/ny1)*nx1 + sign(ny1)*j*dy*nx1/ny1)/dx ));
					
					if(CrntY>NumYPoints-1) break;
					
					for(k=0; ;k--)
					{
						if( CrntX+k>NumXPoints-1 ) continue;
						if( CrntX+k<0 ) break;
						if( (xs[CrntX+k][CrntY]-px)*ny2-(ys[CrntX+k][CrntY]-py)*nx2 < 0 ) break;
						
						double dst = sqrt( (xs[CrntX+k][CrntY]-px)*(xs[CrntX+k][CrntY]-px)+(ys[CrntX+k][CrntY]-py)*(ys[CrntX+k][CrntY]-py) );
						if( dst*dst<Dist[CrntX+k][CrntY]*Dist[CrntX+k][CrntY] ) { 
							Dist[CrntX+k][CrntY] = -dst;
							d_x[CrntX+k][CrntY] = -(xs[CrntX+k][CrntY]-px)/dst;
							d_y[CrntX+k][CrntY] = -(ys[CrntX+k][CrntY]-py)/dst;
							ClsdPnt[CrntX+k][CrntY] = i;
							ClsdPntId[CrntX+k][CrntY] = -1;
						}
					}
				}
				if(ny2<0.0)
				{
					eps = -( floor(py/dy)*dy-py );
					
					for(j=0; ;j++)
					{
						int CrntY = int(round( (py1+sign(ny2)*(eps/ny2)*ny2 + sign(ny2)*j*dy)/dy ));
						int CrntX = int(floor( (px1+sign(ny2)*(eps/ny2)*nx2 + sign(ny2)*j*dy*nx2/ny2)/dx ));
						
						if(CrntY<0 || CrntY>NumYPoints-1) break;
						
						for(k=0; ;k--)
						{
							if( CrntX+k>NumXPoints-1 ) continue;
							if( CrntX+k<0 ) break;
							if( (xs[CrntX+k][CrntY]-px)*ny2-(ys[CrntX+k][CrntY]-py)*nx2 < 0 ) break;
							
							double dst = sqrt( (xs[CrntX+k][CrntY]-px)*(xs[CrntX+k][CrntY]-px)+(ys[CrntX+k][CrntY]-py)*(ys[CrntX+k][CrntY]-py) );
							if( dst*dst<Dist[CrntX+k][CrntY]*Dist[CrntX+k][CrntY] ) { 
								Dist[CrntX+k][CrntY] = -dst;
								d_x[CrntX+k][CrntY] = -(xs[CrntX+k][CrntY]-px)/dst;
								d_y[CrntX+k][CrntY] = -(ys[CrntX+k][CrntY]-py)/dst;
								ClsdPnt[CrntX+k][CrntY] = i;
								ClsdPntId[CrntX+k][CrntY] = -1;
							}
						}
					}
				}
			}
			else	// IV
			{
				double eps = -( floor(py/dy)*dy-py );
				
				for(j=0; ;j++)
				{
					int CrntY = int(round( (py1+sign(ny1)*(eps/ny1)*ny1 + sign(ny1)*j*dy)/dy ));
					int CrntX = int(ceil( (px1+sign(ny1)*(eps/ny1)*nx1 + sign(ny1)*j*dy*nx1/ny1)/dx ));
					
					if(CrntY<0) break;
					
					for(k=0; ;k++)
					{
						if( CrntX+k<0 ) continue;
						if(  CrntX+k>NumXPoints-1) break;
						if( (xs[CrntX+k][CrntY]-px)*ny2-(ys[CrntX+k][CrntY]-py)*nx2 < 0 ) break;
						
						double dst = sqrt( (xs[CrntX+k][CrntY]-px)*(xs[CrntX+k][CrntY]-px)+(ys[CrntX+k][CrntY]-py)*(ys[CrntX+k][CrntY]-py) );
						if( dst*dst<Dist[CrntX+k][CrntY]*Dist[CrntX+k][CrntY] ) { 
							Dist[CrntX+k][CrntY] = -dst;
							d_x[CrntX+k][CrntY] = -(xs[CrntX+k][CrntY]-px)/dst;
							d_y[CrntX+k][CrntY] = -(ys[CrntX+k][CrntY]-py)/dst;
							ClsdPnt[CrntX+k][CrntY] = i;
							ClsdPntId[CrntX+k][CrntY] = -1;
						}
					}
				}
				if( ny2>0.0 )
				{
					eps = ( ceil(py/dy)*dy-py );

					
					for(j=0; ;j++)
					{
						int CrntY = int(round( (py1+sign(ny2)*(eps/ny2)*ny2 + sign(ny2)*j*dy)/dy ));
						int CrntX = int(ceil( (px1+sign(ny2)*(eps/ny2)*nx2 + sign(ny2)*j*dy*nx2/ny2)/dx ));
						
						if(CrntY<0 || CrntY>NumYPoints-1) break;
						
						for(k=0; ;k++)
						{
							if( CrntX+k<0 ) continue;
							if( CrntX+k>NumXPoints-1 ) break;
							if( (xs[CrntX+k][CrntY]-px)*ny2-(ys[CrntX+k][CrntY]-py)*nx2 < 0 ) break;
							
							double dst = sqrt( (xs[CrntX+k][CrntY]-px)*(xs[CrntX+k][CrntY]-px)+(ys[CrntX+k][CrntY]-py)*(ys[CrntX+k][CrntY]-py) );
							if( dst*dst<Dist[CrntX+k][CrntY]*Dist[CrntX+k][CrntY] ) { 
								Dist[CrntX+k][CrntY] = -dst;
								d_x[CrntX+k][CrntY] = -(xs[CrntX+k][CrntY]-px)/dst;
								d_y[CrntX+k][CrntY] = -(ys[CrntX+k][CrntY]-py)/dst;
								ClsdPnt[CrntX+k][CrntY] = i;
								ClsdPntId[CrntX+k][CrntY] = -1;
							}
						}
					}
				}
			}
		}
	}

/////////////////////////////////////////////////////////////////////////////////////
	for(i=0; i<NumBndryPoints; i++)
	{
		int NxtPnt = i+1;
		if( NxtPnt>(NumBndryPoints-1) ) NxtPnt = 0;

		double px = BndryX[i];
		double py = BndryY[i];

		double px1 = BndryX[i]-xs[0][0];
		double py1 = BndryY[i]-ys[0][0];

		double vx = BndryX[NxtPnt]-BndryX[i];
		double vy = BndryY[NxtPnt]-BndryY[i];
		
		double v = sqrt( vx*vx+vy*vy );
		
		vx /= v;
		vy /= v;
		
		double nx = vy;
		double ny = -vx;
		
		Vxs[i] = vx;
		Vys[i] = vy;
		Vs[i] = v;
		
		Nxs[i] = nx;
		Nys[i] = ny;
		
		if( vy>=1.0/sqrt2-1.0e-8 )
		{
			if(vx>0)	// I
			{
				double eps = ceil(px/dx)*dx-px;
				
				for(j=0; ;j++)
				{
					int CrntX = int(round( (px1 + sign(nx)*(eps/nx)*nx + sign(nx)*j*dx)/dx ));
					int CrntY = int(ceil( (py1 + sign(nx)*(eps/nx)*ny + sign(nx)*j*dx*ny/nx)/dy ));
					
					if(CrntX>NumXPoints-1) break;
					
					for(k=0; ;k++)
					{
						if(CrntY+k<0) continue;
						if(CrntY+k>NumYPoints-1) break;
						if( vx*(xs[CrntX][CrntY+k]-px)+vy*(ys[CrntX][CrntY+k]-py) > v ) break;
						
						double dst = nx*(xs[CrntX][CrntY+k]-px)+ny*(ys[CrntX][CrntY+k]-py);
						if( dst*dst<Dist[CrntX][CrntY+k]*Dist[CrntX][CrntY+k] ) { 
							Dist[CrntX][CrntY+k] = dst;
							d_y[CrntX][CrntY+k] = -vx/(vx*vx+vy*vy);
							d_x[CrntX][CrntY+k] =  vy/(vx*vx+vy*vy);
							ClsdPnt[CrntX][CrntY+k] = i;
						}
					}
				}
				for(j=-1; ;j--)
				{
					int CrntX = int(round( (px1 + sign(nx)*(eps/nx)*nx + sign(nx)*j*dx)/dx ));
					int CrntY = int(ceil( (py1 + sign(nx)*(eps/nx)*ny + sign(nx)*j*dx*ny/nx)/dy ));
					
					if(CrntX<0) break;
					
					for(k=0; ;k++)
					{
						if(CrntY+k<0) continue;
						if(CrntY+k>NumYPoints-1) break;
						if( vx*(xs[CrntX][CrntY+k]-px)+vy*(ys[CrntX][CrntY+k]-py) > v ) break;
						
						double dst = nx*(xs[CrntX][CrntY+k]-px)+ny*(ys[CrntX][CrntY+k]-py);
						if( dst*dst<Dist[CrntX][CrntY+k]*Dist[CrntX][CrntY+k] ) {
							Dist[CrntX][CrntY+k] = dst;
							d_y[CrntX][CrntY+k] = -vx/(vx*vx+vy*vy);
							d_x[CrntX][CrntY+k] =  vy/(vx*vx+vy*vy);
							ClsdPnt[CrntX][CrntY+k] = i;
						}
					}
				}
			}
			else	// II
			{
				double eps = ( ceil(px/dx)*dx-px );
				
				for(j=0; ;j++)
				{
					int CrntX = int(round( (px1 + sign(nx)*(eps/nx)*nx + sign(nx)*j*dx)/dx ));
					int CrntY = int(ceil( (py1 + sign(nx)*(eps/nx)*ny + sign(nx)*j*dx*ny/nx)/dy ));
					
					if(CrntX>NumXPoints-1) break;
					
					for(k=0; ;k++)
					{
						if( CrntY+k<0 ) continue;
						if( CrntY+k>NumYPoints-1 ) break;
						if( vx*(xs[CrntX][CrntY+k]-px)+vy*(ys[CrntX][CrntY+k]-py) > v ) break;
						
						double dst = nx*(xs[CrntX][CrntY+k]-px)+ny*(ys[CrntX][CrntY+k]-py);
						if( dst*dst<Dist[CrntX][CrntY+k]*Dist[CrntX][CrntY+k] ) { 
							Dist[CrntX][CrntY+k] = dst;
							d_y[CrntX][CrntY+k] = -vx/(vx*vx+vy*vy);
							d_x[CrntX][CrntY+k] =  vy/(vx*vx+vy*vy);
							ClsdPnt[CrntX][CrntY+k] = i;
						}
					}
				}
				for(j=-1; ;j--)
				{
					int CrntX = int(round( (px1 + sign(nx)*(eps/nx)*nx + sign(nx)*j*dx)/dx ));
					int CrntY = int(ceil( (py1 + sign(nx)*(eps/nx)*ny + sign(nx)*j*dx*ny/nx)/dy ));
					
					if(CrntX<0) break;
					
					for(k=0; ;k++)
					{
						if( CrntY+k<0 ) continue;
						if( CrntY+k>NumYPoints-1 ) break;
						if( vx*(xs[CrntX][CrntY+k]-px)+vy*(ys[CrntX][CrntY+k]-py) > v ) break;
						
						double dst = nx*(xs[CrntX][CrntY+k]-px)+ny*(ys[CrntX][CrntY+k]-py);
						if( dst*dst<Dist[CrntX][CrntY+k]*Dist[CrntX][CrntY+k] ) { 
							Dist[CrntX][CrntY+k] = dst;
							d_y[CrntX][CrntY+k] = -vx/(vx*vx+vy*vy);
							d_x[CrntX][CrntY+k] =  vy/(vx*vx+vy*vy);
							ClsdPnt[CrntX][CrntY+k] = i;
						}
					}
				}
			}
		}
		else if( vy<-1.0/sqrt2+1.0e-8 )
		{
			if(vx<0)	// III
			{
				double eps = -( floor(px/dx)*dx-px );
				
				for(j=0; ;j++)
				{
					int CrntX = int(round( (px1 + sign(nx)*(eps/nx)*nx + sign(nx)*j*dx)/dx ));
					int CrntY = int(floor( (py1 + sign(nx)*(eps/nx)*ny + sign(nx)*j*dx*ny/nx)/dy ));
					
					if(CrntX<0) break;
					
					for(k=0; ;k--)
					{
						if( CrntY+k > (NumYPoints-1) ) continue;
						if( CrntY+k<0 ) break;
						if( vx*(xs[CrntX][CrntY+k]-px)+vy*(ys[CrntX][CrntY+k]-py) > v ) break;
						
						double dst = nx*(xs[CrntX][CrntY+k]-px)+ny*(ys[CrntX][CrntY+k]-py);
						if( dst*dst<Dist[CrntX][CrntY+k]*Dist[CrntX][CrntY+k] ) {
							Dist[CrntX][CrntY+k] = dst;
							d_y[CrntX][CrntY+k] = -vx/(vx*vx+vy*vy);
							d_x[CrntX][CrntY+k] =  vy/(vx*vx+vy*vy);
							ClsdPnt[CrntX][CrntY+k] = i;
						}
					}
				}
				for(j=-1; ;j--)
				{
					int CrntX = int(round( (px1 + sign(nx)*(eps/nx)*nx + sign(nx)*j*dx)/dx ));
					int CrntY = int(floor( (py1 + sign(nx)*(eps/nx)*ny + sign(nx)*j*dx*ny/nx)/dy ));
					
					if(CrntX>NumXPoints-1) break;
					
					for(k=0; ;k--)
					{
						if(CrntY+k>NumYPoints-1) continue;
						if(CrntY+k<0) break;
						if( vx*(xs[CrntX][CrntY+k]-px)+vy*(ys[CrntX][CrntY+k]-py) > v ) break;
						
						double dst = nx*(xs[CrntX][CrntY+k]-px)+ny*(ys[CrntX][CrntY+k]-py);
						if( dst*dst<Dist[CrntX][CrntY+k]*Dist[CrntX][CrntY+k] ) { 
							Dist[CrntX][CrntY+k] = dst;
							d_y[CrntX][CrntY+k] = -vx/(vx*vx+vy*vy);
							d_x[CrntX][CrntY+k] =  vy/(vx*vx+vy*vy);
							ClsdPnt[CrntX][CrntY+k] = i;
						}
					}
				}
			}
			else	// IV
			{
				double eps = -( floor(px/dx)*dx-px );
				
				for(j=0; ;j++)
				{
					int CrntX = int(round( (px1 + sign(nx)*(eps/nx)*nx + sign(nx)*j*dx)/dx ));
					int CrntY = int(floor( (py1 + sign(nx)*(eps/nx)*ny + sign(nx)*j*dx*ny/nx)/dy ));
					
					if(CrntX<0) break;
					
					for(k=0; ;k--)
					{
						if(CrntY+k>NumYPoints-1) continue;
						if(CrntY+k<0) break;
						if( vx*(xs[CrntX][CrntY+k]-px)+vy*(ys[CrntX][CrntY+k]-py) > v ) break;
						
						double dst = nx*(xs[CrntX][CrntY+k]-px)+ny*(ys[CrntX][CrntY+k]-py);
						if( dst*dst<Dist[CrntX][CrntY+k]*Dist[CrntX][CrntY+k] ) { 
							Dist[CrntX][CrntY+k] = dst;
							d_y[CrntX][CrntY+k] = -vx/(vx*vx+vy*vy);
							d_x[CrntX][CrntY+k] =  vy/(vx*vx+vy*vy);
							ClsdPnt[CrntX][CrntY+k] = i;
						}
					}
				}
				for(j=-1; ;j--)
				{
					int CrntX = int(round( (px1 + sign(nx)*(eps/nx)*nx + sign(nx)*j*dx)/dx ));
					int CrntY = int(floor( (py1 + sign(nx)*(eps/nx)*ny + sign(nx)*j*dx*ny/nx)/dy ));
					
					if(CrntX>NumXPoints-1) break;
					
					for(k=0; ;k--)
					{
						if(CrntY+k>NumYPoints-1) continue;
						if(CrntY+k<0) break;
						if( vx*(xs[CrntX][CrntY+k]-px)+vy*(ys[CrntX][CrntY+k]-py) > v ) break;
						
						double dst = nx*(xs[CrntX][CrntY+k]-px)+ny*(ys[CrntX][CrntY+k]-py);
						if( dst*dst<Dist[CrntX][CrntY+k]*Dist[CrntX][CrntY+k] ) { 
							Dist[CrntX][CrntY+k] = dst;
							d_y[CrntX][CrntY+k] = -vx/(vx*vx+vy*vy);
							d_x[CrntX][CrntY+k] =  vy/(vx*vx+vy*vy);
							ClsdPnt[CrntX][CrntY+k] = i;
						}
					}
				}
			}
		}
		else if( vx>=1.0/sqrt2 )
		{
			if(vy>0)	// I'
			{
				double eps = -( floor(py/dy)*dy-py );
				
				for(j=0; ;j++)
				{
					int CrntY = int(round( (py1+sign(ny)*(eps/ny)*ny + sign(ny)*j*dy)/dy ));
					int CrntX = int(ceil( (px1+sign(ny)*(eps/ny)*nx + sign(ny)*j*dy*nx/ny)/dx ));
					
					if(CrntY<0) break;
					
					for(k=0; ;k++)
					{
						if(CrntX+k<0) continue;
						if(CrntX+k>NumXPoints-1) break;
						if( vx*(xs[CrntX+k][CrntY]-px)+vy*(ys[CrntX+k][CrntY]-py) > v ) break;
						
						double dst = nx*(xs[CrntX+k][CrntY]-px)+ny*(ys[CrntX+k][CrntY]-py);
						if( dst*dst<Dist[CrntX+k][CrntY]*Dist[CrntX+k][CrntY] ) { 
							Dist[CrntX+k][CrntY] = dst;
							d_y[CrntX+k][CrntY] = -vx/(vx*vx+vy*vy);
							d_x[CrntX+k][CrntY] =  vy/(vx*vx+vy*vy);
							ClsdPnt[CrntX+k][CrntY] = i;
						}
					}
				}
				for(j=-1; ;j--)
				{
					int CrntY = int(round( (py1+sign(ny)*(eps/ny)*ny + sign(ny)*j*dy)/dy ));
					int CrntX = int(ceil( (px1+sign(ny)*(eps/ny)*nx + sign(ny)*j*dy*nx/ny)/dx ));
					
					if(CrntY>NumYPoints-1) break;
					
					for(k=0; ;k++)
					{
						if(CrntX+k<0) continue;
						if(CrntX+k>NumXPoints-1) break;
						if( vx*(xs[CrntX+k][CrntY]-px)+vy*(ys[CrntX+k][CrntY]-py) > v ) break;
						
						double dst = nx*(xs[CrntX+k][CrntY]-px)+ny*(ys[CrntX+k][CrntY]-py);
						if( dst*dst<Dist[CrntX+k][CrntY]*Dist[CrntX+k][CrntY] ) { 
							Dist[CrntX+k][CrntY] = dst;
							d_y[CrntX+k][CrntY] = -vx/(vx*vx+vy*vy);
							d_x[CrntX+k][CrntY] =  vy/(vx*vx+vy*vy);
							ClsdPnt[CrntX+k][CrntY] = i;
						}
					}
				}
			}
			else	// II'
			{
				double eps = -( floor(py/dy)*dy-py );
				
				for(j=0; ;j++)
				{
					int CrntY = int(round( (py1+sign(ny)*(eps/ny)*ny + sign(ny)*j*dy)/dy ));
					int CrntX = int(ceil( (px1+sign(ny)*(eps/ny)*nx + sign(ny)*j*dy*nx/ny)/dx ));
					
					if(CrntY<0) break;
					
					for(k=0; ;k++)
					{
						if(CrntX+k<0) continue;
						if(CrntX+k>NumXPoints-1) break;
						if( vx*(xs[CrntX+k][CrntY]-px)+vy*(ys[CrntX+k][CrntY]-py) > v ) break;
						
						double dst = nx*(xs[CrntX+k][CrntY]-px)+ny*(ys[CrntX+k][CrntY]-py);
						if( dst*dst<Dist[CrntX+k][CrntY]*Dist[CrntX+k][CrntY] ) { 
							Dist[CrntX+k][CrntY] = dst;
							d_y[CrntX+k][CrntY] = -vx/(vx*vx+vy*vy);
							d_x[CrntX+k][CrntY] =  vy/(vx*vx+vy*vy);
							ClsdPnt[CrntX+k][CrntY] = i;
						}
					}
				}
				for(j=-1; ;j--)
				{
					int CrntY = int(round( (py1+sign(ny)*(eps/ny)*ny + sign(ny)*j*dy)/dy ));
					int CrntX = int(ceil( (px1+sign(ny)*(eps/ny)*nx + sign(ny)*j*dy*nx/ny)/dx ));
					
					if(CrntY>NumYPoints-1) break;
					
					for(k=0; ;k++)
					{
						if(CrntX+k<0) continue;
						if(CrntX+k>NumXPoints-1) break;
						if( vx*(xs[CrntX+k][CrntY]-px)+vy*(ys[CrntX+k][CrntY]-py) > v ) break;
						
						double dst = nx*(xs[CrntX+k][CrntY]-px)+ny*(ys[CrntX+k][CrntY]-py);
						if( dst*dst<Dist[CrntX+k][CrntY]*Dist[CrntX+k][CrntY] ) { 
							Dist[CrntX+k][CrntY] = dst;
							d_y[CrntX+k][CrntY] = -vx/(vx*vx+vy*vy);
							d_x[CrntX+k][CrntY] =  vy/(vx*vx+vy*vy);
							ClsdPnt[CrntX+k][CrntY] = i;
						}
					}
				}
			}
		}
		else
		{
			if(vy>0)	// III'
			{
				double eps = ( ceil(py/dy)*dy-py );
				
				for(j=0; ;j++)
				{
					int CrntY = int(round( (py1+sign(ny)*(eps/ny)*ny + sign(ny)*j*dy)/dy ));
					int CrntX = int(floor( (px1+sign(ny)*(eps/ny)*nx + sign(ny)*j*dy*nx/ny)/dx ));
					
					if(CrntY>NumYPoints-1) break;
					
					for(k=0; ;k--)
					{
						if(CrntX+k>NumXPoints-1) continue;
						if(CrntX+k<0) break;
						if( vx*(xs[CrntX+k][CrntY]-px)+vy*(ys[CrntX+k][CrntY]-py) > v ) break;
						
						double dst = nx*(xs[CrntX+k][CrntY]-px)+ny*(ys[CrntX+k][CrntY]-py);
						if( dst*dst<Dist[CrntX+k][CrntY]*Dist[CrntX+k][CrntY] ) { 
							Dist[CrntX+k][CrntY] = dst;
							d_y[CrntX+k][CrntY] = -vx/(vx*vx+vy*vy);
							d_x[CrntX+k][CrntY] =  vy/(vx*vx+vy*vy);
							ClsdPnt[CrntX+k][CrntY] = i;
						}
					}
				}
				for(j=-1; ;j--)
				{
					int CrntY = int(round( (py1+sign(ny)*(eps/ny)*ny + sign(ny)*j*dy)/dy ));
					int CrntX = int(floor( (px1+sign(ny)*(eps/ny)*nx + sign(ny)*j*dy*nx/ny)/dx ));
					
					if(CrntY<0) break;
					
					for(k=0; ;k--)
					{
						if(CrntX+k>NumXPoints-1) continue;
						if(CrntX+k<0) break;
						if( vx*(xs[CrntX+k][CrntY]-px)+vy*(ys[CrntX+k][CrntY]-py) > v ) break;
						
						double dst = nx*(xs[CrntX+k][CrntY]-px)+ny*(ys[CrntX+k][CrntY]-py);
						if( dst*dst<Dist[CrntX+k][CrntY]*Dist[CrntX+k][CrntY] ) {
							Dist[CrntX+k][CrntY] = dst;
							d_y[CrntX+k][CrntY] = -vx/(vx*vx+vy*vy);
							d_x[CrntX+k][CrntY] =  vy/(vx*vx+vy*vy);
							ClsdPnt[CrntX+k][CrntY] = i;
						}
					}
				}
			}
			else	// IV'
			{
				double eps = ( ceil(py/dy)*dy-py );
				
				for(j=0; ;j++)
				{
					int CrntY = int(round( (py1+sign(ny)*(eps/ny)*ny + sign(ny)*j*dy)/dy ));
					int CrntX = int(floor( (px1+sign(ny)*(eps/ny)*nx + sign(ny)*j*dy*nx/ny)/dx ));
					
					if(CrntY>NumYPoints-1) break;
					
					for(k=0; ;k--)
					{
						if(CrntX+k>NumXPoints-1) continue;
						if(CrntX+k<0) break;
						if( vx*(xs[CrntX+k][CrntY]-px)+vy*(ys[CrntX+k][CrntY]-py) > v ) break;
						
						double dst = nx*(xs[CrntX+k][CrntY]-px)+ny*(ys[CrntX+k][CrntY]-py);
						if( dst*dst<Dist[CrntX+k][CrntY]*Dist[CrntX+k][CrntY] ) {
							Dist[CrntX+k][CrntY] = dst;
							d_y[CrntX+k][CrntY] = -vx/(vx*vx+vy*vy);
							d_x[CrntX+k][CrntY] =  vy/(vx*vx+vy*vy);
							ClsdPnt[CrntX+k][CrntY] = i;
						}
					}
				}
				for(j=-1; ;j--)
				{
					int CrntY = int(round( (py1+sign(ny)*(eps/ny)*ny + sign(ny)*j*dy)/dy ));
					int CrntX = int(floor( (px1+sign(ny)*(eps/ny)*nx + sign(ny)*j*dy*nx/ny)/dx ));
					
					if(CrntY<0) break;
					
					for(k=0; ;k--)
					{
						if(CrntX+k>NumXPoints-1) continue;
						if(CrntX+k<0) break;
						if( vx*(xs[CrntX+k][CrntY]-px)+vy*(ys[CrntX+k][CrntY]-py) > v ) break;
						
						double dst = nx*(xs[CrntX+k][CrntY]-px)+ny*(ys[CrntX+k][CrntY]-py);
						if( dst*dst<Dist[CrntX+k][CrntY]*Dist[CrntX+k][CrntY] ) { 
							Dist[CrntX+k][CrntY] = dst;
							d_y[CrntX+k][CrntY] = -vx/(vx*vx+vy*vy);
							d_x[CrntX+k][CrntY] =  vy/(vx*vx+vy*vy);
							ClsdPnt[CrntX+k][CrntY] = i;
						}
					}
				}
			}
		}
		
	}
}
