/*
Read input files
*/

ifstream* rd1;
ifstream* rd2;

char str1[100];
char str2[100];
char str3[100];

int NumBndryPoints = 0;

rd1 = new ifstream("./meshfiles/x.dat");
while(rd1->getline(str1,100))
{
	NumBndryPoints++;
}
rd1->close();
delete rd1;

// cout<<NumBndryPoints;

double* BndryX = new double[NumBndryPoints*20];
double* BndryY = new double[NumBndryPoints*20];

rd1 = new ifstream("./meshfiles/x.dat");
rd2 = new ifstream("./meshfiles/y.dat");

for(i=0; i<NumBndryPoints; i++)
{
	rd1->getline(str1,100);
	rd2->getline(str2,100);

	BndryX[i] = atof(str1);
	BndryY[i] = atof(str2);
}
rd1->close();
delete rd1;

rd2->close();
delete rd2;

int NumEdges = 0;
rd1 = new ifstream("./meshfiles/edge_1.dat");
while( rd1->getline(str1,100) )
{
	NumEdges++;
}
rd1->close();
delete rd1;

int* Edge1 = new int[NumEdges*20];
int* Edge2 = new int[NumEdges*20];

rd1 = new ifstream("./meshfiles/edge_1.dat");
rd2 = new ifstream("./meshfiles/edge_2.dat");

for(i=0; i<NumEdges; i++)
{
	rd1->getline(str1,100);
	rd2->getline(str2,100);

	Edge1[i] = int( atof(str1) )-1;
	Edge2[i] = int( atof(str2) )-1;
}
rd1->close();
delete rd1;

rd2->close();
delete rd2;

int NumApicalPulled = 0;

rd1 = new ifstream("./meshfiles/apical_pulled.dat");
while( rd1->getline(str1,100) )
{
	NumApicalPulled++;
}
rd1->close();
delete rd1;

int* ApicalPulled = new int[NumApicalPulled];

rd1 = new ifstream("./meshfiles/apical_pulled.dat");

for(i=0; i<NumApicalPulled; i++)
{
	rd1->getline(str1,100);

	ApicalPulled[i] = int( atof(str1) )-1;
}
rd1->close();
delete rd1;

int NumLateralPulled = 0;

rd1 = new ifstream("./meshfiles/lateral_pulled.dat");
while( rd1->getline(str1,100) )
{
	NumLateralPulled++;
}
rd1->close();
delete rd1;

int* LateralPulled = new int[NumLateralPulled];

rd1 = new ifstream("./meshfiles/lateral_pulled.dat");

for(i=0; i<NumLateralPulled; i++)
{
	rd1->getline(str1,100);

	LateralPulled[i] = int( atof(str1) )-1;
}
rd1->close();
delete rd1;

int NumApicalPoints = 0;

rd1 = new ifstream("./meshfiles/apical_points.dat");
while( rd1->getline(str1,100) )
{
	NumApicalPoints++;
}
rd1->close();
delete rd1;

int* ApicalPoints = new int[NumApicalPoints];

rd1 = new ifstream("./meshfiles/apical_points.dat");

for(i=0; i<NumApicalPoints; i++)
{
	rd1->getline(str1,100);

	ApicalPoints[i] = int( atof(str1) )-1;
}
rd1->close();
delete rd1;

int NumBasalEdges = 0;

rd1 = new ifstream("./meshfiles/basal_edges.dat");
while( rd1->getline(str1,100) )
{
	NumBasalEdges++;
}
rd1->close();
delete rd1;

int* BasalEdges = new int[NumBasalEdges];

rd1 = new ifstream("./meshfiles/basal_edges.dat");

for(i=0; i<NumBasalEdges; i++)
{
	rd1->getline(str1,100);

	BasalEdges[i] = int( atof(str1) )-1;
}
rd1->close();
delete rd1;

/*
cout<<"NumBndryPoints="<<NumBndryPoints<<endl;
cout<<"NumEdges="<<NumEdges<<endl;
cout<<"NumApicalPulled="<<NumApicalPulled<<endl;
cout<<"NumLateralPulled="<<NumLateralPulled<<endl;
cout<<"NumApicalPoints="<<NumApicalPoints<<endl;
cout<<"NumBasalEdges="<<NumBasalEdges<<endl;

return 0;
*/

double* L0s = new double[NumEdges*20];
double* Es = new double[NumEdges*20];

for(i=0; i<NumEdges; i++)
{
	int idx1 = Edge1[i];
	int idx2 = Edge2[i];

	x1 = BndryX[idx1];
	y1 = BndryY[idx1];

	x2 = BndryX[idx2];
	y2 = BndryY[idx2];

	rx = x2-x1;
	ry = y2-y1;

	r2 = rx*rx + ry*ry;

	L0s[i] = sqrt(r2);

	Es[i] = E;
}

int* EdgeId = new int[NumEdges*20];

for(i=0; i<NumEdges; i++)
{
	EdgeId[i] = 0;
}

for(i=0; i<NumApicalPulled; i++)
{
	EdgeId[ApicalPulled[i]] = 1;
}

for(i=0; i<NumLateralPulled; i++)
{
	EdgeId[LateralPulled[i]] = 2;
}

for(i=0; i<NumBasalEdges; i++)
{
	EdgeId[BasalEdges[i]] = 3;
}
