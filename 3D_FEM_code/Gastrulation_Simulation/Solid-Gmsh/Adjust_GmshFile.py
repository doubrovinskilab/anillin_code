import os
os.environ["HOME"] = "/mnt/d"
Location = os.getenv('HOME')+ "/Dropbox/python/Tools/"
exec(compile(open(Location+"Functions.py", "rb").read(), Location+"Functions.py", 'exec')) # For python 3
#execfile(Location+"Functions.py") # For python 2.7
np.set_printoptions(suppress=True)

Remove_Half = True
Folders   =["basal_closed","basal_open"]
VTK_Files =["Embryo_slice_closed","Embryo_slice_open"]

R_Cell_Membrane = 80.
Peri_Fluid_Thickness = 4.
R = R_Cell_Membrane - Peri_Fluid_Thickness
Cell_Height = 35. # um

for f in range(len(Folders)):
  folder = Folders[f]
  VTK = readVTKFile("%s.vtk" % VTK_Files[f])

  P, Cl, Ct = VTK['P'], VTK['Cl'], VTK['Ct']

  if len(Cl)>0:
    raise NameError("Cl expected to be 0")

  # Get all the edges
  Cl = GetTriangleEdges(Ct)

  # Rounnd up
  P = P.round(decimals=5)

  # Remove lines that are zero in length (Excess edges)
  Cl = Cl[np.where(CalcLinesLength(P,Cl)!=0)]

  # Remove lines that are very long (Incorrect edges)
  Cl = Cl[np.where(CalcLinesLength(P,Cl)<3.0)]

  # Remove repeated lines
  Cl = np.unique(Cl,axis=0)

  # Circle Edges
  #===========================================
  P_Cl_center = (P[Cl[:,0]] + P[Cl[:,1]])/2.0
  R_p = np.sqrt(P_Cl_center[:,0]**2 + P_Cl_center[:,1]**2)
  I_outer_circle = np.where( (R_p >= R-0.1)*(P_Cl_center[:,2] < np.min(P_Cl_center[:,2])+0.01) )[0]
  CE = np.ones(len(Cl),dtype=float)*-1
  CE[I_outer_circle] = 1.0
  Cp_Gr = np.ones(len(P), dtype=int)*-1 # Point Cells
  Cp_Gr[np.unique(Cl[I_outer_circle].flatten())] = 1

  if Remove_Half:
    tol = 0.1
    # Remove Half of domain
    P_Cl_center = (P[Cl[:,0]] + P[Cl[:,1]])/2.0
    Remo_Cl  = np.where(P_Cl_center[:,0]>tol)[0]  # lines to remove
    Keep_Cl  = np.where(P_Cl_center[:,0]<=tol)[0] # lines to keep
    Cl = Cl[Keep_Cl]
    P_Cl_counter = np.zeros(len(P))
    np.add.at(P_Cl_counter,Cl[:,0],1)
    np.add.at(P_Cl_counter,Cl[:,1],1)

    Remo_P = np.where(P_Cl_counter==0)[0] # Points to remove
    Keep_P = np.where((P_Cl_counter!=0)+(Cp_Gr==1))[0] # Points to stay

    # Traingle to remove
    Remo_Ct = np.in1d(Ct[:,0],Remo_P)+np.in1d(Ct[:,1],Remo_P)+np.in1d(Ct[:,2],Remo_P)
    Keep_Ct = ~Remo_Ct
    Ct = Ct[Keep_Ct]

    # Update P index since points are removed
    New_index_P = np.ones(len(P),dtype=int)*-1
    New_index_P[Keep_P] = np.arange(len(Keep_P))

    P  = P[Keep_P]
    Cp_Gr = Cp_Gr[Keep_P]

    Cl[:,0] = New_index_P[Cl[:,0]] 
    Cl[:,1] = New_index_P[Cl[:,1]] 

    Ct[:,0] = New_index_P[Ct[:,0]] 
    Ct[:,1] = New_index_P[Ct[:,1]] 
    Ct[:,2] = New_index_P[Ct[:,2]] 

    I = np.where((P[:,0]>0.0)*(Cp_Gr!=1))[0]
    P[I,0] = 0.0

    folder = "half_%s" % folder

  if not os.path.exists(folder):
    os.makedirs(folder)

  # Choosing the Cells for the contractil domain
  #=============================================
  # Creating the angle of each line
  P_Cl_center = (P[Cl[:,0]] + P[Cl[:,1]])/2.0
  Theta = np.arctan(np.divide(P_Cl_center[:,1],P_Cl_center[:,0]))
  I = np.where((P_Cl_center[:,0]>=0)*(P_Cl_center[:,1]>=0))[0]
  Theta[I] = Theta[I]
  I = np.where((P_Cl_center[:,0]<0)*(P_Cl_center[:,1]>=0))[0]
  Theta[I] = Theta[I]+np.pi
  I = np.where((P_Cl_center[:,0]>=0)*(P_Cl_center[:,1]<0))[0]
  Theta[I] = Theta[I]+2*np.pi
  I = np.where((P_Cl_center[:,0]<0)*(P_Cl_center[:,1]<0))[0]
  Theta[I] = Theta[I]+np.pi

  # Creatign the Radius of each line
  R_p = np.sqrt(P_Cl_center[:,0]**2 + P_Cl_center[:,1]**2)

  Left_Angle, Right_Angle = np.pi/2+0.52, np.pi/2-0.52
  # Finding the contactile domain

  # Elastic Stresses
  #===========================================
  M_apical = 5.  
  M_lateral = 18.
  M_basal = 18. 

  I_apical = np.where( (R_p >= R-0.1) )
  I_lateral = np.where( (R_p < R-0.1) & (R_p>R-Cell_Height+0.1))
  I_basal = np.where( (R_p < R-Cell_Height+0.1) )

  LE = np.ones(len(Cl),dtype=float)
  LE[I_apical] = M_apical
  LE[I_lateral] = M_lateral
  LE[I_basal] = M_basal #  Contractial Domain

  np.savetxt('%s/LE_Edges_Multi_A%.1lf_L%.1lf_B%.1lf.txt' % (folder, M_apical, M_lateral, M_basal),LE,fmt="%lf") # linear 
  # Surface Tension
  #===========================================
  M_apical = 160.   # Original 12.
  M_lateral = 16.   # Original 1.
  M_basal = 2.     # Original 2.
  M_general = 0.5  # Original 0.5

  I_apical = np.where( (Theta>Right_Angle) & (Theta<Left_Angle) & (R_p >= R-0.1) )
  I_lateral = np.where( (Theta>Right_Angle) & (Theta<Left_Angle) & (R_p < R-0.1) & (R_p>R-Cell_Height+0.1))
  I_basal = np.where( (Theta>Right_Angle) & (Theta<Left_Angle) & (R_p < R-Cell_Height+0.1) )

  ST_Option   = 'Constant'
  ST_Option2  = 'Only_Latteral_Radial_Direction' # Only_Latteral_Radial_Direction or Random

  if ST_Option == 'Constant':
    ST = np.ones(len(Cl),dtype=float)*M_general
    ST[I_apical] = M_apical
    ST[I_lateral] = M_lateral
    ST[I_basal] = M_basal # Contractial Domain
    # Note Cl_Gr in VTK will be converted to int so the 0.5-->0.0

    if ST_Option2 == 'Only_Latteral_Radial_Direction':
      NV = CalcLinesNormalizedVector(P,Cl) # Normalized vector    
      # Remove ST from lines parallel to z-axis
      I_Z_lines = np.where((np.abs(NV[:,2]) > np.abs(NV[:,0]))*(np.abs(NV[:,2])>np.abs(NV[:,1])))[0]
      ST[I_Z_lines] = 0.0
      # Remove ST from line not parallel to radial vector
      P_Cl_center2 = np.zeros([len(P_Cl_center)+1,3])
      P_Cl_center2[:-1] = P_Cl_center
      P_Cl_center2[-1] = [0,0,0]
      Cl2 = np.zeros([len(Cl),2],dtype=int)
      Cl2[:,0] = np.arange(len(Cl))
      Cl2[:,1] = len(P_Cl_center2)-1
      NRV = CalcLinesNormalizedVector(P_Cl_center2,Cl2) # Normalized vector from origin to center of line
      for i in range(len(NV)):
        F =NRV[i,:2].dot(NV[i,:2])
        if np.abs(F) < 0.8:
          ST[i] = 0.0
      ST[I_apical] = M_apical
      if VTK_Files[f] == "Embryo_slice_closed":
        ST[I_basal] = M_basal 
      np.savetxt('%s/ST_Edges_Multi_A%.1lf_L%.1lf_B%.1lf_G%.1lf_LatteralRadialOnly.txt' % (folder, M_apical,M_lateral,M_basal,M_general), ST,fmt="%lf") 
    else:
      np.savetxt('%s/ST_Edges_Multi.txt' % (folder), ST,fmt="%lf") 

  # Final output
  #===========================================
  writeVTKFile("%s/Embryo_Slice.vtk" % folder, P=P, Cl=Cl,  Cp_Gr=Cp_Gr) # VTK only to see the 
  np.savetxt('%s/Triangles.txt' % (folder), Ct, fmt="%d") 


'''
Test Code:

import os
import numpy as np
import matplotlib.pyplot as plt


x = np.genfromtxt('solid-gmsh/basal_closed/x.dat')
y = np.genfromtxt('solid-gmsh/basal_closed/y.dat')
z = np.genfromtxt('solid-gmsh/basal_closed/z.dat')

l0 = np.genfromtxt('solid-gmsh/basal_closed/edge1.dat',dtype='int') 
l1 = np.genfromtxt('solid-gmsh/basal_closed/edge2.dat',dtype='int') 

I = np.genfromtxt('solid-gmsh/basal_closed/vc-cylin-lines.dat',dtype='int') 

A = 0
for i in I:
  p0, p1 = l0[i], l1[i]
  A += x[p0]*y[p1]-y[p0]*x[p1]
A = A/2.0
print(A)

plt.plot([x[l0[I[:10]]],x[l1[I[:10]]]],[y[l0[I[:10]]],y[l1[I[:10]]]])
plt.show()

'''