import numpy as np
import DEMutilities.postprocessing.h5storage as h5s
import DEMutilities.postprocessing.operators.operators as op
import DEMutilities.numeric_functions as nf
import mpacts.io.datasave as ds
import time

import os

# Select cell
cell_id = 0

###

# Unpack archive (.a) file
currFolder = os.path.split(os.getcwd())[1]
folderName = os.environ['VSC_DATA'] + '/006a/' + currFolder
h5FileName = os.environ['VSC_SCRATCH'] + '/006a/' + currFolder + '/simulation_cell_' + str(cell_id) + '.h5'

###

def makeNumpyArray( x_elem, y_elem, z_elem):
    return np.array([x_elem, y_elem, z_elem])

def inTriangle(vtx_idx_sb, vtx_pos_sb, node_pos_pl):
    
    inTriArr = np.zeros((len(vtx_idx_sb)))

    for i in range(len(vtx_idx_sb)):
        for j in range(len(node_pos_pl)):

            pt = [node_pos_pl[j][0],0.0,node_pos_pl[j][2]]
            v0 = vtx_pos_sb[vtx_idx_sb[i][2]] - vtx_pos_sb[vtx_idx_sb[i][0]]
            v1 = vtx_pos_sb[vtx_idx_sb[i][1]] - vtx_pos_sb[vtx_idx_sb[i][0]]
            v2 = pt - vtx_pos_sb[vtx_idx_sb[i][0]]

            dot00 = np.dot(v0,v0);
            dot01 = np.dot(v0,v1);
            dot02 = np.dot(v0,v2);
            dot11 = np.dot(v1,v1);
            dot12 = np.dot(v1,v2);

            invDenom = 1.0/(dot00*dot11 - dot01*dot01);
            a = (dot11*dot02 - dot01*dot12)*invDenom;
            b = (dot00*dot12 - dot01*dot02)*invDenom;

            nodeInTri = ( a>0.0 and b>0.0 and (a+b)<1 )

            if nodeInTri:
                inTriArr[i] = True
                continue
    
    return inTriArr

def compute_eigvals(Mxx,Mxz,Mzx,Mzz):
    return np.linalg.eigvals( [[Mxx,Mxz],[Mzx,Mzz]] )

def get_major_dipole(dipole):
    return dipole[np.argmax(np.abs(dipole))]

def get_minor_dipole(dipole):
    return dipole[np.argmin(np.abs(dipole))]

###

# Make analysis container with sets of commands
a = op.AnalysisContainer()

# Read mpacts simulation frames
dr = ds.DataReader( "simulation", folder=folderName )
#print dr[0] # See list of available data

# Get 'time from a simulation frame'
ti = op.GetData( a, 'time')
t0 = ti.execute( dr[0] )()

# Make a mask (selecting time)
final_frame = dr[-1].index
frame = op.GetData(a, 'DataFrameIndex')
#time_mask = frame > final_frame-5 #Average over the last 5 output steps
#time_mask = ti > -1.0
time_mask = ti < 14400.0

# Print frame index during analysis
op.FrameIndexPrinter(a, final_frame)

# Select nodes per cell
ppidx_node = op.GetData(a, 'cell/nodes/parentIndex').execute(dr[0])()
ppidx_tri = op.GetData(a, 'cell/triangles/parentIndex').execute(dr[0])()

node_ids = np.where(ppidx_node==cell_id)
tri_ids = np.where(ppidx_tri==cell_id)

# Get metrics:
# Center of mass (CoM) at interface
pos = op.GetData(a, 'cell/nodes/x')  
pos_cell = pos[node_ids]
Lrup = op.GetData( a, 'params/Lo_max_FA').execute(dr[0])()
node_pos_pl = pos_cell[pos_cell[:,1]<Lrup]
node_pos_pl_tissue = pos[pos[:,1]<Lrup]
mean_pos_pl_cell = op.Function(a, np.mean, node_pos_pl , 0 )
mean_pos_pl_tissue = op.Function(a, np.mean, node_pos_pl_tissue , 0 )

##

# get nodes with AJs within a distance of substrate (5.5X Lrup)
node_AJ_bool = op.GetData(a, 'cell/nodes/atAJ')
pos_iF = pos[node_AJ_bool==1]                       # postion of nodes with AJs (iF = interFace)
pos_iF_subs = pos_iF[pos_iF[:,1] < 5.5*Lrup]        # consider only those very close to substrate
# select those within a distance (1.5e-6) of z=0
subCnt_low_bool = pos_iF_subs[:,2] <= 1.5e-6        # near substrace central relative to z-axis (subCnt)
subCnt_high_bool = pos_iF_subs[:,2] >= -1.5e-6
subCnt_bool = op.Function(a, np.logical_and, subCnt_low_bool, subCnt_high_bool)
pos_iF_subCnt = pos_iF_subs[subCnt_bool==1]         # position of nodes with AJs considered central
# Average X-position
xyz_iF = op.Function(a, np.mean, pos_iF_subCnt, 0)                                      # average position of selected nodes
mean_pos_pl = op.Function(a, makeNumpyArray, xyz_iF[0], mean_pos_pl_tissue[1], 0.0 )    # calculate origin of cell pair based on interface
#CoM_pl = op.Recorder(a, mean_pos_pl)

##

# Cell contact area (Ac)
Ac_tri = op.GetData(a, 'cell/triangles/contact_area')
Ac_cell = op.Function(a, np.sum, Ac_tri[tri_ids], 0)
Ac_t = op.Recorder(a, Ac_cell)

# Number focal adhesions (FA)
node_FA_bool = op.GetData(a, 'cell/nodes/atFA')
nFA = op.Function(a, np.sum, node_FA_bool[node_ids], 0)
nFA_t = op.Recorder(a, nFA)

# Find triangles in plane in contact with cell
vtx_idx_sb = op.GetData(a, 'ground/triangles/vertexIndices').execute(dr[0])()
vtx_pos_sb = op.GetData(a, 'ground/controlPoints/x').execute(dr[0])()

inTriArr = op.Function(a, inTriangle, vtx_idx_sb, vtx_pos_sb, node_pos_pl)

# Total in-plane traction
F_tri = op.GetData(a, 'ground/triangles/FprimOut')
A_tri = op.GetData(a, 'ground/triangles/area')

Fmag_tri = op.Function(a, np.sqrt , (F_tri[:,0]**2 + F_tri[:,2]**2) ) # sum of X and Z components (in plane)
Fmag_plane = op.Function(a, np.sum, Fmag_tri[inTriArr], 0)
tot_F_cell = op.Recorder(a, Fmag_plane)

Tmag_tri = op.Function(a, np.divide, Fmag_tri, A_tri)
Tmag_plane = op.Function(a, np.sum, Tmag_tri[inTriArr], 0)
tot_T_cell = op.Recorder(a, Tmag_plane)

# Vector sum of cell-substrate (cs) force
#Fcs_cell = op.Function(a, np.sum, F_tri[inTriArr], 0)
#Fcs_t = op.Recorder(a, Fcs_cell)

# Force dipole
#x_tri = op.GetData(a, 'ground/triangles/x')
#xR_tri = op.Function(a, np.subtract,  x_tri , mean_pos_pl_tissue)  # distance relative to CoM on plane

#xR_tri_cell = xR_tri[inTriArr]
#F_tri_cell = F_tri[inTriArr]

#Mxx = op.Function(a, np.dot, xR_tri_cell[:,0] , F_tri_cell[:,0] )
#Mxz = op.Function(a, np.dot, xR_tri_cell[:,0] , F_tri_cell[:,2] )
#Mzx = op.Function(a, np.dot, xR_tri_cell[:,2] , F_tri_cell[:,0] )
#Mzz = op.Function(a, np.dot, xR_tri_cell[:,2] , F_tri_cell[:,2] )

#dipole = op.Function(a, compute_eigvals , Mxx,Mxz,Mzx,Mzz )

#maj_dip = op.Function(a, get_major_dipole, dipole)

#min_dip = op.Function(a, get_minor_dipole, dipole)

#Rdip = op.Function(a, np.divide, maj_dip, min_dip)
#Rdip_t = op.Recorder(a, Rdip)
#Mxx_t = op.Recorder(a, Mxx)

# Cell-substrate force
Ffa = op.GetData(a, 'cell/nodes/Fcs')
Ffa_cell = Ffa[node_ids]

atFA_cell = node_FA_bool[node_ids]
Ffa_cell = Ffa_cell[atFA_cell==1]                                                      

summed_Ffa_cell = op.Function(a, np.sum, Ffa_cell , 0 )
Ffa_t = op.Recorder(a, summed_Ffa_cell)

Ffa_mag = op.Function(a, np.sqrt , (Ffa_cell[:,0]**2 + Ffa_cell[:,1]**2 + Ffa_cell[:,2]**2) )
summed_Ffa_mag = op.Function(a, np.sum, Ffa_mag, 0 )
tot_Ffa_mag = op.Recorder(a, summed_Ffa_mag)

# Cell-substrate force -> Mxx
F_FA_cell = Ffa_cell

pos_FA_cell = pos_cell[atFA_cell==1]
r_FA_cell = op.Function(a, np.subtract, pos_FA_cell , mean_pos_pl ) 

Mxx = op.Function(a, np.dot, r_FA_cell[:,0] , -1.0*F_FA_cell[:,0] )
Mxx_t = op.Recorder(a, Mxx)

# Length cell pair                                                                                                  
max_x_pos_FA = op.Function(a, np.amax, pos_FA_cell[:,0] )
min_x_pos_FA = mean_pos_pl[0]
L_cell_FA = max_x_pos_FA - min_x_pos_FA
L_cell_t = op.Recorder(a, L_cell_FA)

# Intercellular force
Fcc = op.GetData(a, 'cell/nodes/Fcc')
Fcc_cell = Fcc[node_ids]

node_AJ_bool = op.GetData(a, 'cell/nodes/atAJ')
atAJ_cell = node_AJ_bool[node_ids]
Fcc_cell = Fcc_cell[atAJ_cell==1]

summed_Fcc_cell = op.Function(a, np.sum, Fcc_cell , 0 )
Fcc_t = op.Recorder(a, summed_Fcc_cell)

Fcc_mag = op.Function(a, np.sqrt , (Fcc_cell[:,0]**2 + Fcc_cell[:,1]**2 + Fcc_cell[:,2]**2) )
summed_Fcc_mag = op.Function(a, np.sum, Fcc_mag, 0 )
tot_Fcc_mag = op.Recorder(a, summed_Fcc_mag)

# Perform analysis; loop over data frames in dr
a.loop( dr, mask_function = time_mask  )

#--------------------------------------------------------------

storage = h5s.H5Storage( h5FileName, 'a' )

results = storage.data_section( "results/geom", overwrite=True)
#results.add_data( "CoM_pl", CoM_pl()
#                , description = 'CoM of points near surface'
#                , label='$CoM on surface$'
#                , unit = 'm')

results.add_data( "L_cell", L_cell_t()
                , description = 'Length of cell in x-axis'
                , label='$L_cell$'
                , unit = 'm')

results.add_data( "Ac_t", Ac_t()
                , description = 'cell contact area'
                , label='$Ac$'
                , unit = 'm*m')

results = storage.data_section( "results/adh", overwrite=True)
results.add_data( "nFA", nFA_t()
                , description ='number of focal adhesions'
                , label='$n_F_A$'
                , unit = '-')

results.add_data( "Ffa", Ffa_t()
                , description = 'Vector sum force along focal adhesions'
                , label='$F_F_A$'
                , unit = 'N')

results.add_data( "force_fa_t", tot_Ffa_mag()
                , description = 'Total force-substrate magnitude'
                , label='$F_F_Amag$'
                , unit = 'N')

results.add_data( "Fcc", Fcc_t()
                , description = 'Vector sum intercellular force'
                , label='$F_c_c$'
                , unit = 'N')

results.add_data( "force_cc_t", tot_Fcc_mag()
                , description = 'Total intercellular force magnitude'
                , label='$F_c_cmag$'
                , unit = 'N')

results = storage.data_section('results/trac',overwrite=True)
results.add_data( 'force_t', tot_F_cell()
                , description = 'total force exerted by cell on the subrate in the interface plane'
                , label = '$F_cell$'
                , unit = 'N' )

#results.add_data( 'Fcs', Fcs_t()
#                , description = 'vector sum of  forces exerted by cell on the subrate'
#                , label = '$F_c_s$'
#                , unit = 'N' )

results.add_data( 'trac_t', tot_T_cell()
                , description = 'total traction exerted by cell oon the subrate in the interface plane'
                , label = '$T_cell$'
                , unit = 'Pa' )

#results.add_data( 'Rdip_t', Rdip_t()
#                , description = 'Dipole ratio'
#                , label = '$Rdip_t$'
#                , unit = '-' )

results.add_data( 'Mxx_t', Mxx_t()
                , description = 'cell moment in x'
                , label = '$Mxx$'
                , unit = 'J' )


# Close hdf5 file
storage.close()

