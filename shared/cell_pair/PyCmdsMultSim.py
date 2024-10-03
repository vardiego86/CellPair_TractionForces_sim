# Python commands used in cell pair simulations.
#
# This file contains the commands written in Python that are used in the simulation script found in multCellSpr.py
# 
# Specific commands:
#
# UpdateStressFibListCmd:   Update SF list by creating new SFs and deleting SFs. Updates all the SF measurements per node.
# AbsoluteValuePCArray:     Get absolute value of particle container array and store it in same array.
# PatternSubstrateRectCmd:  Determine triangles in which FA can be formed (rectangle pattern).
#
###
#
# Group:        MAtrix - Tissue Engineering and Mechanobiology
# Department:   Mechanical Engineering (Biomechanics Division)
# Institution:  KU Leuven
#
# Diego A. Vargas
# Jan 24, 2020

# mpacts imports
import mpacts.core.command as cmd   

# Other imports
import numpy as np
import random
import os

###

# Python command definitions

class UpdateStressFibListCmd( cmd.PythonCommand):
  def __init__(self, name, sim, x, atFA, atAJ, minLFib, orFib , k_on, k_off, atFib, L_fib, dir_fib, typ_fib, mult_F_am, ppidx, SF_ContactDetector, fibers, t_int, gate):
    cmd.PythonCommand.__init__(self, name, sim("loop_cmds/contact_cmds"))
    self.sim=sim.root_node()                                                                                                                                                       
    self.x = x                      # Cell nodes positions
    self.atFA = atFA                # Boolean value per node indicating if FA is present (i.e. adhered to substrate)
    self.atAJ = atAJ                # Boolean value per node indicating if AJ is present (i.e. adhered to substrate)
    self.minLFib = minLFib          # Minimum SF length for formation
    self.orFib = orFib              # Angle with respect to x-axis above which a SF cannot form
    self.k_on = k_on                # Probability of creating a SF
    self.k_off = k_off              # Probability of deleting a SF
    self.atFib = atFib              # Boolean value per node indicating if SF is present
    self.L_fib = L_fib              # Length of SF
    self.dir_fib = dir_fib          # Unit vector describing direction of SF; direction is node1 -> node2)
    self.typ_fib = typ_fib          # Value per node defining type of connection with SF: (0) node-node (1) node-FA-FA (2) node- FA-AJ (3) node-FA-node(no adh) (4) node-AJ-FA (5) node-AJ-node(no adh)
    self.mult_F_am = mult_F_am      # Strengthening factor of SF
    self.ppidx = ppidx              # Cell (of cell pair) to which the node belongs to (0 or 1)
    self.SF_cd = SF_ContactDetector # SF contact detector: To add and remove contacts (i.e. node pairs connected by SFs)
    self.fibers = fibers            # SF list solely for display/visualization purposes (through VTP files)
    self.t_int = t_int              # time step (UNUSED in this simulation)
    self.set(gate = gate)           # gate dictates in what time steps or how often is python command implemented

  def execute(self):
    # timestep to scale SF formation/deletion probabilities (UNUSED in this simulation)
    dt = self.t_int.get_value()

    # Find nodes that could potentially be connected by a SF (i.e. nodes attached to a FA or AJ but not yet a SF)
    adhPerCellLists = [[] for _ in range(self.x.parent_node().parent_node().size())]        # Makes as many empty lists per cell
    for i in range(len(self.atFA)):                                                         # Loop though all FAs
        if ( (self.atFA[i] > 0.0 or self.atAJ[i] > 0.0 ) and self.atFib[i] == False ):      # If FA meets criteria
            myppidx = self.ppidx[i]
            adhPerCellLists[myppidx].append(i)                                              # Populates appropriate list (i.e. for each cell described by parent index)

    # Create SF based on k_on (set to 1), length, orientation
    for adhNodeList in adhPerCellLists:                     # For adhered nodes (i.e. FAs) in each cell:
        while len(adhNodeList) > 1:
            nodeIdx1 = random.randint(0,len(adhNodeList)-1) # Pick node1 at random 
            node1 = adhNodeList[nodeIdx1]
            del adhNodeList[nodeIdx1]

            nodeIdx2 = random.randint(0,len(adhNodeList)-1) # Pick node2 at random
            node2 = adhNodeList[nodeIdx2]
            del adhNodeList[nodeIdx2]

            # Calculate unit vector (node1->node2), length, and orientation of potential SF
            v = ( self.x[node1][0]-self.x[node2][0] , self.x[node1][1]-self.x[node2][1] , self.x[node1][2]-self.x[node2][2] )
            lenFib = np.sqrt( v[0]**2 + v[1]**2 + v[2]**2  )
            lenFibProj = np.sqrt( v[0]**2 + v[2]**2  )
            angFibAx = np.arccos(v[0]/lenFibProj)

            if ( (angFibAx <= self.orFib or angFibAx >= np.pi-self.orFib) and lenFib >= self.minLFib and random.random() < 1.0-(1.0-self.k_on)**dt ):   # If potential SF meets criteria
                self.atFib[node1] = True
                self.atFib[node2] = True

                self.L_fib[node1] = lenFib
                self.L_fib[node2] = lenFib

                self.dir_fib[node1] = tuple(x/-lenFib for x in v)
                self.dir_fib[node2] = tuple(x/lenFib for x in v)

                if (self.atFA[node1]==True):
                    if (self.atFA[node2]==True):
                        self.typ_fib[node1]=1.0
                        self.typ_fib[node2]=1.0
                    elif (self.atAJ[node2]==True):
                        self.typ_fib[node1]=2.0
                        self.typ_fib[node2]=4.0
                elif (self.atAJ[node1]==True):
                    self.typ_fib[node1]=4.0
                    self.typ_fib[node2]=2.0
                else:
                    print ("Forgot some case when CREATING a fiber.")

                self.SF_cd.addContact(node1,node2)

    # Delete some of existing SFs:
    fiberList = self.SF_cd.get_contact_data().get_contact_list()
    for idx in range( len(fiberList) ):     # Loop though existing SFs 
        select_fiber = fiberList[idx]
        node1 = select_fiber[0]
        node2 = select_fiber[1]
        # - delete if fiber not attached to at least 1 node also attached to a FA or AJ
        if ( ((self.atFA[select_fiber[0]]==0.0 and self.atAJ[select_fiber[0]]==0.0) or (self.atFA[select_fiber[1]]==0.0 and self.atAJ[select_fiber[1]]==0.0)) ):
            self.atFib[node1] = False
            self.atFib[node2] = False
            
            self.L_fib[node1] = 0
            self.L_fib[node2] = 0
            
            self.dir_fib[node1] = (0.,0.,0.)
            self.dir_fib[node2] = (0.,0.,0.)

            self.mult_F_am[node1] = 1.0
            self.mult_F_am[node2] = 1.0

            self.typ_fib[node1] = 0.0
            self.typ_fib[node2] = 0.0

            self.SF_cd.deleteContact( *(select_fiber) )
        # - else update length, unit vector, and classification of node type (this last for analysis purposes)
        else:
            v = ( self.x[node1][0]-self.x[node2][0] , self.x[node1][1]-self.x[node2][1] , self.x[node1][2]-self.x[node2][2] )
            lenFib = np.sqrt( v[0]**2 + v[1]**2 + v[2]**2  )

            self.L_fib[node1] = lenFib
            self.L_fib[node2] = lenFib

            self.dir_fib[node1] = tuple(x/-lenFib for x in v)
            self.dir_fib[node2] = tuple(x/lenFib for x in v)

            if (self.atFA[node1]==True and self.atFA[node2]==True):
                self.typ_fib[node1] = 1.0
                self.typ_fib[node2] = 1.0
            elif (self.atFA[node1]==True and self.atAJ[node2]==True):
                self.typ_fib[node1] = 2.0
                self.typ_fib[node2] = 4.0
            elif (self.atFA[node1]==True and (self.atFA[node2]==False and self.atAJ[node2]==False) ):
                self.typ_fib[node1] = 3.0
                self.typ_fib[node2] = 0.0
            elif (self.atAJ[node1]==True and self.atFA[node2]==True):
                self.typ_fib[node1] = 4.0
                self.typ_fib[node2] = 2.0
            elif (self.atAJ[node1]==True and (self.atFA[node2]==False and self.atAJ[node2]==False) ):
                self.typ_fib[node1] = 5.0
                self.typ_fib[node2] = 0.0
            else:
                print ("Forgot some case when UPDATING a fiber.")

    # Update SF list (visualization purposes)
    fiberList = self.SF_cd.get_contact_data().get_contact_list()
    self.fibers.resize(0)
    if fiberList:
        self.fibers('vertexIndices').set_array(fiberList)

###

class AbsoluteValuePCArray( cmd.PythonCommand):
  def __init__(self, name, sim, Ldiff, gate):
    cmd.PythonCommand.__init__(self, name, sim)
    self.sim=sim.root_node()
    self.Ldiff = Ldiff
    self.set(gate = gate)

  def execute(self):
    for idx in range( len(self.Ldiff) ):
        self.Ldiff[idx] = abs(self.Ldiff[idx])

###

class PatternSubstrateRectCmd( cmd.PythonCommand):
  def __init__(self, name, sim, ratio, ldim,  p_trng, p_vrt, lig, gate):
    cmd.PythonCommand.__init__(self, name, sim("loop_cmds/pre_body_force_cmds"))
    self.sim=sim.root_node()                                                                                                                                                       
    self.ratio = ratio          # Ratio (Length / Width) of rectangular pattern
    self.ldim = ldim            # Length of rectangular pattern
    self.p_trng = p_trng        # Substrate triangles
    self.p_vrt = p_vrt          # Substrate triangle vertices
    self.lig = lig              # Boolean array dictating whether triangle is in pattern (1) or not (0)
    self.set(gate = gate)       # gate dictates in what time steps or how often is python command implemented

  def execute(self):
    l = self.ldim.get_value()
    R = self.ratio.get_value()
    
    w = l/R
    lpat = l/2.0
    wpat = w/2.0

    for tri_it in range(len(self.p_trng)):  # Loop through all substrate triangles
        # Get maximum X and Z (axes that lie on substrate plane) positions from any of the vertices defining the triangle
        maxXvt = max( abs(self.p_vrt[self.p_trng[tri_it][0]][0]) , abs(self.p_vrt[self.p_trng[tri_it][1]][0]) , abs(self.p_vrt[self.p_trng[tri_it][2]][0])  )
        maxZvt = max( abs(self.p_vrt[self.p_trng[tri_it][0]][2]) , abs(self.p_vrt[self.p_trng[tri_it][1]][2]) , abs(self.p_vrt[self.p_trng[tri_it][2]][2])  )
        # If X-position is within length range AND Z-position within width range
        if ( maxXvt <= lpat and maxZvt <= wpat):
            # Make triangle part of pattern (i.e. set lig to 1)
            self.lig[tri_it] = 1.0

