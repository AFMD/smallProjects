#!/usr/bin/env python  2
''''' 
    File name: MNR_Model_NWN.py 
    Authors: C. G. Rocha, C. O'Callaghan, H. G. Manning 
    Date created: 01/10/2014 
    Date last modified: 10/12/2018 
    Python Version: 2.7 
     
    This program is free software: you can redistribute it and/or modify 
    it under the terms of the GNU General Public License as published by 
    the Free Software Foundation, either version 3 of the License, or 
    (at your option) any later version. 
 
    This program is distributed in the hope that it will be useful, 
    but WITHOUT ANY WARRANTY; without even the implied warranty of 
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
    GNU General Public License for more details. 
 
    You should have received a copy of the GNU General Public License 
    along with this program. If not, see <https://www.gnu.org/licenses/>. 
 
    Required packages: numpy, matplotlib, scipy, math, itertools, shapely, networkx, cvxopt, collections. 
    Additionally, this code uses a sparse linear algebra package named symmlq.py developed by the Systems Optimisation Laboratory (SOL) 
    at Stanford University, Department of Management Science and Engineering (MS&E), which is available at 
    <http://stanford.edu/group/SOL/software.html> under the terms of The MIT License (MIT). In some cases, minres.py (also developed 
    by the SOL) provides better results.  
 
Citation: 
    If you use this code in academic publications, please, cite our work appropriately. 
'''  
  

from shapely.geometry import LineString, MultiLineString, MultiPoint, Point  
from shapely.ops import cascaded_union  
from scipy.misc import comb  
from itertools import product  
import scipy.stats as stats  
import matplotlib.pyplot as plt  
import matplotlib.patches as patches  
import math  
import numpy as np  
from itertools import islice  
from cvxopt import matrix, solvers  
from cvxopt.base import sparse  
from cvxopt.base import matrix as m  
from cvxopt.lapack import *  
from cvxopt.blas import *  
import cvxopt.misc as misc  
from pykrylov.symmlq import symmlq
from pykrylov.symmlq import *
from symmlq import *  
#import symmlq
import networkx as nx  
from itertools import islice, combinations  
from collections import Counter, defaultdict  
 

from pykrylov.linop import PysparseLinearOperator
from pykrylov.symmlq import *
import scipy
from scipy.sparse.linalg import *
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import minres

import os.path
import time 
import os
import matplotlib.pyplot as plt
import random
#def Gfun(x,y,trans='N'):
   # '''''Function that passes matrix A to the symmlq routine which solves Ax=B.'''  
    #gemv(Amatrix,x,y,trans)  
 
  
''''' 
Parameters: 
----------- 
   R_junc: float 
           Junction resistance value (in Ohms). 
 
   rho0: float 
         Wire resistivity (in Ohms-um). 
 
  wire_diameter: float 
                  Wire diameter (in nm) 
 
   wire_length: float 
                Wire length (in um). 
 
   box_length: float 
               Squared area of box_length x box_length where wires will be randomly spread (in um).  
               Electrodes are placed on the left and right side of the box and have same dimensions as the size of the box (elec_length). 
 
   samples: int 
            Number of representative NWN samples of the ensemble. 
 
   n_initial: float 
              A scan in wire densities is performed. This is the initial wire density value (in #wires/um^2). 
 
   n_final: float 
            A scan in wire densities is performed. This is the final wire density value (in #wires/um^2). 
 
   nstep: float 
          A scan in wire densities is performed. This is the wire density step (in #wires/um^2). 
 
   distl: boolean  
         Stablishes if a dispersion in wire lengths is to be considered. For a network made of wires all with the same length wire_length, 
         distl = False. Otherwise, distl = True which enables a truncated normal distribution of lengths with average length of lmean, standard deviation 
         of sigmal, truncated at [lower_l,upper_l]. 
          
   A0: float 
       Cross-sectional area of a cylinder (in um^2). 
 
Returns: 'output.txt' (text file) 
         four columns output file containing: wire density (density), average sheet resistance (resAvg), standard deviation in the sheet resistance (resStd), junction density (junc_avg).  
 
'''  
  
# Start ---------- Parameters block -------------  
R_junc = 300100.0   
rho0 = 90 #0.0790 #0.0226   
wire_diameter = 50.0 #30.0  
wire_length= 6.5  #6.0
extinction_coeff = 0.2 #0.1 #0.25 #1 #4.0
box_length = 20.0 #15.0   5x wire length gives good results independent of tol for e-9 to e-15
samples = 1
elec_length = box_length    
box_y = box_length  
lead_sep = box_length  
nstep = 0.16411
n_initial =  7*nstep #1.90079+30*nstep  #0.16411  
n_final =  17*nstep #1.90079+31*nstep  
  
distl = False  
lower_l = 2.2  
upper_l = np.inf  
sigmal = 2.0  
lmean = wire_length  
A0 = math.pi*((wire_diameter*0.001)/2)**2   
# End ---------- Parameters block -------------  


### Parameter to decrease number of nodes, trying to simulate effect of wrapping polymer
percentage_chance = 0.0 ## chance for any occuring node to be deleted from list



# ---------- Parameters for symmlq routine -------------  
tol=1e-10
show=False  
maxit=None   

#-----------  Parameters for Calculation time display --------
start_time = time.clock()
  
  
# ---------- Output file -------------  
res_file = "output.txt"  
if os.path.exists(res_file)==False:
    open(res_file, "w").write("Density AF Transmittance Average_resistance resStdev Junct_density R_junc rho0 wire_diameter box_length samples nstep n_initial n_final tolerance_minres distl lower_l upper_l sigmal calctime\n")
res_dist = open(res_file,"a")  
  
 
  

  
# ---------- Auxiliary lists for ensemble calculation -------------  
res_list=[]  
short_sep_list=[]  
junc_dens=[]  
dens_temp=[]  
avg_res_temp=[]  
st_dev_temp=[]  
  
for density in np.arange(n_initial,n_final,2*nstep):  
    for sample in range(samples):     
        while True:
            try:
                area = box_length**2 # box area (in um^2)  
                box_x = box_length # box width (in um)  
                box_y = box_length # box length (in um)  
                num_junc = 0 # junction counter  
                nwires = area*density # total number of nanowires  
                          
                # Start ---------- Creation of random stick coordinates and electrodes -------------  
                        # a single wire is represented by a set of initial and final coordinates as [(x1,y1),(x2,y2)].  
                x1 = np.random.rand(int(nwires))*box_x  
                y1 = np.random.rand(int(nwires))*box_y  
                length_array = np.zeros(int(nwires))        
          
                if distl == True:  
                    lengths = stats.truncnorm((lower_l - lmean) / sigmal, (upper_l - lmean) / sigmal, loc=lmean, scale=sigmal)  
                    length_array = lengths.rvs(size=nwires)  
                else:  
                    length_array.fill(wire_length)  
                  
                # Sorting the angles that define the wire orientation (in radians from 0 to 2 *pi).  
                theta1 = np.random.rand(int(nwires))*2.0*math.pi  
          
                x2 = length_array * np.cos(theta1) + x1  
                y2 = length_array * np.sin(theta1) + y1  
          
                # Adding to the coordinate list (x1,y1) the points corresponding to the contact leads.  
                x1 = np.insert(x1, 0, 0.0)  
                x1 = np.insert(x1, 0,0)  
          
          
                # Adding to the coordinate list (x2,y2) the points corresponding to the contact leads.  
                x2 = np.insert(x2, 0, 0.0)  
                x2 = np.insert(x2, 0,0)  
                ypostop = box_y/2 + elec_length/2  
                yposbot = box_y/2 - elec_length/2  
                y1 = np.insert(y1, 0,ypostop)  
                y1 = np.insert(y1, 0,ypostop)  
                y2 = np.insert(y2, 0,yposbot)  
                y2 = np.insert(y2, 0, yposbot)  
          
                xposleft = box_x/2-lead_sep/2  
                xposright = box_x/2+lead_sep/2  
          
                x1[0]= xposleft  
                x2[0] = xposleft      
                x1[1] = xposright  
                x2[1] = xposright  
          
                # Merging [(x1,y1),(x2,y2)] in accordance to shapely format.  
                coords1 = zip(x1,y1)  
                coords2 = zip(x2,y2)  
                coords = zip(coords1,coords2)  
                mlines = MultiLineString(coords)  
                nwires_plus_leads = int(nwires+2)  
                        # End ---------- Creation of random stick coordinates and electrodes -------------  
          
          
                        # Start ---------- Identifying intersections between wires -------------  
                # all pair wire combination  
                lines_comb = combinations(mlines, 2)  
                          
                # list storing True or False for pair intersection  
                intersection_check = [pair[0].intersects(pair[1]) for pair in lines_comb]  
          
                # list storing the indexes of intersection_check where the intersection between two wires is TRUE  
                #intersections = [i for i, x in enumerate(intersection_check) if x and random.random()>percentage_chance]  ###intervention ### and random.random() > percentage_chance
                intersections = [i for i, x in enumerate(intersection_check) if x ] #original
                
                #full list containing all non-repeated combinations of wires  
                combination_index = list((i,j) for ((i,_),(j,_)) in combinations(enumerate(mlines), 2))  
          
                # list storing the connection (wire_i, wire_j)   
                intersection_index = [combination_index[intersections[i]] for i in range(len(intersections))]  
          
                # checking the coordinates for interesection points  
                inter_point_coll = [pair[0].intersection(pair[1]) for pair in combinations(mlines, 2)]  
          
                # eliminating empty shapely points from the previous list  
                no_empty_inter_point_coll = [inter_point_coll[intersections[i]] for i in range(len(intersections))]  
          
                # total number of intersections  
                nintersections = len(intersection_index)  
                        # End ---------- Identifying intersections between wires -------------  
          
          
                        # Start ---------- MNR nodal mapping -------------  
                # dictionary containing wire index: [list of wires connected to a given wire]  
                wire_touch_list = defaultdict(list)  
                for k, v in intersection_index:  
                    wire_touch_list[k].append(v)  
                    wire_touch_list[v].append(k)  
          
                # dictionary containing wire index: [label nodes following MNR mapping]  
                wire_touch_label_list = defaultdict(list)  
                each_wire_inter_point_storage = defaultdict(list)  
                label = 2  
          
                # Assigning new node labelling according to MNR mapping  
                for i in iter(wire_touch_list.viewitems()):  
                    for j in range(len(i[1])):  
                        cpoint = mlines[i[0]].intersection(mlines[i[1][j]])  
                        npoint = (cpoint.x,cpoint.y)  
                        each_wire_inter_point_storage[i[0]].append(npoint)  
                      
                        if i[0] > 1:  
                            wire_touch_label_list[i[0]].append(label)  
                            label += 1  
                        else:  
                            wire_touch_label_list[i[0]].append(i[0])  
          
                        maxl = label # dimension of the resistance matrix  
                              
                # flattening intersection_index for counting the amount of occurances of wire i  
                flat = list(sum(intersection_index, ()))   
                conn_per_wire = Counter(flat)  
          
                # checking for isolated wires  
                complete_list = range(nwires_plus_leads)  
                isolated_wires = [x for x in complete_list if not x in flat]  
          
                # list containing the length segments of each wire (if it has a junction)  
                each_wire_length_storage = [[] for _ in range(nwires_plus_leads)]    
          
                # Routine that obtains the segment lengths on each wire  
                for i in each_wire_inter_point_storage:  
              
                    point_ini = Point(mlines[i].coords[0])    
                    point_fin = Point(mlines[i].coords[1])    
                    wlength = point_ini.distance(point_fin)   
                    wire_points = each_wire_inter_point_storage[i]  
          
                    dist = [0.0]*(len(wire_points)+1)  
                    for j in range(len(wire_points)):  
                        point = Point(wire_points[j])  
                        dist[j] = point_ini.distance(point)  
          
                    dist[-1] = wlength    
                    dist.sort()   
          
                    dist_sep = [0.0]*len(dist)  
                    dist_sep[0] = dist[0]  
                    dist_sep[1:len(dist)] = [dist[k]-dist[k-1] for k in range(1,len(dist))]   
                    each_wire_length_storage[i].append(dist_sep)  
                        # End ---------- MNR nodal mapping -------------  
          
                          
                        # The MNR mapping associated to the NWN is also converted into a mathematical graph given by G.  
                        # G contains 2*nintersections nodes and we conventioned that left and right electrodes are labelled as node 0 and 1, respectively.  
                G = nx.Graph()  
                G.add_nodes_from(range(2*nintersections))  
                mr_matrix_plus = np.zeros((2*nintersections,2*nintersections))  
                inner_count = 0  
                inter_count = 0  
                #nx.draw(G)
                #nx.draw_random(G)
                #nx.draw_circular(G)
                nx.draw_spectral(G, node_size= 10)
                ##nx.draw_networkx_nodes(G)
                plt.show()
          
                        # Start ---------- Building resistance matrix -------------  
                for iwire in xrange(nwires_plus_leads):  
                    if each_wire_inter_point_storage[iwire]:  
                        for j, pointj in enumerate(each_wire_inter_point_storage[iwire]):  
                            point = Point(pointj)  
                            for i, pointw in enumerate(each_wire_inter_point_storage[iwire]):  
                                comp_pointw = Point(pointw)  
                                inter_dist = point.distance(comp_pointw)  
                                round_inter_dist = round(inter_dist, 4)  
                                for il in each_wire_length_storage[iwire][0]:  
                                    value = float(il)  
                                    value = round(value,4)  
                                    if value == round_inter_dist and value != 0:  
                                        inner_resis = (float(value) * rho0 / A0)  
                                          
                                        if iwire != 0 and iwire != 1 and mr_matrix_plus[wire_touch_label_list[iwire][i], wire_touch_label_list[iwire][j]] == 0.0:                                  
                                            mr_matrix_plus[wire_touch_label_list[iwire][i], wire_touch_label_list[iwire][j]] = -1.0/inner_resis  
                                            mr_matrix_plus[wire_touch_label_list[iwire][j], wire_touch_label_list[iwire][i]] = -1.0/inner_resis  
                                            G.add_edge(wire_touch_label_list[iwire][i],wire_touch_label_list[iwire][j])  
                                            inner_count += 1  
                            for k, label in enumerate(wire_touch_list[iwire]):  
                                for kk, pointk in enumerate(each_wire_inter_point_storage[label]):  
                                    pointk = Point(pointk)  
                                    inter_dist = point.distance(pointk)  
                                    round_inter_dist = round(inter_dist, 4)  
                                    if round_inter_dist == 0 and mr_matrix_plus[wire_touch_label_list[iwire][j], wire_touch_label_list[label][kk]] == 0:  
                                        G.add_edge(wire_touch_label_list[label][kk],wire_touch_label_list[iwire][j])  
                                        r0 = -1/R_junc  
                                        mr_matrix_plus[wire_touch_label_list[iwire][j], wire_touch_label_list[label][kk]] = r0  
                                        mr_matrix_plus[wire_touch_label_list[label][kk], wire_touch_label_list[iwire][j]] = r0   
                                   
                                                
                sum_rows_mr_plus = mr_matrix_plus.sum(1)  
                np.fill_diagonal(mr_matrix_plus, abs(sum_rows_mr_plus))  
                mr_nozero_rows_plus = mr_matrix_plus[~(mr_matrix_plus==0).all(1),:]  
                          
                        # nonconnected wires are eliminated from the resistance matrix  
                mr_nonconnected_plus = mr_nozero_rows_plus[:,~(mr_nozero_rows_plus==0).all(0)]  
                # End ---------- Building resistance matrix -------------  
          
                        # input current vector  
                i0 = 1.0 # absolute value of the current (in Amp)  
                ic = np.zeros(mr_nonconnected_plus.shape[0])  
                ic[0] = +i0  
                ic[1] = -i0  
                Imatrix = m(ic)  
          
                        # Solving Ohm's law in matrix form, R^(-1)V = I. Resulting voltages are in Volts.  
                #Amatrix = m(mr_nonconnected_plus)  
                
                #Amatrix = np.array(mr_nonconnected_plus) 
                
                #ks = Symmlq(Imatrix)
                #elec_pot_mr = ks.solve(Gfun)
                  
                #print Gfun
                #print Imatrix
                #or
                #ks = Symmlq(Gfun)
                #print Amatrix
                #elec_pot_mr = ks.solve(Imatrix)
                Amatrix = csc_matrix(mr_nonconnected_plus)
                elec_pot_mr = minres(Amatrix, Imatrix, tol=tol)
                
                #elec_pot_mr = Symmlq(Imatrix, Gfun, show=show, rtol=tol, maxit=maxit) 
                
                #elec_pot_mr = minres(Imatrix, Amatrix)  
         
                        # Sheet resistance  
                resistance = ((elec_pot_mr[0][0] - elec_pot_mr[0][1]))/i0  
          
                        # Checking if there is a path connecting electrodes at nodes 0 and 1  
                if nx.has_path(G,0,1):        
                    separation_short = nx.shortest_path_length(G,0,1)     
                    res_list.append(resistance)  
                    short_sep_list.append(separation_short)  
                    junc_dens.append(float(nintersections)/area)  
            except IndexError:
                    continue
            break                      
    AF = density*wire_diameter*wire_length*0.001
    transmittance =    round(math.exp(-AF*extinction_coeff), 4)                    
    junc_avg = np.mean(junc_dens)  
    resAvg = np.mean(res_list)  
    resStd = np.std(res_list)  
    short = np.mean(short_sep_list)   
    dens_temp.append(junc_avg)  
    avg_res_temp.append(resAvg)  
    st_dev_temp.append(resStd)     
    res_dist.write("%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n" %(density,AF,transmittance,resAvg,resStd,junc_avg,R_junc,rho0,wire_diameter,box_length,samples,nstep,n_initial,n_final,tol,distl,lower_l,upper_l,sigmal,round(time.clock() - start_time, 5)))
    print("Density: %s, Transmittance: %s, Average resistance: %s, Standard deviation: %s, Junction density: %s" %(density,transmittance,resAvg,resStd,junc_avg))
    print"runtime was", round(time.clock() - start_time, 5), "seconds"
    res_list=[]  
    short_sep_list=[]  
    junc_dens=[]  
                    
  
res_dist.close() 
duration = 0.1
freq = 1100
#os.system('play -nq -t alsa synth {} sine {}'.format(duration, freq))