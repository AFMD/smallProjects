##read data from simulation output files to plot and analyse


import numpy as np
import matplotlib.pyplot as plt
from statistics import mean


#User enters file name and varied parameter, and adjusts plot/axis labels


#****************************************
##Enter file name ###
with open('output_vary_wire_lengths_0.8to1.2um_0.1um_steps_d1nm_bl5um_uniforml_15samples_final copy.txt','r') as f: 
    data = f.readlines()
    #print(data)

#arrays to store data points #list entry 
network_dens_list=[] #0
AF=[] #1
transmittancelist=[] #1
resistancelist=[] #3
resistance_stdev_list=[] #4
jcn_dens_list=[] #5

Rjcn=[] #6
rho0=[] #7
wire_diameter=[] #8
wire_length=[] #9
box_length=[] #10
samples=[] #11

nstep=[]  #12
n_initial=[] #13
n_final=[] #14
tol_minres=[] #15
disl=[] #16
lower_l=[] #17
sigmal=[] #18
jcn_removal=[] #19

calctime=[] #20


i=0
j=-1
d_prev=100
for line in data:
    
    words = line.split()
    
    if i>0:
        dens=float(words[0])
        #print(dens)
        if dens<d_prev:
             network_dens_list.append([])
             AF.append([])
             transmittancelist.append([])
             resistancelist.append([])
             resistance_stdev_list.append([])
             jcn_dens_list.append([])
             calctime.append([])
             #parameters below constant for each density loop
             Rjcn.append(words[6])
             rho0.append(words[7])
             wire_diameter.append(words[8])
             wire_length.append(words[9])
             box_length.append(words[10])
             samples.append(words[11])
             nstep.append(words[12])
             n_initial.append(words[13])
             n_final.append(words[14])
             tol_minres.append(words[15])
             disl.append(words[16])
             lower_l.append(words[17])
             sigmal.append(words[18])
             jcn_removal.append(words[19])
             j=j+1
             
        resAvg=float(words[3])
        

        if np.isnan(resAvg) != True:
            network_dens_list[j].append(words[0])
            AF[j].append(words[1])
            transmittancelist[j].append(words[2])
            resistancelist[j].append(words[3])
            resistance_stdev_list[j].append(words[4])
            jcn_dens_list[j].append(words[5])
            calctime[j].append(words[20])
            
        d_prev=dens
    i=i+1
   
    

#convert to arrays 
for i in np.arange(0,j+1,1):
    network_dens_list[i]=np.array(network_dens_list[i], dtype=np.float64)
    AF[i]=np.array(AF[i], dtype=np.float64)
    transmittancelist[i]=np.array(transmittancelist[i], dtype=np.float64)
    resistancelist[i]=np.array(resistancelist[i], dtype=np.float64)
    resistance_stdev_list[i]=np.array(resistance_stdev_list[i], dtype=np.float64)
    jcn_dens_list[i]=np.array(jcn_dens_list[i], dtype=np.float64)

Rjcn=np.array(Rjcn, dtype=np.float64)
rho0=np.array(rho0, dtype=np.float64)
wire_diameter=np.array(wire_diameter, dtype=np.float64)
wire_length=np.array(wire_length, dtype=np.float64)
box_length=np.array(box_length, dtype=np.float64)
samples=np.array(samples, dtype=np.float64)
    


#*******************************************
### Set varied parameter as legend label ###
data_label=[]
for i in np.arange(0,j+1,1):
    data_label.append('length=%s um' %(wire_length[i]))

print(data_label)


#T vs Rs plot and fit
from scipy.optimize import curve_fit
Z0=377
def T_perc_func(r, p, n):
    return (1+(1/p)*((Z0/r)**(1/(1+n))))**(-2)
def T_bulk_func(r,sratio):
    return (1+(Z0/(2*sratio*r)))**(-2)

#may need to adjust if further colors necessary 
colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']




for i in np.arange(0,j+1,1):
    popt_perc, pcov_perc = curve_fit(T_perc_func, resistancelist[i], transmittancelist[i])
    popt_perc
    popt_bulk, pcov_bulk = curve_fit(T_bulk_func, resistancelist[i], transmittancelist[i])
    popt_bulk

    resistancelist_srt=np.sort(resistancelist[i])
    #print(resistancelist_srt)
    reslistlength=len(resistancelist_srt)
    res_start=resistancelist_srt[0]
    res_end=resistancelist_srt[reslistlength-1]
    res_step= (res_end - res_start)/25
    
    #print(res_start, res_end, res_step)
    resfitlist=[]

    
    for k in np.arange(res_start,res_end + (res_step/2),res_step):
        resfitlist.append(k)
    
    #print(resfitlist)
    resfitlist=np.array(resfitlist, dtype=np.float64)

    plotcolor=colors[i]
    
    plt.plot(resfitlist, T_perc_func(resfitlist, *popt_perc), plotcolor, linestyle='-', label='Percolative fit: \u03A0 =%5.3f, n=%5.3f' % tuple(popt_perc))
    plt.plot(resfitlist, T_bulk_func(resfitlist, *popt_bulk), plotcolor, linestyle='--', label='Bulk fit: \u03C3 ratio=%5.3f' % tuple(popt_bulk))
    plt.plot(resistancelist[i], transmittancelist[i], plotcolor, marker='o', linestyle='None', label= data_label[i])
   
    
    
plt.title('T vs Rs')
plt.ylabel('T')
plt.xlabel('Rs (Ohm/sq) - Log scale')
plt.xscale('log')
plt.legend(loc='upper left', bbox_to_anchor=(1.0, 1.0))
plt.tight_layout()
#plt.close()
plt.show()



#Plot resistanceStdev vs resistance
for i in np.arange(0,j+1,1):
    plotcolor=colors[i]
    plt.plot(resistancelist[i], resistance_stdev_list[i], plotcolor, marker='o', linestyle='-',label= data_label[i])

plt.title('St.Dev. in Rs vs Rs')
plt.xlabel('Rs (Ohm/sq) - Log scale')
plt.ylabel('Standard Deviation in Rs (Ohm/sq)')
plt.xscale('log')
plt.legend(loc='upper left', bbox_to_anchor=(1.0, 1.0))
plt.tight_layout()
plt.show()

#plot StDevRs/Rs vs network density
for i in np.arange(0,j+1,1):
    ratio_list=[]
    ratio_list = [k / j for k, j in zip(resistance_stdev_list, resistancelist)]
    plotcolor=colors[i]
    plt.plot(network_dens_list[i], ratio_list[i], plotcolor,marker='o', linestyle='-', label= data_label[i])

plt.title('St.Dev.Rs/Rs vs Network Density')
plt.ylabel('St.Dev.Rs/Rs')
plt.xlabel('Network Density (wires/um^2)')
plt.legend(loc='upper left', bbox_to_anchor=(1.0, 1.0))
plt.tight_layout()
plt.show()

#Plot jcn density vs network density
for i in np.arange(0,j+1, 1):
    plotcolor=colors[i]
    plt.plot(network_dens_list[i], jcn_dens_list[i],plotcolor, marker='o', linestyle='-',label=  data_label[i])

plt.title('Junction Density vs Network Density')
plt.xlabel('Network Density (wires/um^2)')
plt.ylabel('Junction Density (junctions/um^2)')
plt.legend(loc='upper left', bbox_to_anchor=(1.0, 1.0))
plt.tight_layout()
plt.show()




#fit for Rs ~l^-1.5
def Rs_l_fit(l,a,b):
    return (a*l**(b))

min_res_list=[]

#Plot Rs with max T value (min Rs) vs wire length
for i in np.arange(0,j+1, 1):
    plotcolor=colors[i]
    resistance_srt=np.sort(resistancelist[i])
    res_length=len(resistance_srt)
    min_res_list.append(resistance_srt[0])

    plt.plot(wire_length[i], resistance_srt[0],plotcolor, marker='o', linestyle='-',label=  data_label[i])
    
#    # excluded for now
#    # shows all Rs points stacked vertically for each wire length
#   # for k in np.arange(1,res_length-2,1):
#       # plt.plot(wire_length[i], resistance_srt[k],plotcolor, marker='o', linestyle='-')

'''
print(min_res_list)
popt, pcov = curve_fit(Rs_l_fit, wire_length, min_res_list)
popt

min_resistance_srt=np.sort(min_res_list)
reslistlength=len(min_resistance_srt)
res_start=min_resistance_srt[0]
res_end=min_resistance_srt[reslistlength-1]
res_step= (res_end - res_start)/25
print(res_start, res_end, res_step)

minresfitlist=[]    
for k in np.arange(res_start,res_end + (res_step/2),res_step):
    minresfitlist.append(k)
    
#print(resfitlist)
minresfitlist=np.array(minresfitlist, dtype=np.float64)

plt.plot(minresfitlist, Rs_l_fit(minresfitlist, *popt), plotcolor, linestyle='-', label='Rs= a*(wire_length)^b : a =%5.3f, b=%5.3f' % tuple(popt))
'''

plt.title('Minimum Rs vs Wire Length')
plt.xlabel('Wire Length (um)')
plt.ylabel('Rs (Ohm)')
plt.legend(loc='upper left', bbox_to_anchor=(1.0, 1.0))
plt.tight_layout()
plt.show()




#Plot max Jcn density vs wire length
for i in np.arange(0,j+1, 1):
    plotcolor=colors[i]
    #find endpoint
    jcnlist_len=len(jcn_dens_list[i])

    plt.plot(wire_length[i], jcn_dens_list[i][jcnlist_len-1],plotcolor, marker='o', linestyle='-',label=  data_label[i])

plt.title('Maximum Junction Density vs Wire Length')
plt.xlabel('Wire Length (um)')
plt.ylabel('Junction Density (junctions/um^2)')
plt.legend(loc='upper left', bbox_to_anchor=(1.0, 1.0))
plt.tight_layout()
plt.show()

#Convert T  to T^(-1/2)-1 and take log of both arrays
for i in np.arange(0,j+1,1):
    transmittancelist[i]=np.log((transmittancelist[i]**(-1/2))-1)
    resistancelist[i]=np.log(resistancelist[i])
    

#function for line fit
def best_fit_slope_and_intercept(xs,ys):
    m = (((mean(xs)*mean(ys)) - mean(xs*ys)) /
         ((mean(xs)*mean(xs)) - mean(xs*xs)))

    b = mean(ys) - m*mean(xs)
    
    return m, b

#line fit and plot data on log scale
for i in np.arange(0,j+1,1):
    m, b = best_fit_slope_and_intercept(resistancelist[i],transmittancelist[i])
    #print(m,b)
    #plot best fit line on graph
    regression_line = [(m*x)+b for x in resistancelist[i]]
    plotcolor=colors[i]
    plt.plot(resistancelist[i],transmittancelist[i], plotcolor, marker= 'o', linestyle='None', label=  data_label[i])
    plt.plot(resistancelist[i], regression_line, plotcolor, linestyle='-', label='Line Fit y = %s x + %s' %(round(m,3),round(b,3)))


    
plt.title('Log(T^(-1/2)-1) vs Log(Rs) with Line Fit')
plt.ylabel('Log(T^(-1/2)-1)')
plt.xlabel('Log(Rs)')
plt.legend(loc='upper left', bbox_to_anchor=(1.0, 1.0))
plt.tight_layout()
plt.show()
