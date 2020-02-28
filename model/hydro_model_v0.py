import sys
import numpy as np
import datetime
sys.path.append('./config/')
import gl
import pandas as pd
from netCDF4 import Dataset

#*****************************************************************************
#definition of paramater and input files 
    
file_modset = "./config/modset.dat"


file_modpar = sys.argv[1] 
file_init = sys.argv[2] 
file_input = sys.argv[3] 
spinup = sys.argv[4]

if file_modpar=="default":
    file_modpar = "./config/modpar.dat"
if file_init=="default":
    file_init="./config/init.dat"


outputfiles=[]
for i in range(5,len(sys.argv)):
    outputfiles=outputfiles+[sys.argv[i]]



#*****************************************************************************
# functions


def read_modset(file_modset):
    print("reading model setup from "+file_modset+"...")
    with open(file_modset, "r") as fmodset:
        #convergence criterion
        gl.convcrit=float(fmodset.readline().strip().split(",")[1])
        #maxiterations
        gl.maxiter=int(fmodset.readline().strip().split(",")[1])
        
        # numbers of cells are read
        gl.nofswcells=int(fmodset.readline().strip().split(",")[1])
        gl.nofgwcells=int(fmodset.readline().strip().split(",")[1])
        gl.nofoutlets=int(fmodset.readline().strip().split(",")[1])
        # links between reservoirs are read
        gl.noflinks=[]
        gl.downcell=[]
        for scell in range(gl.nofswcells):
            nlinks=int(fmodset.readline().strip().split(",")[1])
            gl.noflinks.append(nlinks)
            temp=[]
            if nlinks>0:
                for link in range(nlinks):
                    temp.append(int(fmodset.readline().strip().split(",")[1]))
            gl.downcell.append(temp)
        # ouput flag for cells
        gl.outputflag=[]
        for scell in range(gl.nofswcells):
            gl.outputflag.append(int(fmodset.readline().strip().split(",")[1]))

        # final sum flag
        gl.finalsumflag=[]
        for scell in range(gl.nofswcells):
            gl.finalsumflag.append(int(fmodset.readline().strip().split(",")[1]))

        gl.swcellname=[]
        for scell in range(gl.nofswcells):
            cellname=fmodset.readline().strip().split(",")[1]
            #print(scell,cellname)
            gl.swcellname.append(cellname)
    print("done")




#*****************************************************************************
def read_input(file_input):
    print ("reading input data from: "+file_input)
    gl.recdate=[]
    gl.sq_in0=[]
    gl.p1=[]
    gl.p2=[]
    gl.evap=[]
    gl.tminmax=[]
    with open(file_input, "r") as finput:
        fdata=finput.readlines()[1:]  
    test=fdata[0].strip().split(",")
    if len(test)==5:
        #evaporation
        for aline in fdata:
            temp=aline.strip().split(",")
            gl.recdate.append(temp[0])
            gl.sq_in0.append(float(temp[1]))
            gl.p1.append([float(temp[2]),float(temp[3])])
            gl.evap.append(float(temp[4]))
        gl.noftsteps=len(gl.recdate)
    else:
        #temperatures 
        for aline in fdata:
            temp=aline.strip().split(",")
            gl.recdate.append(temp[0])
            gl.sq_in0.append(float(temp[1]))
            gl.p1.append([float(temp[2]),float(temp[3])])
            gl.tminmax.append([float(temp[4]),float(temp[5])])
        gl.noftsteps=len(gl.recdate)
        evap_calc()
    print (str(gl.noftsteps) + " time steps read")


#*****************************************************************************
def read_init(file_init):
    print ("reading initial condition from: "+file_init)
    with open(file_init, "r") as finit:

        # initial storage of surface cells
        gl.sv_start=[]
        for scell in range(gl.nofswcells):
            gl.sv_start.append(float(finit.readline().strip().split(",")[1]))

        #initial storage of groundwater cells
        gl.fv_start=[]
        gl.iv_start=[]
        for scell in range(gl.nofswcells):
            temp=finit.readline().strip().split(",")

            gl.fv_start.append([float(temp[1])]*gl.nofgwcells)
            gl.iv_start.append([float(temp[2])]*gl.nofgwcells)

    gl.iv_start=np.array(gl.iv_start)
    gl.fv_start=np.array(gl.fv_start)
    print("done")






#*****************************************************************************
def read_modpar(file_modpar):
    print ("reading model parameters from: "+file_modpar)
    with open(file_modpar, "r") as fmodpar:
        # spatially constant parameters
        gl.fdet=float(fmodpar.readline().strip().split(",")[1])
        gl.idet=float(fmodpar.readline().strip().split(",")[1])
        gl.fpor=float(fmodpar.readline().strip().split(",")[1])
        gl.ipor=float(fmodpar.readline().strip().split(",")[1])
    
        # volume-area parameters
        gl.bpar=[]
        gl.exponent=[]
        for scell in range(gl.nofswcells):
            temp=fmodpar.readline().strip().split(",")
            gl.exponent.append(float(temp[1]))
            gl.bpar.append(float(temp[2]))

        # outlet parameters
        gl.k=[]
        gl.V=[]
        for scell in range(gl.nofswcells):
            temp2=[]
            temp3=[]
            if gl.noflinks[scell] > 0:
                for link in range(gl.noflinks[scell]):
                    temp=fmodpar.readline().strip().split(",")
                    temp2.append(float(temp[1]))
                    temp3.append(float(temp[2]))                
            gl.k.append(temp2)
            gl.V.append(temp3)
            

        # delay parameter for units
        gl.delay=[]
        for scell in range(gl.nofswcells):
            gl.delay.append(int(fmodpar.readline().strip().split(",")[1]))

        # maun/shakawe rainfall ratio parameters
        gl.statratio=[]
        for scell in range(gl.nofswcells):
            gl.statratio.append(float(fmodpar.readline().strip().split(",")[1]))

        # groundwater reservoir areas and "transmissivity"
        gl.fa=[]
        gl.ia=[]
        gl.kgw=[]
        gl.fa_max=[]
        for scell in range(gl.nofswcells):
            temp=fmodpar.readline().strip().split(",")
            gl.fa.append(float(temp[1]))
            gl.ia.append(float(temp[2]))
            gl.kgw.append(float(temp[3]))
            gl.fa_max.append(float(temp[1])*gl.nofgwcells)
    print ("done")





#*****************************************************************************
def evap_calc():
    r0 = [16.35261505, 14.95509782, 12.8087226, 10.86376736, 9.847079426, 10.22676382, 11.84785549, 14.00041471, 15.76601788, 16.82545576, 17.20206337, 17.09344496]
    kc= [0.95, 0.9, 0.8, 0.7, 0.63, 0.6, 0.6, 0.63, 0.7, 0.8, 0.9, 0.95]
    gl.evap=[]
    for ts in range(gl.noftsteps):
        curmonth=datetime.datetime.strptime(gl.recdate[ts], "%b-%Y").month
        temp= 31 * kc[curmonth-1] * 0.0023 * r0[curmonth-1] * (gl.tminmax[ts][1] - gl.tminmax[ts][0]) ** 0.5 * (((gl.tminmax[ts][1] + gl.tminmax[ts][0]) / 2) + 17.8)
        gl.evap.append(temp)
    print ("calculated evap...")
    
    
    




#***************************************************************************
def model_calc():
    print ("running model..")
    # this is main calculation program
    #prepare some lists to store data
    gl.sq_in=np.array([[0]*gl.noftsteps]*gl.nofswcells)
    gl.sq_in[0,:]=gl.sq_in0
    gl.sq_out=np.array([[0]*gl.noftsteps]*gl.nofswcells)
    gl.fin_sa_end=np.array([[0]*gl.noftsteps]*gl.nofswcells)
    gl.fin_sv_end=np.array([[0]*gl.noftsteps]*gl.nofswcells)
    gl.fin_iv_end=np.array([[0]*gl.noftsteps]*gl.nofswcells)
    gl.fin_fv_end=np.array([[0]*gl.noftsteps]*gl.nofswcells)    
    gl.fin_sev_end=np.array([[0]*gl.noftsteps]*gl.nofswcells)
    gl.fin_iev_end=np.array([[0]*gl.noftsteps]*gl.nofswcells)
    gl.fin_fev_end=np.array([[0]*gl.noftsteps]*gl.nofswcells)    
    gl.sv_end=0
    gl.fa_frac_start=np.transpose(np.array([[0]*gl.nofswcells]*gl.nofgwcells))
    gl.fa_frac_avg=[0]*gl.nofgwcells
    gl.fv_finish=[0]*gl.nofgwcells
    gl.iv_finish=np.transpose(np.array([[0]*gl.nofswcells]*gl.nofgwcells))
    gl.fa_frac_finish=np.transpose(np.array([[0]*gl.nofswcells]*gl.nofgwcells))
    gl.sv_model=0    
    gl.intcptv_model=0
    gl.sa_model=0
    gl.siv_model=0
    gl.sev_model=0
    gl.fev_model=0
    gl.iev_model=0
    gl.fq_model=0
    
    #iterate through time steps
#    for ts in range(50):
    for ts in range(gl.noftsteps):
        if ts%50==0:
            print ("ts="+str(ts))
        #iterate through model cells
        for scell in range(gl.nofswcells):
            calc_res_iter(ts, scell)
            #get starting volume for the next time step for this cell
            gl.sv_start[scell] = gl.sv_end
            #get inflows to the downstream cells
            for link in range(gl.noflinks[scell]):
                cellno=gl.downcell[scell][link]
                gl.sq_in[cellno][ts] = gl.sq_in[cellno][ts] + gl.sq[link]
                    
            if gl.finalsumflag[scell]== 1:
                get_compositefluxes()
    print("done")



#***************************************************************************    
def get_compositefluxes():
    a=1
    gl.sv_model = gl.sv_model + gl.sv_end
#    gl.intcptv_model = gl.intcptv_model + gl.intcptv
    gl.sa_model = gl.sa_model + gl.sa_end
    gl.siv_model = gl.siv_model + gl.siv_unit
    gl.sev_model = gl.sev_model + gl.sev
#    gl.fev_model = gl.fev_model + gl.fev_unit
#    gl.iev_model = gl.iev_model + gl.iev_unit
#    gl.fq_model = gl.fq_model + gl.fq_unit


def mergecells(glvar):
    #write areas for each cell
    selcells=glvar[np.array(gl.outputflag)==1,:]
    cellnames=np.array(gl.swcellname)
    tooutput=np.array(gl.outputflag)
    selcellnames=cellnames[tooutput==1]
    merged=[[0],[1,2],[3],[4],[5],[6,7],[8],[9],[10]]
    outputtable=[]
    outputcellnames=[]
    for m in merged:
        outputcellnames=outputcellnames+[selcellnames[m[0]]]
        current=0
        for i in m:
            current=current+selcells[i,:]
        outputtable=outputtable+[current]
    outputtable=np.array(outputtable)
    outputcellnames=np.array(outputcellnames)
    return outputcellnames,outputtable



#***************************************************************************
def write_output_cellinundation(file_output):
    print ("writing surface area output file...")
    cellnames,celldata=mergecells(gl.fin_sa_end)
    with open(file_output, "w") as foutput:
#        foutput.write("Inundated surface area [km2]");
        foutput.write("Date,");
        for scell in range(len(cellnames)):
            foutput.write(cellnames[scell]+",");
        foutput.write("\n");
        for ts in range(gl.noftsteps):
            foutput.write(gl.recdate[ts]+",");
            for scell in range(len(cellnames)):
                foutput.write(str(int(celldata[scell][ts]))+",");
            foutput.write("\n");
    print ("done")

def write_output_totalinundation(file_output):
    print ("writing surface area output file...")
    cellnames,celldata=mergecells(gl.fin_sa_end)
    with open(file_output, "w") as foutput:
#        foutput.write("Inundated surface area [km2]");
        foutput.write("Date,");
        foutput.write("Total inundation,");
        foutput.write("\n");
        for ts in range(gl.noftsteps):
            foutput.write(gl.recdate[ts]+",");
            foutput.write(str(int(np.sum(celldata,0)[ts]))+",");
            foutput.write("\n");
    print ("done")


def write_output_cellq(file_output):
#write discharges for each cell
    print ("writing discharge output file...")
    cellnames,celldata=mergecells(gl.sq_out)
    with open(file_output, "w") as foutput:
        foutput.write("Inundated surface area [km2]");
        foutput.write("Date,");
        for scell in range(len(cellnames)):
            foutput.write(cellnames[scell]+",");
        foutput.write("\n");
        for ts in range(gl.noftsteps):
            foutput.write(gl.recdate[ts]+",");
            for scell in range(len(cellnames)):
                foutput.write(str(int(celldata[scell][ts]))+",");
            foutput.write("\n");
    print ("done")



def write_output_totalecoregions(file_output):
    print ("writing ecoregions output file...")
    ecoclasses=["Aquatic","Sedges","Grassland","Savanna"]
    with open(file_output, "w") as foutput:
        foutput.write("Year,");
        for unit in range(4):
            foutput.write(ecoclasses[unit]+",");
        foutput.write("\n");
        nts=len(gl.ecototal[0])
        for ts in range(0, nts):
            foutput.write(str(gl.firstyear+ts)+",")
            for unit in range(4):
                foutput.write(str(int(gl.ecototal[unit][ts]))+",");
            foutput.write("\n");
    print ("output has "+str(nts)+" annual time steps")
    print("done\n")






#***************************************************************************    

def calc_res_iter(ts, scell):
    # this calculates one reservoir
    # prepare data
    #reset outlets outflow
    gl.sq=[0]*gl.nofoutlets

    # read inputs and declare initial STATIC variables
    sqin = gl.sq_in[scell][max(0, ts - gl.delay[scell])]
    gl.sv_beg = gl.sv_start[scell]
    gl.sv_end = gl.sv_beg + 0.5 * sqin
    s_iter_flag = 0
    sitern = 0
    sv_endmin = 0
    sv_endmax = ""
    
    s_iter_flag=0
    while s_iter_flag<1:
        sqout = 0
        sv_av = (gl.sv_end + gl.sv_beg) / 2
           
        if scell==0:
            #Case 1
            if gl.sv_end < 1500:
                gl.sa_end = (0.006 * gl.sv_end)**3
            else:
                gl.sa_end = (gl.bpar[scell] * gl.sv_end)**gl.exponent[scell]
        
            if gl.sv_beg < 1500:
                gl.sa_beg = (0.006 * gl.sv_beg)**3
            else:
                gl.sa_beg = (gl.bpar[scell] * gl.sv_beg)**gl.exponent[scell]
        
        elif scell==1 or scell==6 or scell==5 or scell==13 or scell==14:
            #Case 2, 7, 6, 14, 15
            gl.sa_end = (gl.bpar[scell] * gl.sv_end)**gl.exponent[scell]
            if gl.sa_end > gl.fa_max[scell]:
                gl.sa_end = gl.fa_max[scell]
        
            gl.sa_beg = (gl.bpar[scell] * gl.sv_beg)**gl.exponent[scell]
            if gl.sa_beg > gl.fa_max[scell]:
                gl.sa_beg = gl.fa_max[scell]
        
        else:
            gl.sa_end = (gl.bpar[scell] * gl.sv_end)**gl.exponent[scell]
            gl.sa_beg = (gl.bpar[scell] * gl.sv_beg)**gl.exponent[scell]

        
        sa_av = (gl.sa_beg + gl.sa_end) / 2
        rain = (gl.statratio[scell] * gl.p1[max(ts - gl.delay[scell],0)][0] + (1 - gl.statratio[scell]) * gl.p1[max(0,ts - gl.delay[scell])][1])
    
        spv = rain * sa_av / 1000
    
    
        if ts - gl.delay[scell] == 0:
            # this is to account for ts-delay evap from sw
            gl.sev = 0
        else:
            gl.sev = gl.evap[ts - gl.delay[scell]] * sa_av / 1000
            #gl.sev = kc(Month(inp(ts - delay(scell)).recdate)) * inp(ts - delay(scell)).evap * sa_av / 1000



        #-------------------------------------------------------------------------
        # calculate surface outflows
        for n in range(gl.noflinks[scell]):
            if sv_av > gl.V[scell][n]:
                gl.sq[n] = gl.k[scell][n] * (sv_av - gl.V[scell][n])
            else:
                gl.sq[n] = 0
                
        #calculate groundwater outflow
        calc_gw_iter(scell, ts)

        #calculate surface water outflows
        for n in range(gl.nofoutlets):
            sqout = sqout + gl.sq[n]
        #calculate end water balance
        sv_endc = gl.sv_beg + spv - gl.sev - sqout + sqin - gl.siv_unit
#        if(scell==4 and ts%10==0):
    
        #-------------------------------------------------------------------------
        #check if convergence is achieved
        if abs(sv_endc - gl.sv_end) < gl.convcrit * gl.sv_end:
            s_iter_flag = 1
        elif gl.sv_end < 0.001:
            gl.sv_end = 0
            s_iter_flag = 1
        else:
            if gl.sv_end > sv_endc:
                sv_endmax = gl.sv_end
            else:
                sv_endmin = gl.sv_end
        
            if sv_endmax=="":
                gl.sv_end = sv_endc
            else:
                gl.sv_end = sv_endmin + 0.5 * (sv_endmax - sv_endmin)
        
        #-------------------------------------------------------------------------
        #advance iteration
    
        sitern = sitern + 1
        
        if sitern > gl.maxiter:
            print (str(gl.maxiter) + " s iterations in cell " + str(scell) + " in " + gl.recdate[ts])
            s_iter_flag=1

    #preparation for next time step
    gl.sq_out[scell][ts] = sqout

    gl.fin_sa_end[scell][ts] = gl.sa_end
    gl.fin_sv_end[scell][ts] = gl.sv_end
    gl.fin_iv_end[scell][ts] = np.sum(gl.iv_finish[scell])
    gl.fin_fv_end[scell][ts] = np.sum(gl.fv_finish)
    gl.fin_sev_end[scell][ts] = gl.sev

    for gwcell in range(gl.nofgwcells):
        gl.fv_start[scell][gwcell] = gl.fv_finish[gwcell]
        gl.iv_start[scell][gwcell] = gl.iv_finish[scell][gwcell]
        gl.fa_frac_start[scell][gwcell] = gl.fa_frac_finish[scell][gwcell]




#*****************************************************************************
def calc_gw_iter(scell, ts):
    # this calculates a coupled groundwater reservoir
    # prepare data
    #-----------------------------------------------------------------------------
    #get initial values of variables
    fa_cum = 0
    gl.siv = 0
    gl.siv_unit = 0
    fev_sum = 0
    iev_sum = 0
    fv_sum = 0
    iv_sum = 0
    fq_sum = 0

    fpv_sum = 0
    ipv_sum = 0
    gl.fpv = 0
    ipv = 0

    siv_frontsum = 0
    fa_incrfrac = 0

    #-----------------------------------------------------------------------------
    #check if the flood is increasing in the sw cell
    if gl.sv_end > gl.sv_beg:
        sv_isincr = 1
    else:
        sv_isincr = 0
    #-----------------------------------------------------------------------------
    # calculate status of gwcells flooding
    fa_frac=[1]*gl.nofgwcells
    
    for gwcell in range(gl.nofgwcells):
        fa_cumprev = fa_cum
        fa_cum = fa_cum + gl.fa[scell]
        
        if fa_cum < gl.sa_end or fa_cum==gl.sa_end:
            fa_frac[gwcell] = 1
        elif fa_cum > gl.sa_end and fa_cumprev < gl.sa_end:
            fa_frac[gwcell] = (gl.sa_end - (fa_cumprev)) / gl.fa[scell]
            if gl.sa_beg > fa_cumprev and gl.sa_end > gl.sa_beg:
                 fa_incrfrac = (gl.sa_end - gl.sa_beg) / (fa_cum - gl.sa_beg)
            else:
                fa_incrfrac = fa_frac[gwcell]
        else:
            fa_frac[gwcell] = 0
        gl.fa_frac_avg[gwcell] = (fa_frac[gwcell] + gl.fa_frac_start[scell, gwcell]) / 2

    #*****************************************************************************
    #calculate groundwater flow between floodplains and islands
    for gwcell in range(gl.nofgwcells):
       #-----------------------------------------------------------------------------
       # if the gw cell is not flooded
        if fa_frac[gwcell]== 0:
            gl.fv_beg = gl.fv_start[scell][gwcell]
            siv_front = 0
            gl.f_isfloodflag = 0

            iteration_1(scell,gwcell,ts)

            siv_adv = 0 #fa_frac(gwcell) * gl.fq
            # fv_end = fv_end + siv_adv
        #-----------------------------------------------------------------------------
        # if the gw cell is partly flooded
        elif fa_frac[gwcell] > 0 and fa_frac[gwcell] < 1:
            if sv_isincr==1:
               siv_front = ((gl.idet * gl.fpor * gl.fa[scell]) - gl.fv_start[scell][gwcell]) * fa_incrfrac
            else:
                siv_front = 0
        
            gl.fv_beg = gl.fv_start[scell][gwcell] + siv_front
            gl.f_isfloodflag = 0
        
            iteration_1(scell, gwcell,ts)
        
            siv_adv = fa_frac[gwcell] * gl.fq
            #gl.fv_end = gl.fv_end + siv_adv    
        #-----------------------------------------------------------------------------
        # if the gw cell is entirely flooded
        else:
            gl.fv_max = gl.idet * gl.fpor * gl.fa[scell]    
            siv_front = gl.fv_max - gl.fv_start[scell][gwcell]
            gl.fv_beg = gl.fv_max
            gl.fv_end = gl.fv_max
            gl.f_isfloodflag = 1
            gl.fev = 0
        
            gl.fpv = ((gl.statratio[scell] * gl.p1[ts][0] + (1 - gl.statratio[scell]) * gl.p1[ts][1]) / 1000 * gl.fa[scell]) * (1 - gl.fa_frac_avg[gwcell])
   
            iteration_2(scell, gwcell, ts)
        
            siv_adv = gl.fq
            
        #-----------------------------------------------------------------------------
        #sets initial values for the next time step
    
        gl.fv_finish[gwcell] = gl.fv_end
        gl.iv_finish[scell][gwcell] = gl.iv_end
        gl.fa_frac_finish[scell][gwcell] = fa_frac[gwcell]
    
    
        #-----------------------------------------------------------------------------
        #calculates composite fluxes
#        siv_advsum = siv_advsum + siv_adv
        siv_frontsum = siv_frontsum + siv_front             # sw cell infiltration
        gl.siv = siv_front + siv_adv
        gl.siv_unit = gl.siv_unit + gl.siv              #sw cell infiltration
        fev_sum = fev_sum + gl.fev
        iev_sum = iev_sum + gl.iev
        fv_sum = fv_sum + gl.fv_end
        iv_sum = iv_sum + gl.iv_end
        fq_sum = fq_sum + gl.fq
        fpv_sum = fpv_sum + gl.fpv
        ipv_sum = ipv_sum + ipv
    gl.fin_iev_end[scell][ts]=iev_sum
    gl.fin_fev_end[scell][ts]=fev_sum






def iteration_1(scell, gwcell, ts):
    # this calculates the floodplain groundwater reservoir
    #**********************************************************************************
    #prepare data
    gl.fv_end = gl.fv_beg       #initial guess of end floodplain groundwater volume
    iv_beg = gl.iv_start[scell, gwcell]
    gl.iv_end = iv_beg       #initial guess of end island groundwater volume
    fiter_endflag = 0
    fitern = 0
    fv_endmin = 0
    fv_endmax = ""

    #**********************************************************************************
    #outer iteration start
    while fiter_endflag==0:
        gl.fv_av = (gl.fv_end + gl.fv_beg) / 2
        if gl.fv_av > ((gl.idet - gl.fdet) * gl.fa[scell] * gl.fpor):
            temp = ((gl.fv_av - ((gl.idet - gl.fdet) * gl.fa[scell] * gl.fpor)) / (gl.fdet * gl.fa[scell] * gl.fpor)) * (1 - gl.fa_frac_avg[gwcell])
            if temp > 1:
                temp = 1
            gl.fev = (gl.evap[ts] / 1000 * gl.fa[scell] * temp)
            gl.fpv = ((gl.statratio[scell] * gl.p1[ts][0] + (1 - gl.statratio[scell]) * gl.p1[ts][1]) / 1000 * gl.fa[scell]) * (1 - gl.fa_frac_avg[gwcell])
        else:
            gl.fev = 0
    
        gl.f_isfloodflag = 0
        iteration_2(scell,gwcell, ts)
        fv_endc = gl.fv_beg + gl.fpv - gl.fq - gl.fev
        #-----------------------------------------------------------------------------
        #check if outer iteration convergence is achieved
        if abs(fv_endc - gl.fv_end) < (gl.convcrit * gl.fv_end):
            fiter_endflag = 1
        elif gl.fv_end < 0.001:
            fiter_endflag = 1
        else:
            if gl.fv_end > fv_endc:
                fv_endmax = gl.fv_end
            else:
                fv_endmin = gl.fv_end

            #----------------------------------------------------
            #KKKKKKKKKKKKKKKKKKKKK
            if fv_endmax=="":
                gl.fv_end = fv_endc
            else:
                gl.fv_end = fv_endmin + 0.5 * (fv_endmax - fv_endmin)
    #-----------------------------------------------------------------------------
    #advance the iteration
        fitern = fitern + 1
        if fitern > gl.maxiter:
            print (str(gl.maxiter) + " f iterations in scell " + str(scell) + ", gwcell"+str(gwcell)+" in " + gl.recdate[ts])
            fiter_endflag=1





def iteration_2(scell,gwcell,ts):
    # this calculates the island groundwater reservoir
    #*********************************************************************************
    #prepare data
    i_itern = 0
    i_iter_endflag = 0
    iv_endmin = 0
    iv_endmax = ""
    #-----------------------------------------------------------------------------
    #initial guess of end island volume
    iv_beg = gl.iv_start[scell, gwcell]
    gl.iv_end = iv_beg
    #**********************************************************************************
    #inner iteration start

    while i_iter_endflag==0: 
        iv_av = (gl.iv_end + iv_beg) / 2
        if gl.f_isfloodflag==1:
            gl.fq = ((gl.fv_max / (gl.fa[scell] * gl.fpor)) - (iv_av / ((gl.ia[scell]) * gl.ipor))) * gl.kgw[scell]
        else:
            gl.fq = ((gl.fv_av / (gl.fa[scell] * gl.fpor)) - (iv_av / ((gl.ia[scell]) * gl.ipor))) * gl.kgw[scell]            

        evapi = iv_av / (gl.ipor * gl.idet * gl.ia[scell])    
        if evapi > 0.6:
            evapi = 0.6
        gl.iev = gl.evap[ts] * (gl.ia[scell]) / 1000 * evapi
        ipv = ((gl.statratio[scell] * gl.p1[ts][0] + (1 - gl.statratio[scell]) * gl.p1[ts][1]) / 1000 * (gl.ia[scell]))
        iv_endc = iv_beg + gl.fq - gl.iev + ipv
    
        #-----------------------------------------------------------------------------
        #check if inner iteration convergence is achieved
        if abs(iv_endc - gl.iv_end) < gl.convcrit * gl.iv_end:
            i_iter_endflag = 1
        elif gl.iv_end < gl.convcrit:
            i_iter_endflag = 1
        else:
            if gl.iv_end > iv_endc:
                iv_endmax = gl.iv_end
            else:
                iv_endmin = gl.iv_end

            if iv_endmax=="":
                gl.iv_end = iv_endc
            else:
                gl.iv_end = iv_endmin + 0.5 * (iv_endmax - iv_endmin)

        #-----------------------------------------------------------------------------
        #advance the iteration
        i_itern = i_itern + 1

def eco_calc():
    print("calculating eco model... ")
    size=[]
    nofyears=int(np.floor(gl.noftsteps)/12)
    nmonths=nofyears*12
    cellnames,celldata=mergecells(gl.fin_sa_end)

    dates=pd.to_datetime(gl.recdate)
    firstmonth=dates[0].month
    print(firstmonth)
    first2read=(12-firstmonth+1)%12
    gl.firstyear=dates[first2read].year
    print("input file has "+str(len(dates))+" monthly time and "+str(celldata.shape[1])+" units")
    print("January at timestep "+str(first2read))

    floodsize=np.sum(celldata[(first2read):],0)
    temp=np.copy(floodsize[0:nmonths])
    temp=temp.reshape(nofyears, 12)
    areayear=np.mean(temp,1)
    sizerange=range(1,12000, 1)
    ecoall=[]
    for size in sizerange:
        d=np.sum(temp>size, 1)
        eco=[2]*nofyears
        for y in range(nofyears):
            if eco[y-1]==1: #A
                # rules for Aquatic
                if d[y] == 0:
                    eco[y] = 2 #RS
                else:
                    if d[y - 1] == 12:
                        eco[y] = 1 #"A"
                    elif d[y - 2] == 12 and d[y - 3] == 12 and d[y - 4] == 12:
                        eco[y] == 1 #"A"
                    else:
                        eco[y] = 2 #"RS"
            elif eco[y-1]==2: #"RS"
            # rules for Sedges
                if d[y]== 0:
                    eco[y] = 3 #"G"
                elif d[y - 1] < 12:
                    eco[y] = 2  #"RS"
                else:
                    eco[y]=1 #A
    #            elif d[y - 1] ==12:
    #                if d[y - 2] ==12:
    #                    eco[y] = 1 #"A"
    #                else:
    #                    eco[y] = 2 #"RS"

            elif eco[y-1]==3: #"G"
            # rules for Grassland
                if d[y]== 0:
                    if d[y - 1]== 0:
                        if d[y - 2] > 0 or d[y - 3] > 0 or d[y - 4] > 0:
                            eco[y] = 3 #"G"
                        else:
                            eco[y] = 4 #"S"
                    elif d[y - 1] > 0:
                        eco[y] = 3 #"G"
                elif d[y] > 0:
                    if d[y - 1]==0:
                        eco[y]== 3 #"G"
                    else:
                        if d[y - 2]==0:
                            eco[y] = 3 #"G"
                        else:
                            eco[y] = 2 #"RS"

            elif eco[y-1]==4: #"S"
                # rules for Savanna
                if d[y]== 0:
                    eco[y] = 4 #"S"
                elif d[y] > 1:
                    eco[y] = 3 #"G"
                else:
                    eco[y] = 4 #"S"
        ecoall=ecoall+[eco]
    ecoall=np.array(ecoall)
    gl.ecototal=[[]]*5
    for j in range(0,4):
        gl.ecototal[j]=np.sum(ecoall==j+1, 0)
    gl.ecototal=np.array(gl.ecototal)
    print("done\n")



def inund_calc(outputfile):
    print("Animating inundation...\n")
    import struct
    import numpy as np
    import PIL
    import sys
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt
    import scipy.stats as st
    import matplotlib.animation as animation

    fps=5

    mapdir="./config/"
    nrow,ncol = 300,303
    nofpix = ncol * nrow

    skip=3
    thresh=0.7

    print("reading inundation map parameters...")

    ncdata=Dataset(mapdir+"./m_arc.nc")
    m=ncdata.variables['Band1'][:]
    m[m<0]=1000

    ncdata=Dataset(mapdir+"./sigma_arc.nc")
    sigma=ncdata.variables['Band1'][:]

    ncdata=Dataset(mapdir+"./units_arc.nc")
    units=ncdata.variables['Band1'][:]
    units=units.astype("float")

    units[m<=0]=np.nan
    units[sigma<=0]=np.nan

    codes={1:"Panhandle", 2:"Thaoge", 3:"Xudum", 4:"Boro", 5:"Khwai", 6:"Nqoga-1a", 7:"Selinda", 8:"Nqoga-2a", 9:"Mboroga"}

    unitsf=units.flatten()
    mf=m.flatten()
    sigmaf=sigma.flatten()

    print("done\n")


    cellnames,cellvalues=mergecells(gl.fin_sa_end)
    dates=pd.to_datetime(gl.recdate)
    
    sa=pd.DataFrame(cellvalues.T,index=dates, columns=cellnames)

    nts=sa.shape[0]

    print("read", nts, "time steps")
    print("will skip", skip, "time steps")

    print("done\n")


    print("processing...")

    fig, pl= plt.subplots(figsize=(5, 5))

    #preparing empty frame
    temp=np.zeros_like(units).astype(float)
    temp[:]=0
    im = pl.imshow(temp, cmap=plt.cm.RdBu)

    #plt.show()

    #preparing canvas
    pl.set_yticks([])
    pl.set_xticks([])
    tx=pl.text(0.7,0.9,sa.index.strftime('%Y %B')[0], transform=pl.transAxes)


    def plotflood(ts): 
        if ts%10==0:
            print ("ts="+str(ts))
        inu=sa.iloc[ts,:]
        amap=np.zeros_like(mf).astype(float)
        amap[:]=np.nan
        for key in codes.keys():
            area=inu[codes[key]]
            sel=np.where(unitsf==key)[0]
            for x in sel:
                prob=st.norm.cdf(area,mf[x],sigmaf[x])
                if thresh>0:
                    if prob>thresh:
                        amap[x]=1
                else:
                        amap[x]=prob
        amap=np.flipud(amap.reshape(nrow,ncol))
        tx.set_text(sa.index.strftime('%Y %B')[ts])
        im.set_array(amap)
        return im,


    ani = animation.FuncAnimation(fig, plotflood, range(skip,nts))
    writer = animation.ImageMagickFileWriter(fps=fps)
    ani.save(outputfile, writer=writer) 
    print ("done")








read_modset(file_modset)                                #reading model configuration
read_input(file_input)                                  #reading inputs
read_init(file_init)                                    #reading initial conditions
read_modpar(file_modpar)                                #reading model parameters

model_calc()                                            #this is when the model is actually run

for outputfile in outputfiles:
    if "allinundation" in outputfile:
        write_output_cellinundation(outputfile)                     #inundation by unit
    if "alloutflows" in outputfile:
        write_output_cellq(outputfile)                       #streamflow/unit discharges
    if "totalinundation" in outputfile:
        write_output_totalinundation(outputfile)                       #total inundation
    if "totalecoregions" in outputfile:
        eco_calc()
        write_output_totalecoregions(outputfile)                       #ecoregions
    if "animatedinundation" in outputfile:
        inund_calc(outputfile)
#        write_output_animatedinundation(outputfile)                       #inundation movie

print("success")
 
