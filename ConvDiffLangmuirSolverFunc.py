import numpy as np
from tridag import *

def ConvDiffLangmuirSolver(numpar, physpar, bc):
    numpar.alfa3 = numpar.dt/(physpar.peclet*(numpar.dx**2))
    
    varspatial  = np.rec.array([(np.zeros([numpar.nx,1])), (np.zeros([numpar.nx,1]))], 
          dtype=[('cmod_old','float64'), ('cmod','float64'), ])
    
    pv = np.arange(0.0, numpar.tstop+numpar.dt, numpar.dt)
    vartemporal = np.rec.array([(pv), (np.zeros([len(pv)])), (np.zeros([len(pv)]))], 
          dtype=[('pv','float64'), ('cout','float64'),('time','float64') ])
                 
    a = np.zeros([numpar.nx,3])
    b = np.zeros([numpar.nx,1])
    alfa4 = np.zeros([numpar.nx,1])

    for j in range(0,numpar.nt):         
        alfa4 = 1 + physpar.k/(1+physpar.beta*varspatial.cmod_old)**2
        a[:,0]=-(numpar.alfa1+numpar.alfa3)
        a[1:numpar.nx-1,1]=alfa4[1:numpar.nx-1,0]+numpar.alfa1+2*numpar.alfa3
        a[0,1]=alfa4[0,0]+numpar.alfa1+3*numpar.alfa3
        a[numpar.nx-1,1]=alfa4[numpar.nx-1,0]+numpar.alfa1+numpar.alfa3
        a[:,2]=-numpar.alfa3
        varspatial.cmod_old[0,0]=varspatial.cmod_old[0,0]*alfa4[0,0]+\
                 (numpar.alfa1+2*numpar.alfa3)*bc.bcL
        varspatial.cmod=tridag(a,varspatial.cmod_old,numpar.nx);
        varspatial.cmod_old=varspatial.cmod
        vartemporal.cout[j+1]=varspatial.cmod[numpar.nx-1,0]
        vartemporal.time[j+1]=vartemporal.time[j]+numpar.dt

    return varspatial, vartemporal
