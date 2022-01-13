#!/usr/bin/env python
"""
Computing Energy and Force Using Models Inside Model Zoo
========================================================
TorchANI has a model zoo trained by NeuroChem. These models are shipped with
TorchANI and can be used directly.
"""

###############################################################################
# To begin with, let's first import the modules we will use:
import torch
import torchani
import numpy as np
import sys
import getopt
import time
from torch_search import LineSearch

###############################
def rotation_BCD_xy(coord,i,j,k):
   
    bc = coord[0,i] - coord[0,j]
    dc = coord[0,k] - coord[0,j]
   
    t_n1 = torch.cross(bc,dc)
    Vac=t_n1/t_n1.norm(p=2)
    a=Vac[0].item()
    b=Vac[1].item()
    c=Vac[2].item()

    
    Rby = torch.tensor([[1,0,0],
                           [0,1,0],
                           [0,0,1]],dtype=torch.float32)
    if a**2+c**2>0.0000000000000001:
       Rby = torch.tensor([[c/np.sqrt(a**2+c**2),0,-a/np.sqrt(a**2+c**2)],
                        [0,1,0],
                        [a/np.sqrt(a**2+c**2),0,c/np.sqrt(a**2+c**2)]],dtype=torch.float32,device=device) 

    Rbx = torch.tensor([[1,0,0],
                        [0,np.sqrt(a**2+c**2),-b],
                        [0,b,np.sqrt(a**2+c**2)],
                        ],dtype=torch.float32, device=device)  
 
    Rbcd=torch.mm(Rbx,Rby)
    
    return Rbcd
#################################doh define constrain###########################
def translation_0(coord,k0):
    
    trans=torch.ones(coord.size(), device=device)
    trans[0,:,0]=coord[0,k0,0]*trans[0,:,0]
    trans[0,:,1]=coord[0,k0,1]*trans[0,:,1]
    trans[0,:,2]=coord[0,k0,2]*trans[0,:,2]    
    return coord-trans
def rotation_x(coord,j):
    bc=coord[0,j]
    Vbc=bc/bc.norm(p=2)
    a=Vbc[0].item()
    b=Vbc[1].item()
    c=Vbc[2].item()
    

    
    Rby = torch.tensor([[1,0,0],
                           [0,1,0],
                           [0,0,1]],dtype=torch.float32, device=device)
    if a**2+c**2>0.0000000000000001:
       Rby = torch.tensor([[a/np.sqrt(a**2+c**2),0,c/np.sqrt(a**2+c**2)],
                        [0,1,0],
                        [-c/np.sqrt(a**2+c**2),0,a/np.sqrt(a**2+c**2)]],dtype=torch.float32, device=device) 

    Rbz = torch.tensor([[np.sqrt(a**2+c**2),b,0],
                        [-b,np.sqrt(a**2+c**2),0],
                        [0,0,1]],dtype=torch.float32, device=device)  
 
    Rb=torch.mm(Rbz,Rby)
    
    return Rb
def rotation_z(coord,j):
    bc=coord[0,j]
    Vbc=bc/bc.norm(p=2)
    a=Vbc[0].item()
    b=Vbc[1].item()
    c=Vbc[2].item()
    
    
    Rbz = torch.tensor([[a,b,0],
                        [-b,a,0],
                        [0,0,1]],dtype=torch.float32, device=device)     
    return Rbz
def rotation_xy(Rb,coord,i):
    coord=torch.mm(Rb,coord.squeeze().T)
    ac=coord.T[i]
    Vac=ac/ac.norm(p=2)
    a=Vac[0].item()
    b=Vac[1].item()
    c=Vac[2].item()

    
    Rax = torch.tensor([[1,0,0],
                        [0,1,0],
                    [0,0,1]],dtype=torch.float32, device=device)  
    if b**2+c**2>0.0000000000000001:
              Rax = torch.tensor([[1,0,0],
                        [0,b/np.sqrt(b**2+c**2),c/np.sqrt(b**2+c**2)],
                    [0,-c/np.sqrt(b**2+c**2),b/np.sqrt(b**2+c**2)]],dtype=torch.float32, device=device)
    
    return Rax


def torision_force_bcd(Rbcd,Ra,Rb,Rd,force,i,j,k,l):  
    Rbcd_f=torch.mm(Rbcd,force.squeeze().T)
    Rb_f=torch.mm(Rb,Rbcd_f)
    Ra_f=torch.mm(Ra,Rb_f)
    Rd_f=torch.mm(Rd,Rb_f)
    Ra_f.T[i,2]=0.0
    Rb_f.T[j,1:3]=0.0
    Rb_f.T[k,1:3]=0.0
    Rd_f.T[l,2]=0.0
    Ra_f=torch.mm(Ra.inverse(),Ra_f)
    Ra_f=torch.mm(Rb.inverse(),Ra_f)
    Ra_f=torch.mm(Rbcd.inverse(),Ra_f)
    Rb_f=torch.mm(Rb.inverse(),Rb_f)
    Rb_f=torch.mm(Rbcd.inverse(),Rb_f)
    Rd_f=torch.mm(Rd.inverse(),Rd_f)
    Rd_f=torch.mm(Rb.inverse(),Rd_f)
    Rd_f=torch.mm(Rbcd.inverse(),Rd_f)
    force[0,i]=Ra_f.T[i]
    force[0,j]=Rb_f.T[j]
    force[0,k]=Rb_f.T[k]
    force[0,l]=Rd_f.T[l]
    return force
  

def angle_force(Ra,Rc,force,i,j,k):  
    Rc_f=torch.mm(Rc,force.squeeze().T)
    Ra_f=torch.mm(Ra,Rc_f)
    Ra_f.T[i,1:3]=0.0
    Rc_f.T[j,0:3]=0.0
    Rc_f.T[k,1:3]=0.0
    Ra_f=torch.mm(Ra.inverse(),Ra_f)
    Ra_f=torch.mm(Rc.inverse(),Ra_f)
    Rc_f=torch.mm(Rc.inverse(),Rc_f)
    force[0,i]=Ra_f.T[i]
    force[0,j]=Rc_f.T[j]
    force[0,k]=Rc_f.T[k]
    return force

def distance_force(Rb,force,i,j):  
    Rb_f=torch.mm(Rb,force.squeeze().T)
    Rb_f.T[i,0:3]=0.0
    Rb_f.T[j,0:3]=0.0
    Rb_f=torch.mm(Rb.inverse(),Rb_f)
    force[0,i]=Rb_f.T[i]
    force[0,j]=Rb_f.T[j]
    return force
def t_ang1(coord,i,j,k,l):
    a = coord[0,j] - coord[0,i]
    b = coord[0,k] - coord[0,j]
    c = coord[0,l] - coord[0,k]
    t_n1 = torch.cross(a,b)
    t_n2 = torch.cross(c,b)
    return np.degrees(torch.acos(torch.dot(t_n1,t_n2)/(t_n1.norm(p=2)*t_n2.norm(p=2))).item())
 
    #Rb_coord=torch.mm(Rb_mat,trc_coord.squeeze().T)
def t_ang0(coord,i,j,k,l):
    f = coord[0,i] - coord[0,j]
    g = coord[0,j] - coord[0,k]
    h = coord[0,l] - coord[0,k]
    
    t_n1 = torch.cross(f,g)
    t_n2 = torch.cross(h,g)
    n_cro = torch.cross(t_n1, t_n2)
    cosn = torch.dot(t_n1,t_n2)
    sinn = torch.dot(n_cro, g)/ g.norm(p=2)
    r = -np.arctan2(sinn.item(), cosn.item())
    return np.degrees(r)  
 ##################   Ra_mat=rotation_xy(Rb_mat,trc_coord,t_tor[i]-1)
 
def b_ang0(coord,i,j,k):  
    s1 = coord[0,i]-coord[0,j]
    s2 = coord[0,k]-coord[0,j]
    return np.degrees(torch.acos(torch.dot(s1,s2)/(s1.norm(p=2)*s2.norm(p=2))).item())
def d_con0(coord,i,j):
    return torch.norm(coord[0,j] - coord[0,i]).item()
######extract input################## 
def extract(text,target):
   linenum = 0
   found = False
   for line in text:
      if not found:
        if (line.find(target)) > -1:
          found=True
        else:
          linenum += 1
   return linenum 
def force_projection(t_bond, t_ang,  t_tor, force, coordinates):
        
        if t_bond.size>=2:
            flag=1
            for i in range(0,int(t_bond.size),2):
                force[0,t_bond[i]]=0.0
                force[0,t_bond[i+1]]=0.0
 ### #########################bond angle constrain############################

        if t_ang.size>=3:
            flag=1
            for i in range(0,int(t_ang.size),3):
                trc_coord=translation_0(coordinates,t_ang[i+1])
                Rc_mat=rotation_x(trc_coord,t_ang[2])
                coord_rc=torch.mm(Rc_mat,trc_coord.squeeze().T)
                coord_rc=coord_rc.T
                coord_rc=coord_rc.reshape([1,coord_rc.size()[0],coord_rc.size()[1]])    
                Ra_mat=rotation_x(coord_rc,t_ang[i])
                force=angle_force(Ra_mat,Rc_mat,force,t_ang[i],t_ang[i+1],t_ang[i+2])
    
    ############################torision constrain###########################



        if t_tor.size==4:
            flag=1
            
            coordinates0=coordinates*1
            for i in range(0,int(t_tor.size),4):
                if  i<=int(t_tor.size)-4:
                    Rbcd=rotation_BCD_xy(coordinates,t_tor[i+1]-1,t_tor[i+2]-1,t_tor[i+3]-1)
                    BCD_coordinates=torch.mm(Rbcd, coordinates.squeeze().T)
                    coordinates0[0]=BCD_coordinates.T 
                      
                trc_coord=translation_0(coordinates0,t_tor[i+2]-1)
 #              print('trc_coord:',trc_coord,t_tor[i+2]-1)
                Rb_mat=rotation_x(trc_coord,t_tor[i+1]-1)
 #              print('Rb_mat:',Rb_mat)
                Rb_coord=torch.mm(Rb_mat,trc_coord.squeeze().T)
 #              print('Rb_coord:',Rb_coord.T)
                Ra_mat=rotation_xy(Rb_mat,trc_coord,t_tor[i]-1)
                Ra_coord=torch.mm(Ra_mat,trc_coord.squeeze().T)
                #print('Ra_coord:',Ra_coord.T)
                Rd_mat=rotation_xy(Rb_mat,trc_coord,t_tor[i+3]-1)
                Rd_coord=torch.mm(Rd_mat,trc_coord.squeeze().T)
                #print('Rd_coord:',Rd_coord.T)
                force=torision_force_bcd(Rbcd,Ra_mat,Rb_mat,Rd_mat,force,t_tor[i]-1,t_tor[i+1]-1,t_tor[i+2]-1,t_tor[i+3]-1)
                
        if t_tor.size>4:
            flag=1
                          
            coordinates0=coordinates*1 
            for i in range(0,int(t_tor.size),4):
                if  i<int(t_tor.size)-4:
                    Rbcd=rotation_BCD_xy(coordinates,t_tor[i+1]-1,t_tor[i+2]-1,t_tor[i+3]-1)
                BCD_coordinates=torch.mm(Rbcd, coordinates.squeeze().T)
                coordinates0[0]=BCD_coordinates.T 
                      
                trc_coord=translation_0(coordinates0,t_tor[i+2]-1)
 #              print('trc_coord:',trc_coord,t_tor[i+2]-1)
                Rb_mat=rotation_x(trc_coord,t_tor[i+1]-1)
 #              print('Rb_mat:',Rb_mat)
                Rb_coord=torch.mm(Rb_mat,trc_coord.squeeze().T)
 #              print('Rb_coord:',Rb_coord.T)
                Ra_mat=rotation_xy(Rb_mat,trc_coord,t_tor[i]-1)
                Ra_coord=torch.mm(Ra_mat,trc_coord.squeeze().T)
                #print('Ra_coord:',Ra_coord.T)
                Rd_mat=rotation_xy(Rb_mat,trc_coord,t_tor[i+3]-1)
                Rd_coord=torch.mm(Rd_mat,trc_coord.squeeze().T)
                #print('Rd_coord:',Rd_coord.T)
                force=torision_force_bcd(Rbcd,Ra_mat,Rb_mat,Rd_mat,force,t_tor[i]-1,t_tor[i+1]-1,t_tor[i+2]-1,t_tor[i+3]-1)
                  
                force[0,t_tor[i+2]-1]= 0             
                  
        return  force 
def convergence_fun(max_f,rms_f,max_dp,rms_dp,f_max_th,f_rms_th,dp_max_th,dp_rms_th):
        f_max_flag='NO'
        f_rms_flag='NO'
        dp_max_flag='NO'
        dp_rms_flag='NO'                  
        if max_f < f_max_th:
                   	f_max_flag='YES' 
        else:
                    	f_max_flag='NO'
        if rms_f < f_rms_th:
                  	f_rms_flag='YES'
        else:
                  	f_rms_flag='NO'
        if max_dp < dp_max_th:
                	dp_max_flag='YES'
        else:
                 	dp_max_flag='NO'
        if rms_dp < dp_rms_th:
                	dp_rms_flag='YES'
        else:
                 	dp_rms_flag='NO'
        return f_max_flag, f_rms_flag,  dp_max_flag, dp_rms_flag
def BFGS(model, species, t_elem,  coord, t_bond,t_ang, t_tor, f_max_th, f_rms_th, dp_max_th, dp_rms_th, use_line_search=False,memory=100,  curvature=70.0,  maxstep=0.2, maxiteration=10000):
    # default parameters
    

        """BFGS optimizer.

        Parameters:

        atoms: Atoms object
            The Atoms object to relax.

        restart: string
            Pickle file used to store hessian matrix. If set, file with
            such a name will be searched and hessian matrix stored will
            be used, if the file exists.

        trajectory: string
            Pickle file used to store trajectory of atomic movement.

        logfile: file object or str
            If *logfile* is a string, a file with that name will be opened.
            Use '-' for stdout.

        maxstep: float
            Used to set the maximum distance an atom can move per
            iteration (default value is 0.2 Å).

        master: boolean
            Defaults to None, which causes only rank 0 to save files.  If
            set to true,  this rank will save files.

        alpha: float
            Initial guess for the Hessian (curvature of energy surface). A
            conservative value of 70.0 is the default, but number of needed
            steps to converge might be less if a lower value is used. However,
            a lower value also means risk of instability.
        """

        # initial hessian
        
        H0 = torch.eye(3 * (species.numpy().size))*curvature

        H = None       
       
        alph=1.0
       
        p = None #search direction of LBFGS
 
        iteration = 0
# Store r-r0, g-g0, and 1/(r-r0)(g-g0) into s[], y[], rho[] list respectively. 

        r = coord*1.0
        r0 = coord*1.0
        dr0 =r-r0
        e = model((species,  coord)).energies
        e0 = e*1.0
        g = torch.autograd.grad(e.sum(),  coord)[0]
        f = -1.0*g
        f = force_projection(t_bond,t_ang, t_tor, f, r)
        f0= f*1.0
        f0 = f0.reshape(-1)
        log=1
        coord_size = species.numpy().size*3.0
###########convergence flag##################
        f_max_flag='NO'
        f_rms_flag='NO'
        dp_max_flag='NO'
        dp_rms_flag='NO'
        con_flag = 'NO'
       
######Initial contrain bond, angle, torsion angle########
 

        while iteration<maxiteration:
        
                f = f.reshape(-1)
                f0 = f0.reshape(-1)
                r1 = r.reshape(-1) 
                r0 = r0.reshape(-1) 
                if H is None:
                     H = H0
                if iteration>0:  
                     dr = r1 - r0          
                     df = f - f0
                     a = dr.dot(df)
                     dg = torch.matmul(H, dr)
                     b = dr.dot(dg)
                     H -= torch.ger(df, df) / a + torch.ger(dg, dg) / b
                omega, V = torch.symeig(H, eigenvectors=True)
                #print(torch.matmul(f, V))
                dr = (torch.matmul(V, torch.matmul(f, V)/ torch.abs(omega))).reshape((-1, 3))
                dr = dr0+dr
                dr = force_projection(t_bond,t_ang, t_tor, dr, r)
                longest_step =torch.norm(dr.data.norm(p=2,dim=2),float('inf'))
        
                if longest_step >= maxstep:
                   dr.data *= maxstep/longest_step*1.0
                r0.data = r.data
                f0 = f
                
                r=r0 + dr
                iteration += 1
                f0 = -1.0*g
                e0 = e
               
                
                e = model((species, r)).energies
                #print(time.time()-start_time)
                
                g = torch.autograd.grad(e.sum(), r)[0]
               
               
                #print(time.time()-start_time)
                f = -g*1.0
                f = force_projection(t_bond,t_ang, t_tor, f, r)
###########convergence_FLAG##############################
               
                max_dp=dr.norm(float('inf'))
                rms_dp=torch.sqrt(dr.norm(p=2)**2/coord_size)
                max_f=torch.norm(f.norm(p=2,dim=2),float('inf'))
                rms_f =torch.sqrt(f.norm(p=2)**2/coord_size)
                print(iteration,e.item(), max_f.item(), max_dp.item(), alph)
                
                f_max_flag, f_rms_flag,dp_max_flag,dp_rms_flag = convergence_fun(max_f,rms_f,max_dp,rms_dp,f_max_th,f_rms_th,dp_max_th,dp_rms_th)
                if  f_max_flag=='YES'  and f_rms_flag=='YES' and dp_max_flag=='YES'and dp_rms_flag=='YES':     
                    return  r, iteration,max_f,rms_f,max_dp,rms_dp,f_max_flag,f_rms_flag,dp_max_flag,dp_rms_flag
        return  r, iteration,max_f,rms_f,max_dp,rms_dp,f_max_flag,f_rms_flag,dp_max_flag,dp_rms_flag           
def LBFGS(model, species, t_elem,  coord, t_bond,t_ang, t_tor, f_max_th, f_rms_th, dp_max_th, dp_rms_th, use_line_search=False,memory=100,  curvature=70.0,  maxstep=0.2, maxiteration=10000):
        """Parameters:
        
        species: Atoms object
            The Atoms object to relax.

        t_elem: Element of Atoms

        The convergence criteria:
        Max force         RMS focre        Max displacement  Rms Displacement
        f_max_th=0.00045, f_rms_th=0.0003, dp_max_th=0.0018, dp_rms_th=0.0012
       
        coord float
        Intial coordinates of atoms.

        logfile: file object or str
            If *logfile* is a string, a file with that name will be opened.
            Use '-' for stdout.

      

        maxstep: float
            How far is a single atom allowed to move. This is useful for DFT
            calculations where wavefunctions can be reused if steps are small.
            Default is 0.2 Angstrom.

        memory: int
            Number of steps to be stored. Default value is 100. Three numpy
            arrays of this length containing floats are stored.


        curvature: float
            Initial guess for the Hessian (curvature of energy surface). A
            conservative value of 70.0 is the default, but number of needed
            steps to converge might be less if a lower value is used. However,
            a lower value also means risk of instability.


        """
       

        

        if maxstep > 1.0:
            raise ValueError('You are using a much too large value for ' +
                             'the maximum step size: %.1f Angstrom' %
                             maxstep)

        
        # Initial approximation of inverse Hessian 1./70. is to emulate the behaviour of BFGS. 
        H0 = 1. / curvature
        
        alph=1.0
       
        p = None #search direction of LBFGS
 
        iteration = 0
# Store r-r0, g-g0, and 1/(r-r0)(g-g0) into s[], y[], rho[] list respectively. 
        s = []
        y = []       
        rho = []
        
        r = coord*1.0
        r0 = coord*1.0
        dr =r-r0
        e = model((species,  coord)).energies
        e0 = e*1.0
        g = torch.autograd.grad(e.sum(),  coord)[0]
        f = -1.0*g
        f = force_projection(t_bond,t_ang, t_tor, f, r)
        f0= f*1.0
        q = -f*1.0
        log=1
        coord_size = species.numpy().size*3.0
###########convergence flag##################
        f_max_flag='NO'
        f_rms_flag='NO'
        dp_max_flag='NO'
        dp_rms_flag='NO'
        con_flag = 'NO'
        a = torch.empty((memory,), dtype=torch.float32)
######Initial contrain bond, angle, torsion angle########
        d0=[]
        b0=[]
        t0=[]
        d0=np.array([ d_con0(coordinates,t_bond[i],t_bond[i+1]) for i in range(0,int(t_bond.size),2)])
     
        b0=np.array([ b_ang0(coordinates,t_ang[i],t_ang[i+1],t_ang[i+2]) for i in range(0,int(t_ang.size),3)])

        t0=np.array([ t_ang0(coordinates,t_tor[i]-1,t_tor[i+1]-1,t_tor[i+2]-1,t_tor[i+3]-1) for i in range(0,int(t_tor.size),4)])
        
###########iteration#########################
        while con_flag == 'NO' and iteration < maxiteration:
                start_time = time.time()
                if iteration > 0:
                   s0 = r - r0
                   s.append(s0)
                   y0 = f0 - f
                   y.append(y0)

                   rho0 = 1.0 / (y0*s0).sum()
                   rho.append(rho0)

                if iteration > memory:
                   s.pop(0)
                   y.pop(0)
                   rho.pop(0)
                   
                loopmax = np.min([memory, iteration])
                
               
               # print('loopmax:',loopmax)
        	# ## The algorithm itself:
                q=g*1.0
                
       	        for i in range(loopmax - 1, -1, -1):
            	    a[i] = rho[i] * (s[i]* g).sum()
            	    q -= a[i] * y[i]
                z = H0 * q

                for i in range(loopmax):
            	    b = rho[i] * (y[i]*z).sum()
            	    z += s[i] * (a[i] - b)
                #print('i:',loopmax)
       	        p = -1.0*z.data
        # ##

                p = force_projection(t_bond,t_ang, t_tor, p, r)
                #print('p:',p.numpy())
                dr.data=p*1.0
                if use_line_search==True:
                   ls = LineSearch()
                   alph, e, e0, no_update = \
                     ls._line_search(model,force_projection, species, t_bond,t_ang, t_tor, r, p, g, e, e0, maxstep=0.2, c1=.23,c2=.46, stpmax=50.)
                   if alph is None:
                      print('LBFGS_WS failed!, try LBFGS,CG_WS or CG_BS')
                      break
                longest_step =torch.norm(dr.norm(p=2,dim=2),float('inf'))
        
                if longest_step >= maxstep:
                   dr.data *= maxstep/longest_step*1.0	 		
                
                r0.data=r.data
                r=r0+alph*dr
               
                iteration += 1
                f0 = -1.0*g
                e0 = e
                start_time=time.time()
                
                e = model((species, r)).energies
                #print(time.time()-start_time)
                
                g = torch.autograd.grad(e.sum(), r)[0]
               
               
                #print(time.time()-start_time)
                f = -g*1.0
                f = force_projection(t_bond,t_ang, t_tor, f, r)
###########convergence_FLAG##############################
               
                max_dp=dr.norm(float('inf'))
                rms_dp=torch.sqrt(dr.norm(p=2)**2/coord_size)
                max_f=torch.norm(f.norm(p=2,dim=2),float('inf'))
                rms_f =torch.sqrt(f.norm(p=2)**2/coord_size)
                print(iteration,e.item(), max_f.item(), max_dp.item(), alph)
                
                f_max_flag, f_rms_flag,dp_max_flag,dp_rms_flag = convergence_fun(max_f,rms_f,max_dp,rms_dp,f_max_th,f_rms_th,dp_max_th,dp_rms_th)
                if  f_max_flag=='YES'  and f_rms_flag=='YES' and dp_max_flag=='YES'and dp_rms_flag=='YES':     
                    return  r, iteration,max_f,rms_f,max_dp,rms_dp,f_max_flag,f_rms_flag,dp_max_flag,dp_rms_flag
        return  r, iteration,max_f,rms_f,max_dp,rms_dp,f_max_flag,f_rms_flag,dp_max_flag,dp_rms_flag             	

                   
                	

def CG_WS(model, species, t_elem,coord, t_bond,t_ang, t_tor, f_max_th, f_rms_th, dp_max_th, dp_rms_th,   maxstep=0.2, maxiteration=10000):
        """Parameters:

        """
        iteration = 0    
        rho = []
        alph = 0.05
        r = coord*1.0
        r0 = coord*1.0
        dr =r-r0
        e = model((species, coord)).energies
        e0 = e*1.0
        g = torch.autograd.grad(e.sum(),  coord)[0]
        force = -1.0*g
        force = force_projection(t_bond,t_ang, t_tor, force, r)
        
        g_force = force*1.0
        h_force = force*1.0
        log=1
        beta =1.0
        coord_size = species.numpy().size*3.0
###########convergence flag##################
        f_max_flag='NO'
        f_rms_flag='NO'
        dp_max_flag='NO'
        dp_rms_flag='NO'
        con_flag = 'NO'
        max_dp=dr.norm(float('inf'))
        rms_dp=torch.sqrt(dr.norm(p=2)**2/species.numel()/3)
        max_f=force.norm(float('inf'))
        rms_f =torch.sqrt(force.norm(p=2)**2/species.numel()/3) 
######Initial contrain bond, angle, torsion angle########
        
                #alph=0.9999*2*torch.abs((g_force*h_force).sum())/h_force.norm(p=2)    
        while con_flag == 'NO' and iteration < maxiteration: 
                e0=e*1.0
                ls = LineSearch()
                alph, e, e0, no_update = \
                  ls._line_search(model,force_projection, species, t_bond,t_ang, t_tor, r, h_force, g, e, e0, maxstep=0.2, c1=.0001,c2=.49, stpmax=10.)
                if alph is None:
                   print('CG_WS failed!, try LBFGS or CG_BS')
                   return r, iteration,max_f,rms_f,max_dp,rms_dp,f_max_flag,f_rms_flag,dp_max_flag,dp_rms_flag
                dr.data=alph*h_force.data
                max_dr=torch.norm(dr.norm(p=2,dim=2),float('inf'))
                if max_dr>maxstep:
                   beta=1/max_dr*maxstep 
     #########################################################################
                r0 = r.data*1.0
                
                
                r = r0+beta*dr
                e = g_au* model((species, r)).energies
                g = torch.autograd.grad(e.sum(), r)[0] 
                force = -g*1.0 
                force = force_projection(t_bond,t_ang, t_tor, force, r)
                
                
#################################################################################
               
                max_dp=dr.norm(float('inf'))
                rms_dp=torch.sqrt(dr.norm(p=2)**2/species.numel()/3)
                max_f=force.norm(float('inf'))
                rms_f =torch.sqrt(force.norm(p=2)**2/species.numel()/3) 
                print(iteration,e.item(), max_f.item(), max_dp.item(), alph)   
    ###########################################################################################   
                f_max_flag, f_rms_flag,dp_max_flag,dp_rms_flag = convergence_fun(max_f,rms_f,max_dp,rms_dp,f_max_th,f_rms_th,dp_max_th,dp_rms_th)       
                if  f_max_flag=='YES'  and f_rms_flag=='YES' and dp_max_flag=='YES'and dp_rms_flag=='YES':
                	con_flag = 'YES'               
                b=((force-g_force)*force).sum()
               
               
                c=g_force.norm(p=2)**2
                if  iteration>=0 and b>0 and c>0.0:
     
             
                    gama = b/c
             
                else:
     
                    gama = 0 
             
     
                #print('gama:', gama)
                g_force = force*1.0
                h_force = g_force+gama*h_force 
                beta=1.0
                iteration+=1
        return  r, iteration,max_f,rms_f,max_dp,rms_dp,f_max_flag,f_rms_flag,dp_max_flag,dp_rms_flag 
def CG_BS(model, species, t_elem,coord, t_bond,t_ang, t_tor, f_max_th, f_rms_th, dp_max_th, dp_rms_th, maxstep=0.2, maxiteration=10000):
        """Parameters:

        """
        iteration = 0    
        rho = []
        alph = 0.2/g_au
        max_alph=0.8/g_au
        r = coord*1.0
        r0 = coord*1.0
        dr =r-r0
        e = model((species, coord)).energies
        e0 = e*1.0
        g = torch.autograd.grad(e.sum(),  coord)[0]
        force = -1.0*g
        force = force_projection(t_bond,t_ang, t_tor, force, r)
        
        g_force = force*1.0
        h_force = force*1.0
        log=1
        beta =1.0
        coord_size = species.numpy().size*3.0
###########convergence flag##################
        f_max_flag='NO'
        f_rms_flag='NO'
        dp_max_flag='NO'
        dp_rms_flag='NO'
        con_flag = 'NO'
        max_dp=dr.norm(float('inf'))
        rms_dp=torch.sqrt(dr.norm(p=2)**2/species.numel()/3)
        max_f=force.norm(float('inf'))
        rms_f =torch.sqrt(force.norm(p=2)**2/species.numel()/3) 
        opt_r = r*1.0
        opt_e = e *1.0
        cg_flag=1
        sr=1
######Initial contrain bond, angle, torsion angle########
        
                #alph=0.9999*2*torch.abs((g_force*h_force).sum())/h_force.norm(p=2)    
        while con_flag == 'NO' and iteration < maxiteration: 
                e0=e*1.0        
                
                dr.data=alph*h_force.data
                max_dr=torch.norm(dr.norm(p=2,dim=2),float('inf'))
                if max_dr>maxstep and beta==1.0:
                   beta=1/max_dr*maxstep 
     #########################################################################
                r0 = r.data*1.0
                
                
                r = r0+beta*dr
                e = g_au* model((species, r)).energies
                g = torch.autograd.grad(e.sum(), r)[0] 
                force = -g*1.0 
                force = force_projection(t_bond,t_ang, t_tor, force, r)
                aa=-g_force*h_force
                if  e-e0>0.001*alph*beta*aa.sum() and beta>0.02/g_au:                    
                    beta=beta*0.5
                    r = r0*1.0
                    e = e0*1.0
                    cg_flag=0
                    #print('beta:', beta)
                    continue
                if  e < opt_e:
                    opt_r = r*1.0
                    opt_e = e *1.0
                    
                    
                if  alph<max_alph:
                    alph = alph*1.01
#################################################################################
               
                max_dp=dr.norm(float('inf'))
                rms_dp=torch.sqrt(dr.norm(p=2)**2/species.numel()/3)
                max_f=force.norm(float('inf'))
                rms_f =torch.sqrt(force.norm(p=2)**2/species.numel()/3) 
                print(iteration,e.item(), max_f.item(), max_dp.item(),alph*beta)   
    ###########################################################################################   
                f_max_flag, f_rms_flag,dp_max_flag,dp_rms_flag = convergence_fun(max_f,rms_f,max_dp,rms_dp,f_max_th,f_rms_th,dp_max_th,dp_rms_th)       
                if  f_max_flag=='YES'  and f_rms_flag=='YES' and dp_max_flag=='YES'and dp_rms_flag=='YES':
                    con_flag = 'YES'               
                        
                b=((force-g_force)*force).sum()      
                c=g_force.norm(p=2)**2
                if  cg_flag>0 and b>0 and c>0.0:            
                    gama = b/c
             
                else:
                    gama = 0  

                #print('gama:', gama)
                g_force = force*1.0
                h_force = g_force+gama*h_force 
                beta=1.0
                cg_flag = 1
                iteration+=1
        if e>opt_e:
              print('No convergence, restart using opt_r')
        return  opt_r, iteration,max_f,rms_f,max_dp,rms_dp,f_max_flag,f_rms_flag,dp_max_flag,dp_rms_flag         
#################################################################################
# Let's now manually specify the device we want TorchANI to run:
if __name__ == '__main__':

###############################################input##########
     input_file_name=''
     output_file_name=''
     output_log_name=''
     energy_force_log =''
     gpuid=''
     Method=''
     max_iteration_default=10000
     Converge_kind=np.array([1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8])
     try:
         myopts, args = getopt.getopt(sys.argv[1:],"i:o:m:g:l:")
     except getopt.GetoptError as err:
         print ('Error: ',"%s"%(str(err)))
         print ('Usage: _*.py  -i input.ani -o output.log -g 0 or 1 or other gpu id and no parameter for cpu -m 1-CG-BS,2-CG-WS,3-LBFGS,4-LBFGS-WS,5-BFGS,0-single point')
         sys.exit(2)
     for o, a in myopts:
         if o == "-i":
             input_file_name = a
         elif o == "-o":
             output_file_name = a
         elif o == "-g":
             gpuid = a
             print ('Usage:',gpuid)
         elif o == "-l":
             output_log_name = a
             energy_force_log = open(output_log_name,'w')
         elif o == "-m":
             Method = a
             
         else:
             print ('Usage: _*.py  -i input.ani -o output.log -g 0 or 1 or other gpu id and no parameter for cpu -m 1-CG-BS,2-CG-WS,3-LBFGS,4-LBFGS-WS,5-BFGS,0-single point')
             sys.exit(0)
     try:
          
          input_file  = open(input_file_name, 'r') # File to be processed
          text = input_file.readlines()
          input_file.close()  
         
     except IOError as err:
         print ('Error: ',"%s"%(str(err)))
         sys.exit(1)
     ###############################################################################
     # Let's now load the built-in ANI-1ccx models. The builtin ANI-1ccx contains 8
     # models trained with diffrent initialization. Predicting the energy and force
     # using the average of the 8 models outperform using a single model, so it is
     # always recommended to use an ensemble, unless the speed of computation is an
     # issue in your application.
     #
     # The ``periodic_table_index`` arguments tells TorchANI to use element index
     # in periodic table to index species. If not specified, you need to use
     # 0, 1, 2, 3, ... to index species
     if gpuid=='':
        device = torch.device('cpu')
     else:
        device = torch.device("cuda:"+gpuid if torch.cuda.is_available() else "cpu")
     
     model = torchani.models.ANI1x(periodic_table_index=True).to(device)
     
     ###############################################################################
     # Now let's define the coordinate and species. If you just want to compute the
     # energy and force for a single structure like in this example, you need to
     # make the coordinate tensor has shape ``(1, Na, 3)`` and species has shape
     # ``(1, Na)``, where ``Na`` is the number of atoms in the molecule, the
     # preceding ``1`` in the shape is here to support batch processing like in
     # training. If you have ``N`` different structures to compute, then make it
     # ``N``.
     #
     
     #########################read coordinates and species#####################################################
     max_iteration = text[extract(text,'#MAXSTEP')+1:]
     max_iteration = max_iteration[:extract(max_iteration,'!')]
     max_iteration = np.array(max_iteration[0].split()).astype(int)[0]
     print ('max_iteration:',max_iteration)
     if max_iteration<=1:
        max_iteration=max_iteration_default
     converge = text[extract(text,'#CONVERGE')+1:]
     converge = converge[:extract(converge,'!')]
     converge = np.array(converge[0].split()).astype(int)[0]
     ######read elem#############
     elem = text[extract(text,'#ELEM')+1:]
     elem = elem[:extract(elem,'#ATOMIC_NUM')]
     ######ATOMIC_NUM#############
     a_num = text[extract(text,'#ATOMIC_NUM')+1:]
     a_num = a_num[:extract(a_num,'#COORD')]
     ###read coordinates########
     coord= text[extract(text,'#COORD')+1:]
     coord= coord[:extract(coord,'#CONS_BOND')]
     #######read constrain#######
     cons_bond = text[extract(text,'#CONS_BOND')+1:]
     cons_bond = cons_bond[:extract(cons_bond,'#CONS_ANGLE')]
     
     cons_ang =text[extract(text,'#CONS_ANGLE')+1:]
     cons_ang =cons_ang[:extract(cons_ang,'#CONS_TOR')]
     
     cons_tor =text[extract(text,'#CONS_TOR')+1:]
     cons_tor =cons_tor[:extract(cons_tor,'#END')]
     #text = text[:extract(text,'--Link1--')-3]
     flag=0
     t_spec = [] 
     t_coord = []
     t_bond = []
     t_ang = []
     t_tor = []
     t_elem = []
     for line in elem:
          column=line.split()
          t_elem.extend(column)
     for line in a_num:
          column=line.split()
          t_spec.extend(column)
     t_spec=np.array(t_spec)
     species = torch.tensor([t_spec.astype(int)], device=device)
     for line in coord:
          column=line.split()
          t_coord.extend(column)
     t_coord=np.array(t_coord)
     t_coord= t_coord.reshape(int(t_coord.size/3),3)
     coordinates=torch.tensor([t_coord.astype(float)],requires_grad=True,dtype=torch.float32, device=device)
     for line in cons_bond:
          column=line.split()
          t_bond.extend(column)
     t_bond=np.array(t_bond).astype(int)
     for line in cons_ang:
          column=line.split()
          t_ang.extend(column)
     t_ang=np.array(t_ang).astype(int)
    
     for line in cons_tor:
          column=line.split()
          t_tor.extend(column)
     t_tor=np.array(t_tor).astype(int)
     
     
     ################################################################################################################
     
     ############################intial############################
     d0=[]
     b0=[]
     t0=[]
     d0=np.array([ d_con0(coordinates,t_bond[i],t_bond[i+1]) for i in range(0,int(t_bond.size),2)])
     
     b0=np.array([ b_ang0(coordinates,t_ang[i],t_ang[i+1],t_ang[i+2]) for i in range(0,int(t_ang.size),3)])

     t0=np.array([ t_ang0(coordinates,t_tor[i]-1,t_tor[i+1]-1,t_tor[i+2]-1,t_tor[i+3]-1) for i in range(0,int(t_tor.size),4)])
     d1=d0
     b1=b0
     t1=t0
     print('gpuid:',gpuid)
     print ('MAXSTEP:',max_iteration)
     print ('Converge:',converge)
     print ('Species:',species)
     print ('Coordinates:',coordinates)
     print ('Bond constrain index:',t_bond)
     print ('Bond constrain number list:',d0)
     print ('Ang constrain index:',t_ang)
     print ('Ang constrain number list:',b0)
     print ('Tor constrain index:',t_tor)
     print ('Tor constrain number list:',t0)
     ###############################################################################
     # Now let's compute energy and force:
     phi_flag=np.zeros(int(t_tor.size))
     coordinates0=coordinates*1
     g_au =  1 #27.21138505
     

     
     for i in range(0,int(t_tor.size),4):
              if i<int(t_tor.size)-4:
                   if t_tor[i+1]-1==t_tor[i+4]-1  and t_tor[i+2]-1==t_tor[i+5]-1   and t_tor[i+3]-1==t_tor[i+6]-1:
                         phi_flag[i]=1
               
     energy = g_au*model((species, coordinates)).energies
     energy0 = energy
     derivative = torch.autograd.grad(energy.sum(), coordinates)[0]
     force = -derivative
     

   
 ############################force projection############################
     force = force_projection(t_bond,t_ang, t_tor, derivative, coordinates)
                       
    
     ########opt coordinates############
 
     ###########################################
     ############################################################################
     
     alph = 0.05/g_au
     max_alph=0.5/g_au
     f_max_th=0.00045*Converge_kind[converge-1]
     
     f_rms_th=0.0003*Converge_kind[converge-1]
     
     dp_max_th=0.0018*Converge_kind[converge-1]
     
     dp_rms_th=0.0012*Converge_kind[converge-1]
    
     f_max_flag='NO'
     f_rms_flag='NO'
     dp_max_flag='NO'
     dp_rms_flag='NO'
     con_flag = 'NO'
     coord_size=species.numpy().size*3.0
     
     
     dp=coordinates-coordinates*1.0
     max_dp=dp.norm(float('inf'))
     rms_dp=torch.sqrt(dp.norm(p=2)**2/coord_size)
     max_f=force.norm(float('inf'))
     rms_f =torch.sqrt(force.norm(p=2)**2/coord_size)
     iteration=0
     if energy_force_log!='':
        log=1
     start_time = time.time()
     if Method=='5':
        coordinates, iteration,max_f,rms_f,max_dp,rms_dp,f_max_flag,f_rms_flag,dp_max_flag,dp_rms_flag=BFGS(model, species, t_elem,coordinates, t_bond,t_ang, t_tor, f_max_th, f_rms_th, dp_max_th, dp_rms_th, use_line_search=False,memory=100,  curvature=70.0,  maxstep=0.2, maxiteration=max_iteration)
     if Method=='4':
        coordinates, iteration,max_f,rms_f,max_dp,rms_dp,f_max_flag,f_rms_flag,dp_max_flag,dp_rms_flag=LBFGS(model, species, t_elem,coordinates, t_bond,t_ang, t_tor, f_max_th, f_rms_th, dp_max_th, dp_rms_th, use_line_search=True,memory=100,  curvature=70.0,  maxstep=0.2, maxiteration=max_iteration)
     if Method=='3':
        coordinates, iteration,max_f,rms_f,max_dp,rms_dp,f_max_flag,f_rms_flag,dp_max_flag,dp_rms_flag=LBFGS(model, species, t_elem,coordinates, t_bond,t_ang, t_tor, f_max_th, f_rms_th, dp_max_th, dp_rms_th, use_line_search=False,memory=100,  curvature=70.0,  maxstep=0.2, maxiteration=max_iteration)
     if Method=='2':
        coordinates, iteration,max_f,rms_f,max_dp,rms_dp,f_max_flag,f_rms_flag,dp_max_flag,dp_rms_flag=CG_WS(model, species, t_elem,coordinates, t_bond,t_ang, t_tor, f_max_th, f_rms_th, dp_max_th, dp_rms_th, maxstep=0.2, maxiteration=max_iteration)
     if Method=='1':
        coordinates, iteration,max_f,rms_f,max_dp,rms_dp,f_max_flag,f_rms_flag,dp_max_flag,dp_rms_flag=CG_BS(model, species, t_elem,coordinates, t_bond,t_ang, t_tor, f_max_th, f_rms_th, dp_max_th, dp_rms_th, maxstep=0.2, maxiteration=max_iteration)
     energy = model((species, coordinates)).energies
     ax0 = model[0]((species, coordinates)).energies
     ax1 = model[1]((species, coordinates)).energies
     ax2 = model[2]((species, coordinates)).energies
     ax3 = model[3]((species, coordinates)).energies
     ax4 = model[4]((species, coordinates)).energies
     ax5 = model[5]((species, coordinates)).energies
     ax6 = model[6]((species, coordinates)).energies
     ax7 = model[7]((species, coordinates)).energies
     ax=torch.tensor([ax0.item(),ax1.item(),ax2.item(),ax3.item(),ax4.item(),ax5.item(),ax6.item(),ax7.item()])
     rho_ani=ax.std()/species.numel()**0.5
     d1=np.array([ d_con0(coordinates,t_bond[i],t_bond[i+1]) for i in range(0,int(t_bond.size),2)])
     
     b1=np.array([ b_ang0(coordinates,t_ang[i],t_ang[i+1],t_ang[i+2]) for i in range(0,int(t_ang.size),3)])

     t1=np.array([ t_ang0(coordinates,t_tor[i]-1,t_tor[i+1]-1,t_tor[i+2]-1,t_tor[i+3]-1) for i in range(0,int(t_tor.size),4)])
     outfile=open(output_file_name,'w')
     print('        Item             Value        Threshold      Converged',file=outfile)
     
     print('Maximum Force         ', "%10.6f   %10.6f       %s"%(max_f,f_max_th,f_max_flag),file=outfile)
     print('RMS     Force         ', "%10.6f   %10.6f       %s"%(rms_f,f_rms_th,f_rms_flag),file=outfile)
     print('Maximum Displacement  ', "%10.6f   %10.6f       %s"%(max_dp,dp_max_th,dp_max_flag),file=outfile)
     print('RMS     Displacement  ', "%10.6f   %10.6f       %s"%(rms_dp,dp_rms_th,dp_rms_flag),file=outfile)
     print('Energy:',energy.item(),file=outfile)
     print('rho:',rho_ani.item(),file=outfile)
     print('bond length:',d0, d1,file=outfile)
     print('bond angle:',b0, b1,file=outfile)
     print('Toristion angle:',t0, t1,file=outfile)
     print('Step:',iteration,file=outfile)
     print('Time:',time.time()-start_time ,file=outfile)
     for i in range(0,int (species.size()[1])):
            print ("%-6s%5d %-4s%4s%6d    %16.10f%16.10f%16.10f%6.2lf%6.2lf%12s"%('ATOM',i+1,t_elem[i]+str(i),'MOL',1, coordinates[0,i,0],coordinates[0,i,1],coordinates[0,i,2],1.0,0.0,t_elem[i]), file=outfile)
     outfile.close() 










        
