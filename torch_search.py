# flake8: noqa
import numpy as np
import torch
import torchani

pymin = min
pymax = max

class LineSearch:
    def __init__(self,  xtol=1e-14):

        self.xtol = xtol
        self.task = 'START'
        self.isave = np.zeros((2,), np.intc)
        self.dsave = np.zeros((13,), float)
        self.fc = 0
        self.gc = 0
        self.case = 0
        self.old_stp = 0
        #self.model = torchani.models.ANI1x(periodic_table_index=True).to('cpu')
        
    def _line_search(self, model, force_projection, species,t_bond,t_ang, t_tor, xk, pk, gfk, old_fval, old_old_fval,
                     maxstep=.2, c1=.1, c2=0.46, xtrapl=1.1, xtrapu=4.,
                     stpmax=10., stpmin=1e-8, args=()):
        self.stpmin = stpmin
        
        self.hk=pk*1.0
        np_pk=pk.data.numpy() 
        self.xk=xk*1.0
        self.pk = np_pk.reshape(-1)
        # ??? p_size = np.sqrt((pk **2).sum())
        self.stpmax = stpmax
        self.xtrapl = xtrapl
        self.xtrapu = xtrapu
        self.maxstep = maxstep
        phi0 = old_fval.item()
        g_pk=(gfk*pk).sum()
        derphi0 = g_pk.item()
        
        self.dim = len(self.pk)
        self.gms = np.sqrt(self.dim) * maxstep
        #alpha1 = pymin(maxstep,1.01*2*(phi0-old_old_fval)/derphi0)
        alpha1 = 0.5
        self.no_update = False

        g_au =  1.0
        dr=self.xk-self.xk*1.0
        fval = old_fval
        gval = gfk
        self.steps=[]

        while True:
            stp = self.step(alpha1, phi0, derphi0, c1, c2,
                                             self.xtol,
                                             self.isave, self.dsave)

            if self.task[:2] == 'FG':
                alpha1 = stp
                dr.data=self.hk.data
                coord= self.xk + stp * dr
                energy =g_au* model((species,  coord)).energies
               
                fval = energy.item()
                self.fc += 1
                #stp = 0.5
                
                #print('coord:',stp)
                derivative = torch.autograd.grad(energy.sum(),  coord)[0]
                gval=derivative
                if t_tor.size>=4:
                   gval = force_projection(t_bond,t_ang, t_tor, derivative, coord)
                self.gc += 1
                
                phi0 = fval
                g_pk=(gval*pk).sum()
                derphi0 = g_pk.item()
                self.old_stp = alpha1
                if self.no_update == True:
                    break
            else:
                break

        if self.task[:5] == 'ERROR' or self.task[1:4] == 'WARN':
            stp = None  # failed
        return stp, fval, old_fval, self.no_update

    def step(self, stp, f, g, c1, c2, xtol, isave, dsave):
        if self.task[:5] == 'START':
            # Check the input arguments for errors.
            if stp < self.stpmin:
                self.task = 'ERROR: STP .LT. minstep'
            if stp > self.stpmax:
                self.task = 'ERROR: STP .GT. maxstep'
            if g >= 0:
                self.task = 'ERROR: INITIAL G >= 0'
            if c1 < 0:
                self.task = 'ERROR: c1 .LT. 0'
            if c2 < 0:
                self.task = 'ERROR: c2 .LT. 0'
            if xtol < 0:
                self.task = 'ERROR: XTOL .LT. 0'
            if self.stpmin < 0:
                self.task = 'ERROR: minstep .LT. 0'
            if self.stpmax < self.stpmin:
                self.task = 'ERROR: maxstep .LT. minstep'
            if self.task[:5] == 'ERROR':
                return stp

            # Initialize local variables.
            self.bracket = False
            stage = 1
            finit = f
            ginit = g
            gtest = c1 * ginit
            width = self.stpmax - self.stpmin
            width1 = width / .5
#           The variables stx, fx, gx contain the values of the step,
#           function, and derivative at the best step.
#           The variables sty, fy, gy contain the values of the step,
#           function, and derivative at sty.
#           The variables stp, f, g contain the values of the step,
#           function, and derivative at stp.
            stx = 0
            fx = finit
            gx = ginit
            sty = 0
            fy = finit
            gy = ginit
            stmin = 0
            stmax = stp + self.xtrapu * stp
            self.task = 'FG'
            self.save((stage, ginit, gtest, gx,
                       gy, finit, fx, fy, stx, sty,
                       stmin, stmax, width, width1))
            stp = self.determine_step(stp)
            #return stp, f, g
            #print('stp:',stp)
            return stp
        else:
            if self.isave[0] == 1:
                self.bracket = True
            else:
                self.bracket = False
            stage = self.isave[1]
            (ginit, gtest, gx, gy, finit, fx, fy, stx, sty, stmin, stmax, \
            width, width1) =self.dsave

#           If psi(stp) <= 0 and f'(stp) >= 0 for some step, then the
#           algorithm enters the second stage.
            ftest = finit + stp * gtest
            if stage == 1 and f < ftest and g >= 0.:
                stage = 2

#           Test for warnings.
            if self.bracket and (stp <= stmin or stp >= stmax):
                self.task = 'WARNING: ROUNDING ERRORS PREVENT PROGRESS'
            if self.bracket and stmax - stmin <= self.xtol * stmax:
                self.task = 'WARNING: XTOL TEST SATISFIED'
            if stp == self.stpmax and f <= ftest and g <= gtest:
                self.task = 'WARNING: STP = maxstep'
            if stp == self.stpmin and (f > ftest or g >= gtest):
                self.task = 'WARNING: STP = minstep'

#           Test for convergence.
            if f <= ftest and abs(g) <= c2 * (- ginit):
                self.task = 'CONVERGENCE'

#           Test for termination.
            if self.task[:4] == 'WARN' or self.task[:4] == 'CONV':
                self.save((stage, ginit, gtest, gx,
                           gy, finit, fx, fy, stx, sty,
                           stmin, stmax, width, width1))
                #return stp, f, g
                return stp


#           Call step to update stx, sty, and to compute the new step.

            stx, sty, stp, gx, fx, gy, fy= self.update(stx, fx, gx, sty,
                                               fy, gy, stp, f, g,
                                               stmin, stmax)


#           Decide if a bisection step is needed.

            if self.bracket:
                if abs(sty-stx) >= .66 * width1:
                    stp = stx + .5 * (sty - stx)
                width1 = width
                width = abs(sty - stx)

#           Set the minimum and maximum steps allowed for stp.

            if self.bracket:
                stmin = min(stx, sty)
                stmax = max(stx, sty)
            else:
                stmin = stp + self.xtrapl * (stp - stx)
                stmax = stp + self.xtrapu * (stp - stx)

#           Force the step to be within the bounds maxstep and minstep.

            stp = max(stp, self.stpmin)
            stp = min(stp, self.stpmax)

            if (stx == stp and stp == self.stpmax and stmin > self.stpmax):
                self.no_update = True
#           If further progress is not possible, let stp be the best
#           point obtained during the search.

            if (self.bracket and stp < stmin or stp >= stmax) \
               or (self.bracket and stmax - stmin < self.xtol * stmax):
                stp = stx

#           Obtain another function and derivative.

            self.task = 'FG'
            self.save((stage, ginit, gtest, gx,
                       gy, finit, fx, fy, stx, sty,
                       stmin, stmax, width, width1))
            return stp

    def update(self, stx, fx, gx, sty, fy, gy, stp, fp, gp,
               stpmin, stpmax):
        sign = gp * (gx / abs(gx))

#       First case: A higher function value. The minimum is bracketed.
#       If the cubic step is closer to stx than the quadratic step, the
#       cubic step is taken, otherwise the average of the cubic and
#       quadratic steps is taken.
        if fp > fx:  #case1
            self.case = 1
            theta = 3. * (fx - fp) / (stp - stx) + gx + gp
            s = max(abs(theta), abs(gx), abs(gp))
            gamma = s * np.sqrt((theta / s) ** 2. - (gx / s) * (gp / s))
            if stp < stx:
                gamma = -gamma
            p = (gamma - gx) + theta
            q = ((gamma - gx) + gamma) + gp
            r = p / q
            stpc = stx + r * (stp - stx)
            stpq = stx + ((gx / ((fx - fp) / (stp-stx) + gx)) / 2.) \
                   * (stp - stx)
            if (abs(stpc - stx) < abs(stpq - stx)):
               stpf = stpc
            else:
               stpf = stpc + (stpq - stpc) / 2.

            self.bracket = True

#       Second case: A lower function value and derivatives of opposite
#       sign. The minimum is bracketed. If the cubic step is farther from
#       stp than the secant step, the cubic step is taken, otherwise the
#       secant step is taken.

        elif sign < 0:  #case2
            self.case = 2
            theta = 3. * (fx - fp) / (stp - stx) + gx + gp
            s = max(abs(theta), abs(gx), abs(gp))
            gamma = s * np.sqrt((theta / s) ** 2 - (gx / s) * (gp / s))
            if stp > stx:
                 gamma = -gamma
            p = (gamma - gp) + theta
            q = ((gamma - gp) + gamma) + gx
            r = p / q
            stpc = stp + r * (stx - stp)
            stpq = stp + (gp / (gp - gx)) * (stx - stp)
            if (abs(stpc - stp) > abs(stpq - stp)):
               stpf = stpc
            else:
               stpf = stpq
            self.bracket = True

#       Third case: A lower function value, derivatives of the same sign,
#       and the magnitude of the derivative decreases.

        elif abs(gp) < abs(gx):  #case3
            self.case = 3
#           The cubic step is computed only if the cubic tends to infinity
#           in the direction of the step or if the minimum of the cubic
#           is beyond stp. Otherwise the cubic step is defined to be the
#           secant step.

            theta = 3. * (fx - fp) / (stp - stx) + gx + gp
            s = max(abs(theta), abs(gx), abs(gp))

#           The case gamma = 0 only arises if the cubic does not tend
#           to infinity in the direction of the step.

            gamma = s * np.sqrt(max(0.,(theta / s) ** 2-(gx / s) * (gp / s)))
            if stp > stx:
                gamma = -gamma
            p = (gamma - gp) + theta
            q = (gamma + (gx - gp)) + gamma
            r = p / q
            if r < 0. and gamma != 0:
               stpc = stp + r * (stx - stp)
            elif stp > stx:
               stpc = stpmax
            else:
               stpc = stpmin
            stpq = stp + (gp / (gp - gx)) * (stx - stp)

            if self.bracket:

#               A minimizer has been bracketed. If the cubic step is
#               closer to stp than the secant step, the cubic step is
#               taken, otherwise the secant step is taken.

                if abs(stpc - stp) < abs(stpq - stp):
                    stpf = stpc
                else:
                    stpf = stpq
                if stp > stx:
                    stpf = min(stp + .66 * (sty - stp), stpf)
                else:
                    stpf = max(stp + .66 * (sty - stp), stpf)
            else:

#               A minimizer has not been bracketed. If the cubic step is
#               farther from stp than the secant step, the cubic step is
#               taken, otherwise the secant step is taken.

                if abs(stpc - stp) > abs(stpq - stp):
                   stpf = stpc
                else:
                   stpf = stpq
                stpf = min(stpmax, stpf)
                stpf = max(stpmin, stpf)

#       Fourth case: A lower function value, derivatives of the same sign,
#       and the magnitude of the derivative does not decrease. If the
#       minimum is not bracketed, the step is either minstep or maxstep,
#       otherwise the cubic step is taken.

        else:  #case4
            self.case = 4
            if self.bracket:
                theta = 3. * (fp - fy) / (sty - stp) + gy + gp
                s = max(abs(theta), abs(gy), abs(gp))
                gamma = s * np.sqrt((theta / s) ** 2 - (gy / s) * (gp / s))
                if stp > sty:
                    gamma = -gamma
                p = (gamma - gp) + theta
                q = ((gamma - gp) + gamma) + gy
                r = p / q
                stpc = stp + r * (sty - stp)
                stpf = stpc
            elif stp > stx:
                stpf = stpmax
            else:
                stpf = stpmin

#       Update the interval which contains a minimizer.

        if fp > fx:
            sty = stp
            fy = fp
            gy = gp
        else:
            if sign < 0:
                sty = stx
                fy = fx
                gy = gx
            stx = stp
            fx = fp
            gx = gp
#       Compute the new step.

        stp = self.determine_step(stpf)

        return stx, sty, stp, gx, fx, gy, fy

    def determine_step(self, stp):
        dr = stp - self.old_stp
        x = np.reshape(self.pk, (-1, 3))
        steplengths = ((dr*x)**2).sum(1)**0.5
        maxsteplength = pymax(steplengths)
        if maxsteplength >= self.maxstep:
            dr *= self.maxstep / maxsteplength
        stp = self.old_stp + dr
        return stp

    def save(self, data):
        if self.bracket:
            self.isave[0] = 1
        else:
            self.isave[0] = 0
        self.isave[1] = data[0]
        self.dsave = data[1:]

###############################
    def rotation_BCD_xy(self, coord,i,j,k):
   
        bc = coord[0,i] - coord[0,j]
        dc = coord[0,k] - coord[0,j]
   
        t_n1 = torch.cross(bc,dc)
        Vac=t_n1/t_n1.norm(p=2)
        a=Vac[0].item()
        b=Vac[1].item()
        c=Vac[2].item()
        if a**2+c**2>0:
           Rby = torch.tensor([[c/np.sqrt(a**2+c**2),0,-a/np.sqrt(a**2+c**2)],
                        [0,1,0],
                        [a/np.sqrt(a**2+c**2),0,c/np.sqrt(a**2+c**2)]],dtype=torch.float32) 

        if a**2+c**2==0:
           Rby = torch.tensor([[1,0,0],
                           [0,1,0],
                           [0,0,1]],dtype=torch.float32)
        Rbx = torch.tensor([[1,0,0],
                        [0,np.sqrt(a**2+c**2),-b],
                        [0,b,np.sqrt(a**2+c**2)],
                        ],dtype=torch.float32)  
 
        Rbcd=torch.mm(Rbx,Rby)
    
        return Rbcd
#################################doh define constrain###########################
    def translation_0(self,coord,k0):
    
        trans=torch.ones(coord.size())
        trans[0,:,0]=coord[0,k0,0]*trans[0,:,0]
        trans[0,:,1]=coord[0,k0,1]*trans[0,:,1]
        trans[0,:,2]=coord[0,k0,2]*trans[0,:,2]    
        return coord-trans
    def rotation_x(self, coord,j):
        bc=coord[0,j]
        Vbc=bc/bc.norm(p=2)
        a=Vbc[0].item()
        b=Vbc[1].item()
        c=Vbc[2].item()
    
        if a**2+c**2>0:
           Rby = torch.tensor([[a/np.sqrt(a**2+c**2),0,c/np.sqrt(a**2+c**2)],
                        [0,1,0],
                        [-c/np.sqrt(a**2+c**2),0,a/np.sqrt(a**2+c**2)]],dtype=torch.float32) 

        if a**2+c**2==0:
           Rby = torch.tensor([[1,0,0],
                           [0,1,0],
                           [0,0,1]],dtype=torch.float32)
        Rbz = torch.tensor([[np.sqrt(a**2+c**2),b,0],
                        [-b,np.sqrt(a**2+c**2),0],
                        [0,0,1]],dtype=torch.float32)  
 
        Rb=torch.mm(Rbz,Rby)
    
        return Rb
    def rotation_z(self, coord,j):
        bc=coord[0,j]
        Vbc=bc/bc.norm(p=2)
        a=Vbc[0].item()
        b=Vbc[1].item()
        c=Vbc[2].item()
    
    
        Rbz = torch.tensor([[a,b,0],
                        [-b,a,0],
                        [0,0,1]],dtype=torch.float32)  
 
    
    
        return Rbz
    def rotation_xy(self, Rb,coord,i):
        coord=torch.mm(Rb,coord.squeeze().T)
        ac=coord.T[i]
        Vac=ac/ac.norm(p=2)
        a=Vac[0].item()
        b=Vac[1].item()
        c=Vac[2].item()
        if b**2+c**2>0:
              Rax = torch.tensor([[1,0,0],
                        [0,b/np.sqrt(b**2+c**2),c/np.sqrt(b**2+c**2)],
                    [0,-c/np.sqrt(b**2+c**2),b/np.sqrt(b**2+c**2)]],dtype=torch.float32)
    
        if b**2+c**2<=0.000000000000001:
             Rax = torch.tensor([[1,0,0],
                        [0,1,0],
                    [0,0,1]],dtype=torch.float32)  
        return Rax


    def torision_force_bcd(self, Rbcd,Ra,Rb,Rd,force,i,j,k,l):  
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
    def torision_force(self, Ra,Rb,Rd,force,i,j,k,l):  
       #Rbcd_f=torch.mm(Rbcd,force.squeeze().T)
        Rb_f=torch.mm(Rb,force.squeeze().T)
        Ra_f=torch.mm(Ra,Rb_f)
        Rd_f=torch.mm(Rd,Rb_f)
        Ra_f.T[i,2]=0.0
        Rb_f.T[j,1:3]=0.0
        Rb_f.T[k,1:3]=0.0
        Rd_f.T[l,2]=0.0
        Ra_f=torch.mm(Ra.inverse(),Ra_f)
        Ra_f=torch.mm(Rb.inverse(),Ra_f)
    #Ra_f=torch.mm(Rbcd.inverse(),Ra_f)
        Rb_f=torch.mm(Rb.inverse(),Rb_f)
    #Rb_f=torch.mm(Rbcd.inverse(),Rb_f)
        Rd_f=torch.mm(Rd.inverse(),Rd_f)
        Rd_f=torch.mm(Rb.inverse(),Rd_f)
   # Rd_f=torch.mm(Rbcd.inverse(),Rd_f)
        force[0,i]=Ra_f.T[i]
        force[0,j]=Rb_f.T[j]
        force[0,k]=Rb_f.T[k]
        force[0,l]=Rd_f.T[l]
        return force

    def t_ang1(self, coord,i,j,k,l):
        a = coord[0,j] - coord[0,i]
        b = coord[0,k] - coord[0,j]
        c = coord[0,l] - coord[0,k]
        t_n1 = torch.cross(a,b)
        t_n2 = torch.cross(c,b)
        return np.degrees(torch.acos(torch.dot(t_n1,t_n2)/(t_n1.norm(p=2)*t_n2.norm(p=2))).item())
    def t_ang0(self, coord,i,j,k,l):
    
        coordinates1=coord*1
        Rbcd=rotation_BCD_xy(coord,j,k,l)
        BCD_coordinates=torch.mm(Rbcd, coord.squeeze().T)
        coordinates1[0]=BCD_coordinates.T 
                      
        trc_coord=translation_0(coordinates1,k)
 
        Rb_m=rotation_x(trc_coord,j)
 
    #Rb_coord=torch.mm(Rb_mat,trc_coord.squeeze().T)

 ##################   Ra_mat=rotation_xy(Rb_mat,trc_coord,t_tor[i]-1)
        coord1=torch.mm(Rb_m,trc_coord.squeeze().T)
        aa=coord1.T[i]
        Vaa=aa/aa.norm(p=2)
        ax=Vaa[0].item()
        ay=Vaa[1].item()
        az=Vaa[2].item()
        dd=coord1.T[l]
        Vdd=dd/dd.norm(p=2)
        dy=Vdd[1].item()
    
    
        if ay**2+az**2>0 and ay*ay==0:    
           if dz*az>=0:
              tor=0
           if dz*az<0:
              tor=180
        if ay**2+az**2>0 and ay*dy>=0:
        
           tor=np.degrees(np.arcsin(-az/np.sqrt(ay**2+az**2))) 
        if ay**2+az**2>0 and ay*dy<0:
        
           tor=np.degrees(np.arcsin(-az/np.sqrt(ay**2+az**2))+3.1415926)   
        if ay**2+az**2==0:
           tor=0
    
        return tor

######extract input################## 

#################################################################################
    def force_projection(self, t_tor, force, coordinates):
        if t_tor.size==4:
            flag=1
                          
            coordinates0=coordinates*1
            for i in range(0,int(t_tor.size),4):
                if  i<=int(t_tor.size)-4:
                    Rbcd=self.rotation_BCD_xy(coordinates,t_tor[i+1]-1,t_tor[i+2]-1,t_tor[i+3]-1)
                    BCD_coordinates=torch.mm(Rbcd, coordinates.squeeze().T)
                    coordinates0[0]=BCD_coordinates.T 
                      
                trc_coord=self.translation_0(coordinates0,t_tor[i+2]-1)
 #              print('trc_coord:',trc_coord,t_tor[i+2]-1)
                Rb_mat=self.rotation_x(trc_coord,t_tor[i+1]-1)
 #              print('Rb_mat:',Rb_mat)
                Rb_coord=torch.mm(Rb_mat,trc_coord.squeeze().T)
 #              print('Rb_coord:',Rb_coord.T)
                Ra_mat=self.rotation_xy(Rb_mat,trc_coord,t_tor[i]-1)
                Ra_coord=torch.mm(Ra_mat,trc_coord.squeeze().T)
                #print('Ra_coord:',Ra_coord.T)
                Rd_mat=self.rotation_xy(Rb_mat,trc_coord,t_tor[i+3]-1)
                Rd_coord=torch.mm(Rd_mat,trc_coord.squeeze().T)
                #print('Rd_coord:',Rd_coord.T)
                force=self.torision_force_bcd(Rbcd,Ra_mat,Rb_mat,Rd_mat,force,t_tor[i]-1,t_tor[i+1]-1,t_tor[i+2]-1,t_tor[i+3]-1)
                
        if t_tor.size>4:
            flag=1
                          
            coordinates0=coordinates*1 
            for i in range(0,int(t_tor.size),4):
                if  i<int(t_tor.size)-4:
                    Rbcd=self.rotation_BCD_xy(coordinates,t_tor[i+1]-1,t_tor[i+2]-1,t_tor[i+3]-1)
                BCD_coordinates=torch.mm(Rbcd, coordinates.squeeze().T)
                coordinates0[0]=BCD_coordinates.T 
                      
                trc_coord=self.translation_0(coordinates0,t_tor[i+2]-1)
 #              print('trc_coord:',trc_coord,t_tor[i+2]-1)
                Rb_mat=self.rotation_x(trc_coord,t_tor[i+1]-1)
 #              print('Rb_mat:',Rb_mat)
                Rb_coord=torch.mm(Rb_mat,trc_coord.squeeze().T)
 #              print('Rb_coord:',Rb_coord.T)
                Ra_mat=self.rotation_xy(Rb_mat,trc_coord,t_tor[i]-1)
                Ra_coord=torch.mm(Ra_mat,trc_coord.squeeze().T)
                #print('Ra_coord:',Ra_coord.T)
                Rd_mat=self.rotation_xy(Rb_mat,trc_coord,t_tor[i+3]-1)
                Rd_coord=torch.mm(Rd_mat,trc_coord.squeeze().T)
                #print('Rd_coord:',Rd_coord.T)
                force=self.torision_force_bcd(Rbcd,Ra_mat,Rb_mat,Rd_mat,force,t_tor[i]-1,t_tor[i+1]-1,t_tor[i+2]-1,t_tor[i+3]-1)
                  
                force[0,t_tor[i+2]-1]= 0             
                  
        return force
