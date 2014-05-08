#do 2nd order advection
import numpy
from matplotlib import pyplot as plt
class advect:
    def __init__(self,n=3000,u=1.0,xmax=1.0,C=0.1):
        self.rho=numpy.zeros(n)
        self.rho[n/3:2*n/3]=1.0
        self.u=u
        self.n=n
        self.dx=(0.0+xmax)/(0.0+n)
        self.C=C
    def get_ts(self):
        return numpy.abs(self.dx/self.u)*self.C

    def get_bc(self):
        #we now need to put in 2 ghost cells.  These are periodic
        self.rho[0:2]=self.rho[-4:-2]
        self.rho[-2:]=self.rho[2:4]
    def get_interfaces(self):        
        #calculate the derivatives for the cell centers
        self.myderiv=0.5*(self.rho[2:]-self.rho[0:-2])/self.dx
        #now calculate the flux at the right and left edges of the cell
        #after taking half a step forward in time.
        self.right=self.rho[1:-1]+0.5*self.dx*(1.0-self.C)*self.myderiv
        self.left=self.rho[1:-1]-0.5*self.dx*(1.0+self.C)*self.myderiv
        
    def update(self):
        dt=self.get_ts()
        self.get_bc()
        self.get_interfaces()
        #Solve Riemann problem, and update.  For advection, if velocity is positive,
        #cell to my left flows into me with the flux at its right edge, and I flow into 
        #the cell to my right with my flux at my right edge.
        if self.u>0:
            self.rho[2:-2]=self.rho[2:-2]+(self.right[0:-2]-self.right[1:-1])*self.C
        else:
            self.rho[2:-2]=self.rho[2:-2]+(self.left[1:]-self.left[0:-2])*self.C
        
if __name__=="__main__":

    fac=1000;
    npt=300
    C=1.0/fac
    fwee=advect(n=npt,C=C)
    fwee_org=advect(n=npt,C=C)
    plt.ion()

    #we have 2 ghost cells on either side, so only have npt-4 true cells
    #if we want to get back to the initial state for periodic BC's.
    for i in range(0,npt-4):
        for ii in range(0,fac):
            fwee.update()
        if i%10==-1:  #to plot very often, makes code much faster.
            plt.clf()
            plt.plot(fwee.rho)
            plt.draw()
    print numpy.mean(numpy.abs(fwee.rho-fwee_org.rho))
plt.clf()
plt.plot(fwee.rho)
plt.plot(fwee_org.rho)
plt.draw()
plt.savefig('advect2_output.png')
