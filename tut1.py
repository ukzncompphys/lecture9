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
    #here we add an option to set the initial conditions to be a gaussian
    #placed at the centre of the interval
    def ic_gaussian(self,sig=0.1):
        xvec=numpy.arange(0,self.rho.size)*self.dx
        xvec=xvec-xvec.mean()
        self.rho=numpy.exp(-0.5*xvec**2/sig**2)
        

#we are going to pull a python trick here, which can be very useful.  A first-order and second-order
#class look very similar in many aspects - for instsance the initialization and possibly the boundary 
#conditions are the same (in particular, the 2nd order boundary conditions will also work for a 1st 
#order solver, though the reverse is not true.  However, some parts look different, like evaluating
#the interface conditions.  In python, we can define subclasses, which inherit the methods associated
#with the parent class.  Basically, when I declare
#class advect1(advect):
#that means that advect1 is a subclass of advect, and all routines (like __init__) that are there
#for advect will also be there for advect1.  So for this case, where advect1 is the 1st order advection
#solver, the only extra thing it needs is the update.  The 2nd order class (advect2) needs to have
#the interface solver as well.

class advect1(advect):
    def update(self):
        self.get_bc()
        df=self.rho[0:-1]-self.rho[1:]
        self.rho[1:]+=df*self.C

    def test_update(self):
        rho_org=self.rho.copy()
        #let's calculate what a single step *should* give us
        rhoft=numpy.fft.fft(rho_org)

        kvec=numpy.arange(0,self.rho.size)
        kvec=kvec/(self.rho.size+0.0)*2*numpy.pi

        i=numpy.complex(0,1.0)
        fac=1+self.C*numpy.exp(-i*kvec)-self.C
        newft=rhoft*fac
        dat_check=numpy.real(numpy.fft.ifft(newft)) #is only non-real due to roundoff error

        self.update()
        mean_err=numpy.mean(numpy.abs(dat_check-self.rho))
        if mean_err/numpy.mean(numpy.abs(self.rho))<1e-13:
            print 'test PASSED with mean error ' + repr(mean_err)
        else:
            print 'test FAILED with mean error ' + repr(mean_err)

        self.rho=rho_org


class advect2(advect):
    def get_interfaces(self):        
        #calculate the derivatives for the cell centers
        self.myderiv=0.5*(self.rho[2:]-self.rho[0:-2])/self.dx
        #now calculate the flux at the right and left edges of the cell
        #after taking half a step forward in time.
        self.right=self.rho[1:-1]+0.5*self.dx*(1.0-self.C)*self.myderiv
        #self.left=self.rho[1:-1]-0.5*self.dx*(1.0+self.C)*self.myderiv
        
    def update(self):
        #dt=self.get_ts()
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

    fac=10;
    npt=200
    C=1.0/fac

    fwee1=advect1(n=npt,C=C)
    fwee2=advect2(n=npt,C=C)


    fwee_org=advect2(n=npt,C=C)


    if False:  #we can switch on and off gaussian IC's
        sig=0.1
        fwee1.ic_gaussian(sig=sig)
        fwee2.ic_gaussian(sig=sig)
        fwee_org.ic_gaussian(sig=sig)
        mytag='gaussian'
    else:
        mytag='rectangle'
    plt.ion()

    fwee1.test_update()

    #we have 2 ghost cells on either side, so only have npt-4 true cells
    #if we want to get back to the initial state for periodic BC's.
    for i in range(0,npt-4):
        for ii in range(0,fac):
            fwee1.update()
            fwee2.update()
        if i%10==-1:  #to plot very often, makes code much faster.
            plt.clf()
            plt.plot(fwee1.rho)
            plt.draw()
    print numpy.mean(numpy.abs(fwee1.rho-fwee_org.rho))
    print numpy.mean(numpy.abs(fwee2.rho-fwee_org.rho))
plt.clf()
plt.plot(fwee1.rho)
plt.plot(fwee2.rho)
plt.plot(fwee_org.rho)
plt.draw()
plt.savefig('tut1_' + mytag + '.png')

