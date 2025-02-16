
import numpy as np
import math
import matplotlib.pyplot as plot

# Diffusion- interior points
def diffusion(c, dx):

    ngrid = len(c)
    dc = np.zeros(ngrid)

    dc[1:ngrid-1] = c[0:ngrid-2] - 2.0*c[1:ngrid-1] + c[2:ngrid]

    dc = dc/dx**2

    return dc

# BCs
def bcDirichlet(c, val0, valL):

    c[0] = val0
    c[len(c)-1] = valL

    return c

def bcNeumann(c):

    c[0] = c[1]
    c[len(c)-1] = c[len(c)-2]

    return c


# Move frame - for fronts
def moveFrame1(c):

    retval = 0
    ngrid = len(c)
    halfngrid = int(ngrid/2)

    if ( c[halfngrid] > 0.5 ):
        c[1:ngrid-2] = c[2:ngrid-1]

        retval = 1

    return retval


def moveFrame2(c):

    retval = 0
    ngrid = len(c)
    halfngrid = int(ngrid/2)

    if ( c[halfngrid] < 0.5 ):
        c[2:ngrid-1] = c[1:ngrid-2] 

        retval = 1

    return retval






# Setup systems
nloops = 1000
ngrid  = 100
Lbox = 10.0

Da, k = 1.0, 1.0
# Initial value
x = np.linspace(0, Lbox, ngrid)    

c = np.sin(math.pi*x/Lbox)

# Derived parameters
dx = Lbox/ngrid
dt = 0.1*dx**2 #divided with diffusion coeff.

val = Da*(math.pi/Lbox)**2
print("Value ", Da )
# Loop
t = 0.0
for n in range(nloops):

    # Get diffusion term
    diff_c = diffusion(c, dx)

    # Define reactions
    r = k*c
    
    # Update with Euler
    a = a + dt*(r + Da*diff_c)
    
    # Boundaries
    a = bcDirichlet(a, 0, 0)
        
    # Increment
    t = t + dt

    # Plot and print (again, why is python smart..?)
    if ( n%100 == 0 ):
        plot.clf()
        plot.plot(x, a)
        plot.ylim([0, 1])
        plot.pause(0.1)
        
