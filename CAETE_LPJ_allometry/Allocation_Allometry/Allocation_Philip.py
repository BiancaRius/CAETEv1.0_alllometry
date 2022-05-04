import numpy as np
import matplotlib.pyplot as plt

class Tree:

    S = 0.0 #sapwood    
    H = 0.0 #heartwood
    L = 0.0 #leaf biomass
    R = 0.0 #root biomass

    def __init__(self, S = 0, H =0, L =0  , R = 0):
        self.S = S
        self.H = H
        self.L = L
        self.R = R


# Allocation class that takes one tree as input
# I put the all the constants in that class, they could also be part of the Tree class

class Allocation:
    kLatoSa = 6000.0 #leaf_area:sapwood_area
    ltor = 0.77302587552347657 #leaf:root 
    dw = 200.0 #wood density (variant trait)
    SLA = 15.365607091853349 #specific leaf area (variant trait)
    k2 = 36.0
    k3 = 0.22

    tol = 0.00000001 #tolerance for iteration

    tree = Tree()

    #cmass_leaf_inc
    DeltaL = 0.0

    def __init__(self, tree):
        self.tree = tree

    def Distribute(self, bminc): #biomass increment : basically NPP // The next tau1,2,3 and SS are fixed, they ensure that the allometric equations are respected
        self.tau1 = np.power(self.k2 , 2.0 / self.k3) * 4.0 / 3.14159 / self.dw
        self.tau2 = 1.0 + 2.0 / self.k3
        self.tau3 = self.kLatoSa / self.dw / self.SLA
        self.SS = self.tree.S + bminc - self.tree.L / self.ltor + self.tree.R

        

        #Main Code: Use the bisection method to solve for the leaf mass increment
        self.DeltaL = self.bisection_method(0.0, 3.0)
        #Once we have the leaf mass increment we can cant get root mass increment based on the LTOR constant
        self.DeltaR = (self.DeltaL + self.tree.L) / self.ltor - self.tree.R
        #Finally using the cmass_increment mass conservation we can calculate sapwood increment
        self.DeltaS = bminc - self.DeltaL - self.DeltaR

        print('delta_leaf', self.DeltaL)
        print('delta_root', self.DeltaR)
        print('delta_sap', self.DeltaS)

    # Minimization funnction with delta_cmass_leaf as x
    # A value of x is searched with f(x) = 0
    def f(self, x):
        return self.tau1 * (self.SS - x - x / self.ltor + self.tree.H) \
               - ((self.SS - x - x / self.ltor) / (self.tree.L + x) * self.tau3) ** self.tau2

    # Numerical root-finding method to find a solution x with f(x) = 0
    # there are still some stop-criteria missing, e.g. stop distribution after n=40
    def bisection_method(self, a, b):
        if self.f(a) * self.f(b) > 0:
            # end function, no root.
            return -2.0
        else:
            nstep = 0
            while (b - a) / 2.0 > self.tol:
                midpoint = (a + b) / 2.0
                if self.f(midpoint) == 0:
                    return (midpoint)  # The midpoint is the x-intercept/root.
                elif self.f(a) * self.f(midpoint) < 0:  # Increasing but below 0 case
                    b = midpoint
                else:
                    a = midpoint
                nstep += 1

            #print(nstep)
            return (midpoint)

def main():

    oldTree = Tree(S= 29.790591253578555, H = 108.91909828977032, L =1.2279169651518438, R = 0.88026193814051623)
    #oldTree = Tree(S= 12.790591253578555, H = 18.91909828977032, L =1.2279169651518438, R = 0.88026193814051623)
    #oldTree = Tree(S= 00, H = 0, L =0.001, R = 00)
    #Set up allocation Class with an old tree Example
    allocation = Allocation(oldTree)

    #Create a range of bminc from 1.55. to 6.5 kgC/m²
    data = np.zeros((100, 4))
    for i in range(0,100):
        bminc = np.interp(i, (0, 99), (0.55, 6.5))
        allocation.Distribute(bminc)
        data[i] = bminc, allocation.DeltaL/bminc, allocation.DeltaR/bminc, allocation.DeltaS/bminc

    fig = plt.figure()
    ax = plt.axes()
    lineDL = plt.plot(data[:,0], data[:,1], label = "$\Delta L$")
    lineDR = plt.plot(data[:,0], data[:,2], label = "$\Delta R$")
    lineDS = plt.plot(data[:,0], data[:,3], label = "$\Delta S$")
    ax.legend(loc='upper right', bbox_to_anchor=(1.00, 1.00), shadow=True)
    ax.set_xlabel(r'BMinc (NPP) in kgC $m^{-2}$', size=14)
    ax.set_ylabel('Fraction of BMinc (NPP)', size=14)
    plt.show()





if __name__ == '__main__':
    main()