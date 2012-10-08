from __future__ import division
import math

class spline:
    def cubic(self):
        print 'cubic spline'
    
class homework6:
    sxV = []
    def printHeader(self):
        print 'Principles of Numerical Computation'
        print 'David Parker'
        print 'Homework 6, November 6, 2011'

    def eqn(self,x):
        return float(math.sqrt(x))

    def sx(self,z,h,t,y,x,i):
        self.sxV[i] = (z[i]/(6.0*h[i]))*pow((t[i+1] - x),3.0) + \
            (z[i+1]/(6.0*h[i]))*pow((x-t[i]),3.0) + \
            ((y[i+1]/h[i]) - (z[i+1]*h[i]/6.0))*(x - t[i]) + \
            ((y[i]/h[i]) - (z[i]*h[i]/6.0))*(t[i+1] - x)
        return self.sxV[i]

    def example653(self):
        print '\nSpline problem page 352-354'
        sp = spline()
        knots = 10.0
        pts = 36.0
        incr = float(2.5 / knots)
        t = []
        y = []
        i = 0
        j = 0
        x = 0
        n = 9
        while i <= 2.25:
            t.append(i) 
            y.append(self.eqn(i))
            i += incr
        print 't',t
        print 'y',y

        # algorithm
        h = [0 for col in range(n+1)]
        b = [0 for col in range(n+1)]
        u = [0 for col in range(n+1)]
        v = [0 for col in range(n+1)]
        z = [0 for col in range(n+1)]
        for i in range(0,n):
            h[i] = t[i+1] - t[i]
            b[i] = 6.0*(y[i+1] - y[i]) / h[i]
        u[1] = 2.0*(h[0] + h[1])
        v[1] = b[1] - b[0]

        for i in range(2,n):
            u[i] = 2.0*(h[i]+h[i-1]) - pow(h[i-1],2.0)/u[i-1]
            v[i] = b[i] - b[i-1] - h[i-1]*v[i-1]/u[i-1]
        z[n] = 0

        for i in range(n-1,0,-1):
            z[i] = (v[i] - h[i]*z[i+1])/u[i]
        z[0] = 0
        print 'z',z
        # end algorithm

        # this problem
        # create fx and sxV
        self.sxV = [0 for col in range(n+1)]
        while x <= 2.25 and j < n:
            fx = self.eqn(x)
            self.sxV[j] = self.sx(z,h,t,y,x,j)
            ex = self.sxV[j] - fx
            print 'x = ',x,'f(x) = ',fx,'s(x)',self.sxV[j],'e(x)',ex

            x += 2.25/pts
            j += 1
            
        for i in range(0,n):
            print 
        
printDetails = False
hw6 = homework6()
hw6.printHeader()
hw6.example653()

