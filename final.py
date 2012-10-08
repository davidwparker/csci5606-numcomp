from __future__ import division
import math

class final:
    epsilon = 0.0
    """ final project class """
    def __init__(self):
        self.epsilon = self.machineEpsilon()
    
    def printHeader(self):
        """ Prints the header for the class """
        print '\nPrinciples of Numerical Computation'
        print 'David Parker'
        print 'Final Project, December 3, 2011'

########################################################
##  Run tests on our algorithms to make sure they're ok
########################################################

    def testFunctions(self):
        """ Tests the functionality of the algorithms """
        print '\nTests the functionality of the algorithms'
        print 'Machine epsilon =',self.epsilon
        # Tests for Newton
        det = self.determinant2x2([[1,1],[2,10]])
        assert det == 8
        assert self.matrixScale([[1,2],[-3,1]],2) == [[2,4],[-6,2]]
        assert self.inverse2x2([[1,1],[2,10]],det) == [[1.25,-0.125],[-0.25,0.125]]
        assert self.jacobi2x2(1,5,self.testF1dx,self.testF1dy,self.testF2dx,self.testF2dy) == [[1,1],[2,10]]
        self.testNewton()
        self.test2Newton()
        # Tests for Conjugate Gradient
        assert self.dot([1,2],[3,-4]) == -5
        assert self.vectorScale([1,2],5) == [5,10]
        assert self.multiplyAx([[1,2],[0,1]],[2,1]) == [4,1]
        assert self.vectorSubtract([2,1],[1,-1]) == [1,2]
        assert self.vectorAddition([2,1],[1,-1]) == [3,0]
        self.testConjugateGradient([[4,1],[1,3]],[1,2],[2,1],50)
        # Tests for Cubic Spline
        self.testCubicSpline()
        # Tests for PQ Decomposition
        p,q = self.testPQ(2,[[3, 15],[-1,-1]])
        assert p == [[0,3],[4,-1]]
        assert q == [[0,1],[1,5]]
        p,q = self.testPQ(3,[[2,6,4],[4,15,5],[5,17,7]])
        assert p == [[0,0,2],[0,3,4],[-1,2,5]]
        assert q == [[0,0,1],[0,1,-1],[1,3,2]]
        assert self.solveAxb([[0,3],[4,-1]],[33,-3]) == [2,11]
        assert self.solvePQxb([[0,3],[4,-1]],[[0,1],[1,5]],[33,-3]) == [1,2]

########################################################
##  Algorithms and Functions
########################################################
    def machineEpsilon(self,func=float):
        """ Returns the machine epsilon """
        mEpsilon = func(1)
        while func(1)+func(mEpsilon) != func(1):
            mEpsilonLast = mEpsilon
            mEpsilon = func(mEpsilon) / func(2)
        return mEpsilonLast        

    def multiplyAx(self,a,x):
        """ Return a matrix-vector product """
        ax = [(sum([a[i][j]*x[j] for j in range(len(x))])) for i in range(len(x))]
        if printDetails:
            print 'Matrix multiplied by a vector'
            print 'A  =',a
            print 'x  =',x
            print 'Ax =',ax
        return ax

    def determinant2x2(self,a):
        """ Return the determinant of a 2x2 matrix """
        det = a[1][1]*a[0][0]-a[0][1]*a[1][0]
        if printDetails:
            print 'Determinant of a 2x2 matrix A = a[1][1]*a[0][0]-a[0][1]*a[1][0]'
            print 'A =',a
            print 'Det A =',det
        return det

    def matrixScale(self,a,s):
        """ Return a matrix-constant product (a scaled matrix) """
        scaled = [[s * a[i][j] for j in range(len(a))] for i in range(len(a))]
        if printDetails:
            print 'Matrix scaled by S'
            print 'A =',a
            print 'S =',s
            print 'Scaled matrix =',scaled
        return scaled

    def inverse2x2(self,a,d):
        """ 
        Return the inverse of a 2x2 matrix 
          A^-1 = (1/Determinant(A))(  D -B )
                                   ( -C  A )
        """
        # change positions
        ainv = [[0 for j in range(0,2)] for i in range(0,2)]
        ainv[0][0] =  a[1][1]
        ainv[0][1] = -a[0][1]
        ainv[1][0] = -a[1][0]
        ainv[1][1] =  a[0][0]
        # scale
        ainv = self.matrixScale(ainv,1.0/d)
        if printDetails:
            print 'Inverse of a 2x2 matrix'
            print 'a     =',a
            print 'det a =',d
            print 'ainv  =',ainv
        return ainv

    def jacobi2x2(self,x,y,f1,f2,f3,f4):
        """ Return the Jacobi of a 2x2 matrix """
        j = [[0 for j in range(0,2)] for i in range(0,2)]
        j[0][0] = f1(x,y)
        j[0][1] = f2(x,y)
        j[1][0] = f3(x,y)
        j[1][1] = f4(x,y)
        return j

    def newton(self,x,y,max,f1,f1dx,f1dy,f2,f2dx,f2dy):
        """ Return the newton value for a given 2x2 matrix """

        # Newton => x^(n+1) = x^n - J(x^n)F(xn)
        for i in range(0,max):
            # Setup Jacobi
            # J = ( dx/df1 dy/df1 )
            #     ( dx/df2 dy/df2 )
            print 'x =',x,'y =',y

            f1v = f1(x,y)
            f2v = f2(x,y)
            f = []
            f.append(f1v)
            f.append(f2v)
            print 'f =',f
            # Do I actually want this check in here?
            if abs(f1v) < self.epsilon or abs(f2v) < self.epsilon:
                print 'f1v or f2v below machine epsilon'
                return
       
            j = self.jacobi2x2(x,y,f1dx,f1dy,f2dx,f2dy)
            if printDetails: print 'j =',j

            # Get inverse of 2x2 Jacobi
            jdet = self.determinant2x2(j)
            jinv = self.inverse2x2(j,jdet)
            if printDetails: print 'jinv =',jinv
            
            negjinv = self.matrixScale(jinv,-1)
            if printDetails: print 'negjinv =',negjinv

            # solve for h
            h = self.multiplyAx(negjinv,f)
            if printDetails: print 'h =',h

            # x^n+1 = x^n - h
            x = x + h[0]
            y = y + h[1]
        
    def dot(self,v1,v2):
        """ Return the dot product of two vectors """
        dotProd = sum([v1[i] * v2[i] for i in range(len(v1))])
        if printDetails:
            print 'Dot product for two vectors is the sum of the product of the items'
            print 'v1        =',v1
            print 'v2        =',v2
            print 'v1 dot v2 =',dotProd
        return dotProd

    def vectorScale(self,v,s):
        """ Return a vector-constant product (a scale vector) """
        scaled = [s * v[i] for i in range(len(v))]
        if printDetails:
            print 'Vector scaled by a constant'
            print 'v   =',v
            print 's   =',s
            print 'v*s =',scaled
        return scaled

    def vectorSubtract(self,v1,v2):
        """ Return the differenc of two vectors """
        vsub = [v1[i] - v2[i] for i in range(len(v1))]
        if printDetails:
            print 'Vector subtraction'
            print 'v1      =',v1
            print 'v2      =',v2
            print 'v1 - v2 =',vsub
        return vsub
        
    def vectorAddition(self,v1,v2):
        """ Return the sum of two vectors """
        vadd = [v1[i] + v2[i] for i in range(len(v1))]
        if printDetails:
            print 'Vector addition'
            print 'v1      =',v1
            print 'v2      =',v2
            print 'v1 + v2 =',vadd
        return vadd

    def conjugateGradient(self,x,a,b,max,delta=0.0):
        """ Solve Ax = b with conjugate gradient method """
        print '\nSolving Ax = b with conjugate gradient method'
        if delta == 0.0:
            delta = self.epsilon

        # r = residual
        # r = b - Ax
        r = self.vectorSubtract(b,self.multiplyAx(a,x))
        if printDetails: print 'r=',r
        v = list(r)
        if printDetails: print 'v=',v
        c = self.dot(r,r)
        if printDetails: print 'c=',c

        for k in range(1,max+1):
            if c < self.epsilon:
                print 'finish due to c < epsilon'
                break
            if math.sqrt(self.dot(v,v)) < delta: 
                print 'finish due to sqrt(dot(v,v)) < delta'
                break

            # z = a*v
            z = self.multiplyAx(a,v)
            if printDetails: print 'z=',z

            t = c / self.dot(v,z)
            if printDetails: print 't=',t

            # x = x + tv
            x = self.vectorAddition(x,self.vectorScale(v,t))
            if printDetails: print 'x=',x

            # r = r - tz
            r = self.vectorSubtract(r,self.vectorScale(z,t))
            if printDetails: print 'r=',r

            d = self.dot(r,r)
            if printDetails: print 'd=',d

            # NOTE: I moved this check for d < self.epsilon to be
            #       at the top of the function so that we still
            #       output the last possible k,x,r
            # if (d < self.epsilon): 
            #     print 'd < epsilon'
            #     break

            # v = r + (d/c)*v
            v = self.vectorAddition(r,self.vectorScale(v,d/c))
            c = d
            print 'k=',k,'\nx=',x,'\nr=',r

    def solveHB(self,n,t,y):
        """ Solving for h and b """
        h = [0.0 for i in range(0,n)]
        b = [0.0 for i in range(0,n)]
        # n... (actually n-1)
        for i in range(0,n):
            h[i] = t[i+1] - t[i]
            b[i] = 6.0*(y[i+1] - y[i]) / h[i]
        return h,b

    def solveUV(self,n,h,b):
        """ Solving for u and v """
        u = [0 for i in range(0,n)]
        v = [0 for i in range(0,n)]
        u[1] = 2.0*(h[0]+h[1])
        v[1] = b[1]-b[0]
        # n... (actually n-1)
        for i in range(2,n):
            u[i] = 2.0*(h[i]+h[i-1]) - pow(h[i-1],2)/u[i-1]
            v[i] = b[i] - b[i-1] - h[i-1]*v[i-1]/u[i-1]
        return u,v

    def solveZ(self,n,h,u,v):
        """ Solving for z in Cubic Spline """
        z = [0 for i in range(0,n+1)]
        z[n] = 0
        # 0... (actually 1)
        for i in range(n-1,0,-1):
            z[i] = (v[i] - h[i]*z[i+1])/u[i]
        z[0] = 0
        return z

    def solveA(self,n,h,z):
        """ Solving for a in Cubic Spline """
        return [((1.0/(6.0*h[i]))*(z[i+1]-z[i])) for i in range(0,n)]

    def solveB(self,n,z):
        """ Solving for b in Cubic Spline """
        return [(z[i]/2.0) for i in range(0,n)]

    def solveC(self,n,h,z,y):
        """ Solving for c in Cubic Spline """
        return [(-h[i]/6.0*z[i+1] - h[i]/3.0*z[i] + 1/h[i]*(y[i+1]-y[i])) for i in range(0,n)]

    def solveS(self,n,xlist,y,t,c,b,a):
        """ Solve for s in a Cubic Spline """
        s = [0 for i in range(len(xlist))]
        for j,x in enumerate(xlist):
            for i in range(0,n):
                if t[i] <= x and x < t[i+1]:
                    s[j] = (y[i] + (x-t[i])*(c[i] + (x-t[i])*(b[i] + (x-t[i])*a[i])))

        return s

    def solveE(self,n,s,y):
        """ Solve for the error in a Cubic Spline """
        return [abs(s[i]-y[i]) for i in range(0,n)]

    def cubicSpline(self,n,t,y,xlist,fx,m):
        """ Cubic Spline Algorithm """
        print '\nPerforming Cubic Spline'

        h,b = self.solveHB(n,t,y)
        if printDetails: print 'h =',h,'\nb =',b
        u,v = self.solveUV(n,h,b)
        if printDetails: print 'u =',u,'\nv =',v
        z = self.solveZ(n,h,u,v)
        if printDetails: print 'z =',z
        a = self.solveA(n,h,z)
        if printDetails: print 'a =',a
        b = self.solveB(n,z)
        if printDetails: print 'b =',b
        c = self.solveC(n,h,z,y)
        if printDetails: print 'c =',c

        if printDetails: print 'x =',xlist

        s = self.solveS(n,xlist,y,t,c,b,a)
        if printDetails: print 's =',s

        e = self.solveE(m,s,fx)
        if printDetails: print 'e =',e

        # PRINT OUTPUT
        for i in range(0,m):
            text = 'i= %2d' % i
            text += ' x[i]= %10f' % xlist[i]
            text += ' s[i]= %10f' % s[i]
            text += ' f(x)= %10f' % fx[i]
            text += ' e[i]= %10f' % e[i]
            print text

########################################################
##  NEWTON METHOD TESTS
########################################################

    def testF1(self,x,y):
        return x + y - 3
    def testF1dx(self,x,y):
        return 1
    def testF1dy(self,x,y):
        return 1

    def testF2(self,x,y):
        return pow(x,2) + pow(y,2) - 9
    def testF2dx(self,x,y):
        return 2*x
    def testF2dy(self,x,y):
        return 2*y

    def testNewton(self):
        """ Tests the Newton Method """
        print "\nTesting the Newton Method (Test 1)"
        x = 1
        y = 5
        max = 3
        self.newton(x,y,max,self.testF1,self.testF1dx,self.testF1dy, \
                        self.testF2,self.testF2dx,self.testF2dy)

    def test2F1(self,x,y):
        return pow(x,2) - 2*x - y + 0.5
    def test2F1dx(self,x,y):
        return 2*x - 2
    def test2F1dy(self,x,y):
        return -1

    def test2F2(self,x,y):
        return pow(x,2) + 4*pow(y,2) - 4
    def test2F2dx(self,x,y):
        return 2*x
    def test2F2dy(self,x,y):
        return 8*y
    
    def test2Newton(self):
        """ 
        Tests the Newton Method 
        http://math.fullerton.edu/mathews/numerical/n2.htm
        Exercise 1
        """
        print "\nTesting the Newton Method (Test 2)"
        x = 2.0
        y = 0.25
        max = 5
        self.newton(x,y,max,self.test2F1,self.test2F1dx,self.test2F1dy, \
                        self.test2F2,self.test2F2dx,self.test2F2dy)

########################################################
##  Page 93, Computer Problem 14 
########################################################

    # PART A
    def function14a1(self,x,y):
        """ 
        f(x,y) = 4y^2 + 4y + 52x - 19
        """
        return 4*pow(y,2) + 4*y + 52*x - 19

    def function14a1dx(self,x,y):
        """ 
        f(x,y) = 4y2 + 4y + 52x -19
        dx/df  = 52
        """
        return 52

    def function14a1dy(self,x,y):
        """ 
        f(x,y) = 4y2 + 4y + 52x - 19
        dy/df  = 8y + 4
        """
        return 8*y + 4

    def function14a2(self,x,y):
        """
        f(x,y) = 169x^2 + 3y^2 + 111x - 10y - 10
        """
        return 169*pow(x,2) + 3*pow(y,2) + 111*x - 10*y - 10

    def function14a2dx(self,x,y):
        """
        f(x,y) = 169x^2 + 3y^2 + 111x - 10y - 10
        dx/df  = 338x + 111
        """
        return 338*x + 111

    def function14a2dy(self,x,y):
        """
        f(x,y) = 6y - 10 - 10
        """
        return 6*y - 10

    # PART B
    def function14b1(self,x,y):
        """
        f(x,y) = x + e^-1x + y^3 = 0
        """
        return x + math.exp(-1*x) + pow(y,3)

    def function14b1dx(self,x,y):
        """
        f(x,y) = x + e^-1x + y^3 = 0
        dx/df  = 1 - e^-1x
        """
        return 1 - math.exp(-1*x)

    def function14b1dy(self,x,y):
        """
        f(x,y) = x + e^-1x + y^3 = 0
        dy/df  = 3y^2
        """
        return 3*pow(y,2)

    def function14b2(self,x,y):
        """
        f(x,y) = x^2 + 2xy - y^2 + tan(x)
        """
        return pow(x,2) + 2*x*y - pow(y,2) + math.tan(x)

    def function14b2dx(self,x,y):
        """
        f(x,y) = x^2 + 2xy - y^2 + tan(x)
        dx/df  = 2x + 2y + sec^2(x) = 2x + 2y + 1 + tan(x)**2
        """
        return 2*x + 2*y + 1 + math.tan(x)**2

    def function14b2dy(self,x,y):
        """
        f(x,y) = x^2 + 2xy - y^2 + tan(x)
        dy/df  = 2x - 2y
        """
        return 2*x - 2*y

    def problem93_14(self):
        """ Solves computer problem 14 on page 93 """
        print '\nPage 93, Computer Problem 14'

        # GUESS FOR STARTING POSITION
        x = 0
        y = 0
        max = 10
        print 'Part a)'
        self.newton(x,y,max,self.function14a1,self.function14a1dx,self.function14a1dy, \
                        self.function14a2,self.function14a2dx,self.function14a2dy)
        # GUESS FOR STARTING POSITION
        x = 1
        y = 1
        max = 10
        print 'Part b)'
        self.newton(x,y,max,self.function14b1,self.function14b1dx,self.function14b1dy, \
                        self.function14b2,self.function14b2dx,self.function14b2dy)

########################################################
##  CONJUGATE GRADIENT METHOD TESTS
########################################################

    def testConjugateGradient(self,a,b,x,max):
        """ Tests the conjugate gradient """
        print 'Testing Conjugate Gradient'
        print 'Matrix A =',a
        print 'Vector b =',b
        print 'Vector x =',x
        print 'Answer should be x = [0.0909, 0.6364] (see below)'
        self.conjugateGradient(x,a,b,max)

########################################################
##  Page 245, Computer Problem 1
########################################################

    def problem245_1(self):
        """ Solves computer problem 1 on page 245 """
        print '\nPage 245, Computer Problem 1'
        
        # setup hilbert matrix
        # n_elements = 5 = true number of elements in A
        # n_range = n_elements+1 = number used for Python's range function
        # aij = (i + j - 1)^-1
        n_elements = 10
        n_range = n_elements+1
        a = [[(1.0/(row+col-1)) for col in range(1,n_range)] for row in range(1,n_range)]

        # Optional print Hilbert matrix
        if printDetails:
            print 'Hilbert',n_elements,'x',n_elements,'matrix A:'
            for i in range(0,n_elements):
                print a[i]

        # setup b-vector
        # bi = 1/3 (SUM j=1 to n of aij)
        b = [0 for row in range(1,n_range)]
        for i in range(1,n_range):
            sum = 0.0
            for j in range(1,n_range):
                sum += a[i-1][j-1]
            b[i-1] = 1.0/3.0*sum

        # Optional print b vector
        if printDetails:
            print 'b vector:\n',b

        # initial x-vector
        # x = [0]
        x = [0 for i in range(1,n_range)]
        if printDetails:
            print 'x vector:\n',x

        # Solve via Conjugate Gradient
        max = 50
        self.conjugateGradient(x,a,b,max)

########################################################
##  CUBIC SPLINE TESTS
########################################################

    def testCubicSpline(self):
        """ Tests the Cubic Spline functions """
        print '\nTesting Cubic Spline'

        # Setup Problem
        n = 10
        nMinus = n - 1
        t = [i/nMinus*2.25 for i in range(0,n+1)]
        y = [math.sqrt(t[i]) for i in range(0,n+1)]

        m = 37
        xlist = [i/(m-1)*2.25 for i in range(0,m+1)]
        fx = [math.sqrt(xlist[i]) for i in range(0,m+1)]
        if printDetails:
            print 't     =',t
            print 'y     =',y
            print 'xlist =',xlist
            print 'fx    =',fx

        # algorithm
        self.cubicSpline(n,t,y,xlist,fx,m)

########################################################
##  Page 365, Computer Problem 4
########################################################

    def problem365_4(self):
        """ Solves computer problem 4 on page 365 """
        print '\nPage 365, Computer Problem 4'

        # Setup Problem
        n = 10
        nMinus = n - 1
        t = [i/nMinus*2.25 for i in range(0,n+1)]
        y = [math.sqrt(t[i]) for i in range(0,n+1)]

        m = 37
        xlist = [i/(m-1)*2.25 for i in range(0,m+1)]
        fx = [math.sqrt(xlist[i]) for i in range(0,m+1)]
        if printDetails:
            print 't     =',t
            print 'y     =',y
            print 'xlist =',xlist
            print 'fx    =',fx

        # algorithm
        self.cubicSpline(n,t,y,xlist,fx,m)

########################################################
##  PQ factorization
########################################################

    def testPQ(self,n,a):
        """ Testing PQ decomposition """
        print '\nTesting PQ factorization algorithm'
        p,q = self.pqdecomposition(n,a)
        return p,q

    def pqdecomposition(self,n,a):
        p = [[0 for j in range(0,n)] for i in range(0,n)]
        q = [[0 for j in range(0,n)] for i in range(0,n)]

        for k in range(0,n):
            # 1's in the Q matrix - diagonal
            for j in range(0,n):
                if k + j == n-1:
                    q[k][j] = 1

            # P matrix first
            for j in range(0,k+1):
                total = sum([p[k][n-s-1]*q[n-s-1][j] for s in range(0,k+1)])
                p[k][n-1-j] = a[k][j] - total
                if printDetails: print 'SET p[',k,'][',n-1-j,'] =',p[k][n-1-j]

            # Q matrix second
            for j in range(n,k+1,-1):
                total = sum([p[k][n-s-1]*q[n-s-1][j-1] for s in range(0,k+1)])
                q[n-k-1][j-1] = (a[k][j-1] - total) / p[k][n-k-1]
                if printDetails: print 'SET q[',n-k-1,'][',j-1,'] =',q[n-k-1][j-1]


        print 'a =',a,'\np =',p,'\nq =',q
        return p,q

    def solvePQxb(self,p,q,b):
        """ Return x in the equation PQx = b """
        n = len(b)
        
        # first solve Pz = b
        z = self.solveAxb(p,b)

        # second solve Qx = z
        x = self.solveAxb(q,z)
        return x

    def solveAxb(self,a,b):
        """ Return x in the equation Ax = b """
        n = len(b)
        x = [0 for i in range(n)]
        for i in range(0,n):
            total = sum([a[i][j]*x[j] for j in range(n)])
            x[n-i-1] = (b[i] - total) / a[i][n-i-1]

        if printDetails:
            print 'Solving Ax = b for x'
            print 'A =',a
            print 'b =',b
            print 'x =',x
        return x

performTests = True
performProblems = True
printDetails = False

final = final()
final.printHeader()
if performTests: final.testFunctions()
if performProblems:
    final.problem93_14()
    final.problem245_1()
    final.problem365_4()
