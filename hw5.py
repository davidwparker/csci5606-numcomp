from __future__ import division

class interpolator:
    n = 0
    x = 0
    c = 0
    def dividedDifferences(self,n,x,c):
        self.n = n
        self.x = x
        self.c = c
        print 'Divided differences'
        # Python doesn't include end for range, so need to +1
        for j in range(1,self.n+1):
            for i in range(0,self.n-j+1):
                self.c[i][j] = (self.c[i+1][j-1] - self.c[i][j-1])/(self.x[i+j] - self.x[i])

                if printDetails:
                    print 'i+1,j-1    =', i+1,j-1
                    print 'x[i+j]-x[i]=', (self.x[i+j] - self.x[i])
                    print 'c[i+1][j-1]=', self.c[i+1][j-1]
                    print 'c[i][j-1]  =', self.c[i][j-1]
                    print 'c[i][j]    =', self.c[i][j]
    
    def printX(self):
        print self.x
    def printC(self):
        print self.c
    def printN(self):
        print self.n

    def printPx(self):
        stmt = ''
        for i in range(0,self.n+1):
            stmt += str(abs(self.c[0][i]))
            if (i > 0):
                for j in range(0,i):
                    stmt += '(x'
                    if (self.x[j] < 0):
                        stmt += ' + ' + str(abs(self.x[j]))
                    else:
                        stmt += ' - ' + str(self.x[j])
                    stmt += ')'
            if (i < self.n):
                if (self.c[0][i+1] < 0):
                    stmt += ' - '
                else:
                    stmt += ' + '
        print 'p(x) = ' + stmt

class homework5:
    def printHeader(self):
        print 'Principles of Numerical Computation'
        print 'David Parker'
        print 'Homework 5, October 30, 2011'

    def example62_2(self):
        print '\nProblem 6.2 example 2, page 331'
        inter = interpolator()
        n = 3
        x = [3,1,5,6]
        c = [[0 for col in range(n+1)] for row in range(n+1)]
        c[0][0] = 1
        c[1][0] = -3
        c[2][0] = 2
        c[3][0] = 4
        inter.dividedDifferences(n,x,c)
        inter.printPx()
        
    def problem61_1(self):
        print '\nProblem 6.1 #1'
        self.problem61_1a()
        self.problem61_1b()
        self.problem61_1c()

    def problem61_1a(self):
        print '\nPart a'
        inter = interpolator()
        n = 1
        x = [3,7]
        c = [[0 for col in range(n+1)] for row in range(n+1)]
        c[0][0] = 5
        c[1][0] = -1
        inter.dividedDifferences(n,x,c)
        inter.printPx()

    def problem61_1b(self):
        print '\nPart b'
        inter = interpolator()
        n = 2
        x = [7,1,2]
        c = [[0 for col in range(n+1)] for row in range(n+1)]
        c[0][0] = 146
        c[1][0] = 2
        c[2][0] = 1
        inter.dividedDifferences(n,x,c)
        inter.printPx()

    def problem61_1c(self):
        print '\nPart c'
        inter = interpolator()
        n = 3
        x = [3,7,1,2]
        c = [[0 for col in range(n+1)] for row in range(n+1)]
        c[0][0] = 10
        c[1][0] = 146
        c[2][0] = 2
        c[3][0] = 1
        inter.dividedDifferences(n,x,c)
        inter.printPx()

    def problem61_21(self):
        print '\nProblem 6.1 #21'
        inter = interpolator()
        n = 2
        x = [2,0,3]
        c = [[0 for col in range(n+1)] for row in range(n+1)]
        c[0][0] = 11
        c[1][0] = 7
        c[2][0] = 28
        inter.dividedDifferences(n,x,c)
        inter.printPx()
        
    def problem62_1(self):
        print '\nComputer Problem 6.2 #1'
        print 'For n = 5, 10, and 15, find the Newton interpolating polynomial pn for the function f(x) = 1/(1+x^2) on the interval [-5,5]. Use equally spaced nodes. In each case, compute f(x) - pn(x) for 30 equally spaced points in [-5,5] in order to see the divergence of pn from f.'
        low = -5
        high = 5
        print 'n = 5'
        self.problem62_1s(low,high,5)
        print 'n = 10'
        self.problem62_1s(low,high,10)
        print 'n = 15'
        self.problem62_1s(low,high,15)

    def problem62_1s(self,low,high,n):
        inter = interpolator()
        step = (high - low) / n
        x = []
        c = [[0 for col in range(n+1)] for row in range(n+1)]

        if printDetails: print 'step',step
        for i in range(0,n+1):
            x.append(low+i*step)
        if printDetails: print 'x',x

        for i in range(0,n+1):
            c[i][0] = (1/(pow(x[i],2)+1))
        if printDetails: print c

        inter.dividedDifferences(n,x,c)
        inter.printPx()

    def final_1(self):
        print '\nfinal problem 1'
        inter = interpolator()
        n = 1
        x = [-2,-1,0,1,2]
        c = [[0 for col in range(5)] for row in range(5)]
        c[0][0] = 0.5
        c[1][0] = 0.5
        c[2][0] = 2.0
        c[3][0] = 3.5
        c[4][0] = 3.5
        inter.dividedDifferences(n,x,c)
        inter.printPx()
        

printDetails = False
hw5 = homework5()
hw5.printHeader()
hw5.example62_2()
hw5.problem61_1()
hw5.problem61_21()
hw5.problem62_1()
hw5.final_1()
