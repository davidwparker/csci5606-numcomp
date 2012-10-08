from __future__ import division
#import math

def f(x):
    return math.atan(x)

s = math.sqrt(2)
h = 1.0
M = 26

for k in range (0,M+1):
    F2 = f(s+h)
    F1 = f(s-h)
    d = F2-F1
    r = d/(2*h)
    print k,h,F2,F1,d,r
    h /= 2
