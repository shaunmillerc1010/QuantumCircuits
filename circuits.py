from projectq.ops import All, CNOT, H, Measure, Rz, X, Z, Swap, C, Toffoli
from projectq import MainEngine
from projectq.meta import Dagger, Control, Compute, Uncompute
from projectq.backends import CircuitDrawer
from projectq.types import Qubit, Qureg, WeakQubitRef
from projectq.ops import BasicMathGate
from projectq.backends import ResourceCounter
from math import ceil, log2
import numpy as np

def get_binrep(intrep, register_width):
    a = np.binary_repr(intrep)
    a_list =[0 for i in range(register_width - len(a))] + [int(_) for _ in a]
    return a_list[::-1]
    #our binary rep is ascending unlike python's default

def get_intrep(binrep):
    binstring = ''
    for i in binrep: binstring = binstring+str(i)
    return int(binstring[::-1],2)
    #our binary rep is ascending unlike python's default

def _2scomplement(x = Qureg(), z = Qureg()):
    ''' computes x in 2's complement.
    Does not include the original sign bit.
    Use control gate with this function
    outputs x, s (nbit 2s complement of x), z (overflow bit)
    x original:  [0, 0, 0, 1, 0, 1, 0, 0]
    x:           [0, 0, 0, 1, 1, 0, 1, 1]'''
    "Enter n+1 0s for carry over and storing n bit register for 1 value"

    X      | z[0]
    All(X) | x
    ADD(z[:-1],x,z[-1])#should be able to get rid of zero bit completely
    #wise to continue to use z[-1] for this part.
    X      | z[0]

def SHIFT(s = Qureg(), x = Qureg()):
    '''shift x to the left by s bits
    https://arxiv.org/pdf/1807.02023.pdf'''
    j = len(s)
    for k in range(j):
        i = 0
        while i + 2**k <= 2**j - 1:
            C(Swap,1) | (s[k], x[i], x[i + (2**k)])
            i += 1

def RSHIFT(s = Qureg(), x = Qureg()):
    '''shift x to the right by s bits
    https://arxiv.org/pdf/1807.02023.pdf'''
    j = len(s)
    for k in range(j):
        i = 0
        while i + 2**k <= 2**j - 1:
            C(Swap,1) | (s[k], x[-1*i], x[-1*(i + (2**k))])
            i += 1

def ADD(x = Qureg(), y = Qureg(), z = Qureg()):
    #### [ x | y | 0 > --> [ x | x+y >
    n = len(x)
    for i in range(1,n): CNOT | (x[i],y[i])#1
    for i in range(n-1,0,-1):
        if i+1 == n : CNOT | (x[i], z)#2
        else: CNOT | (x[i], x[i+1])
    for i in range(n):
        if i == n-1: Toffoli | (y[i], x[i], z)#3
        else: Toffoli | (y[i], x[i], x[i+1])
    for i in range(n-1,0,-1): #4
        CNOT | (x[i],y[i])
        Toffoli | (y[i-1], x[i-1], x[i])
    for i in range(1, n-1): CNOT | (x[i], x[i+1])#5
    for i in range(n): CNOT | (x[i], y[i])#6
    return 0

def fixed_MUL(b = Qureg() ,a = Qureg(), z = Qureg(), eng = MainEngine()):
    '''see https://arxiv.org/pdf/0910.2530.pdf'''
    ''' unsigned multiplication is done shifted controlled addition.'''
        ##input\\ qubit registers a, b, z --- z twice the length of a
        ## suppose len(a) = n.
        ##output\\ qubit registers b, a, a*b
        #### [ b | a | z > --> [ b | a | a*b >
    n = len(b)
    for _ in range(n):
        with Control(eng, b[_]): ADD( a, z[_:n+_], z[n+_])
    return 0

def F8(f = Qureg(), p = Qureg(), x = Qureg()):
    '''circuit for finding the first one in the bit representation
    of x. The flag f is toggled to 0 as soon as the first 1 has been
    found. The position of the first one is store in the p-register,
    consisting of 3 bits in this example. For more details see
    https://arxiv.org/pdf/1807.02023.pdf'''
    X       |   f                   #1
    CNOT    |   (x[0],f)            #2
    C(X,2)  |   (f,x[1],p[0])       #3
    CNOT    |   (p[0],f)            #4
    C(X,2)  |   (f,x[2],p[1])       #5
    CNOT    |   (p[1],f)            #6
    C(X,2)  |   (f,x[3],p[0])       #7
    C(X,2)  |   (f,x[3],p[1])       #8
    C(X,2)  |   (p[0],p[1],f)       #9
    C(X,2)  |   (f,x[4],p[2])       #10
    CNOT    |   (p[2],f)            #11
    C(X,2)  |   (f,x[5],p[0])       #12
    C(X,2)  |   (f,x[5],p[2])       #13
    C(X,2)  |   (p[0],p[2],f)       #14
    C(X,2)  |   (f,x[6],p[1])       #15
    C(X,2)  |   (f,x[6],p[2])       #16
    C(X,2)  |   (p[1],p[2], f)      #17
    C(X,2)  |   (f,x[7],p[0])       #18
    C(X,2)  |   (f,x[7],p[1])       #19
    C(X,2)  |   (f,x[7],p[2])       #20
    C(X,3)  |   (p,f)               #21
#####################End 8 bits, start 9###########################

def Fn(f = Qureg(), p = Qureg(), x = Qureg()):
    n = len(x)
    X       |   f                   #1
    CNOT    |   (x[0],f)            #2
    C(X,2)  |   (f,x[1],p[0])       #3
    CNOT    |   (p[0],f)            #4
    for i in range(1, len(p)):
        if n < 2**(i+1): max = n
        else: max = 2**(i+1)
        for j in range(2**i, max):
            a = get_binrep(j, 2**(i))
            ones_place = []
            for k in range(len(a)):
                if a[k]==1:
                    C(X,2)  |   (f,x[j],p[k])
                    ones_place = ones_place + [k]
            C(X,len(ones_place))  | ([p[k] for k in ones_place],f)
