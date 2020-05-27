from projectq.ops import All, CNOT, H, Measure, Rz, X, Z, Swap, Toffoli, C
from projectq import MainEngine
from projectq.meta import Dagger, Control, Compute, Uncompute
from projectq.backends import CircuitDrawer, ResourceCounter
from projectq.types import Qubit, Qureg, WeakQubitRef
from projectq.ops import BasicMathGate
from projectq.backends import ResourceCounter
from circuits import *
from math import log2, ceil, floor
from circuits import _2scomplement

def list_to_reg(state_list,eng):
    starting_state = state_list
    nqubits = len(starting_state)
    reg = eng.allocate_qureg(nqubits)
    for _ in range(nqubits):
        if starting_state[_] == 1: X | reg[_]
    return reg

def get_int(a = Qureg()):
    f = [int(_) for _ in a]
    res = int("".join(str(x) for x in f[::-1]), 2)#https://www.geeksforgeeks.org/python-binary-list-to-integer/
    return res

def ADD_resources(number_of_bits = 16):
    print('getting costs of ' + str(number_of_bits) + ' bit addition...')
    #___________________initial__________________
    resource_counter = ResourceCounter()
    eng = MainEngine(resource_counter)
    a = [0 for i in range(number_of_bits)]
    b = [0 for i in range(number_of_bits)]
    a = list_to_reg( a, eng) #one extra bit (in most significant) needed for addition
    z = list_to_reg([0],eng)
    b = list_to_reg(b, eng)
    #___________________circuit___________________
    ADD(a,b,z)
    #___________________measure___________________
    eng.flush()
    print(resource_counter)

def ADD_test(x = [0,1,0,1,1,1,1,1], y = [1,0,1,1,0,0,1,1], z = [0]):
	print('Addition of two 8 bit numbers' )
	eng = MainEngine()
	#<least ---significance---most>
	print('x: ',x)
	print('y: ',y,' (205 in binary)')
	x = list_to_reg(x, eng) #one extra bit needed for addition
	y = list_to_reg(y,eng)
	z = list_to_reg(z, eng)
	ADD(x,y,z)
	All(Measure) | x+z+y
	eng.flush()
	print('x value: ',get_int(x))
	print('(sum) value: ',get_int(y+z))
	print('')

def MUL_test():
    print('testing multiplication circuit...')
    #___________________initial__________________
    eng = MainEngine()
    #<least ---significance---most>
    a = list_to_reg([0,1,1],eng) #7
    n = len(a)
    b = list_to_reg([1,1,1],eng)#7
    z = list_to_reg([0 for _ in range(2*n)],eng)
    fixed_MUL(a,b,z,eng = eng)
    All(Measure) | a+b+z
    eng.flush()
    print('a value: ',get_int(a))
    print('prod ',get_int(z))
    print([int(i) for i in a+b+z])
    print('')


def MUL_resources(number_of_bits = 16):
    print('getting costs of ' + str(number_of_bits) + ' unsigned multiplication...')
    #___________________initial__________________
    resource_counter = ResourceCounter()
    eng = MainEngine(resource_counter)
    a = [0 for i in range(number_of_bits)]
    b = [0 for i in range(number_of_bits)]
    n = len(a)
    z = list_to_reg([0 for _ in range(2*n)],eng)
    a = list_to_reg( a, eng) #one extra bit (in most significant) needed for addition
    b = list_to_reg(b, eng)
    fixed_MUL(a,b,z, eng = eng)
    eng.flush()
    print(resource_counter)

def Fn_resources(n):
    #___________________initial__________________
    resource_counter = ResourceCounter()
    eng = MainEngine(resource_counter)
    f = list_to_reg([0],eng)
    p = list_to_reg([0 for i in range(ceil(log2(n)))],eng)
    x = list_to_reg([0 for i in range(n)],eng)
    print('Getting resources for ' + str(len(x)) + 'bit first ones circuit.')
    Fn(f,p,x)
    eng.flush()
    print(resource_counter)

def RSHIFT_resources(n):
    #___________________initial__________________
    resource_counter = ResourceCounter()
    eng = MainEngine(resource_counter)
    p = list_to_reg([0 for i in range(floor(log2(n)))],eng)
    x = list_to_reg([0 for i in range(n)],eng)
    print('Getting resources for ' + str(len(x)) + ' right shift circuit.')
    RSHIFT(p,x)
    eng.flush()
    print(resource_counter)

def CRightShift_resources(k):
    #___________________initial__________________
    resource_counter = ResourceCounter()
    eng = MainEngine(resource_counter)
    control = list_to_reg([0],eng)
    m = k #case when mantissa is 13 bits
    x = list_to_reg([0 for i in range(m)],eng)
    x = x[::-1]#for right swap
    print('Getting resources for ' + str(len(x)) + ' controlled shift by 1 circuit.')
    for i in range(0,m-1): C(Swap,1)|(control,x[i],x[i+1] )
    Measure | control
    Measure | x #for tex purposes
    eng.flush()
    print(resource_counter)

def _2sresources(number_of_bits = 16):
    print('getting costs of ' + str(number_of_bits) + ' bit 2s complement...')
    #___________________initial__________________
    resource_counter = ResourceCounter()
    eng = MainEngine(resource_counter)
    a = [0 for i in range(number_of_bits)]
    a = list_to_reg( a, eng) #one extra bit (in most significant) needed for addition
    z = list_to_reg([0 for _ in a]+[0],eng)
    _2scomplement(a,z)
    All(Measure) | a+z
    eng.flush()
    print(resource_counter)

if __name__ == '__main__':
    ADD_resources()
    ADD_test()

