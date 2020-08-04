from modsim import *
from sympy import *
x=symbols('x')
psi=Function('psi')
def Coeff_maker(num):
    '''
    This function creates a coefficient matrix
    recieves the number of terms
    returns the coeffient matrix
    '''
    d=symbols('d')
    da=type(d)
    A=np.empty(num,dtype=da)
    for i in range(num):
        A[i]=symbols(f'a_{i}')
    return A
def func_maker(num,var):
    '''
    the function recieves the number of terms startting from
    zeroth power of the variable var
    var should be a symbol object
    returns the final function and coefficient matrix
    '''
    A=Coeff_maker(num)
    expr=0
    for i in range(num):#making the function
        expr=expr+A[i]*var**i
    return expr,A  #returns the final polynomial and coefficeint matrix
def coeffsolver(f,co):
    '''
    solves the function and gets the value of coeffecient
    '''
    for i in range(len(co)):
        red=

func,coeff=func_maker(5,x)#function call
x=symbols('x')
f=Function('f')
expr=diff(func,x)-func
print(expr)
