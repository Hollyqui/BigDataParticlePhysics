import numpy as np
f = 500
m_one = 10
m_two = 10
A = 356*24*60*60
M = m_one+m_two
S_O = 6*10**(-49)
f_s = 20
f_0 = 70
t = 100
psi = 100
x = f/f_0
nu = m_one*m_two/M**2
params = [psi, t, M, nu]
alpha = []
alpha.append(1)
alpha.append(0)

def get_h_f(f, params):
    psi, t, M, nu = params
    v = (np.pi*M*f)**(-1/3)
    sum_k = 0
    for k in range(len(alpha)):
        sum_k += alpha[k]*v**2
    phi_f = 2*np.pi*f*t-psi-np.pi/4+3/(128*nu*v**5)*sum_k
    h_f = A*f**(-7/6)*np.exp(1j * phi_f)
    # print(h_f)
    return h_f

def s_h(x):
    return (4.49*x)**-56+(0.16*x)**(-4.52)+0.52+0.32*x**2

def get_h_f_derivative(i):
    delta_x = params[i]/10
    new_params =  params.copy()
    new_params[i] = new_params[i]+params[i]/10
    delta_y = get_h_f(f, new_params) - get_h_f(f, params)
    return delta_y/delta_x
low = 1
high = 10

def create_fisher_matrix(increment):
    fisher_matrix = []
    for i in range(len(params)):
        fisher_matrix.append([])
        for j in range(len(params)):
            integral = 0
            for k in range(int(low*(1/increment)), int(high*(1/increment))):
                integral += np.real(get_h_f_derivative(i)*get_h_f_derivative(j)/s_h(k*increment))
            print(integral)
            fisher_matrix[-1].append(integral)
    return(fisher_matrix)
matrix = create_fisher_matrix(0.001)
matrix
