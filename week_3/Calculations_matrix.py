################### IMPORTS ####################################
import numpy as np
import matplotlib.pyplot as plt
import pylab
import math


################# DECLARE GLOBAL VARIABLES ###################

# convert 10**30 kg (normal black hole mass) to natural numbers
m_one = 10**30*9.109*10**(-31)
m_two = 10**30*9.109*10**(-31)

M = m_one+m_two
t = 100
psi = 100

nu = (m_one*m_two)/M**2
alpha = []
alpha.append(1)
alpha.append(0)
SNR = 10


# modify the variables in params to whichever ones you want to be taken into
# account for the fishermatrix (from the waveform) If a parameter is taken out
# the first line of the wavefunction has to be adjusted to not unpack that
# nonexistent parameter
params = [psi, t, M, nu]

# change this function to adjust the wavefunction to be something different
def phi_f(f, params):
    return phi_f1(f, params)

def phi_f1(f, params):
    psi, t, M, nu = params
    v = (np.pi*M*f)**(-1/3)
    sum_k = 0
    for k in range(len(alpha)):
        sum_k += alpha[k]*v**k
    return 2*np.pi*f*t-psi-np.pi/4+3/(128*nu*v**5)*sum_k

def phi_f2(f, params):
    pass # here you can add another wave function if wanted

###################### FUNCTIONS #######################

''' Finds A via intergation '''
#first find A (integral from low to high f^-7/6)/S(f)df
def get_A(freq_list, pds, SNR):
    integral = 0
    for f, s_h in zip(freq_list, pds):
        integral+=(f)**(-7/6)/s_h
    return SNR/np.sqrt(abs(integral))


''' Calculates h(f)'''
def get_h_f(A, f, params):
    phi = phi_f(f,params)
    h_f = A*f**(-7/6)*np.exp(1j * phi)
    return h_f

'''load Einstein data'''
def load_einstein():
    file = open('EinsteinPDS')
    string_read = file.read()
    arr = list(map(float, string_read.split()))
    x_arr = []
    y_arr1 = []
    y_arr2 = []
    y_arr3 = []
    for i in range(len(arr)):
        if i%4 == 0:
            x_arr.append(arr[i])
        if i%4 == 1:
            y_arr1.append(arr[i])
        if i%4 == 2:
            y_arr2.append(arr[i])
        if i%4 == 3:
            y_arr3.append(arr[i])
    plt.plot(x_arr, y_arr1)
    plt.plot(x_arr, y_arr2)
    plt.plot(x_arr, y_arr3)
    plt.yscale('log')
    plt.xscale('log')
    return x_arr, y_arr3

'''load Advanced Ligo data'''
def load_adv_ligo():
    file = open('LigoPDS')
    string_read = file.read()
    arr = list(map(float, string_read.split()))
    x_arr = []
    y_arr1 = []
    y_arr2 = []
    y_arr3 = []
    y_arr4 = []
    y_arr5 = []
    y_arr6 = []
    for i in range(len(arr)):
        if i%7 == 0:
            x_arr.append(arr[i])
        if i%7 == 1:
            y_arr1.append(arr[i])
        if i%7 == 2:
            y_arr2.append(arr[i])
        if i%7 == 3:
            y_arr3.append(arr[i])
        if i%7 == 4:
            y_arr4.append(arr[i])
        if i%7 == 5:
            y_arr5.append(arr[i])
        if i%7 == 6:
            y_arr6.append(arr[i])
    plt.plot(x_arr, y_arr1)
    plt.plot(x_arr, y_arr2)
    plt.plot(x_arr, y_arr3)
    plt.plot(x_arr, y_arr4)
    plt.plot(x_arr, y_arr5)
    plt.plot(x_arr, y_arr6)
    plt.yscale('log')
    plt.xscale('log')
    return x_arr, y_arr6


''' Calculates S(h)'''
def s_h_ligo(f):
    S_0 = 9*10**(-46)
    f_s = 40
    f_0 = 150
    x = f/f_0
    if(f>=f_s):
        r = S_0*(((4.49*x)**(-56))+0.16*(x**(-4.52))+0.52+(0.32*(x**2)))
    else:
        r = 10000000 # if infinity leads to numerical errors
    return(r)

def s_h_adv_ligo(f):
    f_0 = 215
    f_s = 20
    S_0 = 6*10**(-49)
    x = f/f_0
    if (f>=f_s):
        r = S_0*(x**(-4.14)-5*x**(-2) + ((111*(1-x**2 + (x**4)/2)/(1 + (x**2)/2))))
    else:
        r = 10000000 # if infinity leads to numerical errors
    return(r)

''' Finds derivative h'(f)'''
def get_h_f_derivative(A, f, i):
    delta_x = params[i]/10
    new_params =  params.copy()
    new_params[i] = new_params[i]+params[i]/10
    delta_y = get_h_f(A, f, new_params) - get_h_f(A, f, params)
    return delta_y/delta_x


''' Creates matrix using previous functions'''
def create_fisher_matrix(A, freq_list, pds):
    fisher_matrix = []
    for i in range(len(params)):
        fisher_matrix.append([])
        for j in range(len(params)):
            integral = 0
            for f, sh in zip(freq_list, pds):
                integral += np.real(get_h_f_derivative(A, f, i)*get_h_f_derivative(A, f, j)/sh)
            print(integral)
            fisher_matrix[-1].append(integral)
    return(fisher_matrix)




############################## MAIN #######################################

# Comment colour coding in the main part:
'''green = physics'''
#  grey = programming

# first, read in the pds of advanced ligo/the Einstein telescope
freq_list, pds = load_adv_ligo()
freq_list_einstein, pds_einstein = load_einstein()

''' As you can see here the PDS of the Einstein telescope is smaller across the
board, indicating that it is capable of detecting much smaller signals. It is
also notable that the x-axis doesn't have the same range, where ligo already
significantly worsens around f<10^2 the Einstein telescope keeps up until
f<10^1'''

#now, a is set via integration
A = get_A(freq_list, pds, SNR)
A_einstein = get_A(freq_list_einstein, pds, SNR)

print(A)
print(A_einstein)


matrix = create_fisher_matrix(A, freq_list, pds)

matrix_einstein = create_fisher_matrix(A, freq_list_einstein, pds_einstein)

# invert fisher_matrix
inverse = np.linalg.inv(matrix)
inverse_einstein = np.linalg.inv(matrix_einstein)

# extract diagonal and square root
diag = np.diagonal(inverse)
diag_einstein = np.diagonal(inverse_einstein)
print(diag)
print(diag_einstein)
std = np.sqrt(diag)
std_einstein = np.sqrt(diag_einstein)
print(std)
print(std_einstein)

'''As clearly visible, the standard deviation of the einstein telescope is almost
3 orders of magnitude smaller. This shows that there is a very clear improvement in
the measurement capabilities, not only in terms of what the smallest signal that
can be detected is but also in terms of how accurate those detections are'''


'''The following section plots Figure 1 of the paper for the initial ligo and
advanced ligo
based on the analytical formula (formula 3.7 and 3.8) given in the arxiv paper as a proof of concept.
This is to prove that our '''

def graph_s_h(telescope):
    graph_array = []
    try:
        if telescope == 'ligo':
            for i in range(0,1000):
                graph_array.append(np.sqrt(s_h_ligo(i)))
        elif telescope == 'adv_ligo':
            for i in range(0,1000):
                graph_array.append(np.sqrt(s_h_adv_ligo(i)))
        plt.figure()
        plt.plot(graph_array)
        plt.yscale('log')
        plt.xscale('log')
        axes = plt.gca()
        axes.set_ylim([10**(-24), 10**(-21)])
        axes.set_xlim([10**1, 10**3])
    except:
        print("no analytical solve for this telescope known")

graph_s_h('ligo')
graph_s_h('adv_ligo')
