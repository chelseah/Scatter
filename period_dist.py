import numpy as np
import random
from scipy import integrate, optimize
import matplotlib.pylab as plt

# distribution in TD12
a0=0.085
def f_a(a):
    if a<=0:
        return 0
    return 0.656*(a/a0)**(3.1)/(1.+(a/a0)**(3.6))/a
F_a = lambda a: integrate.quad(lambda x: f_a(x), 0, a)

def semi_major(P_cut_min, P_cut_max, Np, N_s):
    i=0
    P_12 = np.zeros(N_s)
    P_23 = np.zeros(N_s)
    while i < N_s : # 1000 systems
        np.random.seed()
        u = np.random.random_sample(Np)
        a_p = np.zeros([1,Np])[0]
        for j in range(Np):
            a_p[j] = optimize.fsolve(lambda a: F_a(a) - u[j], 0.03)
        ap = np.sort(a_p)
        P_12_aux = (ap[1]/ap[0])**(3./2.)
        P_23_aux = (ap[2]/ap[1])**(3./2.)
        if min(P_12_aux, P_23_aux) > P_cut_min and max(P_12_aux, P_23_aux) < P_cut_max:
            P_12[i] = P_12_aux
            P_23[i] = P_23_aux
            i = i+1
    if __name__ == "__main__":
        return ap, P_12, P_23
    else:
        return ap

def main():
    ap, P_12, P_23 = semi_major(P_cut_min=1.4, P_cut_max=5., Np=3, N_s=1e3)
    plt.figure()
    plt.hist(P_12, 20, color='white', edgecolor='black', label='$P_2/P_1$')
    plt.hist(P_23, 20, color='white', edgecolor='red', label='$P_3/P_2$')
    plt.xlabel('period ratio:$P_{i+1}/P_i$')
    plt.legend(loc='best')
    plt.show()

if __name__ == "__main__":
    main()
