import math
import matplotlib.pyplot as plt
import numpy as np
plt.rcParams.update({'figure.max_open_warning': 0})
filename = 'a.txt'#0INPme.txt'
f1 = open(filename, 'w')
filename = 'a.txt'#0INPem.txt'
f2 = open(filename, 'w')
filename = 'a.txt'#0DEPme.txt'
f3 = open(filename, 'w')
filename = 'a.txt'#0DEPme.txt'
f4 = open(filename, 'w')
filename = 'a.txt'#0CPme.txt'
f5 = open(filename, 'w')
filename = 'a.txt'#0CPem.txt'
f6 = open(filename, 'w')
def Unitary_matrix(energy, flavors, potential, baseline):    
    from sympy import sin, cos, symbols
    from numpy.linalg import eig
    from scipy.linalg import expm
    #variable setup
    flavors = 3
    A = potential
    #L, dmSq, EE, Theta = symbols('L dmSq EE Theta')
    L = baseline
    d_CP = 0.0#9000*math.pi/180#197.00*math.pi/180#232.00*math.pi/180#0.0#3.1415926535897932384#/2.0#0
    dmSq21 = 7.41e-5#0.001 #0.0547723
    dmSq31 = 2.511e-3
    dmSq32 = dmSq31 - dmSq21
    EE = energy
    Theta12 = 33.44*math.pi/180#0.58
    Theta13 = 8.57*math.pi/180#0.58
    Theta23 = 49.2*math.pi/180#0.58
    s12 = math.sin(Theta12)
    c12 = math.cos(Theta12)
    s13 = math.sin(Theta13)
    c13 = math.cos(Theta13)
    s23 = math.sin(Theta23)
    c23 = math.cos(Theta23)
    proj = [0.0] * flavors
    ve = [0.0] * flavors
    vm = [0.0] * flavors
    vt = [0.0] * flavors
    Ident = [[0 for x in range(flavors)] for y in range(flavors)]
    U_s = [0.0] * flavors
    U_f = 0.0

    t = 0.0
    di = [[0 for x in range(flavors)] for y in range(flavors)]
    u_dtest = [[0 for x in range(flavors)] for y in range(flavors)]
    u_test = [[0 for x in range(flavors)] for y in range(flavors)]


    if(energy==0.5):
        print(energy)
        print("c12",c12)
        print("Theta12",Theta12)
        print("Theta23",Theta23)
        print("Theta13",Theta13)
        #print(s23c13)


    #Alpha = (5.076)*(dmSq)/(2.0*EE)
    #Beta = Alpha*math.sin(Theta)*math.sin(Theta)
    #Gamma = Alpha*math.cos(Theta)*math.sin(Theta)
    #Eta = Alpha*math.cos(Theta)*math.cos(Theta)
    M23 = np.array([[1, 0,0],
                  [0, c23,s23],
                  [0, -s23, c23]])
    M13 = np.array([[c13, 0,s13*np.exp(complex(0,-d_CP))],
                  [0, 1,0],
                  [-s13*np.exp(complex(0,d_CP)), 0, c13]])
    M12 = np.array([[c12, s12,0],
                  [-s12, c12,0],
                  [0, 0, 1]])
    #matrix setup + eigensystem
  #  AA = np.array([[Beta + A, Gamma], 
   #                [Gamma, Eta]])
    PMNS = ((np.dot(M23,np.dot(M13,M12))))#array([[Beta, Gamma,0],
                  #[Gamma, Eta,0],
                 #[0, 0, 1]])
    #val,vec=eig(PMNS)
    #PMNS = AA
#    print('E-value:', val)
#   print('E-vector', vec)
    print('PMNS')
    print(PMNS)


    M_0 = np.array([[0,0,0],[0,dmSq21,0],[0,0,dmSq31]]) # matrix with mass squared diffs to use in hamiltonian
    UM_0UT = np.matmul(PMNS,np.matmul(M_0,np.transpose(PMNS)))
    #print(' ')
#(5.076)*
    AA = (5.076)*np.array([[A , 0, 0],
                    [0, 0, 0],
                    [0,0,0]])
    
    H_0 =(5.076)*( (1/(2*EE))*(UM_0UT))

    HT = H_0+AA

#    sz = complex(0,-L)
    
    t =  expm(-1j*HT*L)

    return(t)
def probEE(matrix, flavors):
    from sympy import sin, cos, symbols
    from numpy.linalg import eig
    U_matrix = matrix
    ve = [0.0] * flavors
    vm = [0.0] * flavors
    Tee = U_matrix[0][0]
    Tee_C = np.conjugate(Tee)
    Prob_ee = Tee_C*Tee    
    return(np.abs(Tee*Tee))

def probME(matrix, flavors):
    from sympy import sin, cos, symbols
    from numpy.linalg import eig
    U_matrix = matrix
    ve = [0.0] * flavors
    vm = [0.0] * flavors
    Tme = U_matrix[1][0]
    Tme_C = np.conjugate(Tme)
    Prob_me = Tme_C*Tme
    return(np.abs(Tme*Tme))


def probEM(matrix, flavors):
    from sympy import sin, cos, symbols
    from numpy.linalg import eig
    U_matrix = matrix
    ve = [0.0] * flavors
    vm = [0.0] * flavors
    Tem = U_matrix[0][1]
    Tem_C = np.conjugate(Tem)
    Prob_em = Tem_C*Tem
    return(np.abs(Tem*Tem))

def probMM(matrix, flavors):
    from sympy import sin, cos, symbols
    from numpy.linalg import eig
    U_matrix = matrix
    ve = [0.0] * flavors
    vm = [0.0] * flavors
    Tmm = U_matrix[1][1]
    Tmm_C = np.conjugate(Tmm)
    Prob_mm = Tmm_C*Tmm
    return(np.abs(Tmm*Tmm))


def probET(matrix, flavors):
    from sympy import sin, cos, symbols
    from numpy.linalg import eig
    U_matrix = matrix
    ve = [0.0] * flavors
    vt = [0.0] * flavors
    Tet = U_matrix[0][2]
    Tet_C = np.conjugate(Tet)
    Prob_et = Tet_C*Tet
    return(np.abs(Tet*Tet))#Prob_et.real)

def probMT(matrix, flavors):
    from sympy import sin, cos, symbols
    from numpy.linalg import eig
    U_matrix = matrix
    vt = [0.0] * flavors
    vm = [0.0] * flavors
    Tmt = U_matrix[1][2]
    Tmt_C = np.conjugate(Tmt)
    Prob_mt = Tmt_C*Tmt
    return(np.abs(Tmt*Tmt))#Prob_mt.real)

def probTE(matrix, flavors):
    from sympy import sin, cos, symbols
    from numpy.linalg import eig
    U_matrix = matrix
    ve = [0.0] * flavors
    vt = [0.0] * flavors
    Tte = U_matrix[2][0]
    Tte_C = np.conjugate(Tte)
    Prob_te = Tte_C*Tte
    return(np.abs(Tte*Tte))#Prob_te.real)

def probTM(matrix, flavors):
    from sympy import sin, cos, symbols
    from numpy.linalg import eig
    U_matrix = matrix
    vm = [0.0] * flavors
    vt = [0.0] * flavors
    Ttm = U_matrix[2][1]
    Ttm_C = np.conjugate(Ttm)
    Prob_tm = Ttm_C*Ttm
    return(np.abs(Ttm*Ttm))#Prob_tm.real)

def probTT(matrix, flavors):
    from sympy import sin, cos, symbols
    from numpy.linalg import eig
    U_matrix = matrix
    #ve = [0.0] * flavors
    #vt = [0.0] * flavors
    Ttt = U_matrix[2][2]
    Ttt_C = np.conjugate(Ttt)
    Prob_tt = Ttt_C*Ttt
    return(np.abs(Ttt*Ttt))#.real)

def U_multiply(matrix1, matrix2, flavors):
    from sympy import sin, cos, symbols
    from numpy.linalg import eig
    U1_matrix = matrix1
    U2_matrix = matrix2

   # print(np.dot(U2_matrix,U1_matrix))
    #print(np.dot(U1_matrix,U2_matrix))
    U_total = np.dot(U2_matrix,U1_matrix)
    return(U_total)
    
N = 1000

EE = [0.0] * N
Pem = [0.0] * N
Pee = [0.0] * N
Pmm = [0.0] * N

Ptt = [0.0] * N
Pet = [0.0] * N
Pmt = [0.0] * N

Ptm = [0.0] * N
Pte = [0.0] * N

Ptt_reversed = [0.0] * N
Ptm_reversed = [0.0] * N
Pmt_reversed = [0.0] * N
Pte_reversed = [0.0] * N
Pet_reversed = [0.0] * N

Pmt_vac = [0.0] * N
Ptm_vac = [0.0] * N
Pte_vac = [0.0] * N
Pet_vac = [0.0] * N
Ptt_vac = [0.0] * N
Pem_single_matter = [0.0] * N
Pem_vac = [0.0] * N
Pme_single_matter = [0.0] * N
Pme_vac = [0.0] * N
Pee_single_matter = [0.0] * N
Pee_vac = [0.0] * N
Pmm_single_matter = [0.0] * N

Pte_single_matter = [0.0] * N
Pet_single_matter = [0.0] * N
Ptm_single_matter = [0.0] * N
Pmt_single_matter = [0.0] * N
Ptt_single_matter = [0.0] * N
Pmm_vac = [0.0] * N
Pee_reversed = [0.0] * N
Pmm_reversed = [0.0] * N
Pem_reversed = [0.0] * N
Pme = [0.0] * N
Pme_reversed = [0.0] * N
Unitary1 = [0.0] * N
Unitary1b = [0.0] * N
LL1 = 800.0
Unitary2 = [0.0] * N
Unitary2b = [0.0] * N
LL2 = 2000-LL1
Unitary_total = [0.0] * N
Unitary_total_reversed = [0.0] * N

Unitary_single_matter = [0.0]*N
Unitary_vac = [0.0]*N
L = [0.0] * N
types = 3
V1 = 0.0001#0.0000107*3.0
V2 = 0.0001*3.0
max_energy = 5.0
Step_size = (max_energy-0.5)/(N)
Step_size_L = (2000.0-0.0)/(N)

for nu in range(N):
    EE[nu] = 0.5 + Step_size*nu
        
    Unitary1[nu] = Unitary_matrix(EE[nu], types, V1, LL1)
    Unitary2[nu] = Unitary_matrix(EE[nu], types, V2, LL2)

    Unitary_total[nu] = U_multiply(Unitary1[nu],Unitary2[nu],types)
    Pee[nu] = probEE(Unitary_total[nu], types)
    Pmm[nu] = probMM(Unitary_total[nu], types)
    Pme[nu] = probME(Unitary_total[nu], types)
    Pem[nu] = probEM(Unitary_total[nu], types)

    Ptt[nu] = probTT(Unitary_total[nu], types)
    Pmt[nu] = probMT(Unitary_total[nu], types)
    Pte[nu] = probTE(Unitary_total[nu], types)
    Ptm[nu] = probTM(Unitary_total[nu], types)
    Pet[nu] = probET(Unitary_total[nu], types)


    
    Unitary_total_reversed[nu] = U_multiply(Unitary2[nu],Unitary1[nu],types)
    Pee_reversed[nu] = probEE(Unitary_total_reversed[nu], types)
    Pmm_reversed[nu] = probMM(Unitary_total_reversed[nu], types)
    Pme_reversed[nu] = probME(Unitary_total_reversed[nu], types)
    Pem_reversed[nu] = probEM(Unitary_total_reversed[nu], types)

    Ptt_reversed[nu] = probTT(Unitary_total_reversed[nu], types)
    Pmt_reversed[nu] = probMT(Unitary_total_reversed[nu], types)
    Pte_reversed[nu] = probTE(Unitary_total_reversed[nu], types)
    Ptm_reversed[nu] = probTM(Unitary_total_reversed[nu], types)
    Pet_reversed[nu] = probET(Unitary_total_reversed[nu], types)
    
    Unitary_single_matter[nu] = Unitary_matrix(EE[nu], types, V1,2000)

    Pee_single_matter[nu] = probEE(Unitary_single_matter[nu], types)
    Pmm_single_matter[nu] = probMM(Unitary_single_matter[nu], types)
    Pme_single_matter[nu] = probME(Unitary_single_matter[nu], types)
    Pem_single_matter[nu] = probEM(Unitary_single_matter[nu], types)

    Pet_single_matter[nu] = probET(Unitary_single_matter[nu], types)
    Pmt_single_matter[nu] = probMT(Unitary_single_matter[nu], types)
    Pte_single_matter[nu] = probTE(Unitary_single_matter[nu], types)
    Ptm_single_matter[nu] = probTM(Unitary_single_matter[nu], types)
    Ptt_single_matter[nu] = probTT(Unitary_single_matter[nu], types)
    
    Unitary_vac[nu] = Unitary_matrix(EE[nu], types, 0.0, 2000)

    Pee_vac[nu] = probEE(Unitary_vac[nu], types)
    Pmm_vac[nu] = probMM(Unitary_vac[nu], types)
    Pme_vac[nu] = probME(Unitary_vac[nu], types)
    Pem_vac[nu] = probEM(Unitary_vac[nu], types)  

    Pet_vac[nu] = probET(Unitary_vac[nu], types)
    Pmt_vac[nu] = probMT(Unitary_vac[nu], types)
    Pte_vac[nu] = probTE(Unitary_vac[nu], types)
    Ptm_vac[nu] = probTM(Unitary_vac[nu], types)
    Ptt_vac[nu] = probTT(Unitary_vac[nu], types)
    print('vac prob check')                                     
    print(Pee_vac[nu]+Pem_vac[nu]+Pet_vac[nu])                          
    print(Pme_vac[nu]+Pmm_vac[nu]+Pmt_vac[nu])
    print(Pte_vac[nu]+Ptm_vac[nu]+Ptt_vac[nu])
    print(Pee_vac[nu]+Pme_vac[nu]+Pte_vac[nu])
    print(Pem_vac[nu]+Pmm_vac[nu]+Ptm_vac[nu])
    print(Pet_vac[nu]+Pmt_vac[nu]+Ptt_vac[nu])
    print(' ')
    print('matter prob check')
    print(Pee_single_matter[nu]+Pem_single_matter[nu]+Pet_single_matter[nu])
    print(Pme_single_matter[nu]+Pmm_single_matter[nu]+Pmt_single_matter[nu])
    print(Pte_single_matter[nu]+Ptm_single_matter[nu]+Ptt_single_matter[nu])
    print(Pee_single_matter[nu]+Pme_single_matter[nu]+Pte_single_matter[nu])
    print(Pem_single_matter[nu]+Pmm_single_matter[nu]+Ptm_single_matter[nu])
    print(Pet_single_matter[nu]+Pmt_single_matter[nu]+Ptt_single_matter[nu])
    print(' ')
    print('increasing matter prob check')
    print(Pee[nu]+Pem[nu]+Pet[nu])
    print(Pme[nu]+Pmm[nu]+Pmt[nu])
    print(Pte[nu]+Ptm[nu]+Ptt[nu])
    print(Pee[nu]+Pme[nu]+Pte[nu])
    print(Pem[nu]+Pmm[nu]+Ptm[nu])
    print(Pet[nu]+Pmt[nu]+Ptt[nu])
    print(' ')
    print('decreasing matter prob check')
    print(Pee_reversed[nu]+Pem_reversed[nu]+Pet_reversed[nu])
    print(Pme_reversed[nu]+Pmm_reversed[nu]+Pmt_reversed[nu])
    print(Pte_reversed[nu]+Ptm_reversed[nu]+Ptt_reversed[nu])
    print(Pee_reversed[nu]+Pme_reversed[nu]+Pte_reversed[nu])
    print(Pem_reversed[nu]+Pmm_reversed[nu]+Ptm_reversed[nu])
    print(Pet_reversed[nu]+Pmt_reversed[nu]+Ptt_reversed[nu])
    print(' ')
"""
fig1 =plt.figure(1)
plt.plot(EE, Pme_vac, color='m',linestyle='dashed', label='P_me vacuum ')
plt.plot(EE, Pee_vac, color='y', label='P_ee vacuum')
plt.plot(EE, Pte_vac, color='b', linestyle='dotted',label='P_te vacuum')
plt.title("Vacuum prob nu_electron ",fontsize=10)
plt.xlabel('                                  Energy (GeV)')
plt.legend(fontsize=10)
fig2 =plt.figure(2)
plt.plot(EE, Pem_vac, color='m',linestyle='dashed', label='P_em vacuum')
plt.plot(EE, Pmm_vac, color='y', label='P_mm vacuum')
plt.plot(EE, Ptm_vac, color='b', linestyle='dotted',label='P_tm vacuum')
plt.title("Vacuum prob nu_muon ",fontsize=10)
plt.xlabel('                                  Energy (GeV)')
plt.legend(fontsize=10)
fig3 =plt.figure(3)
plt.plot(EE, Pet_vac, color='m',linestyle='dashed', label='P_et vacuum')
plt.plot(EE, Ptt_vac, color='y', label='P_tt vacuum')
plt.plot(EE, Pmt_vac, color='b', linestyle='dotted',label='P_mt vacuum')
plt.title("Vacuum prob nu_tau ",fontsize=10)
plt.xlabel('                                  Energy (GeV)')
plt.legend(fontsize=10)

fig4 =plt.figure(4)
plt.plot(EE, Pme_single_matter, color='m',linestyle='dashed', label='P_me single_matter ')
plt.plot(EE, Pee_single_matter, color='y', label='P_ee single_matter')
plt.plot(EE, Pte_single_matter, color='b', linestyle='dotted',label='P_te single_matter')
plt.title("single_matter prob nu_electron ",fontsize=10)
plt.xlabel('                                  Energy (GeV)')
plt.legend(fontsize=10)
fig5 =plt.figure(5)
plt.plot(EE, Pem_single_matter, color='m',linestyle='dashed', label='P_em single_matter ')
plt.plot(EE, Pmm_single_matter, color='y', label='P_mm single_matter')
plt.plot(EE, Ptm_single_matter, color='b', linestyle='dotted',label='P_tm single_matter')
plt.title("single_matter prob nu_muon ",fontsize=10)
plt.xlabel('                                  Energy (GeV)')
plt.legend(fontsize=10)
fig6 =plt.figure(6)
plt.plot(EE, Pet_single_matter, color='m',linestyle='dashed', label='P_et single_matter ')
plt.plot(EE, Ptt_single_matter, color='y', label='P_tt single_matter')
plt.plot(EE, Pmt_single_matter, color='b', linestyle='dotted',label='P_mt single_matter')
plt.title("single_matter prob nu_tau ",fontsize=10)
plt.xlabel('                                  Energy (GeV)')
plt.legend(fontsize=10)




fig17 =plt.figure(17)
plt.plot(EE, Pee, color='g', label='P_ee constant increase piecewise matter ')
plt.plot(EE, Pme, color='y', linestyle='dashed',label='P_me constant increase piecewise matter ')
plt.plot(EE, Pte, color='b', linestyle='dotted',label='P_te constant increase piecewise matter ')
plt.title("nu_electron constant increase piecewise matter prob",fontsize=10)
plt.xlabel('                                  Energy (GeV)')
plt.legend(fontsize=10)


fig18 =plt.figure(18)
plt.plot(EE, Pmm, color='g', label='P_mm constant increase piecewise matter ')
plt.plot(EE, Pem, color='y', linestyle='dashed',label='P_em constant increase piecewise matter ')
plt.plot(EE, Ptm, color='b', linestyle='dotted',label='P_tm constant increase piecewise matter ')
plt.title("nu_muon constant increase piecewise matter prob",fontsize=10)
plt.xlabel('                                  Energy (GeV)')
plt.legend(fontsize=10)

fig19 =plt.figure(19)
plt.plot(EE, Ptt, color='g', label='P_tt constant increase piecewise matter ')
plt.plot(EE, Pet, color='y', linestyle='dashed',label='P_et constant increase piecewise matter ')
plt.plot(EE, Pmt, color='b', linestyle='dotted',label='P_mt constant increase piecewise matter ')
plt.title("nu_tau constant increase piecewise matter prob",fontsize=10)
plt.xlabel('                                  Energy (GeV)')
plt.legend(fontsize=10)


fig20 =plt.figure(20)
plt.plot(EE, Pee_reversed, color='g', label='P_ee constant decrease piecewise matter ')
plt.plot(EE, Pme_reversed, color='y', linestyle='dashed',label='P_me constant decrease piecewise matter ')
plt.plot(EE, Pte_reversed, color='b', linestyle='dotted',label='P_te constant decrease piecewise matter ')
plt.title("nu_electron constant decrease piecewise matter prob",fontsize=10)
plt.xlabel('                                  Energy (GeV)')
plt.legend(fontsize=10)

fig21 =plt.figure(21)
plt.plot(EE, Pmm_reversed, color='g', label='P_mm constant decrease piecewise matter ')
plt.plot(EE, Pem_reversed, color='y', linestyle='dashed',label='P_em constant decrease piecewise matter ')
plt.plot(EE, Pmt_reversed, color='b', linestyle='dotted',label='P_tm constant decrease piecewise matter ')
plt.title("nu_muon constant decrease piecewise matter prob",fontsize=10)
plt.xlabel('                                  Energy (GeV)')
plt.legend(fontsize=10)

fig22 =plt.figure(22)
plt.plot(EE, Ptt_reversed, color='g', label='P_tt constant decrease piecewise matter ')
plt.plot(EE, Pmt_reversed, color='y', linestyle='dashed',label='P_mt constant decrease piecewise matter ')
plt.plot(EE, Pet_reversed, color='b', linestyle='dotted',label='P_et constant decrease piecewise matter ')
plt.title("nu_tau constant decrease piecewise matter prob",fontsize=10)
plt.xlabel('                                  Energy (GeV)')
plt.legend(fontsize=10)





####################################################################################


fig45 =plt.figure(45)
plt.plot(EE, Pme_vac, color='b', label='P_me vacuum ')
plt.plot(EE, Pem_vac, color='y',linestyle='dashed', label='P_em vacuum')
plt.title("CP= 0: Vacuum prob: nu_e and nu_mu",fontsize=10)
plt.xlabel('                                  Energy (GeV)')
plt.legend(fontsize=10)


fig46 =plt.figure(46)
plt.plot(EE, Pmt_vac, color='b', label='P_mt vacuum ')
plt.plot(EE, Ptm_vac, color='y',linestyle='dashed', label='P_tm vacuum')
plt.title("Vacuum prob: nu_tau and nu_mu ",fontsize=10)
plt.xlabel('                                  Energy (GeV)')
plt.legend(fontsize=10)


fig47 =plt.figure(47)
plt.plot(EE, Pet_vac, color='b', label='P_et vacuum ')
plt.plot(EE, Pte_vac, color='y',linestyle='dashed', label='P_te vacuum')
plt.title("Vacuum prob: nu_tau and nu_e ",fontsize=10)
plt.xlabel('                                  Energy (GeV)')
plt.legend(fontsize=10)

fig48 =plt.figure(48)
plt.plot(EE, Pme_single_matter, color='b', label='P_me single matter ')
plt.plot(EE, Pem_single_matter, color='y',linestyle='dashed', label='P_em single matter')
plt.title("CP=90: single matter prob: nu_e and nu_mu",fontsize=10)
plt.xlabel('                                  Energy (GeV)')
plt.legend(fontsize=10)


fig49 =plt.figure(49)
plt.plot(EE, Pmt_single_matter, color='b', label='P_mt single matter ')
plt.plot(EE, Ptm_single_matter, color='y',linestyle='dashed', label='P_tm single matter')
plt.title("single matter prob: nu_tau and nu_mu ",fontsize=10)
plt.xlabel('                                  Energy (GeV)')
plt.legend(fontsize=10)


fig50 =plt.figure(50)
plt.plot(EE, Pet_single_matter, color='b', label='P_et single matter')
plt.plot(EE, Pte_single_matter, color='y',linestyle='dashed', label='P_te single matter')
plt.title("single matter prob: nu_tau and nu_e ",fontsize=10)
plt.xlabel('                                  Energy (GeV)')
plt.legend(fontsize=10)

fig51 =plt.figure(51)
plt.plot(EE, Pme_reversed, color='b', label='P_me decreasing constant matter ')
plt.plot(EE, Pem_reversed, color='y',linestyle='dashed', label='P_em decreasing constant matter')
plt.title("CP=0: decreasing constant matter prob: nu_e and nu_mu",fontsize=10)
plt.xlabel('                                  Energy (GeV)')
plt.legend(fontsize=10)


fig52 =plt.figure(52)
plt.plot(EE, Pmt_reversed, color='b', label='P_mt decreasing constant matter ')
plt.plot(EE, Ptm_reversed, color='y',linestyle='dashed', label='P_tm decreasing constant matter')
plt.title("decreasing constant matter prob: nu_tau and nu_mu ",fontsize=10)
plt.xlabel('                                  Energy (GeV)')
plt.legend(fontsize=10)


fig53 =plt.figure(53)
plt.plot(EE, Pet_reversed, color='b', label='P_et decreasing constant matter ')
plt.plot(EE, Pte_reversed, color='y',linestyle='dashed', label='P_te decreasing constant matter')
plt.title("decreasing constant matter prob: nu_tau and nu_e ",fontsize=10)
plt.xlabel('                                  Energy (GeV)')
plt.legend(fontsize=10)


fig54 =plt.figure(54)
plt.plot(EE, Pme, color='b', label='P_me increasing constant matter ')
plt.plot(EE, Pem, color='y',linestyle='dashed', label='P_em increasing constant matter')
plt.title("CP=0:increasing constant matter prob: nu_e and nu_mu",fontsize=10)
plt.xlabel('                                  Energy (GeV)')
plt.legend(fontsize=10)

fig55 =plt.figure(55)
plt.plot(EE, Pmt, color='b', label='P_mt increasing constant matter ')
plt.plot(EE, Ptm, color='y',linestyle='dashed', label='P_tm increasing constant matter')
plt.title("increasing constant matter prob: nu_tau and nu_mu ",fontsize=10)
plt.xlabel('                                  Energy (GeV)')
plt.legend(fontsize=10)


fig56 =plt.figure(56)
plt.plot(EE, Pet, color='b', label='P_et increasing constant matter ')
plt.plot(EE, Pte, color='y',linestyle='dashed', label='P_te increasing constant matter')
plt.title("increasing constant matter prob: nu_tau and nu_e ",fontsize=10)
plt.xlabel('                                  Energy (GeV)')
plt.legend(fontsize=10)

fig60 =plt.figure(60)
plt.plot(EE, Pme, color='b', label='P_me increasing constant matter ')
plt.plot(EE, Pme_reversed, color='r',linestyle='dashed', label='P_me decreasing constant matter ')
plt.plot(EE, Pem, color='g', linestyle='dotted',label='P_em increasing constant matter ')
plt.plot(EE, Pem_reversed, color='y',linestyle='dashdot', label='P_em decreasing constant matter ')
plt.title("piecewise constant matter prob: nu_e and nu_mu ",fontsize=10)
plt.xlabel('                                  Energy (GeV)')
plt.legend(fontsize=10)

fig61 =plt.figure(61)
plt.plot(EE, Aem, color='b', label='Asym. factor: increasing ')
plt.plot(EE, Aem_r, color='m', linestyle='dashed',label='Asym. factor: decreasing')
plt.plot(EE, Ame_proper, color='g', label='Asym. factor: P_nu_mu->nu_e increasing first')
plt.plot(EE, Aem_proper, color='y', linestyle='dashed',label='Asym. factor: P_nu_mu->nu_e decreasing first')
plt.title("Asym. factor ",fontsize=10)
plt.xlabel('                                  Energy (GeV)')
plt.legend(fontsize=10)


plt.show()

fig1.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23b/3_flavor_oscillations_10^-4V_2steps_CP90_Method2/Vacuum_NuElectron.jpg')
fig2.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23b/3_flavor_oscillations_10^-4V_2steps_CP90_Method2/Vacuum_NuMuon.jpg')
fig3.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23b/3_flavor_oscillations_10^-4V_2steps_CP90_Method2/Vacuum_NuTau.jpg')
fig4.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23b/3_flavor_oscillations_10^-4V_2steps_CP90_Method2/Single_matter_NuElectron.jpg')
fig5.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23b/3_flavor_oscillations_10^-4V_2steps_CP90_Method2/Single_matter_NuMuon.jpg')
fig6.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23b/3_flavor_oscillations_10^-4V_2steps_CP90_Method2/Single_matter_NuTau.jpg')

fig7.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23b/3_flavor_oscillations_10^-4V_2steps_CP90_Method2/diff_single_Vac_prob_nuElectron.jpg')
fig8.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23b/3_flavor_oscillations_10^-4V_2steps_CP90_Method2/diff_single_Vac_prob_nuMuon.jpg')
fig9.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23b/3_flavor_oscillations_10^-4V_2steps_CP90_Method2/diff_single_Vac_prob_nuTau.jpg')

fig10.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23b/3_flavor_oscillations_10^-4V_2steps_CP90_Method2/Vac_prob_appear.jpg')
fig11.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23b/3_flavor_oscillations_10^-4V_2steps_CP90_Method2/Vac_prob_disappear.jpg')
fig12.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23b/3_flavor_oscillations_10^-4V_2steps_CP90_Method2/Single_matter_prob_appear.jpg')
fig13.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23b/3_flavor_oscillations_10^-4V_2steps_CP90_Method2/Single_matter_prob_disappear.jpg')

fig14.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23b/3_flavor_oscillations_10^-4V_2steps_CP90_Method2/diff_single_Vac_e.jpg')
fig15.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23b/3_flavor_oscillations_10^-4V_2steps_CP90_Method2/diff_single_Vac_m.jpg')
fig16.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23b/3_flavor_oscillations_10^-4V_2steps_CP90_Method2/diff_single_Vac_t.jpg')

fig17.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23b/3_flavor_oscillations_10^-4V_2steps_CP90_Method2/nuE_constant_increase_piecewise_matter_prob.jpg')
fig18.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23b/3_flavor_oscillations_10^-4V_2steps_CP90_Method2/nuM_constant_increase_piecewise_matter_prob.jpg')
fig19.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23b/3_flavor_oscillations_10^-4V_2steps_CP90_Method2/nuT_constant_increase_piecewise_matter_prob.jpg')

fig20.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23b/3_flavor_oscillations_10^-4V_2steps_CP90_Method2/nuE_constant_decrease_piecewise_matter_prob.jpg')
fig21.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23b/3_flavor_oscillations_10^-4V_2steps_CP90_Method2/nuM_constant_decrease_piecewise_matter_prob.jpg')
fig22.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23b/3_flavor_oscillations_10^-4V_2steps_CP90_Method2/nuT_constant_decrease_piecewise_matter_prob.jpg')

fig23.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23b/3_flavor_oscillations_10^-4V_2steps_CP90_Method2/diff_increase_piecewise_Vac_prob_nuE.jpg')
fig24.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23b/3_flavor_oscillations_10^-4V_2steps_CP90_Method2/diff_increase_piecewise_Vac_prob_nuM.jpg')
fig25.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23b/3_flavor_oscillations_10^-4V_2steps_CP90_Method2/diff_increase_piecewise_Vac_prob_nuT.jpg')

fig26.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23b/3_flavor_oscillations_10^-4V_2steps_CP90_Method2/diff_decrease_piecewise_Vac_prob_nuE.jpg')
fig27.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23b/3_flavor_oscillations_10^-4V_2steps_CP90_Method2/diff_decrease_piecewise_Vac_prob_nuM.jpg')
fig28.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23b/3_flavor_oscillations_10^-4V_2steps_CP90_Method2/diff_decrease_piecewise_Vac_prob_nuT.jpg')

fig29.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23b/3_flavor_oscillations_10^-4V_2steps_CP90_Method2/increase_piecewise_prob_appear.jpg')
fig30.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23b/3_flavor_oscillations_10^-4V_2steps_CP90_Method2/increase_piecewise_prob_disappear.jpg')
fig31.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23b/3_flavor_oscillations_10^-4V_2steps_CP90_Method2/decrease_piecewise_prob_appear.jpg')
fig32.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23b/3_flavor_oscillations_10^-4V_2steps_CP90_Method2/decrease_piecewise_prob_disappear.jpg')

fig33.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23b/3_flavor_oscillations_10^-4V_2steps_CP90_Method2/diff_decrease_piecewise_Vac_e.jpg')
fig34.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23b/3_flavor_oscillations_10^-4V_2steps_CP90_Method2/diff_decrease_piecewise_Vac_m.jpg')
fig35.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23b/3_flavor_oscillations_10^-4V_2steps_CP90_Method2/diff_decrease_piecewise_Vac_t.jpg')

fig36.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23b/3_flavor_oscillations_10^-4V_2steps_CP90_Method2/diff_increase_piecewise_Vac_e.jpg')
fig37.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23b/3_flavor_oscillations_10^-4V_2steps_CP90_Method2/diff_increase_piecewise_Vac_m.jpg')
fig38.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23b/3_flavor_oscillations_10^-4V_2steps_CP90_Method2/diff_increase_piecewise_Vac_t.jpg')

fig39.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23b/3_flavor_oscillations_10^-4V_2steps_CP90_Method2/diff_piecewise_appearance.jpg')
fig40.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23b/3_flavor_oscillations_10^-4V_2steps_CP90_Method2/diff_piecewise_disappearance.jpg')

fig41.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23b/3_flavor_oscillations_10^-4V_2steps_CP90_Method2/diff_piecewise_increase.jpg')
fig42.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23b/3_flavor_oscillations_10^-4V_2steps_CP90_Method2/diff_piecewise_decrease.jpg')

fig43.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23b/3_flavor_oscillations_10^-4V_2steps_CP90_Method2/diff_single.jpg')
fig44.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23b/3_flavor_oscillations_10^-4V_2steps_CP90_Method2/diff_vac.jpg')

fig45.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23b/3_flavor_oscillations_10^-4V_2steps_CP90_Method2/Vacuum_compar_nuE_muMU.jpg')
fig46.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23b/3_flavor_oscillations_10^-4V_2steps_CP90_Method2/Vacuum_compar_nuMU_muTAU.jpg')
fig47.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23b/3_flavor_oscillations_10^-4V_2steps_CP90_Method2/Vacuum_compar_nuE_muTAU.jpg')

fig48.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23b/3_flavor_oscillations_10^-4V_2steps_CP90_Method2/Single_matter_compar_nuE_muMU.jpg')
fig49.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23b/3_flavor_oscillations_10^-4V_2steps_CP90_Method2/Single_matter_compar_nuMU_muTAU.jpg')
fig50.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23b/3_flavor_oscillations_10^-4V_2steps_CP90_Method2/Single_matter_compar_nuE_muTAU.jpg')


fig51.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23b/3_flavor_oscillations_10^-4V_2steps_CP90_Method2/Decreasing_piecewise_matter_compar_nuE_muMU.jpg')
fig52.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23b/3_flavor_oscillations_10^-4V_2steps_CP90_Method2/Decreasing_piecewise_matter_compar_nuMU_muTAU.jpg')
fig53.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23b/3_flavor_oscillations_10^-4V_2steps_CP90_Method2/Decreasing_piecewise_matter_compar_nuE_muTAU.jpg')


fig54.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23b/3_flavor_oscillations_10^-4V_2steps_CP90_Method2/Increasing_piecewise_matter_compar_nuE_muMU.jpg')
fig55.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23b/3_flavor_oscillations_10^-4V_2steps_CP90_Method2/Increasing_piecewise_matter_compar_nuMU_muTAU.jpg')
fig56.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23b/3_flavor_oscillations_10^-4V_2steps_CP90_Method2/Increasing_piecewise_matter_compar_nuE_muTAU.jpg')

fig100.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23b/3_flavor_oscillations_10^-4V_2steps_CP90_Method2/Pem_Pme_propVSimprop.jpg')
fig110.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23b/3_flavor_oscillations_10^-4V_2steps_CP90_Method2/Ptm_Pmt_propVSimprop.jpg')
fig120.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23b/3_flavor_oscillations_10^-4V_2steps_CP90_Method2/Pet_Pte_propVSimprop.jpg')
fig130.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23b/3_flavor_oscillations_10^-4V_2steps_CP90_Method2/Pem_Pme_propANDimprop.jpg')

fig140.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23b/3_flavor_oscillations_10^-4V_2steps_CP90_Method2/Ptm_Pmt_propANDimprop.jpg')
fig150.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23b/3_flavor_oscillations_10^-4V_2steps_CP90_Method2/Pet_Pte_propANDimprop.jpg')

fig60.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23b/3_flavor_oscillations_10^-4V_2steps_CP90_Method2/Pem_Pem_Im_P_COMPARE.jpg')



"""
