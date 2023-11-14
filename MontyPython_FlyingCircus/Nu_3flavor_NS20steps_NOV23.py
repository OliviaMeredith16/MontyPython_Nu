import math
import matplotlib.pyplot as plt
import numpy as np
plt.rcParams.update({'figure.max_open_warning': 0})


def Unitary_matrix(energy, flavors, potential, baseline):    
    from sympy import sin, cos, symbols
    from numpy.linalg import eig
    from scipy.linalg import expm
    #variable setup
    flavors = 3
    A = potential
    #L, dmSq, EE, Theta = symbols('L dmSq EE Theta')
    L = baseline
    d_CP = 0.0#90.00*math.pi/180#197.00*math.pi/180#232.00*math.pi/180#0.0#3.1415926535897932384#/2.0#0
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

    M23 = np.array([[1, 0,0],
                  [0, c23,s23],
                  [0, -s23, c23]])
    M13 = np.array([[c13, 0,s13*np.exp(complex(0.0,-d_CP))],
                  [0, 1,0],
                  [-s13*np.exp(complex(0.0,d_CP)), 0, c13]])
    M12 = np.array([[c12, s12,0],
                  [-s12, c12,0],
                  [0, 0, 1]])
    PMNS = ((np.dot(M23,np.dot(M13,M12))))
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

def U_multiply(matrix1, matrix2, matrix3, matrix4,matrix5, matrix6,matrix7, matrix8,matrix9, matrix10,matrix11, matrix12,matrix13, matrix14,matrix15, matrix16,matrix17, matrix18,matrix19, matrix20,flavors):
    from sympy import sin, cos, symbols
    from numpy.linalg import eig
    U1_matrix = matrix1
    U2_matrix = matrix2
    U3_matrix = matrix3
    U4_matrix = matrix4
    U5_matrix = matrix5
    U6_matrix = matrix6
    U7_matrix = matrix7
    U8_matrix = matrix8
    U9_matrix = matrix9
    U10_matrix = matrix10
    U11_matrix = matrix11
    U12_matrix = matrix12
    U13_matrix = matrix13
    U14_matrix = matrix14
    U15_matrix = matrix15
    U16_matrix = matrix16
    U17_matrix = matrix17
    U18_matrix = matrix18
    U19_matrix = matrix19
    U20_matrix = matrix20

   # print(np.dot(U2_matrix,U1_matrix))
    #print(np.dot(U1_matrix,U2_matrix))
#t = np.dot(U2_matrix,U1_matrix)
 #   U_total = np.dot(U3_matrix,t)
    k1= np.dot(U2_matrix,U1_matrix)
    k2= np.dot(U3_matrix,k1)
    k3= np.dot(U4_matrix,k2)
    k4= np.dot(U5_matrix,k3)
    k5= np.dot(U6_matrix,k4)
    k6= np.dot(U7_matrix,k5)
    k7= np.dot(U8_matrix,k6)
    k8= np.dot(U9_matrix,k7)
    k9= np.dot(U10_matrix,k8)
    k10= np.dot(U11_matrix,k9)
    k11= np.dot(U12_matrix,k10)
    k12= np.dot(U13_matrix,k11)
    k13= np.dot(U14_matrix,k12)
    k14= np.dot(U15_matrix,k13)
    k15= np.dot(U16_matrix,k14)
    k16= np.dot(U17_matrix,k15)
    k17= np.dot(U18_matrix,k16)
    k18= np.dot(U19_matrix,k17)
    k19= np.dot(U20_matrix,k18)
    

    U_total = k19#np.dot(U2_matrix,U1_matrix)
    return(U_total)



N = 1000
U1_modsqd01 = [0.0] * N
U2_modsqd01 = [0.0] * N
U1_modsqd10 = [0.0] * N
U2_modsqd10 = [0.0] * N

U1_modsqd00 = [0.0] * N
U2_modsqd00 = [0.0] * N

U1U2_modsqd01 = [0.0] * N
U2U1_modsqd01 = [0.0] * N

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

Aem = [0.0] * N
Aem_r = [0.0] * N

Aem_proper = [0.0] * N
Ame_proper = [0.0] * N

Pmm_vac = [0.0] * N
Pee_reversed = [0.0] * N
Pmm_reversed = [0.0] * N
Pem_reversed = [0.0] * N
Pme = [0.0] * N
Pme_reversed = [0.0] * N
Unitary1 = [0.0] * N
Unitary3 = [0.0] * N
Unitary4 = [0.0] * N
Unitary5 = [0.0] * N
Unitary6 = [0.0] * N
Unitary7 = [0.0] * N
Unitary8 = [0.0] * N
Unitary9 = [0.0] * N
Unitary10 = [0.0] * N
Unitary11 = [0.0] * N
Unitary12 = [0.0] * N
Unitary13 = [0.0] * N
Unitary14 = [0.0] * N
Unitary15 = [0.0] * N
Unitary16 = [0.0] * N
Unitary17 = [0.0] * N
Unitary18 = [0.0] * N
Unitary19 = [0.0] * N
Unitary20 = [0.0] * N
Unitary1b = [0.0] * N
LL1 = 800.0
Unitary2 = [0.0] * N
Unitary2b = [0.0] * N
LL2 = 2000-LL1
Unitary_total = [0.0] * N
Unitary_total_reversed = [0.0] * N
diff_em = [0.0] * N
diff_me = [0.0] * N
diff = [0.0] * N
diff_b = [0.0] * N
diff_em_PiecewiseSingle = [0.0] * N
diff_me_PiecewiseSingle = [0.0] * N
diff_em_PiecewiseVac = [0.0] * N
diff_me_PiecewiseVac = [0.0] * N
diff_em_PiecewiseSingle_r = [0.0] * N
diff_me_PiecewiseSingle_r = [0.0] * N
diff_em_PiecewiseVac_r = [0.0] * N
diff_me_PiecewiseVac_r = [0.0] * N
av_em = 0.0
av_me =	0.0
to_em =	0.0
to_me =	0.0
av_me_r =	0.0
to_me_r =	0.0
av_em_r =	0.0
to_em_r =	0.0

av_tm = 0.0
av_mt = 0.0
to_tm = 0.0
to_mt = 0.0
av_mt_r =       0.0
to_mt_r =       0.0
av_tm_r =       0.0
to_tm_r =       0.0

av_et = 0.0
av_te = 0.0
to_et = 0.0
to_te = 0.0
av_te_r =       0.0
to_te_r =       0.0
av_et_r =       0.0
to_et_r =       0.0



av_em_Vac =	0.0
av_me_Vac =	0.0
to_em_Vac =	0.0
to_me_Vac =	0.0

av_mm_Vac =     0.0
av_ee_Vac =     0.0
to_mm_Vac =     0.0
to_ee_Vac =     0.0

av_et_Vac =     0.0
av_te_Vac =     0.0
to_et_Vac =     0.0
to_te_Vac =     0.0

to_tt_Vac =     0.0
av_tt_Vac =     0.0
av_tm_Vac =     0.0
av_mt_Vac =     0.0
to_tm_Vac =     0.0
to_mt_Vac =     0.0
av_em_Single =	0.0
av_me_Single =	0.0
to_em_Single =	0.0
to_me_Single =	0.0


diff_me_single_Vac = [0.0] * N
diff_em_single_Vac = [0.0] * N
diff_ee_single_Vac = [0.0] * N
diff_mm_single_Vac = [0.0] * N
diff_et_single_Vac = [0.0] * N
diff_te_single_Vac = [0.0] * N
diff_tm_single_Vac = [0.0] * N
diff_mt_single_Vac = [0.0] * N
diff_tt_single_Vac = [0.0] * N

Unitary_single_matter = [0.0]*N
Unitary_vac = [0.0]*N
L = [0.0] * N
types = 3
V1 = 0.0001#0.0000107*3.0
V2 = 0.00015
V3 = 0.0002
V4 = 0.00025
V5 = 0.0003
V6 = 0.00035
V7 = 0.0004
V8 = 0.00045
V9 = 0.0005
V10 = 0.00055
V11 = 0.0006
V12 = 0.00065
V13 = 0.0007
V14 = 0.00075
V15 = 0.0008
V16 = 0.00085
V17 = 0.0009
V18 = 0.00095
V19 = 0.001
V20 = 0.00105

max_energy = 5.0
Step_size = (max_energy-0.5)/(N)
Step_size_L = (2000.0-0.0)/(N)

diff_ee_increase_piecewise_Vac = [0.0] * N
diff_me_increase_piecewise_Vac = [0.0] * N
diff_te_increase_piecewise_Vac = [0.0] * N

diff_mm_increase_piecewise_Vac = [0.0] * N
diff_em_increase_piecewise_Vac = [0.0] * N
diff_tm_increase_piecewise_Vac = [0.0] * N

diff_tt_increase_piecewise_Vac = [0.0] * N
diff_et_increase_piecewise_Vac = [0.0] * N
diff_mt_increase_piecewise_Vac = [0.0] * N

diff_ee_decrease_piecewise_Vac = [0.0] * N
diff_me_decrease_piecewise_Vac = [0.0] * N
diff_te_decrease_piecewise_Vac = [0.0] * N

diff_mm_decrease_piecewise_Vac = [0.0] * N
diff_em_decrease_piecewise_Vac = [0.0] * N
diff_tm_decrease_piecewise_Vac = [0.0] * N

diff_tt_decrease_piecewise_Vac = [0.0] * N
diff_et_decrease_piecewise_Vac = [0.0] * N
diff_mt_decrease_piecewise_Vac = [0.0] * N


diff_piecewise_e_decrease = [0.0] * N
diff_piecewise_m_decrease = [0.0] * N
diff_piecewise_t_decrease = [0.0] * N


diff_piecewise_e_increase = [0.0] * N
diff_piecewise_m_increase = [0.0] * N
diff_piecewise_t_increase = [0.0] * N


diff_single_e = [0.0] * N
diff_single_m = [0.0] * N
diff_single_t = [0.0] * N

diff_vac_e = [0.0] * N
diff_vac_m = [0.0] * N
diff_vac_t = [0.0] * N



diff_Pme_in_de_piecewise_matter = [0.0] * N
diff_Pem_in_de_piecewise_matter = [0.0] * N
diff_PmePem_increase_piecewise_matter = [0.0] * N
diff_PmePem_decrease_piecewise_matter = [0.0] * N

diff_Pmt_in_de_piecewise_matter = [0.0] * N
diff_Ptm_in_de_piecewise_matter = [0.0] * N
diff_PmtPtm_increase_piecewise_matter = [0.0] * N
diff_PmtPtm_decrease_piecewise_matter = [0.0] * N

diff_Pte_in_de_piecewise_matter = [0.0] * N
diff_Pet_in_de_piecewise_matter = [0.0] * N
diff_PetPte_increase_piecewise_matter = [0.0] * N
diff_PetPte_decrease_piecewise_matter = [0.0] * N

#doing 20 small increases for modeling a linear matter potential with equal baseline spacings of 100
vv=[0.0] *20
ll=[0.0] *20
for q in range(20):
    vv[q] = 0.0001+0.00001*q
    ll[q] = 100+100*q


for nu in range(N):
    EE[nu] = 0.5 + Step_size*nu

    
    Unitary1[nu] = Unitary_matrix(EE[nu], types, vv[0], 100)     
    Unitary2[nu] = Unitary_matrix(EE[nu], types, vv[1], 100)
    Unitary3[nu] = Unitary_matrix(EE[nu], types, vv[2], 100)
    Unitary4[nu] = Unitary_matrix(EE[nu], types, vv[3], 100)
    Unitary5[nu] = Unitary_matrix(EE[nu], types, vv[4], 100)
    Unitary6[nu] = Unitary_matrix(EE[nu], types, vv[5], 100)
    Unitary7[nu] = Unitary_matrix(EE[nu], types, vv[6], 100)
    Unitary8[nu] = Unitary_matrix(EE[nu], types, vv[7], 100)
    Unitary9[nu] = Unitary_matrix(EE[nu], types, vv[8], 100)
    Unitary10[nu] = Unitary_matrix(EE[nu], types, vv[9], 100)
    Unitary11[nu] = Unitary_matrix(EE[nu], types, vv[10], 100)
    Unitary12[nu] = Unitary_matrix(EE[nu], types, vv[11], 100)
    Unitary13[nu] = Unitary_matrix(EE[nu], types, vv[12], 100)
    Unitary14[nu] = Unitary_matrix(EE[nu], types, vv[13], 100)
    Unitary15[nu] = Unitary_matrix(EE[nu], types, vv[14], 100)
    Unitary16[nu] = Unitary_matrix(EE[nu], types, vv[15], 100)
    Unitary17[nu] = Unitary_matrix(EE[nu], types, vv[16], 100)
    Unitary18[nu] = Unitary_matrix(EE[nu], types, vv[17], 100)
    Unitary19[nu] = Unitary_matrix(EE[nu], types, vv[18], 100)
    Unitary20[nu] = Unitary_matrix(EE[nu], types, vv[19], 100)
    
  #  for numU in range(10):
 #       Unitary1[nu][numU] = Unitary_matrix(EE[nu], types, V[numU], LL[numU])
    #Unitary2[nu] = Unitary_matrix(EE[nu], types, V2, LL2)

    Unitary_total[nu] = U_multiply(Unitary1[nu],Unitary2[nu],Unitary3[nu],Unitary4[nu],Unitary5[nu],Unitary6[nu],Unitary7[nu],Unitary8[nu],Unitary9[nu],Unitary10[nu],Unitary11[nu],Unitary12[nu],Unitary13[nu],Unitary14[nu],Unitary15[nu],Unitary16[nu],Unitary17[nu],Unitary18[nu],Unitary19[nu],Unitary20[nu],types)
    #k2 = U_multiply(Unitary18[nu],k1,types)
#    Unitary_total[nu][numU] = U_multiply(Unitary1[nu][numU+1],Unitary2[nu][numU],types)
    

    Pee[nu] = probEE(Unitary_total[nu], types)
    Pmm[nu] = probMM(Unitary_total[nu], types)
    Pme[nu] = probME(Unitary_total[nu], types)
    Pem[nu] = probEM(Unitary_total[nu], types)

    Ptt[nu] = probTT(Unitary_total[nu], types)
    Pmt[nu] = probMT(Unitary_total[nu], types)
    Pte[nu] = probTE(Unitary_total[nu], types)
    Ptm[nu] = probTM(Unitary_total[nu], types)
    Pet[nu] = probET(Unitary_total[nu], types)

    
    Unitary_total_reversed[nu] = U_multiply(Unitary20[nu],Unitary19[nu],Unitary18[nu],Unitary17[nu],Unitary16[nu],Unitary15[nu],Unitary14[nu],Unitary13[nu],Unitary12[nu],Unitary11[nu],Unitary10[nu],Unitary9[nu],Unitary8[nu],Unitary7[nu],Unitary6[nu],Unitary5[nu],Unitary4[nu],Unitary3[nu],Unitary2[nu],Unitary1[nu],types)#U_multiply(Unitary2[nu],Unitary1[nu],types)

    print('test')
    print(Unitary1[nu])
    print(Unitary2[nu])
    #print(Unitary3[nu])
    print(Unitary_total_reversed[nu])
    print(Unitary_total[nu])

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


    
    Aem[nu] = (Pme[nu]-Pem[nu])/((Pme[nu]+Pem[nu]))
    Aem_r[nu] = (Pme_reversed[nu]-Pem_reversed[nu])/((Pme_reversed[nu]+Pem_reversed[nu]))
    Ame_proper[nu] = (Pme[nu]-Pem_reversed[nu])/((Pme[nu]+Pem_reversed[nu]))
    Aem_proper[nu] = (Pme_reversed[nu]-Pem[nu])/((Pme_reversed[nu]+Pem[nu]))

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

#plt.xlabel('                                  Energy (GeV)')

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




fig45 =plt.figure(45)
plt.plot(EE, Pme_vac, color='b', label='P_me vacuum ')
plt.plot(EE, Pem_vac, color='y',linestyle='dashed', label='P_em vacuum')
plt.title("Vacuum prob: nu_e and nu_mu",fontsize=10)
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
plt.title("single matter prob: nu_e and nu_mu",fontsize=10)
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
plt.title("decreasing constant matter prob: nu_e and nu_mu",fontsize=10)
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
plt.title("increasing constant matter prob: nu_e and nu_mu",fontsize=10)
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
"""
fig60 =plt.figure(60)
plt.plot(EE, Pme, color='b', label='P_me increasing constant matter ')
plt.plot(EE, Pme_reversed, color='r',linestyle='dashed', label='P_me decreasing constant matter ')
plt.plot(EE, Pem, color='g', linestyle='dotted',label='P_em increasing constant matter ')
plt.plot(EE, Pem_reversed, color='y',linestyle='dashdot', label='P_em decreasing constant matter ')
plt.title("linear matter prob: nu_e and nu_mu ",fontsize=10)
plt.xlabel('                                  Energy (GeV)')
plt.legend(fontsize=10)

fig61 =plt.figure(61)
plt.plot(ll, vv, color='b',marker="o")#, label='Asym. factor: increasing ')
plt.vlines(ll, ymin = 0.0001, ymax = vv,colors = 'purple')

#plt.axhline(y=vv[0],xmin = 100, xmax = 200)
#plt.axhline(y=vv[1],xmin = 200, xmax = 300)
#plt.plot(EE, Aem_r, color='m', linestyle='dashed',label='Asym. factor: decreasing')
#plt.plot(EE, Ame_proper, color='g', label='Asym. factor: P_nu_mu->nu_e increasing first')
#plt.plot(EE, Aem_proper, color='y', linestyle='dashed',label='Asym. factor: P_nu_mu->nu_e decreasing first')
plt.title("A(L)= 0.0001+0.00001*L ",fontsize=10)
plt.xlabel('                                  L (km)')
plt.ylabel('                                  A(L)')
#plt.legend(fontsize=10)

#print(vv,ll)   

plt.show()

