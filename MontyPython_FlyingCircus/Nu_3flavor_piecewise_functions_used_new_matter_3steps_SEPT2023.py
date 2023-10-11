import math
import matplotlib.pyplot as plt
import numpy as np
plt.rcParams.update({'figure.max_open_warning': 0})


def Unitary_matrix(energy, flavors, potential, baseline):    
    from sympy import sin, cos, symbols
    from numpy.linalg import eig

    #variable setup
    flavors = 3
    A = potential
    #L, dmSq, EE, Theta = symbols('L dmSq EE Theta')
    L = baseline
    d_CP = 90.00*math.pi/180#197.00*math.pi/180#232.00*math.pi/180#0.0#3.1415926535897932384#/2.0#0
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


    #Alpha = (5.076)*(dmSq)/(2.0*EE)
    #Beta = Alpha*math.sin(Theta)*math.sin(Theta)
    #Gamma = Alpha*math.cos(Theta)*math.sin(Theta)
    #Eta = Alpha*math.cos(Theta)*math.cos(Theta)
    M23 = np.array([[1, 0,0],
                  [0, c23,s23],
                  [0, -s23, c23]])
    M13 = np.array([[c13, 0,s13*np.exp(complex(0.0,-d_CP))],
                  [0, 1,0],
                  [-s13*np.exp(complex(0.0,d_CP)), 0, c13]])
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

    H_0 =(5.067)*(1/(2*EE))*(UM_0UT)
    
    val,vec=eig(H_0)

    print(val[0],vec[0])
    print(' ')
    print(val[1],vec[1])
    print(' ')
    print(val[2],vec[2])
    
    print('ON check')
    for a in range(flavors):
        for b in range(flavors):
             print(a,b,np.dot(vec[a],vec[b]))
    print(' ')

    #account for -i
    arg_exp = -1j*L#complex(0.0,-L)

    #unit vector for electron flavor
    for c in range(flavors):
        if c == 0:
            ve[c] = 1
        if c != 0:
            ve[c] = 0
            #ve = np.array([1,0,0])

    for c2 in range(flavors):
        if c2 == 1:
            vm[c2] = 1
        if c2 != 1:
            vm[c2] = 0

    for c3 in range(flavors):
        if c3 == 2:
            vt[c3] = 1
        if c3 != 1:
            vt[c3] = 0        
            
        #2x2 identity matrix
        for d in range(flavors):
            for e in range(flavors):
                if d == e:
                    Ident[d][e] = 1
                if d != e:
                    Ident[d][e] = 0
                    #np.array([[1.0, 0.0], 
                    #          [0.0, 1.0]])

        Id = np.array(Ident)        

        for a3 in range(flavors):
            proj[a3] = Id

        #for row construction and U_dagger
        for t1 in range(flavors):
            for t2 in range(flavors):
                u_dtest[t1][t2] = vec[t2][t1]

        print(u_dtest)
        print('U dag ')
        u_dtest2 = np.array(u_dtest)
        u_test = np.conjugate(np.transpose(u_dtest2))
        print(u_dtest2)
        print('corr. U ')
        print(u_test)
        print(' ')
        #print(u_test)
        #print(' ')

        for d2 in range(flavors):
            for e2 in range(flavors):
                if d2 == e2:
                    di[d2][e2] = np.exp(arg_exp*val[d2])#val[d2]
                if d2 != e2:
                    di[d2][e2] = 0.0

        Dig = np.array(di)
        print('di')
        print(Dig)
        print('Unitary')

#        t = np.dot(,np.dot(Dig,u_test))
        t = np.dot(u_dtest2,np.dot(Dig,u_test))
        print(t)
        #print(' ')
        #print(AA)
        #print(' ')

        
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
    return(Prob_ee.real)

def probME(matrix, flavors):
    from sympy import sin, cos, symbols
    from numpy.linalg import eig
    U_matrix = matrix
    ve = [0.0] * flavors
    vm = [0.0] * flavors
    Tme = U_matrix[1][0]
    Tme_C = np.conjugate(Tme)
    Prob_me = Tme_C*Tme
    return(Prob_me.real)


def probEM(matrix, flavors):
    from sympy import sin, cos, symbols
    from numpy.linalg import eig
    U_matrix = matrix
    ve = [0.0] * flavors
    vm = [0.0] * flavors
    Tem = U_matrix[0][1]
    Tem_C = np.conjugate(Tem)
    Prob_em = Tem_C*Tem
    return(Prob_em.real)

def probMM(matrix, flavors):
    from sympy import sin, cos, symbols
    from numpy.linalg import eig
    U_matrix = matrix
    ve = [0.0] * flavors
    vm = [0.0] * flavors
    Tmm = U_matrix[1][1]
    Tmm_C = np.conjugate(Tmm)
    Prob_mm = Tmm_C*Tmm
    return(Prob_mm.real)


def probET(matrix, flavors):
    from sympy import sin, cos, symbols
    from numpy.linalg import eig
    U_matrix = matrix
    ve = [0.0] * flavors
    vt = [0.0] * flavors
    Tet = U_matrix[0][2]
    Tet_C = np.conjugate(Tet)
    Prob_et = Tet_C*Tet
    return(Prob_et.real)

def probMT(matrix, flavors):
    from sympy import sin, cos, symbols
    from numpy.linalg import eig
    U_matrix = matrix
    vt = [0.0] * flavors
    vm = [0.0] * flavors
    Tmt = U_matrix[1][2]
    Tmt_C = np.conjugate(Tmt)
    Prob_mt = Tmt_C*Tmt
    return(Prob_mt.real)

def probTE(matrix, flavors):
    from sympy import sin, cos, symbols
    from numpy.linalg import eig
    U_matrix = matrix
    ve = [0.0] * flavors
    vt = [0.0] * flavors
    Tte = U_matrix[2][0]
    Tte_C = np.conjugate(Tte)
    Prob_te = Tte_C*Tte
    return(Prob_te.real)

def probTM(matrix, flavors):
    from sympy import sin, cos, symbols
    from numpy.linalg import eig
    U_matrix = matrix
    vm = [0.0] * flavors
    vt = [0.0] * flavors
    Ttm = U_matrix[2][1]
    Ttm_C = np.conjugate(Ttm)
    Prob_tm = Ttm_C*Ttm
    return(Prob_tm.real)

def probTT(matrix, flavors):
    from sympy import sin, cos, symbols
    from numpy.linalg import eig
    U_matrix = matrix
    #ve = [0.0] * flavors
    #vt = [0.0] * flavors
    Ttt = U_matrix[2][2]
    Ttt_C = np.conjugate(Ttt)
    Prob_tt = Ttt_C*Ttt
    return(Prob_tt.real)

def U_multiply(matrix1, matrix2, matrix3,flavors):
    from sympy import sin, cos, symbols
    from numpy.linalg import eig
    U1_matrix = matrix1
    U2_matrix = matrix2
    U3_matrix = matrix3
    t = np.dot(U2_matrix,U1_matrix)
    U_total = np.dot(U3_matrix,t)
#    U_total = np.dot(U3_matrix,U2_matrix,U1_matrix)
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
Unitary3 = [0.0] * N
LL2 = 300.0
LL3 = 2000.0-LL1-LL2
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

U1_modsqd01 = [0.0] * N
U2_modsqd01 = [0.0] * N
U1_modsqd10 = [0.0] * N
U2_modsqd10 = [0.0] * N

U1_modsqd00 = [0.0] * N
U2_modsqd00 = [0.0] * N

U3_modsqd00 = [0.0] * N
U3_modsqd10 = [0.0] * N
U3_modsqd01 = [0.0] * N
U1U2U3_modsqd01 = [0.0] * N
U3U2U1_modsqd01 = [0.0] * N
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
V1 = 0.000107#0.0000107*3.0
V2 = 0.000107*3.0
V3 = 0.000107#*6.0
max_energy = 5.0
Step_size = (max_energy-0.5)/(N)
Step_size_L = (2000.0-0.0)/(N)

e_appearance_increase_piecewise_matter = [0.0] * N
e_disappearance_increase_piecewise_matter = [0.0] * N
e_appearance_decrease_piecewise_matter = [0.0] * N
e_disappearance_decrease_piecewise_matter = [0.0] * N

m_appearance_increase_piecewise_matter = [0.0] * N
m_disappearance_increase_piecewise_matter = [0.0] * N
m_appearance_decrease_piecewise_matter = [0.0] * N
m_disappearance_decrease_piecewise_matter = [0.0] * N

t_appearance_increase_piecewise_matter = [0.0] * N
t_disappearance_increase_piecewise_matter = [0.0] * N
t_appearance_decrease_piecewise_matter = [0.0] * N
t_disappearance_decrease_piecewise_matter = [0.0] * N

e_appearance_vac = [0.0] * N
m_appearance_vac = [0.0] * N
t_appearance_vac = [0.0] * N

e_disappearance_vac = [0.0] * N
m_disappearance_vac = [0.0] * N
t_disappearance_vac = [0.0] * N

e_appearance_single_matter = [0.0] * N
m_appearance_single_matter = [0.0] * N
t_appearance_single_matter = [0.0] * N

e_disappearance_single_matter = [0.0] * N
m_disappearance_single_matter = [0.0] * N
t_disappearance_single_matter = [0.0] * N

diff_single_Vac_e_appearance = [0.0] * N
diff_single_Vac_m_appearance = [0.0] * N
diff_single_Vac_t_appearance = [0.0] * N

diff_single_Vac_e_disappearance = [0.0] * N
diff_single_Vac_m_disappearance = [0.0] * N
diff_single_Vac_t_disappearance = [0.0] * N

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

diff_increase_piecewise_Vac_e_appearance = [0.0] * N
diff_decrease_piecewise_Vac_e_appearance = [0.0] * N
diff_increase_piecewise_Vac_e_disappearance = [0.0] * N
diff_decrease_piecewise_Vac_e_disappearance = [0.0] * N

diff_increase_piecewise_Vac_m_appearance = [0.0] * N
diff_decrease_piecewise_Vac_m_appearance = [0.0] * N
diff_increase_piecewise_Vac_m_disappearance = [0.0] * N
diff_decrease_piecewise_Vac_m_disappearance = [0.0] * N

diff_increase_piecewise_Vac_t_appearance = [0.0] * N
diff_decrease_piecewise_Vac_t_appearance = [0.0] * N
diff_increase_piecewise_Vac_t_disappearance = [0.0] * N
diff_decrease_piecewise_Vac_t_disappearance = [0.0] * N


diff_piecewise_e_appearance = [0.0] * N
diff_piecewise_m_appearance = [0.0] * N
diff_piecewise_t_appearance = [0.0] * N

diff_piecewise_e_disappearance = [0.0] * N
diff_piecewise_m_disappearance = [0.0] * N
diff_piecewise_t_disappearance = [0.0] * N

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




#for val in range(N):
 #   L[val] = Step_size_L*val
#print(L)

for nu in range(N):
    EE[nu] = 0.5 + Step_size*nu
        
    Unitary1[nu] = Unitary_matrix(EE[nu], types, V1, LL1)
    Unitary2[nu] = Unitary_matrix(EE[nu], types, V2, LL2)
    Unitary3[nu] = Unitary_matrix(EE[nu], types, V3, LL3)
    Unitary_total[nu] = U_multiply(Unitary1[nu],Unitary2[nu],Unitary3[nu],types)

    #increasing constant matter
    Pee[nu] = probEE(Unitary_total[nu], types)
    Pmm[nu] = probMM(Unitary_total[nu], types)
    Pme[nu] = probME(Unitary_total[nu], types)
    Pem[nu] = probEM(Unitary_total[nu], types)

    Ptt[nu] = probTT(Unitary_total[nu], types)
    Pmt[nu] = probMT(Unitary_total[nu], types)
    Pte[nu] = probTE(Unitary_total[nu], types)
    Ptm[nu] = probTM(Unitary_total[nu], types)
    Pet[nu] = probET(Unitary_total[nu], types)



    
    
    Unitary_total_reversed[nu] = U_multiply(Unitary3[nu],Unitary2[nu],Unitary1[nu],types)
    #decreasing constant matter    

    
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
    #single matter
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
    #vacuum
    
    Pee_vac[nu] = probEE(Unitary_vac[nu], types)
    Pmm_vac[nu] = probMM(Unitary_vac[nu], types)
    Pme_vac[nu] = probME(Unitary_vac[nu], types)
    Pem_vac[nu] = probEM(Unitary_vac[nu], types)  

    Pet_vac[nu] = probET(Unitary_vac[nu], types)
    Pmt_vac[nu] = probMT(Unitary_vac[nu], types)
    Pte_vac[nu] = probTE(Unitary_vac[nu], types)
    Ptm_vac[nu] = probTM(Unitary_vac[nu], types)
    Ptt_vac[nu] = probTT(Unitary_vac[nu], types)



    U3U2U1_modsqd01[nu] = Unitary_total[nu][0][1]*np.conjugate(Unitary_total[nu][0][1])
    U1U2U3_modsqd01[nu] = Unitary_total_reversed[nu][0][1]*np.conjugate(Unitary_total_reversed[nu][0][1])


    U3U2U1_modsqd01[nu]= U3U2U1_modsqd01[nu].real
    U1U2U3_modsqd01[nu] = U1U2U3_modsqd01[nu].real

    print('check U3U2U1_modsqd01 == U1U2U3_modsqd01')
    print(U3U2U1_modsqd01[nu])
    print(U1U2U3_modsqd01[nu])

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
    to_em_Vac += Pem_vac[nu]
    to_me_Vac += Pme_vac[nu]
    to_mm_Vac += Pmm_vac[nu]
    to_ee_Vac += Pee_vac[nu]

    to_et_Vac += Pet_vac[nu]
    to_te_Vac += Pte_vac[nu]
    to_tm_Vac += Ptm_vac[nu]
    to_mt_Vac += Pmt_vac[nu]
    to_tt_Vac += Ptt_vac[nu]


    to_em_r += Pem_reversed[nu]
    to_me_r += Pme_reversed[nu]
    to_em += Pem[nu]
    to_me += Pme[nu]

    to_tm_r += Ptm_reversed[nu]
    to_mt_r += Pmt_reversed[nu]
    to_tm += Ptm[nu]
    to_mt += Pmt[nu]


    to_et_r += Pet_reversed[nu]
    to_te_r += Pte_reversed[nu]
    to_et += Pet[nu]
    to_te += Pte[nu]

    av_em_r =to_em_r/N
    av_me_r =to_me_r/N
    av_em =to_em/N
    av_me =to_me/N

    av_tm_r =to_tm_r/N
    av_mt_r =to_mt_r/N
    av_tm =to_tm/N
    av_mt =to_mt/N

    av_et_r =to_et_r/N
    av_te_r =to_te_r/N
    av_et =to_et/N
    av_te =to_te/N

av_em_Vac =to_em_Vac/N
av_me_Vac =to_me_Vac/N 

av_ee_Vac =to_ee_Vac/N
av_mm_Vac =to_mm_Vac/N
av_et_Vac =to_et_Vac/N
av_te_Vac =to_te_Vac/N
av_tm_Vac =to_tm_Vac/N
av_mt_Vac =to_mt_Vac/N
av_tt_Vac =to_tt_Vac/N

for nu2 in range(N):
    diff_ee_single_Vac[nu2] =(Pee_vac[nu2]-Pee_single_matter[nu2])/av_ee_Vac
    diff_me_single_Vac[nu2] =(Pme_vac[nu2]-Pme_single_matter[nu2])/av_me_Vac
    diff_te_single_Vac[nu2] =(Pte_vac[nu2]-Pte_single_matter[nu2])/av_te_Vac

    diff_mm_single_Vac[nu2] =(Pmm_vac[nu2]-Pmm_single_matter[nu2])/av_mm_Vac
    diff_em_single_Vac[nu2] =(Pem_vac[nu2]-Pem_single_matter[nu2])/av_em_Vac
    diff_tm_single_Vac[nu2] =(Ptm_vac[nu2]-Ptm_single_matter[nu2])/av_tm_Vac

    diff_tt_single_Vac[nu2] =(Ptt_vac[nu2]-Ptt_single_matter[nu2])/av_tt_Vac
    diff_et_single_Vac[nu2] =(Pet_vac[nu2]-Pet_single_matter[nu2])/av_et_Vac
    diff_mt_single_Vac[nu2] =(Pmt_vac[nu2]-Pmt_single_matter[nu2])/av_mt_Vac


    diff_ee_increase_piecewise_Vac[nu2] =(Pee_vac[nu2]-Pee[nu2])/av_ee_Vac
    diff_me_increase_piecewise_Vac[nu2] =(Pme_vac[nu2]-Pme[nu2])/av_me_Vac
    diff_te_increase_piecewise_Vac[nu2] =(Pte_vac[nu2]-Pte[nu2])/av_te_Vac

    diff_mm_increase_piecewise_Vac[nu2] =(Pmm_vac[nu2]-Pmm[nu2])/av_mm_Vac
    diff_em_increase_piecewise_Vac[nu2] =(Pem_vac[nu2]-Pem[nu2])/av_em_Vac
    diff_tm_increase_piecewise_Vac[nu2] =(Ptm_vac[nu2]-Ptm[nu2])/av_tm_Vac

    diff_tt_increase_piecewise_Vac[nu2] =(Ptt_vac[nu2]-Ptt[nu2])/av_tt_Vac
    diff_et_increase_piecewise_Vac[nu2] =(Pet_vac[nu2]-Pet[nu2])/av_et_Vac
    diff_mt_increase_piecewise_Vac[nu2] =(Pmt_vac[nu2]-Pmt[nu2])/av_mt_Vac



    diff_ee_decrease_piecewise_Vac[nu2] =(Pee_vac[nu2]-Pee_reversed[nu2])/av_ee_Vac
    diff_me_decrease_piecewise_Vac[nu2] =(Pme_vac[nu2]-Pme_reversed[nu2])/av_me_Vac
    diff_te_decrease_piecewise_Vac[nu2] =(Pte_vac[nu2]-Pte_reversed[nu2])/av_te_Vac

    diff_mm_decrease_piecewise_Vac[nu2] =(Pmm_vac[nu2]-Pmm_reversed[nu2])/av_mm_Vac
    diff_em_decrease_piecewise_Vac[nu2] =(Pem_vac[nu2]-Pem_reversed[nu2])/av_em_Vac
    diff_tm_decrease_piecewise_Vac[nu2] =(Ptm_vac[nu2]-Ptm_reversed[nu2])/av_tm_Vac

    diff_tt_decrease_piecewise_Vac[nu2] =(Ptt_vac[nu2]-Ptt_reversed[nu2])/av_tt_Vac
    diff_et_decrease_piecewise_Vac[nu2] =(Pet_vac[nu2]-Pet_reversed[nu2])/av_et_Vac
    diff_mt_decrease_piecewise_Vac[nu2] =(Pmt_vac[nu2]-Pmt_reversed[nu2])/av_mt_Vac




    diff_Pme_in_de_piecewise_matter[nu2] = (Pme_reversed[nu2]-Pme[nu2])/av_me
    diff_Pem_in_de_piecewise_matter[nu2] = (Pem_reversed[nu2]-Pem[nu2])/av_em
    diff_PmePem_increase_piecewise_matter[nu2] = (Pme[nu2]-Pem[nu2])/av_em
    diff_PmePem_decrease_piecewise_matter[nu2] = (Pme_reversed[nu2]-Pem_reversed[nu2])/av_em_r

    diff_Pmt_in_de_piecewise_matter[nu2] = (Pmt_reversed[nu2]-Pmt[nu2])/av_mt
    diff_Ptm_in_de_piecewise_matter[nu2] = (Ptm_reversed[nu2]-Ptm[nu2])/av_tm
    diff_PmtPtm_increase_piecewise_matter[nu2] = (Pmt[nu2]-Ptm[nu2])/av_mt
    diff_PmtPtm_decrease_piecewise_matter[nu2] = (Pmt_reversed[nu2]-Ptm_reversed[nu2])/av_mt_r

    diff_Pte_in_de_piecewise_matter[nu2] = (Pte_reversed[nu2]-Pte[nu2])/av_te
    diff_Pet_in_de_piecewise_matter[nu2] = (Pet_reversed[nu2]-Pet[nu2])/av_et
    diff_PetPte_increase_piecewise_matter[nu2] = (Pte[nu2]-Pet[nu2])/av_te
    diff_PetPte_decrease_piecewise_matter[nu2] = (Pte_reversed[nu2]-Pet_reversed[nu2])/av_te_r

    

    diff_piecewise_e_decrease[nu2] = e_appearance_decrease_piecewise_matter[nu2]-e_disappearance_decrease_piecewise_matter[nu2]
    diff_piecewise_m_decrease[nu2] = m_appearance_decrease_piecewise_matter[nu2]-m_disappearance_decrease_piecewise_matter[nu2]
    diff_piecewise_t_decrease[nu2] = t_appearance_decrease_piecewise_matter[nu2]-t_disappearance_decrease_piecewise_matter[nu2]

    diff_piecewise_e_increase[nu2] = e_appearance_increase_piecewise_matter[nu2]-e_disappearance_increase_piecewise_matter[nu2]
    diff_piecewise_m_increase[nu2] = m_appearance_increase_piecewise_matter[nu2]-m_disappearance_increase_piecewise_matter[nu2]
    diff_piecewise_t_increase[nu2] = t_appearance_increase_piecewise_matter[nu2]-t_disappearance_increase_piecewise_matter[nu2]



    diff_single_e[nu2] = e_appearance_single_matter[nu2]-e_disappearance_single_matter[nu2]
    diff_single_m[nu2] = m_appearance_single_matter[nu2]-m_disappearance_single_matter[nu2]
    diff_single_t[nu2] = t_appearance_single_matter[nu2]-t_disappearance_single_matter[nu2]

    diff_vac_e[nu2] = e_appearance_vac[nu2]-e_disappearance_vac[nu2]
    diff_vac_m[nu2] = m_appearance_vac[nu2]-m_disappearance_vac[nu2]
    diff_vac_t[nu2] = t_appearance_vac[nu2]-t_disappearance_vac[nu2]

fig16 =plt.figure(16)
plt.plot(EE, U3U2U1_modsqd01, color='b', label='(0,1) of |U3U2U1|^2 ')
plt.plot(EE, U1U2U3_modsqd01, color='m',linestyle='dashed', label='(0,1) of |U1U2U3|^2 ')
                       
plt.legend(fontsize=10)
          
plt.xlabel('                                  Energy (GeV)')

#fig16.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop/3_flavor_oscillations_10^-4V_3steps_CP90/off_dig_mod_sq.jpg')


#plt.show()

fig100 =plt.figure(100)
plt.plot(EE, diff_Pme_in_de_piecewise_matter, color='m',linestyle='dashed', label='diff_Pme_propVSimprop')
plt.plot(EE, diff_Pem_in_de_piecewise_matter, color='y', label='diff_Pem_propVSimprop')
plt.title("nu e and mu propVSimprop",fontsize=10)
plt.xlabel('                                  Energy (GeV)')
plt.legend(fontsize=10)


fig110 =plt.figure(110)
plt.plot(EE, diff_Pmt_in_de_piecewise_matter, color='m',linestyle='dashed', label='diff_Pmt_propVSimprop')
plt.plot(EE, diff_Ptm_in_de_piecewise_matter, color='y', label='diff_Ptm_propVSimprop')
plt.title("nu tau and mu propVSimprop",fontsize=10)
plt.xlabel('                                  Energy (GeV)')
plt.legend(fontsize=10)


fig120 =plt.figure(120)
plt.plot(EE, diff_Pte_in_de_piecewise_matter, color='m',linestyle='dashed', label='diff_Pte_propVSimprop')
plt.plot(EE, diff_Pet_in_de_piecewise_matter, color='y', label='diff_Ptm_propVSimprop')
plt.title("nu e and tau propVSimprop",fontsize=10)
plt.xlabel('                                  Energy (GeV)')
plt.legend(fontsize=10)

fig130 =plt.figure(130)
plt.plot(EE, diff_PmePem_increase_piecewise_matter, color='m',linestyle='dashed', label='diff e mu improp')
plt.plot(EE, diff_PmePem_decrease_piecewise_matter, color='y', label='diff e mu prop')
plt.title("nu e and mu prop and improp",fontsize=10)
plt.xlabel('                                  Energy (GeV)')
plt.legend(fontsize=10)


fig140 =plt.figure(140)
plt.plot(EE, diff_PmtPtm_increase_piecewise_matter, color='m',linestyle='dashed', label='diff t mu improp')
plt.plot(EE, diff_PmtPtm_decrease_piecewise_matter, color='y', label='diff t mu prop')
plt.title("nu tau and mu prop and improp",fontsize=10)
plt.xlabel('                                  Energy (GeV)')
plt.legend(fontsize=10)

fig150 =plt.figure(150)
plt.plot(EE, diff_PetPte_increase_piecewise_matter, color='m',linestyle='dashed', label='diff t e improp')
plt.plot(EE, diff_PetPte_decrease_piecewise_matter, color='y', label='diff t e prop')
plt.title("nu tau and e prop and improp",fontsize=10)
plt.xlabel('                                  Energy (GeV)')
plt.legend(fontsize=10)


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


fig7 =plt.figure(7)
plt.plot(EE, diff_me_single_Vac, color='m',linestyle='dashed', label='P_me diff_me_single_Vac ')
plt.plot(EE, diff_ee_single_Vac, color='y', label='P_ee diff_ee_single_Vac')
plt.plot(EE, diff_te_single_Vac, color='b', linestyle='dotted',label='P_te diff_te_single_Vac')
plt.title("diff_single_Vac prob nu_electron ",fontsize=10)
plt.xlabel('                                  Energy (GeV)')
plt.legend(fontsize=10)
fig8 =plt.figure(8)
plt.plot(EE, diff_em_single_Vac, color='m',linestyle='dashed', label='P_em diff_em_single_Vac ')
plt.plot(EE, diff_mm_single_Vac, color='y', label='P_mm diff_mm_single_Vac')
plt.plot(EE, diff_tm_single_Vac, color='b', linestyle='dotted',label='P_tm diff_tm_single_Vac')
plt.title("diff_single_Vac prob nu_muon ",fontsize=10)
plt.xlabel('                                  Energy (GeV)')
plt.legend(fontsize=10)
fig9 =plt.figure(9)
plt.plot(EE, diff_et_single_Vac, color='m',linestyle='dashed', label='P_et diff_et_single_Vac ')
plt.plot(EE, diff_tt_single_Vac, color='y', label='P_tt diff_tt_single_Vac')
plt.plot(EE, diff_mt_single_Vac, color='b', linestyle='dotted',label='P_mt diff_mt_single_Vac')
plt.title("diff_single_Vac prob nu_tau ",fontsize=10)
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



fig23 =plt.figure(23)
plt.plot(EE, diff_me_increase_piecewise_Vac, color='m',linestyle='dashed', label='P_me diff_me_increase_piecewise_Vac ')
plt.plot(EE, diff_ee_increase_piecewise_Vac, color='y', label='P_ee diff_ee_increase_piecewise_Vac')
plt.plot(EE, diff_te_increase_piecewise_Vac, color='b', linestyle='dotted',label='P_te diff_te_increase_piecewise_Vac')
plt.title("diff_increase_piecewise_Vac prob nu_electron ",fontsize=10)
plt.xlabel('                                  Energy (GeV)')
plt.legend(fontsize=10)
fig24 =plt.figure(24)
plt.plot(EE, diff_em_increase_piecewise_Vac, color='m',linestyle='dashed', label='P_em diff_em_increase_piecewise_Vac ')
plt.plot(EE, diff_mm_increase_piecewise_Vac, color='y', label='P_mm diff_mm_increase_piecewise_Vac')
plt.plot(EE, diff_tm_increase_piecewise_Vac, color='b', linestyle='dotted',label='P_tm diff_tm_increase_piecewise_Vac')
plt.title("diff_increase_piecewise_Vac prob nu_muon ",fontsize=10)
plt.xlabel('                                  Energy (GeV)')
plt.legend(fontsize=10)
fig25 =plt.figure(25)
plt.plot(EE,diff_et_increase_piecewise_Vac, color='m',linestyle='dashed', label='P_et diff_et_increase_piecewise_Vac ')
plt.plot(EE, diff_tt_increase_piecewise_Vac, color='y', label='P_tt diff_tt_increase_piecewise_Vac')
plt.plot(EE, diff_mt_increase_piecewise_Vac, color='b', linestyle='dotted',label='P_mt diff_mt_increase_piecewise_Vac')
plt.title("diff_increase_piecewise_Vac prob nu_tau ",fontsize=10)
plt.xlabel('                                  Energy (GeV)')
plt.legend(fontsize=10)


fig26 =plt.figure(26)
plt.plot(EE, diff_me_decrease_piecewise_Vac, color='m',linestyle='dashed', label='P_me diff_me_decrease_piecewise_Vac ')
plt.plot(EE, diff_ee_decrease_piecewise_Vac, color='y', label='P_ee diff_ee_decrease_piecewise_Vac')
plt.plot(EE, diff_te_decrease_piecewise_Vac, color='b', linestyle='dotted',label='P_te diff_te_decrease_piecewise_Vac')
plt.title("diff_decrease_piecewise_Vac prob nu_electron ",fontsize=10)
plt.xlabel('                                  Energy (GeV)')
plt.legend(fontsize=10)
fig27 =plt.figure(27)
plt.plot(EE, diff_em_decrease_piecewise_Vac, color='m',linestyle='dashed', label='P_em diff_em_decrease_piecewise_Vac ')
plt.plot(EE, diff_mm_decrease_piecewise_Vac, color='y', label='P_mm diff_mm_decrease_piecewise_Vac')
plt.plot(EE, diff_tm_decrease_piecewise_Vac, color='b', linestyle='dotted',label='P_tm diff_tm_decrease_piecewise_Vac')
plt.title("diff_decrease_piecewise_Vac prob nu_muon ",fontsize=10)
plt.xlabel('                                  Energy (GeV)')
plt.legend(fontsize=10)
fig28 =plt.figure(28)
plt.plot(EE,diff_et_decrease_piecewise_Vac, color='m',linestyle='dashed', label='P_et diff_et_decrease_piecewise_Vac ')
plt.plot(EE, diff_tt_decrease_piecewise_Vac, color='y', label='P_tt diff_tt_decrease_piecewise_Vac')
plt.plot(EE, diff_mt_decrease_piecewise_Vac, color='b', linestyle='dotted',label='P_mt diff_mt_decrease_piecewise_Vac')
plt.title("diff_decrease_piecewise_Vac prob nu_tau ",fontsize=10)
plt.xlabel('                                  Energy (GeV)')
plt.legend(fontsize=10)


fig29 =plt.figure(29)
plt.plot(EE, e_appearance_increase_piecewise_matter, color='m',linestyle='dashed', label='e_appearance_increase_piecewise_matter')
plt.plot(EE, m_appearance_increase_piecewise_matter, color='y', label='m_appearance_increase_piecewise_matter')
plt.plot(EE, t_appearance_increase_piecewise_matter, color='b', label='t_appearance_increase_piecewise_matter',linestyle='dotted')
plt.title("increase_piecewise prob appear",fontsize=10)
plt.xlabel('                                  Energy (GeV)')
plt.legend(fontsize=10)

fig30 =plt.figure(30)
plt.plot(EE,m_disappearance_increase_piecewise_matter,color='m',linestyle='dashed', label='m_disappearance_increase_piecewise_matter')
plt.plot(EE, e_disappearance_increase_piecewise_matter, color='y', label='e_disappearance_increase_piecewise_matter')
plt.plot(EE, t_disappearance_increase_piecewise_matter, color='b', linestyle='dotted',label='t_disappearance_increase_piecewise_matter')
plt.title("increase_piecewise prob disappear",fontsize=10)
plt.xlabel('                                  Energy (GeV)')
plt.legend(fontsize=10)

fig31 =plt.figure(31)
plt.plot(EE, e_appearance_decrease_piecewise_matter, color='m',linestyle='dashed', label='e_appearance_decrease_piecewise_matter')
plt.plot(EE, m_appearance_decrease_piecewise_matter, color='y', label='m_appearance_decrease_piecewise_matter')
plt.plot(EE, t_appearance_decrease_piecewise_matter, color='b', label='t_appearance_decrease_piecewise_matter',linestyle='dotted')
plt.title("decrease_piecewise prob appear",fontsize=10)
plt.xlabel('                                  Energy (GeV)')
plt.legend(fontsize=10)
fig32 =plt.figure(32)
plt.plot(EE,m_disappearance_decrease_piecewise_matter,color='m',linestyle='dashed', label='m_disappearance_decrease_piecewise_matter')
plt.plot(EE, e_disappearance_decrease_piecewise_matter, color='y', label='e_disappearance_decrease_piecewise_matter')
plt.plot(EE, t_disappearance_decrease_piecewise_matter, color='b', linestyle='dotted',label='t_disappearance_decrease_piecewise_matter')
plt.title("decrease_piecewise prob disappear",fontsize=10)
plt.xlabel('                                  Energy (GeV)')
plt.legend(fontsize=10)



fig33 =plt.figure(33)
plt.plot(EE,diff_decrease_piecewise_Vac_e_appearance, color='y', label='diff_decrease_piecewise_Vac_e_appearance')
plt.plot(EE,diff_decrease_piecewise_Vac_e_disappearance, color='b', linestyle='dotted',label='diff_decrease_piecewise_Vac_e_disappearance')
plt.title("diff_decrease_piecewise_Vac_e",fontsize=10)
plt.xlabel('                                  Energy (GeV)')
plt.legend(fontsize=10)
fig34 =plt.figure(34)
plt.plot(EE,diff_decrease_piecewise_Vac_m_appearance, color='y', label='diff_decrease_piecewise_Vac_m_appearance')
plt.plot(EE,diff_decrease_piecewise_Vac_m_disappearance, color='b', linestyle='dotted',label='diff_decrease_piecewise_Vac_m_disappearance')
plt.title("diff_decrease_piecewise_Vac_m",fontsize=10)
plt.xlabel('                                  Energy (GeV)')
plt.legend(fontsize=10)
fig35 =plt.figure(35)
plt.plot(EE,diff_decrease_piecewise_Vac_t_appearance, color='y', label='diff_decrease_piecewise_Vac_t_appearance')
plt.plot(EE,diff_decrease_piecewise_Vac_t_disappearance, color='b', linestyle='dotted',label='diff_decrease_piecewise_Vac_t_disappearance')
plt.title("diff_decrease_piecewise_Vac_t",fontsize=10)
plt.xlabel('                                  Energy (GeV)')
plt.legend(fontsize=10)


fig36 =plt.figure(36)
plt.plot(EE,diff_increase_piecewise_Vac_e_appearance, color='y', label='diff_increase_piecewise_Vac_e_appearance')
plt.plot(EE,diff_increase_piecewise_Vac_e_disappearance, color='b', linestyle='dotted',label='diff_increase_piecewise_Vac_e_disappearance')
plt.title("diff_increase_piecewise_Vac_e",fontsize=10)
plt.xlabel('                                  Energy (GeV)')
plt.legend(fontsize=10)
fig37 =plt.figure(37)
plt.plot(EE,diff_increase_piecewise_Vac_m_appearance, color='y', label='diff_increase_piecewise_Vac_m_appearance')
plt.plot(EE,diff_increase_piecewise_Vac_m_disappearance, color='b', linestyle='dotted',label='diff_increase_piecewise_Vac_m_disappearance')
plt.title("diff_increase_piecewise_Vac_m",fontsize=10)
plt.xlabel('                                  Energy (GeV)')
plt.legend(fontsize=10)
fig38 =plt.figure(38)
plt.plot(EE,diff_increase_piecewise_Vac_t_appearance, color='y', label='diff_increase_piecewise_Vac_t_appearance')
plt.plot(EE,diff_increase_piecewise_Vac_t_disappearance, color='b', linestyle='dotted',label='diff_increase_piecewise_Vac_t_disappearance')
plt.title("diff_increase_piecewise_Vac_t",fontsize=10)
plt.xlabel('                                  Energy (GeV)')
plt.legend(fontsize=10)

fig39 =plt.figure(39)
plt.plot(EE,diff_piecewise_e_appearance, color='y', label='diff_piecewise_e_appearance',linestyle='dashed')
plt.plot(EE,diff_piecewise_m_appearance,	color='b', label='diff_piecewise_m_appearance')
plt.plot(EE,diff_piecewise_t_appearance,	color='g', label='diff_piecewise_t_appearance',linestyle='dotted')
plt.title("diff_piecewise_appearance",fontsize=10)
plt.xlabel('                                  Energy (GeV)')
plt.legend(fontsize=10)

fig40 =plt.figure(40)
plt.plot(EE,diff_piecewise_e_disappearance,	color='y', label='diff_piecewise_e_disappearance',linestyle='dashed')
plt.plot(EE,diff_piecewise_m_disappearance, color='b', label='diff_piecewise_m_disappearance',linestyle='dotted')
plt.plot(EE,diff_piecewise_t_disappearance, color='g', label='diff_piecewise_t_disappearance')
plt.title("diff_piecewise_disappearance",fontsize=10)
plt.xlabel('                                  Energy (GeV)')
plt.legend(fontsize=10)


fig41 =plt.figure(41)
plt.plot(EE,diff_piecewise_e_increase,color='y', label='diff_piecewise_e_increase',linestyle='dashed')
plt.plot(EE,diff_piecewise_m_increase, color='b', label='diff_piecewise_m_increase',linestyle='dotted')
plt.plot(EE,diff_piecewise_t_increase, color='g', label='diff_piecewise_t_increase')
plt.title("diff_piecewise_increase",fontsize=10)
plt.xlabel('                                  Energy (GeV)')
plt.legend(fontsize=10)


fig42 =plt.figure(42)
plt.plot(EE,diff_piecewise_e_decrease,color='y', label='diff_piecewise_e_decrease',linestyle='dashed')
plt.plot(EE,diff_piecewise_m_decrease, color='b', label='diff_piecewise_m_decrease',linestyle='dotted')
plt.plot(EE,diff_piecewise_t_decrease, color='g', label='diff_piecewise_t_decrease')
plt.title("diff_piecewise_decrease",fontsize=10)
plt.xlabel('                                  Energy (GeV)')
plt.legend(fontsize=10)


fig43 =plt.figure(43)
plt.plot(EE,diff_single_e,color='y', label='diff_single_e',linestyle='dashed')
plt.plot(EE,diff_single_m,color='b', label='diff_single_m',linestyle='dotted')
plt.plot(EE,diff_single_t,color='g', label='diff_single_t')
plt.title("diff_single",fontsize=10)
plt.xlabel('                                  Energy (GeV)')
plt.legend(fontsize=10)


fig44 =plt.figure(44)
plt.plot(EE,diff_vac_e,color='y', label='diff_vac_e',linestyle='dashed')
plt.plot(EE,diff_vac_m,color='b', label='diff_vac_m',linestyle='dotted')
plt.plot(EE,diff_vac_t,color='g', label='diff_vac_t')
plt.title("diff_vac",fontsize=10)
plt.xlabel('                                  Energy (GeV)')
plt.legend(fontsize=10)




####################################################################################


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


fig60 =plt.figure(60)
plt.plot(EE, Pme, color='b', label='P_me increasing constant matter ')
plt.plot(EE, Pme_reversed, color='r',linestyle='dashed', label='P_me decreasing constant matter ')
plt.plot(EE, Pem, color='g', linestyle='dotted',label='P_em increasing constant matter ')
plt.plot(EE, Pem_reversed, color='y',linestyle='dashdot', label='P_em decreasing constant matter ')
plt.title("piecewise constant matter prob: nu_e and nu_mu ",fontsize=10)
plt.xlabel('                                  Energy (GeV)')
plt.legend(fontsize=10)


plt.show()
#
#time_invar_propVSimprop/3_flavor_oscillations_10^-4V_CP90

fig1.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23/3_flavor_oscillations_10^-4V_3steps_CP90/Vacuum_NuElectron.jpg')
fig2.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23/3_flavor_oscillations_10^-4V_3steps_CP90/Vacuum_NuMuon.jpg')
fig3.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23/3_flavor_oscillations_10^-4V_3steps_CP90/Vacuum_NuTau.jpg')
fig4.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23/3_flavor_oscillations_10^-4V_3steps_CP90/Single_matter_NuElectron.jpg')
fig5.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23/3_flavor_oscillations_10^-4V_3steps_CP90/Single_matter_NuMuon.jpg')
fig6.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23/3_flavor_oscillations_10^-4V_3steps_CP90/Single_matter_NuTau.jpg')

fig7.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23/3_flavor_oscillations_10^-4V_3steps_CP90/diff_single_Vac_prob_nuElectron.jpg')
fig8.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23/3_flavor_oscillations_10^-4V_3steps_CP90/diff_single_Vac_prob_nuMuon.jpg')
fig9.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23/3_flavor_oscillations_10^-4V_3steps_CP90/diff_single_Vac_prob_nuTau.jpg')

fig17.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23/3_flavor_oscillations_10^-4V_3steps_CP90/nuE_constant_increase_piecewise_matter_prob.jpg')
fig18.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23/3_flavor_oscillations_10^-4V_3steps_CP90/nuM_constant_increase_piecewise_matter_prob.jpg')
fig19.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23/3_flavor_oscillations_10^-4V_3steps_CP90/nuT_constant_increase_piecewise_matter_prob.jpg')

fig20.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23/3_flavor_oscillations_10^-4V_3steps_CP90/nuE_constant_decrease_piecewise_matter_prob.jpg')
fig21.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23/3_flavor_oscillations_10^-4V_3steps_CP90/nuM_constant_decrease_piecewise_matter_prob.jpg')
fig22.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23/3_flavor_oscillations_10^-4V_3steps_CP90/nuT_constant_decrease_piecewise_matter_prob.jpg')

fig23.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23/3_flavor_oscillations_10^-4V_3steps_CP90/diff_increase_piecewise_Vac_prob_nuE.jpg')
fig24.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23/3_flavor_oscillations_10^-4V_3steps_CP90/diff_increase_piecewise_Vac_prob_nuM.jpg')
fig25.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23/3_flavor_oscillations_10^-4V_3steps_CP90/diff_increase_piecewise_Vac_prob_nuT.jpg')

fig26.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23/3_flavor_oscillations_10^-4V_3steps_CP90/diff_decrease_piecewise_Vac_prob_nuE.jpg')
fig27.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23/3_flavor_oscillations_10^-4V_3steps_CP90/diff_decrease_piecewise_Vac_prob_nuM.jpg')
fig28.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23/3_flavor_oscillations_10^-4V_3steps_CP90/diff_decrease_piecewise_Vac_prob_nuT.jpg')

fig36.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23/3_flavor_oscillations_10^-4V_3steps_CP90/diff_increase_piecewise_Vac_e.jpg')
fig37.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23/3_flavor_oscillations_10^-4V_3steps_CP90/diff_increase_piecewise_Vac_m.jpg')
fig38.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23/3_flavor_oscillations_10^-4V_3steps_CP90/diff_increase_piecewise_Vac_t.jpg')


fig41.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23/3_flavor_oscillations_10^-4V_3steps_CP90/diff_piecewise_increase.jpg')
fig42.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23/3_flavor_oscillations_10^-4V_3steps_CP90/diff_piecewise_decrease.jpg')

fig43.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23/3_flavor_oscillations_10^-4V_3steps_CP90/diff_single.jpg')
fig44.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23/3_flavor_oscillations_10^-4V_3steps_CP90/diff_vac.jpg')

fig45.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23/3_flavor_oscillations_10^-4V_3steps_CP90/Vacuum_compar_nuE_muMU.jpg')
fig46.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23/3_flavor_oscillations_10^-4V_3steps_CP90/Vacuum_compar_nuMU_muTAU.jpg')
fig47.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23/3_flavor_oscillations_10^-4V_3steps_CP90/Vacuum_compar_nuE_muTAU.jpg')

fig48.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23/3_flavor_oscillations_10^-4V_3steps_CP90/Single_matter_compar_nuE_muMU.jpg')
fig49.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23/3_flavor_oscillations_10^-4V_3steps_CP90/Single_matter_compar_nuMU_muTAU.jpg')
fig50.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23/3_flavor_oscillations_10^-4V_3steps_CP90/Single_matter_compar_nuE_muTAU.jpg')


fig51.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23/3_flavor_oscillations_10^-4V_3steps_CP90/Decreasing_piecewise_matter_compar_nuE_muMU.jpg')
fig52.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23/3_flavor_oscillations_10^-4V_3steps_CP90/Decreasing_piecewise_matter_compar_nuMU_muTAU.jpg')
fig53.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23/3_flavor_oscillations_10^-4V_3steps_CP90/Decreasing_piecewise_matter_compar_nuE_muTAU.jpg')


fig54.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23/3_flavor_oscillations_10^-4V_3steps_CP90/Increasing_piecewise_matter_compar_nuE_muMU.jpg')
fig55.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23/3_flavor_oscillations_10^-4V_3steps_CP90/Increasing_piecewise_matter_compar_nuMU_muTAU.jpg')
fig56.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23/3_flavor_oscillations_10^-4V_3steps_CP90/Increasing_piecewise_matter_compar_nuE_muTAU.jpg')

fig100.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23/3_flavor_oscillations_10^-4V_3steps_CP90/Pem_Pme_propVSimprop.jpg')
fig110.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23/3_flavor_oscillations_10^-4V_3steps_CP90/Ptm_Pmt_propVSimprop.jpg')
fig120.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23/3_flavor_oscillations_10^-4V_3steps_CP90/Pet_Pte_propVSimprop.jpg')
fig130.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23/3_flavor_oscillations_10^-4V_3steps_CP90/Pem_Pme_propANDimprop.jpg')

fig140.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23/3_flavor_oscillations_10^-4V_3steps_CP90/Ptm_Pmt_propANDimprop.jpg')
fig150.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23/3_flavor_oscillations_10^-4V_3steps_CP90/Pet_Pte_propANDimprop.jpg')

fig60.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop_OCT23/3_flavor_oscillations_10^-4V_3steps_CP90/Pem_Pem_Im_P_COMPARE.jpg') 
