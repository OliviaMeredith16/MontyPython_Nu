import math
import matplotlib.pyplot as plt
import numpy as np

filename = 'output_2F_10-4V_U1_U2_transpose.txt'
f = open(filename, 'w')
filename = 'output_2F_10-4V_U1_U2_conjugate.txt'
fc = open(filename, 'w')
filename = 'output_2F_10-4V_U1_U2_daggered.txt'
fd = open(filename, 'w')
filename = 'output_2F_10-4V_U1U2_U2U1_neg.txt'
fw = open(filename, 'w')
filename = 'output_2F_10-4V_U1U2_U2U1_flipped.txt'
fr = open(filename, 'w')

def Unitary_matrix(energy, flavors, potential, baseline):    
    from sympy import sin, cos, symbols
    from numpy.linalg import eig

    #variable setup
    #from sin2thetaSquared of 0.8 we get Theta value
    #flavors = 2
    A = potential
    #L, dmSq, EE, Theta = symbols('L dmSq EE Theta')
    L = baseline
    dmSq = 0.001 #0.0547723
    EE = energy
    Theta = 0.58
    proj = [0.0] * flavors
    ve = [0.0] * flavors
    vm = [0.0] * flavors
    Ident = [[0 for x in range(flavors)] for y in range(flavors)]
    U_s = [0.0] * flavors
    U_f = 0.0

    t = 0.0
    di = [[0 for x in range(flavors)] for y in range(flavors)]
    u_dtest = [[0 for x in range(flavors)] for y in range(flavors)]
    u_test = [[0 for x in range(flavors)] for y in range(flavors)]


    Alpha = (5.076)*(dmSq)/(2.0*EE)
    Beta = Alpha*math.sin(Theta)*math.sin(Theta)
    Gamma = Alpha*math.cos(Theta)*math.sin(Theta)
    Eta = Alpha*math.cos(Theta)*math.cos(Theta)

    AA2 = np.array([[(5.076)*A , 0],
                   [0, 0]])
    #matrix setup + eigensystem
    AA1 = np.array([[Beta , Gamma], 
                   [Gamma, Eta]])
    AA = AA1+AA2#np.array([[Beta, Gamma,0],
    #              [Gamma, Eta,0],
    #             [0, 0, 1]])
    val,vec=eig(AA)

#    print('E-value:', val)
 #   print('E-vector', vec)
  #  print(' ')

    #ON check
#    for a in range(flavors):
 #       for b in range(flavors):
    #         print(np.dot(vec[a],vec[b]))
    #print(' ')

    #account for -i
    arg_exp = complex(0.0,-L)

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

        #3x3 identity matrix                                                      
            #Ident3 = np.array([[1.0, 0.0, 0.0],
            #                  [0.0, 1.0, 0.0],
            #                   [0.0, 0.0, 1.0]])

        #projction operators
            #for f in range(flavors-1):
            #    proj[f] = 1
            #   for k in range(flavors):
            #        if k == f:
            #           continue
            #      if k != f:
        #         proj[f] =*(AA - val[k]*Id)/(val[f]-val[k])

        #print(proj[f])        
        #print(' ')
            
        for a1 in range(flavors):
            for a2 in range(flavors):
                if a1 != a2:
                    proj[a1] = np.dot(proj[a1],( (AA - val[a2]*Id)/(val[a1]-val[a2]) ) ) 
                if a1 == a2:
                    continue 
     #   print(proj[a1])
      #  print(" ")
    
        #proj[0] = (AA - val[1]*Id)/(val[0]-val[1])
        #proj[1] = (AA - val[0]*Id)/(val[1]-val[0])

        #for f in range(flavors-1):
        #   for k in range(flavors-1):
        #      proj[f] *= (AA - val[k+1]*Id)/(val[f]-val[k+1])#*(AA - val[2]*Ident)/(val[0]-val[2])
        # proj2 *= (AA - val[f]*Id)/(val[1]-val[f])#*(AA - val[2]*Ident3)/(val[1]-val[2])
        #if f == 1:
        #   proj2 *= (AA - val[f+1]*Id)/(val[1]-val[f+1])
        #proj[2] = (AA - val[0]*Ident3)/(val[2]-val[0])*(AA - val[1]*Ident3)/(val[2]-val[1])

        #full unitary evolution matrix
        for g in range(flavors):   
            U_s[g] = np.exp(arg_exp*val[g])*proj[g]# + np.exp(arg_exp*val[1])*proj[1]
            #+ np.exp(arg_exp*val[2])*proj[2]
            #print(U_s)
        #print(' ')

        for h in range(flavors):
            U_f += U_s[h]

           # u_z1 = np.array([[0 , 0],
            #                [0, 0]])
  #      tt = 0.0#[[0 for x in range(flavors)] for y in range(flavors)]
        
        print(' from function: unitary evolution matrix ')
        print(U_f)
   #     for g1 in range(flavors):
    #        for g2 in range(flavors):
     #           if g1 != g2:
      #              print(g1,g2,file=f)
       #             print(U_f[g1][g2],file=f)

                    
        #            tt += U_f[g1][g2]*complex(U_f[g1][g2])
                    
          #          print(g1,g2,file=fm)
         #           print(U_f[g1][g2]*complex(U_f[g1][g2]),file=fm)
        #print(' ',file=f)
        #print(tt,file=fm)
        #print(' ',file=fm)

        
        #check everything is working right

        #for row construction and U_dagger
        for t1 in range(flavors):
            for t2 in range(flavors):
                u_dtest[t1][t2] = vec[t1][t2]

                #print(vec[0])
                #print(vec[1])
        #print('matrix for U_dagger')
        #print(u_dtest)  
        
        u_dtest2 = np.array(u_dtest)
        u_test = np.conjugate(np.transpose(u_dtest2))
        #print(u_test)
        #print(' ')

        for d2 in range(flavors):
            for e2 in range(flavors):
                if d2 == e2:
                    di[d2][e2] = val[d2]
                if d2 != e2:
                    di[d2][e2] = 0

        Dig = np.array(di)
        #print(Dig)
        #print(' ')

        t = np.dot(u_dtest,np.dot(Dig,u_test))
        #print(t)
        #print(' ')
        #print(AA)
        #print(' ')


        return(U_f)

    


    #####
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
Pem_single_matter = [0.0] * N
Pem_vac = [0.0] * N
Pme_single_matter = [0.0] * N
Pme_vac = [0.0] * N
Pee_single_matter = [0.0] * N
Pee_vac = [0.0] * N
Pmm_single_matter = [0.0] * N
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
LL2 = 2000.0-LL1
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
diff_me_constVac = [0.0] * N
diff_em_constVac = [0.0] * N
diff_single = [0.0] * N
diff_vac = [0.0] * N
av_em = 0.0
av_me =	0.0
to_em =	0.0
to_me =	0.0
av_me_r =	0.0
to_me_r =	0.0
av_em_r =	0.0
to_em_r =	0.0
av_em_Vac =	0.0
av_me_Vac =	0.0
to_em_Vac =	0.0
to_me_Vac =	0.0
av_em_Single =	0.0
av_me_Single =	0.0
to_em_Single =	0.0
to_me_Single =	0.0

Unitary_single_matter = [0.0]*N
Unitary_vac = [0.0]*N
L = [0.0] * N
types = 2
V1 = 0.000107
V2 = 0.000107*3.0
max_energy = 5.0
Step_size = (max_energy-0.5)/(N)
Step_size_L = (2000.0-0.0)/(N)

rUT00_con = [0.0] * N
rUT10_con = [0.0] * N
rUT01_con = [0.0] * N
rUT11_con = [0.0] * N
rUT00R = [0.0] * N
rUT10R = [0.0] * N
rUT01R = [0.0] * N
rUT11R = [0.0] * N

iUT00_con = [0.0] * N
iUT10_con = [0.0] * N
iUT01_con = [0.0] * N
iUT11_con = [0.0] * N
iUT00R = [0.0] * N
iUT10R = [0.0] * N
iUT01R = [0.0] * N
iUT11R = [0.0] * N
iUT00 = [0.0] * N
iUT10 = [0.0] * N
iUT01 = [0.0] * N
iUT11 = [0.0] * N
rUT00 = [0.0] * N
rUT10 = [0.0] * N
rUT01 = [0.0] * N
rUT11 = [0.0] * N

U1_modsqd01 = [0.0] * N
U2_modsqd01 = [0.0] * N
U1_modsqd10 = [0.0] * N
U2_modsqd10 = [0.0] * N

U1_modsqd00 = [0.0] * N
U2_modsqd00 = [0.0] * N

U1U2_modsqd01 = [0.0] * N
U2U1_modsqd01 = [0.0] * N
#for val in range(N):
 #   L[val] = Step_size_L*val
#print(L)

for nu in range(N):
    print(nu,file=f)
    print(nu,file=fc)
    print(nu,file=fd)
    EE[nu] = 0.5 + Step_size*nu
    Unitary1[nu] = Unitary_matrix(EE[nu], types, V1, LL1)
    Unitary2[nu] = Unitary_matrix(EE[nu], types, V2, LL2)

    print('U1')
    print(Unitary1[nu])
    print('U2')
    print(Unitary2[nu])
    U1_modsqd00[nu] = (Unitary1[nu][0][0]*np.conjugate(Unitary1[nu][0][0])).real
    U2_modsqd00[nu] = (Unitary1[nu][0][0]*np.conjugate(Unitary2[nu][0][0])).real

    U1_modsqd01[nu] = (Unitary1[nu][0][1]*np.conjugate(Unitary1[nu][0][1])).real
    U2_modsqd01[nu] = (Unitary1[nu][0][1]*np.conjugate(Unitary2[nu][0][1])).real
    U1_modsqd10[nu] = (Unitary1[nu][1][0]*np.conjugate(Unitary1[nu][1][0])).real
    U2_modsqd10[nu] = (Unitary1[nu][1][0]*np.conjugate(Unitary2[nu][1][0])).real

    print('testing U1 diag mod terms')
    print(U1_modsqd01[nu])
    print(U1_modsqd10[nu])
    print('unitarity check')
    print(U1_modsqd00[nu]+U1_modsqd01[nu])
    print(U2_modsqd00[nu]+U2_modsqd01[nu])

    
    u1 = Unitary1[nu]
    u2 = Unitary2[nu]
    u_z = np.array([[0 , 0],
                   [0, 0]])
   # for s in range(types):
    #    for q in range(types):
     #       if(s!=q)
    # print('sum of off diagonals')
   # print(sum1)
   # print(sum2)
    
    #    Unitary3[nu] = Unitary_matrix(EE[nu], types, V1*3.0, 1200)
    Unitary_total[nu] = U_multiply(Unitary1[nu],Unitary2[nu],types)
    print('U2U1')
    print(Unitary_total[nu])

    Pee[nu] = probEE(Unitary_total[nu], types)
    Pmm[nu] = probMM(Unitary_total[nu], types)
    Pme[nu] = probME(Unitary_total[nu], types)
    Pem[nu] = probEM(Unitary_total[nu], types)


    #print('prob check')
    #print(Pee[nu]+Pem[nu])
    #print(Pem[nu]+Pmm[nu])
    #print(' ')
    Unitary_total_reversed[nu] = U_multiply(Unitary2[nu],Unitary1[nu],types)

    print('U1U2')
    print(Unitary_total_reversed[nu])
    Pee_reversed[nu] = probEE(Unitary_total_reversed[nu], types)
    Pmm_reversed[nu] = probMM(Unitary_total_reversed[nu], types)
    Pme_reversed[nu] = probME(Unitary_total_reversed[nu], types)
    Pem_reversed[nu] = probEM(Unitary_total_reversed[nu], types)

    rUT00_con[nu] = (np.conjugate(Unitary_total[nu][0][0])).real
    rUT01_con[nu] = (np.conjugate(Unitary_total[nu][0][1])).real
    rUT10_con[nu] = (np.conjugate(Unitary_total[nu][1][0])).real
    rUT11_con[nu] = (np.conjugate(Unitary_total[nu][1][1])).real

    rUT00R[nu] = (Unitary_total_reversed[nu][0][0]).real
    rUT01R[nu] = (Unitary_total_reversed[nu][0][1]).real
    rUT10R[nu] = (Unitary_total_reversed[nu][1][0]).real
    rUT11R[nu] = (Unitary_total_reversed[nu][1][1]).real
    
    rUT00[nu] = (Unitary_total[nu][0][0]).real
    rUT01[nu] = (Unitary_total[nu][0][1]).real
    rUT10[nu] = (Unitary_total[nu][1][0]).real
    rUT11[nu] = (Unitary_total[nu][1][1]).real

    iUT00_con[nu] = (np.conjugate(Unitary_total[nu][0][0])).imag
    iUT01_con[nu] = (np.conjugate(Unitary_total[nu][0][1])).imag
    iUT10_con[nu] = (np.conjugate(Unitary_total[nu][1][0])).imag
    iUT11_con[nu] = (np.conjugate(Unitary_total[nu][1][1])).imag

    iUT00R[nu] = (Unitary_total_reversed[nu][0][0]).imag
    iUT01R[nu] = (Unitary_total_reversed[nu][0][1]).imag
    iUT10R[nu] = (Unitary_total_reversed[nu][1][0]).imag
    iUT11R[nu] = (Unitary_total_reversed[nu][1][1]).imag

    iUT00[nu] = (Unitary_total[nu][0][0]).imag
    iUT01[nu] = (Unitary_total[nu][0][1]).imag
    iUT10[nu] = (Unitary_total[nu][1][0]).imag
    iUT11[nu] = (Unitary_total[nu][1][1]).imag


    U2U1_modsqd01[nu] = Unitary_total[nu][0][1]*np.conjugate(Unitary_total[nu][0][1])
    U1U2_modsqd01[nu] = Unitary_total_reversed[nu][0][1]*np.conjugate(Unitary_total_reversed[nu][0][1])

    print('check U2U1_modsqd01 == U1U2_modsqd01')
    print(U2U1_modsqd01[nu])
    print(U1U2_modsqd01[nu])
    
#    print(rUT00R[nu])
 #   print(rUT01R[nu])
  #  print(rUT10R[nu])
   # print(rUT11R[nu])
    
    if abs(Unitary_total[nu][0][0] - np.conjugate(Unitary_total_reversed[nu][0][0])) <0.000001:
        if abs(Unitary_total[nu][1][1] - np.conjugate(Unitary_total_reversed[nu][1][1])) <0.000001:
            if abs(Unitary_total[nu][1][0] - np.conjugate(Unitary_total_reversed[nu][0][1])) <0.000001:
                if abs(Unitary_total[nu][0][1] - np.conjugate(Unitary_total_reversed[nu][1][0])) <0.000001:
                    print('U1U2 and U2U1 are complex conjuagates of each other at ', nu,EE[nu],file=fd)

    if abs(Unitary_total[nu][0][0] - np.conjugate(Unitary_total_reversed[nu][0][0])) <0.000001:
        if abs(Unitary_total[nu][1][1] - np.conjugate(Unitary_total_reversed[nu][1][1])) <0.000001:
            if abs(Unitary_total[nu][1][0] - np.conjugate(Unitary_total_reversed[nu][1][0])) <0.000001:
                if abs(Unitary_total[nu][0][1] - np.conjugate(Unitary_total_reversed[nu][0][1])) <0.000001:
                    print('U1U2 and U2U1 are conjuagates of each other at ', nu,EE[nu],file=fc)

    if abs(Unitary_total[nu][0][0] - Unitary_total_reversed[nu][0][0]) <0.000001:
        if abs(Unitary_total[nu][1][1] - Unitary_total_reversed[nu][1][1]) <0.000001:
            if abs(Unitary_total[nu][1][0] - Unitary_total_reversed[nu][0][1]) <0.000001:
                if abs(Unitary_total[nu][0][1] - Unitary_total_reversed[nu][1][0]) <0.000001:
                    print('U1U2 and U2U1 are transposes of each other at ', nu,EE[nu],file=f)    

    v01 = 0.0
    v10  = 0.0 
    v01 = Unitary_total[nu][0][1] - Unitary_total_reversed[nu][0][1]
    v10 = Unitary_total[nu][0][1] - Unitary_total_reversed[nu][0][1]

    n01 = Unitary_total[nu][0][1] - Unitary_total_reversed[nu][1][0]
    n10 = Unitary_total[nu][0][1] - Unitary_total_reversed[nu][1][0]
    if v01 +v10 <0.01:
        print('U2U1 - U1U2: 0,1 is the negative of 1,0 at ', nu,EE[nu],file=fw)

    if n10<0.01:
        if n01<0.01:
            print('U2U1(0,1) = U1U2(1,0) and U2U1(1,0) = U1U2(0,1) at ', nu,EE[nu],file=fr)    
        
    print(Unitary_total[nu] - Unitary_total_reversed[nu])



    
   # print('prob check')
    #print(Pee_reversed[nu]+Pem_reversed[nu])
    #print(Pem_reversed[nu]+Pmm_reversed[nu])
    #print(' ')

    #print('ignore after here ')    
    Unitary_single_matter[nu] = Unitary_matrix(EE[nu], types, V1,2000)

    Pee_single_matter[nu] = probEE(Unitary_single_matter[nu], types)
    Pmm_single_matter[nu] = probMM(Unitary_single_matter[nu], types)
    Pme_single_matter[nu] = probME(Unitary_single_matter[nu], types)
    Pem_single_matter[nu] = probEM(Unitary_single_matter[nu], types)

    #print('prob check')
    #print(Pee_single_matter[nu]+Pem_single_matter[nu])
    #print(Pem_single_matter[nu]+Pmm_single_matter[nu])
    #print(' ')
    Unitary_vac[nu] = Unitary_matrix(EE[nu], types, 0.0, 2000)

    Pee_vac[nu] = probEE(Unitary_vac[nu], types)
    Pmm_vac[nu] = probMM(Unitary_vac[nu], types)
    Pme_vac[nu] = probME(Unitary_vac[nu], types)
    Pem_vac[nu] = probEM(Unitary_vac[nu], types)  

    #print('prob check')
    #print(Pee_vac[nu]+Pem_vac[nu])
    #print(Pem_vac[nu]+Pmm_vac[nu])
    #print(' ')
    to_em_r += Pem_reversed[nu]
    to_me_r += Pme_reversed[nu]
    to_em += Pem[nu]
    to_me += Pme[nu]
    to_em_Single += Pem_single_matter[nu]
    to_me_Single += Pme_single_matter[nu]
    to_em_Vac += Pem_vac[nu]
    to_me_Vac += Pme_vac[nu]


#    diff_me[nu] =(Pme[nu]-Pme_reversed[nu])
 #   diff_em[nu] =(Pem[nu]-Pem_reversed[nu])
av_em_r =to_em_r/N
av_me_r =to_me_r/N
av_em =to_em/N
av_me =to_me/N
av_em_Single =to_em_Single/N
av_me_Single =to_me_Single/N
av_em_Vac =to_em_Vac/N
av_me_Vac =to_me_Vac/N 

for nu2 in range(N):  
    diff_me[nu2] =(Pme[nu2]-Pme_reversed[nu2])/av_me
    diff_em[nu2] =(Pem[nu2]-Pem_reversed[nu2])/av_em
    diff[nu2] =(Pme[nu2]-Pem[nu2])/av_me
    diff_b[nu2] =(Pme_reversed[nu2]-Pem_reversed[nu2])/av_me_r
    print(EE[nu2])
    print(Pem[nu2])
    print(Pem_reversed[nu2])
    print(Pem[nu2]-Pem_reversed[nu2])
    print(' ')
    
    diff_me_PiecewiseVac_r[nu2] = (Pme_reversed[nu2]-Pme_vac[nu2])/av_me_r
    diff_em_PiecewiseVac_r[nu2] = (Pem_reversed[nu2]-Pem_vac[nu2])/av_em_r
    diff_me_PiecewiseSingle_r[nu2] = (Pme_reversed[nu2]-Pme_single_matter[nu2])/av_me_r
    diff_em_PiecewiseSingle_r[nu2] = (Pem_reversed[nu2]-Pem_single_matter[nu2])/av_em_r

    diff_me_PiecewiseVac[nu2] = (Pme[nu2]-Pme_vac[nu2])/av_me
    diff_em_PiecewiseVac[nu2] = (Pem[nu2]-Pem_vac[nu2])/av_em
    diff_me_PiecewiseSingle[nu2] = (Pme[nu2]-Pme_single_matter[nu2])/av_me
    diff_em_PiecewiseSingle[nu2] = (Pem[nu2]-Pem_single_matter[nu2])/av_em

    diff_me_constVac[nu2] = (Pme_single_matter[nu2]-Pme_vac[nu2])/av_me_Vac
    diff_em_constVac[nu2] = ( Pem_single_matter[nu2]-Pem_vac[nu2])/av_em_Vac

    diff_single[nu2] = (Pme_single_matter[nu2]-Pem_single_matter[nu2])/av_me_Single

    diff_vac[nu2] = (Pme_vac[nu2]-Pem_vac[nu2])/av_me_Vac
#figure, axis1 = plt.subplots(1, 1)
 
    
fig1 =plt.figure(1)
plt.plot(EE, Pme_vac, color='b', label='P_me vacuum ')        
plt.plot(EE, Pem_vac, color='m',linestyle='dashed', label='P_em vacuum ')  
#plt.plot(EE, Pee_vac, color='y', label='P_ee vacuum')          
#plt.plot(EE, Pmm_vac, color='g', linestyle='dashed',label='P_mm vacuum') 


plt.title("Vacuum",fontsize=10)
plt.xlabel('                                  Energy (GeV)')
plt.legend(fontsize=10)
fig2= plt.figure(2)
plt.plot(EE, Pme_single_matter, color='g', label='P_me constant matter ')
plt.plot(EE, Pem_single_matter, color='y', linestyle='dashed',label='P_em constant matter')
#plt.plot(EE, Pee, color='m', label='P_ee ')
#plt.plot(EE, Pmm, color='b', linestyle='dashed',label='P_mm ')

plt.title("Constant Matter",fontsize=10)
plt.legend(fontsize=10)
plt.xlabel('                                  Energy (GeV)')
fig3 =plt.figure(3)
plt.plot(EE, Pme, color='g', label='P_me constant piecewise matter ')
plt.plot(EE, Pem, color='y', linestyle='dashed',label='P_em constant piecewise matter')
#plt.plot(EE, Pee, color='m', label='P_ee ')                        
#plt.plot(EE, Pmm, color='b', linestyle='dashed',label='P_mm ')     

plt.title("Increasing Piecewise Constant Matter",fontsize=10)
plt.legend(fontsize=10)
plt.xlabel('                                  Energy (GeV)')
fig4=plt.figure(4)

plt.plot(EE, Pme_reversed, color='g', label='P_me Reversed ')  
plt.plot(EE, Pem_reversed, color='y', linestyle='dashed',label='P_em Reversed') 
#plt.plot(EE, Pee_reversed, color='y', label='P_ee Reversed ')  
#plt.plot(EE, Pmm_reversed, color='g', linestyle='dashed',label='P_mm Reversed ')  

plt.title("Decreasing Piecewise Constant Matter",fontsize=10)
plt.legend(fontsize=10)
plt.xlabel('                                  Energy (GeV)')
fig5=plt.figure(5)
plt.plot(EE, diff_me_constVac,color='y',label='P_me difference between constant and vacuum ')
plt.plot(EE, diff_em_constVac,color='b',linestyle='dashed',label='P_em difference between constant and vacuum ')
plt.title("difference betweeen Vacuum and Constant Matter",fontsize=10)
plt.legend(fontsize=10)
plt.xlabel('                                  Energy (GeV)')
fig6 =plt.figure(6)
plt.plot(EE, diff, color='b',label='P_me - P_em Increasing piecewise(improper)')
plt.title("difference of P_me - P_em increasing piecewise(improper)",fontsize=10)
plt.legend(fontsize=10)
plt.xlabel('                                  Energy (GeV)')
fig7=plt.figure(7)
plt.plot(EE, diff_b, color='b',label='P_me - P_em Decreasing piecewise(proper)')
#axis6[0, 0].set_title("difference of P_me - P_em increasing piecewise",fontsize=5)
plt.title("difference of P_me - P_em decreasing piecewise(proper)",fontsize=10)
plt.legend(fontsize=10)
plt.xlabel('                                  Energy (GeV)')
fig8=plt.figure(8)
plt.plot(EE, diff_em,color='b', label='P_em difference in piecewise order')
plt.plot(EE, diff_me,color='m', label='P_me difference in piecewise order',linestyle='dashdot')
plt.title("difference of P_me and P_em for increasing to decreasing piecewise(improper vs proper)",fontsize=10)
plt.legend(fontsize=10)
plt.xlabel('                                  Energy (GeV)')
fig9=plt.figure(9)
#plt.plot(EE, diff_b,color='y',label='P_em to P_me difference in increasing piecewise potential')#,linestyle='dashed') 
#plt.plot(EE, diff,color='g',label='P_em to P_me difference in decreasing piecewise potential',linestyle='dashdot')

plt.plot(EE, diff_me_PiecewiseSingle,color='b', label='P_me difference between Increasing Piecewise and single potential')
plt.plot(EE, diff_em_PiecewiseSingle,color='y', label='P_em difference between Increasing Piecewise and single potential',linestyle='dotted')
plt.title("difference between Increasing Piecewise potential and single constant potential",fontsize=10)
plt.legend(fontsize=10)
plt.xlabel('                                  Energy (GeV)')
fig10=plt.figure(10)

plt.plot(EE, diff_me_PiecewiseSingle_r,color='b', label='P_me difference between Decreasing Piecewise and single potential')
plt.plot(EE, diff_em_PiecewiseSingle_r,color='y', label='P_em difference between Decreasing Piecewise and single potential',linestyle='dotted')

plt.title("difference between Decreasing Piecewise potential and single constant potential",fontsize=10)
plt.legend(fontsize=10)
plt.xlabel('                                  Energy (GeV)')
fig11=plt.figure(11)
plt.plot(EE, diff_me_PiecewiseVac,color='m', label='P_me difference between Increasing Piecewise potential and vacuum')
plt.plot(EE, diff_em_PiecewiseVac,color='g', label='P_em difference between Increasing Piecewise potential and vacuum',linestyle='dotted')

plt.title("difference between Increasing Piecewise potential and vacuum",fontsize=10)
plt.legend(fontsize=10)
plt.xlabel('                                  Energy (GeV)')
fig12=plt.figure(12)
#plt.plot(EE, diff_me_PiecewiseSingle_r,color='b', label='P_me difference between Piecewise and single potential Reversed')              
#plt.plot(EE, diff_em_PiecewiseSingle_r,color='y', label='P_em difference between Piecewise and single potential Reversed',linestyle='dotted')       
plt.plot(EE, diff_me_PiecewiseVac_r,color='m', label='P_me difference between Decreasing Piecewise potential and vacuum ')                 
plt.plot(EE, diff_em_PiecewiseVac_r,color='g', label='P_em difference between Decreasing Piecewise potential and vacuum',linestyle='dotted') 

plt.title("difference between Decreasing Piecewise potential and vacuum",fontsize=10)
plt.legend(fontsize=10)
plt.xlabel('                                  Energy (GeV)')
fig13 =plt.figure(13)
plt.plot(EE, Pme_vac, color='b', label='P_me vacuum ')
plt.plot(EE, Pem_vac, color='m',linestyle='dashed', label='P_em vacuum ')
plt.plot(EE, Pme_single_matter, color='g', label='P_me constant matter ')
plt.plot(EE, Pem_single_matter, color='y', linestyle='dashed',label='P_em constant matter')
plt.title("single constant matter potential and vacuum",fontsize=10)
plt.legend(fontsize=10)
plt.xlabel('                                  Energy (GeV)')
fig14 =plt.figure(14)
plt.plot(EE, Pme_vac, color='b', label='P_me vacuum ')
plt.plot(EE, Pem_vac, color='m',linestyle='dashed', label='P_em vacuum ')
plt.plot(EE, Pme, color='g', label='P_me constant piecewise matter ')
plt.plot(EE, Pem, color='y', linestyle='dashed',label='P_em constant piecewise matter')
plt.title("increasing piecewise matter potential and vacuum",fontsize=10)
plt.legend(fontsize=10)
plt.xlabel('                                  Energy (GeV)')
fig15 =plt.figure(15)
plt.plot(EE, Pme_vac, color='b', label='P_me vacuum ')
plt.plot(EE, Pem_vac, color='m',linestyle='dashed', label='P_em vacuum ')
plt.plot(EE, Pme_reversed, color='g', label='P_me Reversed')
plt.plot(EE, Pem_reversed, color='y', linestyle='dashed',label='P_em Reversed')
plt.title("decreasing piecewise matter potential and vacuum",fontsize=10)
plt.legend(fontsize=10)
#axis[0, 1].legend(fontsize=5)
#axis[0, 2].legend(fontsize=5)
#axis[1, 0].legend(fontsize=5)
#axis[1, 1].legend(fontsize=5)
#axis[1, 2].legend(fontsize=5)
#axis[2, 0].legend(fontsize=5)
#axis[2, 1].legend(fontsize=5)
#axis[2, 2].legend(fontsize=5)
#axis[0, 3].legend(fontsize=5)
#axis[1, 3].legend(fontsize=5)
#axis[2, 3].legend(fontsize=5)




#txt="V1 = 1.07*10^-5 eV^2/GeV, V2 = 3.57*10^-6 eV^2/GeV, L1 = 800 km, L2 = 1200 km         "#)# step-down ')
#plt.title('Piecewise')#Piecewise Decreasing Potential')
plt.xlabel('                                  Energy (GeV)')
#axis[0, 1].xlabel('                                  Energy (GeV)')
#axis[0, 2].xlabel('                                  Energy (GeV)')
#axis[1, 0].xlabel('                                  Energy (GeV)')
#axis[1, 1].xlabel('                                  Energy (GeV)')
#axis[1, 2].xlabel('                                  Energy (GeV)')
#plt.text(0.6,1,txt,fontsize = 7)
plt.show()

fig1.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop/2_flavor_oscillations_10^-4V/Vacuum.jpg')
fig2.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop/2_flavor_oscillations_10^-4V/Constant_matter.jpg')
fig3.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop/2_flavor_oscillations_10^-4V/Increasing_Piecewise_Constant_Matter.jpg')
fig4.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop/2_flavor_oscillations_10^-4V/Decreasing_Piecewise_Constant_Matter.jpg')
fig5.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop/2_flavor_oscillations_10^-4V/difference_Vacuum_Constant_Matter.jpg')
fig6.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop/2_flavor_oscillations_10^-4V/P_me_minus_P_em_increasing_piecewise.jpg')
fig7.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop/2_flavor_oscillations_10^-4V/P_me_minus_P_em_decreasing_piecewise.jpg')
fig8.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop/2_flavor_oscillations_10^-4V/P_me_minus_P_em_increasing_to_decreasing_piecewise.jpg')
fig9.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop/2_flavor_oscillations_10^-4V/difference_Increasing_Piecewise_single_constant_potential.jpg')
fig10.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop/2_flavor_oscillations_10^-4V/difference_Decreasing_Piecewise_single_constant_potential.jpg')
fig11.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop/2_flavor_oscillations_10^-4V/difference_Increasing_Piecewise_vacuum.jpg')
fig12.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop/2_flavor_oscillations_10^-4V/difference_Decreasing_Piecewise_vacuum.jpg')
fig13.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop/2_flavor_oscillations_10^-4V/Vac_single_matter.jpg')
fig14.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop/2_flavor_oscillations_10^-4V/Vac_increasing_piecewise.jpg')
fig15.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop/2_flavor_oscillations_10^-4V/Vac_decreasing_piecewise.jpg')
              

"""
fig16 =plt.figure(16)                                                                                                         
plt.plot(EE, U1_modsqd01, color='b', label='|U1(0,1)|^2 ')
plt.plot(EE, U1_modsqd10, color='m',linestyle='dashed', label='|U1(1,0)|^2 ')                                                  
plt.plot(EE, U2_modsqd01, color='g', label='|U2(0,1)|^2 ')                                                 
plt.plot(EE, U2_modsqd10, color='y', linestyle='dashed',label='|U2(1,0)|^2')                               
plt.title("mod squared terms: 10-5V",fontsize=10)                                                      
plt.legend(fontsize=10)                                                                                                        
plt.xlabel('                                  Energy (GeV)')


fig17 =plt.figure(17)
plt.plot(EE, rUT00_con, color='b', label='U2U1(0,0) conjugate ')
plt.plot(EE, rUT00R, color='r', linestyle='dashed',label='U1U2(0,0)  ')

plt.plot(EE, rUT01_con, color='k', label='U2U1(0,1) conjugate ')
plt.plot(EE, rUT01R, color='y', linestyle='dashdot',label='U1U2(0,1)  ')

plt.plot(EE, rUT10_con, color='g', label='U2U1(1,0) conjugate ')
plt.plot(EE, rUT10R, color='b', linestyle='dashed',label='U1U2(1,0)  ')

plt.plot(EE, rUT11_con, color='m', label='U2U1(1,1) conjugate ')
plt.plot(EE, rUT11R, color='c', linestyle='dashdot',label='U1U2(1,1)  ')

plt.title("conjugate check(real part): 10-5V",fontsize=10)
plt.legend(fontsize=10)
plt.xlabel('                                  Energy (GeV)')

fig18 =plt.figure(18)
plt.plot(EE, iUT00_con, color='b', label='U2U1(0,0) conjugate ')
plt.plot(EE, iUT00R, color='r', linestyle='dashed',label='U1U2(0,0)  ')

plt.plot(EE, iUT01_con, color='k', label='U2U1(0,1) conjugate ')
plt.plot(EE, iUT01R, color='y', linestyle='dashdot',label='U1U2(0,1)  ')

plt.plot(EE, iUT10_con, color='g', label='U2U1(1,0) conjugate ')
plt.plot(EE, iUT10R, color='b', linestyle='dashed',label='U1U2(1,0)  ')

plt.plot(EE, iUT11_con, color='m', label='U2U1(1,1) conjugate ')
plt.plot(EE, iUT11R, color='c', linestyle='dashdot',label='U1U2(1,1)  ')

plt.title("conjugate check(imag part): 10-5V",fontsize=10)
plt.legend(fontsize=10)
plt.xlabel('                                  Energy (GeV)')



fig19 =plt.figure(19)
plt.plot(EE, rUT00, color='b', label='U2U1(0,0)  ')
plt.plot(EE, rUT00R, color='r', linestyle='dashdot',label='U1U2(0,0)  ')

plt.plot(EE, rUT01, color='y', label='U2U1(0,1)  ')
plt.plot(EE, rUT10R, color='b', linestyle='dashed',label='U1U2(1,0)  ')

plt.plot(EE, rUT10, color='b', label='U2U1(1,0)  ')
plt.plot(EE, rUT01R, color='r', linestyle='dashed',label='U1U2(0,1)  ')

plt.plot(EE, rUT11, color='y', label='U2U1(1,1)  ')
plt.plot(EE, rUT11R, color='b', linestyle='dashdot',label='U1U2(1,1)  ')

plt.title("transpose check(real part): 10-5V",fontsize=10)
plt.legend(fontsize=10)
plt.xlabel('                                  Energy (GeV)')


fig20 =plt.figure(20)
plt.plot(EE, iUT00, color='b', label='U2U1(0,0)  ')
plt.plot(EE, iUT00R, color='r', linestyle='dashdot',label='U1U2(0,0)  ')

plt.plot(EE, iUT01, color='y', label='U2U1(0,1)  ')
plt.plot(EE, iUT10R, color='b', linestyle='dashed',label='U1U2(1,0)  ')

plt.plot(EE, iUT10, color='b', label='U2U1(1,0)  ')
plt.plot(EE, iUT01R, color='r', linestyle='dashed',label='U1U2(0,1)  ')

plt.plot(EE, iUT11, color='y', label='U2U1(1,1)  ')
plt.plot(EE, iUT11R, color='b', linestyle='dashdot',label='U1U2(1,1)  ')

plt.title("transpose check(imag part): 10-5V",fontsize=10)
plt.legend(fontsize=10)
plt.xlabel('                                  Energy (GeV)')


#plt.plot(EE, UT10_con, color='r',linestyle='dashed', label='U2U1(1,0) conjugate ')
#plt.plot(EE, UT01_con, color='g'linestyle='dotted', label='U2U1(0,1) conjugate ')
#plt.plot(EE, UT11_con, color='y',linestyle='dashdot', label='U2U1(1,1) conjugate ')
plt.show()  

fig16.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop/2_flavor_oscillations_10^-4V/MOD_check.jpg')
fig17.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop/2_flavor_oscillations_10^-4V/CONreal_check.jpg')
fig18.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop/2_flavor_oscillations_10^-4V/CONimag_check.jpg')
fig19.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop/2_flavor_oscillations_10^-4V/TRAreal_check.jpg')
fig20.savefig('/Users/oliviabitter/Desktop/Nu_plots/time_invar_propVSimprop/2_flavor_oscillations_10^-4V/Traimage_check.jpg')
"""
