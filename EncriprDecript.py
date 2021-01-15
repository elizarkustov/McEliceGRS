import numpy as np
import time
import sys
import random 
def random_nonsingular_matrix(size,base_ring=GF(2)):
    V = base_ring**size
    vectors = []
    for i in range(size):
        v = V.random_element()
        while v in V.span(vectors):
            v = V.random_element()
        vectors.append(v)
    return(matrix(vectors))
    
    def McEliceGRS(m):
    start_time = time.time()
    m = m
    Limd = m // 2
    F = GF(2^m, 'a') 
    E = GF(2)
    n, k = 2^m, (2^m) // 3
    C = codes.GeneralizedReedSolomonCode(F.list()[:n], k)
    D = codes.decoders.GRSBerlekampWelchDecoder(C)
    H = C.parity_check_matrix()
    Srand = []
    l = list(range(m-Limd))
    random.shuffle(l)
    for p in range(n):
        Srand.append(l)
    Hisom =[]
    a = F('a')
    pkH = (sys.getsizeof(H))
    for i in range(0,H.nrows()):
        Hrows = []
        for j in range(0,H.ncols()):
            for z in range(0,m):
                if z not in Srand[j]:
                    Hrows.append(  (a^z) * H[i][j] )
        Hisom.append(Hrows)
    HisomMat = matrix(F,Hisom)
    Hisom =[]
    HisomCheck =[]
    V, fr, to  = F.vector_space(GF(2), map=True)
    for i in range(0,HisomMat.ncols()):
        Hcols = []
        HcolsCheck = []
        for j in range(0,HisomMat.nrows()):
        
            for z in range(0,m):
                HcolsCheck.append(to(HisomMat[j][i])[z])
        HisomCheck.append(HcolsCheck)
    HisomMatCheck2 = matrix(GF(2),HisomCheck)
    HisomMatCheck2 = HisomMatCheck2.transpose()
    HGRE = []
    for i in range(n):
       HGRE.append(random_nonsingular_matrix(Limd))
    T = block_diagonal_matrix(HGRE)
    T = matrix(GF(2),T)
    Sigma = matrix(GF(2),n,[random.random()<0.5 for _ in range(n^2)])
    while (rank(Sigma) < n) :
        Sigma[floor(n*random.random()),floor(n*random.random())] +=1
    iden =  matrix.identity(Limd)
    Psigma = np.kron(Sigma, iden)
    Psigma = matrix(GF(2),Psigma)
    Q =  T * Psigma
    Q = matrix(GF(2),Q)   
    Hfinal=HisomMatCheck2 * Q
    ifs=HisomMatCheck2 * Q
    Hfinal = matrix(GF(2),HisomMatCheck2)
    ifs=HisomMatCheck2 * Q
    wy = [0] * (n*Limd)
    wy[8] = 1
    wy[18] = 1
    word = vector( F, wy)
    word = vector( E,  wy )    
    cip = Hfinal * word
    i=0
    isocip = []
    while i < len(cip):
        isocip.append(fr(cip[0+i:m+i]))
        i = i +m
    C2=[0] * k   
    C3=[]
    C3=isocip
    C3.extend(C2)
    cod = vector(F, C3)
    mes = D.decode_to_message(cod)
    koef = mes.list()
    koef.extend([0]*(n-k))   
    isokoef = []
    for i in range(0, len(koef)):
        isokoef.extend(to(koef[i]))
    Sranddel = []    
    for i in range(0,n):
        for j in range(0,m-Limd):
            Sranddel.append(Srand[i][j] + i*(m-1))
            
    isokoef[:] = [x for i,x in enumerate(isokoef) if i not in Sranddel]
    mesag = matrix(GF(2),isokoef)
    Qinv = Q.inverse()
    final = mesag * Qinv   


