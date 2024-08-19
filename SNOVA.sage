

S = {}

K.<x> = GF(16)

two = x
three = x + 1
four = x*x
five = x*x + 1
six = x*x + x
seven = x*x + x + 1
eight = x*x*x

S[2] = Matrix(K, [[eight, seven],
                  [seven, six]])

S[3] = Matrix(K, [[eight, seven, six],
                  [seven, six  , five],
                  [six,   five, four]])

S[4] = Matrix(K, [[eight, seven, six,   five],
                  [seven, six,   five,  four],
                  [six,   five,  four,  three],
                  [five,  four,  three, two]])

def sample_from_FS(l):
    Out = zero_matrix(K,l,l)
    while Out == 0:
        for _ in range(l):
            Out *= S[l]
            Out += K.random_element()*identity_matrix(K,l)
    
    return Out

def FS_coeffs(M, l):
    MM = Matrix(K, [ (S[l]**i).row(0) for i in range(l) ] ).transpose()
    return MM.inverse()*M.row(0)

def coeffs_to_FS(cfs,l):
    return sum([cfs[i]*S[l]**i for i in range(l)])

def sample_block(l):
    return Matrix(K, [[ K.random_element() for _ in range(l) ] for _ in range(l)] )

def sample_block_invertible(l):
    B = sample_block(l)
    while not B.is_invertible():
        B += K.random_element()*S[l]
    return B

def KeyGen(params):
    v, o, q , l = params

    T12 = block_matrix(K, [ [sample_from_FS(l) for _ in range(o)] for _ in range(v) ] )

    P11 = [  block_matrix(K, [ [sample_block(l) for _ in range(v)] for _ in range(v) ] ) for _ in range(m) ]
    P12 = [  block_matrix(K, [ [sample_block(l) for _ in range(o)] for _ in range(v) ] ) for _ in range(m) ]
    P21 = [  block_matrix(K, [ [sample_block(l) for _ in range(v)] for _ in range(o) ] ) for _ in range(m) ]
    P22 = [ T12.transpose()*(p11*T12 + p12) + p21*T12  for (p11,p12,p21) in zip(P11,P12,P21) ]

    P = [ block_matrix([[P11[i],P12[i]],[P21[i],P22[i]]]) for i in range(m) ]
    T = block_matrix([[T12],[identity_matrix(K,o*l)]])

    A = [ sample_block_invertible(l) for _ in range(l*l) ] 
    B = [ sample_block_invertible(l) for _ in range(l*l) ] 
    Q1 = [ sample_from_FS(l) for _ in range(l*l) ] 
    Q2 = [ sample_from_FS(l) for _ in range(l*l) ] 

    randomness = (A,B,Q1,Q2)

    return ((P, randomness), T)

def evaluate_public_map( params, pk, U ):
    v, o, q , l = params
    n = v + o
    P, randomness = pk
    A,B,Q1,Q2 = randomness

    out = [ zero_matrix(K,l,l) for _ in range(o)]

    for alpha in range(l*l):
        QQ1 = identity_matrix(n).tensor_product(Q1[alpha])
        QQ2 = identity_matrix(n).tensor_product(Q2[alpha])
        Left = A[alpha]*U.transpose()*QQ1
        Right = QQ2*U*B[alpha]
        for i in range(o):
            out[i] += Left*P[i]*Right
        
    out = block_matrix(m,1, out)
    out = out.list()
    out = vector(K,out)

    return out

    
def evaluate_public_map_derivative( params, pk, U, scalars):
    v, o, q , l = params
    n = v + o
    P, randomness = pk
    A,B,Q1,Q2 = randomness

    PR2 = PolynomialRing(K, n*l ,'t')

    out = [ zero_matrix(PR2,l,l) for _ in range(o)]

    vars = PR2.gens()
    U2 = [ vector(PR2, vars) ]

    for j in range(1,l):
        a = scalars[l*j:l*(j+1)]
        Scalar = identity_matrix(n).tensor_product(coeffs_to_FS( a ,l))
        U2.append(Scalar*U2[0])

    U2 = Matrix(PR2, U2).transpose()


    for alpha in range(l*l):
        QQ1 = identity_matrix(n).tensor_product(Q1[alpha])
        QQ2 = identity_matrix(n).tensor_product(Q2[alpha])
        Left = A[alpha]*U.transpose()*QQ1
        Right = QQ2*U*B[alpha]
        Left2 = A[alpha]*U2.transpose()*QQ1
        Right2 = QQ2*U2*B[alpha]
        for i in range(o):
            out[i] += Left*P[i]*Right2 + Left2*P[i]*Right
        
    out = block_matrix(m,1, out)
    out = out.list()
    out = vector(PR2,out)

    D = { x:i for i,x in enumerate(vars) }

    OUT = zero_matrix(K, len(out), n*l)
    for i in range(len(out)):
        for x,y in out[i]:
            OUT[i,D[y]] = x

    return OUT

## reading pk from file ##

field_dict = {};
field_dict[0] = K(0)
field_dict[1] = K(1)
field_dict[2] = K(x)
field_dict[3] = K(x+1)
field_dict[4] = K(x*x)
field_dict[5] = K(x*x+1)
field_dict[6] = K(x*x+x)
field_dict[7] = K(x*x+x+1)
field_dict[8] = K(x**3)
field_dict[9] = K(x**3 + 1)
field_dict[10] = K(x**3 + x)
field_dict[11] = K(x**3 + x+1)
field_dict[12] = K(x**3 + x*x)
field_dict[13] = K(x**3 + x*x+1)
field_dict[14] = K(x**3 + x*x+x)
field_dict[15] = K(x**3 + x*x+x+1)
invese_field_dict = { field_dict[x]:x for x in field_dict }

def read_sequence_of_matrices(string, l):
    L = list(map(lambda x : field_dict[int(x)] , string.split()))
    Matrices = []

    if len(L) != l**4:
        print("wrong length")
        exit()

    return [ Matrix(K, l,l, L[alpha*l*l:(alpha+1)*l*l]) for alpha in range(l*l) ]

def read_bilinear_form(string, v,o,l):

    L = list(map(lambda x : field_dict[int(x)] , string.split()))

    if len(L) != ((v+o)**2)*(l**2):
        print("wrong length")
        exit()

    ML = [ Matrix(K,l,l, L[i*l*l:(i+1)*l*l]) for i in range((v+o)**2) ]

    P11 = ML[:v*v]
    P12 = ML[v*v:v*v+v*o]
    P21 = ML[v*v+v*o:v*v+2*v*o]
    P22 = ML[v*v+2*v*o:]

    P11 = block_matrix(v,v, P11)
    P12 = block_matrix(v,o, P12)
    P21 = block_matrix(o,v, P21)
    P22 = block_matrix(o,o, P22)

    P = block_matrix([[P11,P12],[P21,P22]])
    return P


def read_pk(filename, params, sanity_check = False):
    v, o, q , l = params

    with open(filename) as file:
        lines = [line.rstrip() for line in file]

    A = read_sequence_of_matrices(lines[1], l)
    B = read_sequence_of_matrices(lines[3], l)
    Q1 = read_sequence_of_matrices(lines[5], l)
    Q2 = read_sequence_of_matrices(lines[7], l)
    
    randomness = (A,B,Q1,Q2)
    P = [ read_bilinear_form(str, v,o,l) for str in lines[9:9+o] ]

    return (P, randomness)


    
