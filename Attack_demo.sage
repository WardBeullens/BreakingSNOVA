import random
import time    
import os

load("SNOVA.sage")

## subroutines ##

# find coefficients of R such that tilde E_R has target rank
def find_minrank_R(params, pk, target_rank):
    v, o, q , l = params
    P, randomness = pk
    A,B,Q1,Q2 = randomness

    # compute quadratic forms in E_tilde 
    matrices = [[ zero_matrix(K,l*l,l*l) for _ in range (l*l) ] for _ in range(l*l) ]

    for alpha in range(l*l):
        for i in range(l*l):
            for j in range(l*l):
                Lc = FS_coeffs(Q1[alpha]*S[l]^(i%l),l)
                Rc = FS_coeffs(Q2[alpha]*S[l]^(j%l),l)

                V = [ lce*rce for lce in Lc for rce in Rc ]
                U = [ lce*rce for lce in A[alpha].column(i//l) for rce in B[alpha].row(j//l)]

                V = Matrix(K, 1, l*l, V)
                U = Matrix(K, l*l, 1, U)

                if i <= j:
                    matrices[i][j] += U*V
                else:
                    matrices[j][i] += U*V

    # coefficients of the Ri
    R = [1,0,0,0]
    for a in K:
        R[2] = a
        for b in K:
            R[3] = b
            E_tilde = zero_matrix(K, 4, 4)
            for i in range(4):
                for j in range(i-1,4):
                    E_tilde += R[i]*R[j]*matrices[i][j]
            if rank(E_tilde) <= target_rank:
                return R, E_tilde

# check if first nonzero entry is equal to 1
def is_normalized(vec):
    for x in vec:
        if x == K(1):
            return True
        if x == K(0):
            continue
        return False
        
# extend linearly independent set to a basis
def completementary_basis(partial_basis, space):
    S = span(partial_basis)
    d = space.dimension()
    tail = []
    for x in space.basis():
        if not x in S:
            tail.append(x)
            if len(partial_basis) + len(tail) == d:
                break
            else:
                S = span(partial_basis+tail)
    return tail

def _Hashimoto(Q):
    B = [ q + q.transpose() for q in Q]
    n = Q[0].ncols()

    VS = VectorSpace(K,n)
    r = VS.random_element()
    found_any_vector = False
    for o_minus_r in VS:
        o = o_minus_r + r
        if not is_normalized(o):
            continue
        
        if all([ o*q*o == K(0) for q in Q ]):
            M = Matrix(K, [ o*b for b in B])
            subspace = M.right_kernel().intersection(span(completementary_basis([o], VS)))
            basis = Matrix(K,subspace.basis()).transpose()

            restricted_Q = [ basis.transpose()*q*basis for q in Q ]

            for O in _Hashimoto(restricted_Q):
                yield block_matrix(1,2, [basis*O, Matrix(K,n,1, o) ])

            if found_any_vector:
                return
            found_any_vector = True
    
    if not found_any_vector:
        yield Matrix(K, n, 0, [ ] )



def Hashimoto(Q, L, C, number_of_extra_guesses = 1, eliminate = 3):
        
    # homogenize
    n = Q[0].ncols()
    Q = [ block_matrix(2,2,[ Matrix(K,1,1,[C[i]]), Matrix(K,1,n,L[i]), zero_matrix(K,n,1), Q[i] ]) for i in range(len(Q)) ]
    n = n+1

    os.system("rm -f system.out")

    attempts = 0
    for O in _Hashimoto(Q[:eliminate]):
    
        attempts += 1
        O_dim = O.ncols()
        PR = PolynomialRing(K, O_dim , "S")
        vars = vector(PR,PR.gens())

        print(f"Attempt {attempts}:")
        print("Write system for m4gb")

        # write file for m4gb
        with open('system.in', 'w') as f:
            f.write('$fieldsize 16\n')
            f.write(f'$vars {O_dim-1} S\n')    
            for q in Q[eliminate:]:
                equation = vars*O.transpose()*q*O*vars      
                line = ""
                for c,M in equation:
                    line += str(invese_field_dict[c])+"*"+str(M)+" + "
                line = line[:-2]
                line = line.replace("*"+str(vars[-1])+"^2","")
                line = line.replace("*"+str(vars[-1]),"")
                f.write(line+"\n")

        print(f"calling m4gb on system of {len(Q[eliminate:])} quadratic equations in {O_dim-1} variables.")
        T = time.time()
        os.system("cd m4gb; ./solver.sh ../system.in ../system.out > log.txt 2>&1 ")
        print(f"m4gb took {(time.time() - T):.2f} seconds.")

        if os.path.isfile("system.out"):
            print("solution found")

            # read solution:
            sol = []

            with open("system.out") as file:
                lines = [line.rstrip() for line in file]

                for line in lines:
                    if "+" in line:
                        sol.append(field_dict[ int(line.split()[-1]) ])
                    else:
                        sol.append(K(0)) 

            sol.reverse()
            sol += [K(1)]
            Osol = O*vector(K,sol)

            if Osol[0] == K(0):
                print("solution was at infinity :(")
                continue
            
            Osol = Osol[1:]/Osol[0]
            return Osol

        else:
            print("no solution\n")

## demonstration of attack

params = (37, 17, 16, 2)
v, o, q , l = params
m = o
n = v + o
target_rank = 1

# sample a random hash to forge a signature for
target = vector(K, [K.random_element() for _ in range(m*l*l)])
print("randomly chosen target (represents the hash of a salted message):")
print(target)

start_time = int(time.time())
total_start_time = start_time

pk = read_pk('pk.txt', params)
P, randomness = pk
A,B,Q1,Q2 = randomness

print(f"reading public key took {(time.time() - start_time):.2f} seconds.\n")
start_time = int(time.time())

R, E_tilde = find_minrank_R(params, pk, target_rank)

if R is None:
    print("Attack does not work on this public key")
    exit(0)

scalars = R
E_tilde
print("Minrank Solution found:", scalars)
print("E_tilde:")
print(E_tilde)

print(f"finding E_tilde of rank {target_rank} took {(time.time() - start_time):.2f} seconds.\n")
start_time = int(time.time())

print("composing system")

# hack to get transformation matrix T to put E_tilde in echelon form
MT = block_matrix([[E_tilde, identity_matrix(l*l)]])
MT.echelonize()
M2 = MT[:,:l*l]
T = MT[:,l*l:]

# find quadratic part of public key restricted to our affine subspace
quadratic_part_1 = [ zero_matrix(K, n*l, n*l) for _ in range(m*target_rank) ]

SS = identity_matrix(n).tensor_product(S[l])
powers_of_SS = [ SS^i for i in range(l) ]

for k,p in enumerate(P):
    twists = [  ]
    for i in range(l):
        for j in range(l):
            twists.append( powers_of_SS[i]*p*powers_of_SS[j] )

    for i in range(target_rank):
        for j in range(l*l):
            quadratic_part_1[k*target_rank + i] += M2[i][j] * twists[j]


V = Matrix(K,n*l,l,[K.random_element() for _ in range(n*l*l)])

TT_1 = identity_matrix(K,o).tensor_product(T[:target_rank,:])
TT_2 = identity_matrix(K,o).tensor_product(T[target_rank:,:])

# constant part 
constant_part = evaluate_public_map(params, pk, V) - target
constant_part_1 = TT_1*constant_part
constant_part_2 = TT_2*constant_part

L = evaluate_public_map_derivative(params, pk, V, scalars)
linear_part_1 = TT_1*L
linear_part_2 = TT_2*L

## specialize to get rid of linear equations

number_of_variables_after_specializtion = n*l - (l*l-target_rank)*m 

Kernel = linear_part_2.right_kernel().basis()
Kernel = Matrix(K, Kernel).transpose()
offset2 = linear_part_2.solve_right(constant_part_2)

# specialize quadratic part
quadratic_part_specialized = [ Kernel.transpose()*Q*Kernel for Q in quadratic_part_1]

#specialize linear part
linear_part_specialized = linear_part_1*Kernel

L = [ vector(K,offset2*(Q+Q.transpose())*Kernel) for Q in quadratic_part_1]
L = Matrix(K,L)

linear_part_specialized = linear_part_specialized + L

# specialize constant part
constant_part_specialized = constant_part_1 + linear_part_1*offset2 + vector(K, [ offset2 * Q* offset2 for Q in quadratic_part_1 ])


print(f'composing system "P**(y) = t**" took {(time.time() - start_time):.2f} seconds.\n')
start_time = int(time.time())

print("start solving system with Hashimoto's algorithm for underdetermined systems.\n")
sol = Hashimoto(quadratic_part_specialized, linear_part_specialized, constant_part_specialized, eliminate = 3)

print(f'finding solution to "P**(y) = t**" took {(time.time() - start_time):.2f} seconds.\n')
start_time = int(time.time())

u_0 = Kernel*sol + offset2

# compute U
R1 = R[2]*identity_matrix(2) + R[3]*S[2]
u_1 = identity_matrix(n).tensor_product(R1)*u_0
U = Matrix(K, [u_0,u_1]).transpose() + V

print(f'finishing signature took {(time.time() - start_time):.2f} seconds.\n')
start_time = int(time.time())

# check if signature is valid
eval = evaluate_public_map(params, pk, U)
print("evaluation of public map at signature:")
print(eval)
if eval == target:
    print("forgery success!")
else:
    print("Oops, forgery is invalid. Something went wrong.")



print(f'total attack time is {(time.time() - total_start_time):.2f} seconds.\n')