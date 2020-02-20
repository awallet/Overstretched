import time
from fpylll import *
from fpylll.algorithms.bkz2 import BKZReduction as BKZ2
from sage.modules.free_module_integer import IntegerLattice

def CenteredMod(a, q):
    a = a.mod(q)
    if a <= floor(q/2):
        return a
    else:
        return a-q

PATH = "SET PATH OF DATA HERE"
    
# parameters
n = 2^10
q = 2^42
# precomputed data; this is for reproducibility 
print("Loading secret key...")
S = load(PATH+"/sk.sobj") # secret key
S_ = load(PATH+"/skinv.sobj") # inverse mod q

print ("Loading gadget...")
G = load(PATH+"/gadget.sobj") # gadget matrix

## If you want to test on your own data, uncomment below
# Zq = ZZ.quo(q)
# print("build secret key...")
# while True:
#     print("Try...")
#     S0 = Matrix(Zq, [ [ randint(0,1) for _ in range(n) ] for _ in range(n) ])
#     if S0.is_invertible():
#         S = S0.change_ring(ZZ)
#         break
# S_ = S0.inverse().change_ring(ZZ) #S^-1 mod q seen in ZZ
# gadget matrix
# G = block_matrix( [ [ 2^i*identity_matrix(n) for i in range(log(q,2)) ] ] )

### BUILD PUBLIC DATA ###
print("Computing cipher of Id...")
TIME = time.time()
E_Id = random_matrix(ZZ, n, log(q,2)*n, x=0, y=2) # E is secret encryption noise
C_Id = S_ * (E_Id[:n,:n] - G[:n, :n]) 
C_Id = Matrix(ZZ, [ [ CenteredMod(e, q) for e in C_Id[i] ] for i in range(n) ] ) #first block of cipher
print("Done. time: ", round(time.time() - TIME)) # takes ~20s

print("Building basis of NTRU lattice...")
TIME = time.time()
B = block_matrix([ [ q*identity_matrix(n), zero_matrix(n) ], [ C_Id, identity_matrix(n) ] ] )
B = IntegerMatrix.from_matrix(B) # type conversion, between 3-5 mins
print("Done. time: ", round(time.time() - TIME))

### PREPROC PHASE FOR FPYLLL ###
k = 280 
print("Building submatrix of size 2*k =", 2*k, "...")
TIME = time.time()
Bk = B.submatrix(range(n-k,n+k), range(n-k,n+k)) # ~30s
FPLLL.set_precision(180) 
GSO_Bk = GSO.Mat(Bk, float_type='mpfr')
print("Preprocessing...")
Bk_BKZ = BKZ2(GSO_Bk) 
print("Done. time: ", round(time.time() - TIME))   

### BKZ PHASE ###
flags = BKZ.AUTO_ABORT|BKZ.MAX_LOOPS|BKZ.GH_BND|BKZ.VERBOSE
beta = 20
par = BKZ.Param(block_size=beta, strategies=BKZ.DEFAULT_STRATEGY, flags=flags, max_loops=4)
print("BKZ reduction with beta =", beta, "...")
TIME = time.time()
DONE = Bk_BKZ(par) #actual BKZ algorithm; updates Bk in place; ~15 hours
# if it fails because infinite loop in babai, set higher precision in FPLLL.set_precision() and rerun without restarting.
print("Done. time: ", round(time.time() - TIME))

### WE WON ###
print( all(Bk[i].norm() < 2^24 for i in range(k)) )


