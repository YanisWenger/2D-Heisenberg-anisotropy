using Distributions, Base.Threads, Random, JLD2   # To have random distrib of spin, parallelize, interpolate, Random, save data in compact file
include("./Heisenberg_Main_functions.jl")

N=20000                         # Number of lattice sweeps
burn = Int(min(N/4,100000))     # Burning period
L=[8,12]                        # Lattice sizes
T=Temperatures(.65,1.15,8)        # Temperatures
d=Float32.([-1,1])              # anisotropic term
PBC   = true                    # Periodic Boundary Conditions
pi32  = Float32(π)
Save  = true                    # to save Data in a folder named L_N (that needs to be created before, i.e.: [8, 12, 20]_1000000

println("  --  $N sweeps  --  ")
starttime = time()
for z in d
    for l in L
        E1, M1 = MH_parallel_tempering(l, N, T, burn, pi32, L, z, Save)
    end
end
t = round(Int, time()-starttime)
if Save == true
    open("Data/$(L)_$N/elapsed_time.txt", "a") do file; write(file, "$t\n$d\n$T\n\n"); end
end

println(t)


# Cerror : error propag from Eerror or binning ?