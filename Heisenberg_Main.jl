using Distributions, Base.Threads, Dierckx, LsqFit, Random, JLD2   # To have random distrib of spin, parallelize, interpolate, Random, save data in compact file
include("./Heisenberg_Main_functions.jl")

N=20000
burn = Int(min(N/4,100000))
L=[8,10]
T=Float32.([.2,.7,1.1,1.6])
d=Float32.([-1000,0,1000])
pIsing = .6f0
pXY = .7f0
burn = Int(min(N/4,100000))     # Burning period
PBC   = true                    # Periodic Boundary Conditions
pi32  = Float32(π)
Save  = true

println("  --  $N sweeps  --  ")
starttime = time()
for z in d
    for l in L
        E1, M1 = MH_parallel_tempering(l, N, T, burn, pi32, L, z, pIsing, pXY, Save)
    end
end
t = round(Int, time()-starttime)
if Save == true
    open("Data/$(L)_$N/elapsed_time.txt", "a") do file; write(file, "$t\n$d\n$T\n\n"); end
end

println(t)




# Katzgraber, MH parallel, with right temperature (and eq19 later)

# C, Susc and Corr max