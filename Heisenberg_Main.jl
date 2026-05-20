using Distributions, Base.Threads, Random, JLD2, StaticArrays, Base.Filesystem   # To have random distrib of spin, parallelize, interpolate, Random, save data in compact file, to use SVector which are optimized, to create folder
include("./Heisenberg_Main_functions.jl")

const N=20000                         # Number of lattice sweeps
const burn = Int(min(N/4,100000))     # Burning period
const L=[8,12]                        # Lattice sizes
const T=Temperatures(.55,1.15,32)     # Temperatures
const D=Float32.([-1,1])              # Anisotropic term
const PBC   = true                    # Periodic Boundary Conditions
const Save  = true                    # To save Data in a folder named L_N
const SaveLattices = false            # To save the (199 by default) last lattices of each chain

println("  --  $N sweeps  --  ")
starttime = time()
if !isdir("Data/$(L)_$N")   # create the folder if it does not exist
    mkpath("Data/$(L)_$N")
end
for d in D
    for l in L
        E1_end = MH_parallel_tempering(l, N, T, burn, L, d, Save, SaveLattices)
    end
end

t = round(Int, time()-starttime)    # record and save time
if Save == true
    open("Data/$(L)_$N/elapsed_time.txt", "a") do file
        write(file, "$t\n$D\n$T\n\n")
    end
end
println(t)