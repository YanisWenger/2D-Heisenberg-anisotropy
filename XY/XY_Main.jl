using Plots, Distributions, Colors, CSV, DataFrames, Base.Threads, Dierckx, LsqFit, Random, JLD2   # To plot, To have random distrib of spin, colors for lattices, write/read csv (CSV & DataFrames), parallelize, interpolate, save data in compact file
include("./XY_functions.jl")

starttime = time()
N = 100_000 # Lattice flip
# T = unique(sort(Float32.(vcat(collect(0.955:.015:1.21),collect(.1:.1:.5),collect(.56:.04:.92),collect(1.24:.04:1.64),collect(1.7:.1:2)))))
T = Float32.(collect(.7:.1:1.4))            # Temperature
L = [8,20]           # Lattice size
# L = [8,12,20,32,50,70,100]
start = Int(min(N/4,100000))    # Burning period
PBC   = true                        # Periodic Boundary Conditions
twopi = Float32(2π)
Save  = true

# --------------------------------------------------------------- #

println("  --  $N sweeps  --  ")
for l in eachindex(L)
    Threads.@threads for t in eachindex(T)
        Energy = MH(rand(Float32, L[l], L[l])*twopi, N, T[t], start, twopi, Save, L, T[1]) # zeros(Float32, L[l],L[l])
        println("Energy at $(T[t]) in $(L[l]) lattice: \t", Energy)
    end
end
t = round(Int, time()-starttime)
if Save == true
    open("Data/$(L)_$(T[1])_$N/elapsed_time.txt", "w") do file; write(file, "$t\n"); end
end
println(t)

# a,b,c,d,e,f,g,h,i,j,k = MHvideo(rand(Float32, 40, 40)*twopi, 2000, .2, Nbin, 200, twopi, PBC) # zeros(L[l],L[l])
# println("    E = ", round(a;digits=4), " ± ", round(b;digits=4), "    Mag = ", round(c;digits=4), " ± ", round(d;digits=4), "    Suscept: ", round(e;digits=5), " ± ", round(f;digits=6), "    C = ", round(g;digits=7), " ± ", round(h;digits=7), "    accept ", round(k;digits=3))
# COLOR MAP !!!


#=
200000   48T   8,10,12,14,16,20    58min
200000   48T   8,12,20             29min   -> 18'

200_000   8T   8                    17 (skip 10), 21 (skip 1)
1_000_000 8T   8                    1'33 (skip=10)
=#