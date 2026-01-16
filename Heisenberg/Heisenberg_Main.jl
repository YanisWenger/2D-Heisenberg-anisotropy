using Plots, Distributions, Colors, Base.Threads, Dierckx, LsqFit, Random, JLD2   # To plot, To have random distrib of spin, colors for lattices, parallelize, interpolate, Random, save data in compact file
include("./Heisenberg_functions.jl")

N = 10_000 # Lattice flip
# T = unique(sort(Float32.(vcat(collect(0.955:.015:1.21),collect(.1:.1:.5),collect(.56:.04:.92),collect(1.24:.04:1.64),collect(1.7:.1:2))))) # XY
T = Float32.(collect(.7:.1:1.4))            # Temperature
# T = unique(sort(Float32.(vcat(collect(.1:.1:1.7),collect(.6:.05:1.3)))))
L = [8]           # Lattice size
# L = [8,12,20,32,50,70,100]
# Δ    = [-1,-.25,0,.25,1]
Δ = [.01]
p = .5
start = Int(min(N/4,100000))    # Burning period
PBC   = true                        # Periodic Boundary Conditions
pi32  = Float32(π)
Save  = false
starttime = time()

# --------------------------------------------------------------- #

println("  --  $N sweeps  --  ")
for z in eachindex(Δ)
    for l in eachindex(L)
        Threads.@threads for t in eachindex(T)
            Energy = MH(Initial_lattice(L[l], pi32), L[l], N, T[t], start, pi32, L, Δ[z], .5, Save) # zeros(Float32, L[l],L[l])
            println("T = $(T[t]) \tL = $(L[l]) Δz = $(Δ[z]) \tE = ", Energy)
        end
    end
end
t = round(Int, time()-starttime)
if Save == true
    open("Data/$(L)_$N/elapsed_time.txt", "w") do file; write(file, "$t\n"); end
end


println(t)

# a,b,c,d,e,f,g,h,i,j,k = MHvideo(rand(Float32, 40, 40)*2*pi32, 40, 2000, .2, Nbin, 200, pi32, PBC) # zeros(L[l],L[l])
# println("    E = ", round(a;digits=4), " ± ", round(b;digits=4), "    Mag = ", round(c;digits=4), " ± ", round(d;digits=4), "    Suscept: ", round(e;digits=5), " ± ", round(f;digits=6), "    C = ", round(g;digits=7), " ± ", round(h;digits=7), "    accept ", round(k;digits=3))
# COLOR MAP !!!


#=
200000   48T   8,10,12,14,16,20    58min
200000   48T   8,12,20             29min   -> 18'

200_000   8T   8                    17 (skip 10), 21 (skip 1)
1_000_000 8T   8                    1'33 (skip=10)
=#



m=12
lattice = Initial_lattice(m, pi32)
# Energy(lattice, m, PBC, Δz)/ m^2
# heatmap(matrixcolor(lattice, m), aspect_ratio = 1, size = (400,400), colormap = :coolwarm, legend = false, framestyle=:box, title = "T")
E = MH(lattice, m, N, .2f0, start, pi32, L, -.01, .5, false)

a,b,c,d#=,e,f,g,h,i,j=#,k = MHvideo(lattice, m, 2000, .1f0, 20, 1600, pi32, PBC, -.9, .5) # zeros(L[l],L[l])
println(k)

mean(cos.(lattice[2,:,:]).^2)






# Pouvoir tune Delta theta indépendamment.
# For low temperature from High temperature to low
# very the limit (Tbkt for xy, Tc for ising (with susceptibility))
# cluster update (wolf algo)

