using Plots, Distributions, Colors, Base.Threads, Dierckx, LsqFit, Random, JLD2   # To plot, To have random distrib of spin, colors for lattices, parallelize, interpolate, Random, save data in compact file
include("./Heisenberg_functions.jl")

N = 50_000_000 # Lattice flip
# T = unique(sort(Float32.(vcat(collect(0.955:.015:1.21),collect(.1:.1:.5),collect(.56:.04:.92),collect(1.24:.04:1.64),collect(1.7:.1:2))))) # XY
# T = unique(sort(Float32.(vcat(collect(2:.1:2.7), collect(2.1:.02:2.44))))) # Ising
# T = Float32.(collect(.1:.25:2))            # Temperature
T = unique(sort(Float32.(vcat(collect(.4:.1:1.7),collect(.6:.05:1.4), collect(.65:.025:.9))))) # XYZ
L = [8, 12, 20]           # Lattice size
# L = [8,12,20,32,50,70,100]
# d    = Float32.([-1000])
d = Float32.([-3,-1,-.3,-.1,0,.1,.5,4,100])
p = .6f0
burn = Int(min(N/4,100000))    # Burning period
PBC   = true                    # Periodic Boundary Conditions
pi32  = Float32(π)
Save  = true

# ---------------------------   Old one   ------------------------------------ #

println("  --  $N sweeps  --  ")
starttime = time()
for z in eachindex(d)
    for l in eachindex(L)
        Threads.@threads for t in eachindex(T)
            Energy, Mag = MH(Initial_lattice(L[l], pi32), L[l], N, T[t], burn, pi32, L, d[z], .5, Save) # zeros(Float32, L[l],L[l])
            println("T = $(T[t]) \tL = $(L[l]) \td = $(d[z]) \tE = ", Energy, "\t Mag : ", Mag)
        end
    end
end
t = round(Int, time()-starttime)
if Save == true
    open("Data/$(L)_$N/elapsed_time.txt", "a") do file; write(file, "$t\n$d\n$T\n\n"); end
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



m=20
ZM = zeros(20)
for i=1:20
    lattice = Initial_lattice(m, pi32)
    E, magz, a, σ,aa, ab = MHXY(lattice, m, 4000, .2f0, 3000, pi32, L, 0.2, false)
    ZM[i] = mean(magz)
    println(ZM[i])
end
count(abs.(ZM) .> .9)

lattice = Initial_lattice(m, pi32)
E, magz, a, σ,sigma1, sigma2, absmagz = MHXY(lattice, m, 20000, .2f0, 0, pi32, L, .2, false)
Energy(lattice, m, PBC, 1.32)/m^2
E
plot(sigma1)
plot(sigma2)
# plot(magz)
plot(absmagz)
println("\n",a)

lattice = Initial_lattice(m, pi32)
E = MH(lattice, m, 4000, .2f0, 200, pi32, L, 0.3, p, false)
heatmap(matrixcolor(lattice, m), aspect_ratio = 1, size = (400,400), colormap = :coolwarm, legend = false, framestyle=:box, title = "T = 0.2")
Energy(lattice, m, PBC, 0.3)/m^2


a,b,c,d#=,e,f,g,h,i,j=#,k = MHvideo(lattice, m, 200, .2f0, 20, 60, pi32, PBC, -.5, p)
println("E = ", a, " ± ", b, "\t M = ", c, " ± ", d, "\t acceptance: ", k)


# The energy is not very accurate, do we need more accuracy ?

# Be able to tune Delta and Theta independently : tunable coef<1 ? theta_new = (theta_new - theta_old)*coef + theta_old       Tune independently theta and phi => no more isotropic, even the isingflip brakes isotropy ?
# ?? Use the same z repartition of the "function Initial_lattice" and no more spherical -> spherical ?          Gaussian multiplied by sqrt(1-z2)


x = collect(-1:.01:1)
y = sqrt.(1 .-x.^2)
z = ℯ.^(-10*(x .- 1/sqrt(2)).^2)
yz = sqrt.(1 .-x.^2) .* ℯ.^(-10*(x .- 1/sqrt(2)).^2)

plot(x,y)
plot!(x,z)
plot!(x,yz)


# Parallel tempering

# Remove Energy2 from MH

x=collect(.001:.001:1)
y=sqrt.(abs.((asin.(2x.-1) .+1 -2*x)/(pi32/2-1))).*(sign.(x .- .5).+1)/2
plot(x, y)


# --- New one with parallel tempering --- #


N=20000
L=[8,12]
T=Float32.([.2,.7,1.1,1.6])
d=Float32.([0,2])

println("  --  $N sweeps  --  ")
starttime = time()
for z in d
    for l in L
        E1, M1 = MH_parallel_tempering(l, N, T, burn, pi32, L, z, p, Save)
    end
end
t = round(Int, time()-starttime)
if Save == true
    open("Data/$(L)_$N/elapsed_time.txt", "a") do file; write(file, "$t\n$d\n$T\n\n"); end
end

println(t)