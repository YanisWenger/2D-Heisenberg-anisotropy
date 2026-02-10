using Plots, Distributions, Colors, CSV, DataFrames, Dierckx, LsqFit, Random, JLD2   # To plot, To have random distrib of spin, colors for lattices, write/read csv (CSV & DataFrames), parallelize, interpolate, save data in compact file
include("./Heisenberg_functions.jl")

N = 100_000 # Lattice flip
T = unique(sort(Float32.(vcat(collect(.1:.1:1.7),collect(.6:.05:1.3)))))
# T = unique(sort(Float32.(vcat(collect(.1:.1:2.4),collect(.28:.04:.96), collect(.62:.02:.78)))))
# T = unique(sort(Float32.(vcat(collect(0.955:.015:1.21),collect(.1:.1:.5),collect(.56:.04:.92),collect(1.24:.04:1.64),collect(1.7:.1:2)))))
# T = Float32.(collect(.7:.1:1.4))            # Temperature
L = [20]           # Lattice size
# L=[8,12,20,32,50]
Nbin = 50                       # Number of bins
PBC=true                        # Periodic Boundary Conditions
Δ = [-1,-.25,0,.25,1]

E    =    zeros(Float32, length(L), length(T), length(Δ))
EΔ   =    zeros(Float32, length(L), length(T), length(Δ))
M    =    zeros(Float32, length(L), length(T), length(Δ))
MΔ   =    zeros(Float32, length(L), length(T), length(Δ))
# vor      =    zeros(Float32, length(L), length(T))
Acceptance =  zeros(3, length(L), length(T), length(Δ))
χ    =    zeros(Float32, length(L), length(T), length(Δ))
χΔ   =    zeros(Float32, length(L), length(T), length(Δ))
C    =    zeros(Float32, length(L), length(T), length(Δ))
CΔ   =    zeros(Float32, length(L), length(T), length(Δ))
Corr = Array{Any}(undef, length(L), length(T), length(Δ))

# ----------------------------------------------- #

for z in eachindex(Δ)
    for l in eachindex(L)
        for t in eachindex(T)
            Data = load("Data/$(L)_$N/$(L[l])_$(T[t])_$(Δ[z]).jld2")
            E[l,t,z], EΔ[l,t,z] = mean(Data[:"Energies"]),                            std(Binor(Data[:"Energies"], Nbin))
            M[l,t,z], MΔ[l,t,z] = mean(Data[:"Mag"]),                                 std(Binor(Data[:"Mag"], Nbin))
            χ[l,t,z], χΔ[l,t,z] = (mean(Data[:"Mag"].^2)-M[l,t,z]^2)/T[t]*L[l]^2,   std(Binor2nd(Data[:"Mag"], Nbin, T[t], L[l]))
            C[l,t,z], CΔ[l,t,z] = (mean(Data[:"Energies"].^2)-E[l,t,z]^2)/T[t]^2*L[l]^2, Errorpropagation(Binor(Data[:"Energies"], Nbin), EΔ[l,t,z])/T[t]^2*L[l]^2
            # vor[l,t]        = Data[:"vor"]
            Corr[l,t,z]       = Data[:"corr"]
            # Acceptance[:,l,t,z] = Data[:"accept"]
            println("N = ", L[l], "\tT = ", T[t], "\tΔ = ", Δ[z], "\t E = ", round(E[l,t,z];digits=3), " ± ", round(EΔ[l,t,z];digits=3), "\tM = ", round(M[l,t,z];digits=3), " ± ", round(MΔ[l,t,z];digits=3), " \t χ = ", round(χ[l,t,z];digits=7), " ± ", round(χΔ[l,t,z];digits=7), "\tC = ", round(C[l,t,z];digits=2), " ± ", round(CΔ[l,t,z];digits=5), "\taccept = ", round.(Acceptance[:,l,t,z];digits=3))
        end
    end
end

basicplot1(L, T, χ, "Title", χΔ)
basicplotΔ(T, Δ, χ, "Title", χΔ)

α, σα = crit(L,T, C, "Capacity")
γ, σγ = crit(L,T, χ, "Susceptibility")

ξ, σξ = critlength(T, Corrmean, 1.0, false)