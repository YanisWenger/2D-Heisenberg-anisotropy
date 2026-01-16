using Plots, Distributions, Colors, CSV, DataFrames, Dierckx, LsqFit, Random, JLD2   # To plot, To have random distrib of spin, colors for lattices, write/read csv (CSV & DataFrames), parallelize, interpolate, save data in compact file
include("./XY_functions.jl")

N = 1000_000 # Lattice flip
# T = Float32.(collect(1.08:.01:1.23))
T = unique(sort(Float32.(vcat(collect(0.955:.015:1.21),collect(.1:.1:.5),collect(.56:.04:.92),collect(1.24:.04:1.64),collect(1.7:.1:2)))))
# T = Float32.(collect(.7:.1:1.4))            # Temperature
# L = [8]           # Lattice size
L=[8,12,20,32,50,70,100]
Nbin = 50                       # Number of bins
PBC=true                        # Periodic Boundary Conditions

E    =    zeros(Float32, length(L), length(T))
EΔ       =    zeros(Float32, length(L), length(T))
M    =    zeros(Float32, length(L), length(T))
MΔ       =    zeros(Float32, length(L), length(T))
vor      =    zeros(Float32, length(L), length(T))
Acceptance =  zeros(Float32, length(L), length(T))
χ    =    zeros(Float32, length(L), length(T))
χΔ       =    zeros(Float32, length(L), length(T))
C    =    zeros(Float32, length(L), length(T))
CΔ       =    zeros(Float32, length(L), length(T))
Corr = Array{Any}(undef, length(L), length(T))

# ----------------------------------------------- #


for l in eachindex(L)
    for t in eachindex(T)
        Data = load("Data/$(L)_$(T[1])_$N/$(L[l])_$(T[t]).jld2")
        E[l,t], EΔ[l,t] = mean(Data[:"Energies"]),                            std(Binor(Data[:"Energies"], Nbin))
        M[l,t], MΔ[l,t] = mean(Data[:"Mag"]),                                 std(Binor(Data[:"Mag"], Nbin))
        χ[l,t], χΔ[l,t] = (mean(Data[:"Mag"].^2)-Mmean[l,t]^2)/T[t]*L[l]^2,   std(Binor2nd(Data[:"Mag"], Nbin, T[t], L[l]))
        C[l,t], CΔ[l,t] = (mean(Data[:"Energies"].^2)-Emean[l,t]^2)/T[t]^2*L[l]^2, Errorpropagation(Binor(Data[:"Energies"], Nbin), EΔ[l,t])/T[t]^2*L[l]^2
        vor[l,t]        = Data[:"vor"]
        Corr[l,t]       = Data[:"corr"]
        Acceptance[l,t]     = Data[:"accept"]
        println("N = ", L[l], "\tT = ", T[t], " \t E = ", round(Emean[l,t];digits=3), " ± ", round(EΔ[l,t];digits=3), "\tM = ", round(Mmean[l,t];digits=3), " ± ", round(MΔ[l,t];digits=3), " \t χ = ", round(χmean[l,t];digits=7), " ± ", round(χΔ[l,t];digits=7), "\tC = ", round(Cmean[l,t];digits=2), " ± ", round(CΔ[l,t];digits=5), "\taccept = ", round(Acceptance[l,t];digits=3))
    end
end

basicplot(L, T, Corr, "Title")#, CΔ)

α, σα = crit(L,T, Cmean, "Capacity")
γ, σγ = crit(L,T, χmean, "Susceptibility")

ξ, σξ = critlength(T, Corrmean, 1.0, false)

