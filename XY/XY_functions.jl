function Energy(Lattice, PBC)    # Calculate the Energy of a configuration
    L=Int(sqrt(length(Lattice)))
    E=0
    n = PBC ? L : L-1
    for i=1:n
        for j=1:n
            E -= cos(Lattice[i,j]-Lattice[i,mod(j,L)+1]) + cos(Lattice[i,j]-Lattice[mod(i,L)+1,j])
        end
    end
    return E
end

function SingleFlip(Lattice, i, j, T, σ, twopi, PBC, L)      # Propose a new configuration and either accept or reject it
    newspin = rand(TaskLocalRNG(), Normal(Lattice[i,j], σ))
    Δ = 0f0
    if PBC==false
        if i!=1
            Δ += cos(newspin-Lattice[i-1,j]) - cos(Lattice[i,j]-Lattice[i-1,j])
        end
        if i != L
            Δ += cos(newspin-Lattice[i+1,j]) - cos(Lattice[i,j]-Lattice[i+1,j])
        end
        if j != 1
            Δ += cos(newspin-Lattice[i,j-1]) - cos(Lattice[i,j]-Lattice[i,j-1])
        end
        if j != L
            Δ += cos(newspin-Lattice[i,j+1]) - cos(Lattice[i,j]-Lattice[i,j+1])
        end
    else
        Δ += cos(newspin-Lattice[mod(i-2,L)+1,j])-cos(Lattice[i,j]-Lattice[mod(i-2,L)+1,j]) + cos(newspin-Lattice[mod(i,L)+1,j])-cos(Lattice[i,j]-Lattice[mod(i,L)+1,j]) + cos(newspin-Lattice[i,mod(j-2,L)+1])-cos(Lattice[i,j]-Lattice[i,mod(j-2,L)+1]) + cos(newspin-Lattice[i,mod(j,L)+1])-cos(Lattice[i,j]-Lattice[i,mod(j,L)+1])
    end
    if exp(Δ/T) > rand(TaskLocalRNG())
        Lattice[i,j] = mod(newspin, twopi)
        acceptance=true
        E=-Δ
    else
        acceptance=false
        E=0
    end
    return Lattice[i,j], acceptance, E
end

# function MultipleFlips(Lattice,T,σ,twopi, PBC)      # Propose a new configuration and either accept or reject it
#     L=Int(sqrt(length(Lattice)))
#     rotation = rand(TaskLocalRNG(), Normal(0, 2*σ))
#     n = round(Int, rand(TaskLocalRNG(), Float32)*20+15)
#     # Δ = (cos(newspin-Lattice[i,j+1]) + cos(newspin-Lattice[i,j-1]) + cos(newspin-Lattice[i+1,j]) + cos(newspin-Lattice[i-1,j]))  -  (cos(Lattice[i,j]-Lattice[i,j+1]) + cos(Lattice[i,j]-Lattice[i,j-1]) + cos(Lattice[i,j]-Lattice[i+1,j]) + cos(Lattice[i,j]-Lattice[i-1,j]))
#     Δ = 0f0
#     if PBC==false
        
#         E = Energy(Lattice, PBC); Δ=0
#     else
#         Δ += cos(newspin-Lattice[mod(i-2,L)+1,j])-cos(Lattice[i,j]-Lattice[mod(i-2,L)+1,j]) + cos(newspin-Lattice[mod(i,L)+1,j])-cos(Lattice[i,j]-Lattice[mod(i,L)+1,j]) + cos(newspin-Lattice[i,mod(j-2,L)+1])-cos(Lattice[i,j]-Lattice[i,mod(j-2,L)+1]) + cos(newspin-Lattice[i,mod(j,L)+1])-cos(Lattice[i,j]-Lattice[i,mod(j,L)+1])
#     end
#     if exp(Δ/T) > rand(TaskLocalRNG(), Float32)
#         Lattice[i,j] = mod(newspin, twopi)
#         acceptance=true
#         E=-Δ
#     else
#         acceptance=false
#         E=0
#     end
#     return Lattice, acceptance, E
# end

function MHvideo(Lattice, N, T, Nbin, start, twopi, PBC)     # Sampler
    L=Int(sqrt(length(Lattice)))
    σ=.2f0+T^1.5f0+T^.3f0/5
    acceptance::Int64=0
    Energies= zeros(Float32, N-start+2)
    Energies[1] = Energy(Lattice, PBC)
    magx = zeros(Float32, N-start+1)
    Mag  = zeros(Float32, N-start+1)
    Energies[1] = Energy(Lattice, PBC)
    corr = zeros(Float32, length(RowCorr(Lattice, PBC, L)))
    Lattices=[]
    vor = 0f0
    
    for i=1:(start-1)                   # Calculate the energy at each flip
        for j=1:L
              for k=1:L
                Lattice[j,k], a, ΔE=SingleFlip(Lattice,j,k,T,σ,twopi, PBC, L)
                acceptance+=a
                Energies[1]+=ΔE
              end
        end
        if mod(i,Int(N/400))==0
            push!(Lattices, copy(Lattice))
        end
    end
    for i=start:N                       # 2nd for loop to optimize the code, to not push useless values in vectors
        E=Energies[i-start+1]
        for j=1:L
            for k=1:L
                Lattice[j,k], a, ΔE=SingleFlip(Lattice,j,k,T,σ,twopi, PBC, L)
                acceptance+=a
                E+=ΔE
            end
        end
        Energies[i-start+2] = E
        magx[i-start+1] = mean([cos(Lattice[j]) for j=1:L^2])
        Mag[i-start+1] = abs(mean([ℯ^(im*Lattice[j]) for j=1:L^2]))
        vor += Vortex(Lattice,twopi, PBC)
        corr += RowCorr(Lattice, PBC, L)
        if mod(i,Int(N/400))==0
            push!(Lattices, copy(Lattice))
        end
    end
    popfirst!(Energies)
    Energies /= !PBC*L*(L-1) + PBC*L^2
    Magbin = Binor(Mag, Nbin)
    # χxxmean   = mean(magx.^2)-mean(magx)^2/T*L^2
    # χxxbin    = Binor2nd(magx, Nbin, T, L)
    χmean     = (mean(Mag.^2)-mean(Mag)^2)/T*L^2
    # χbin    = Binor2nd(Mag, Nbin, T, L)
    Cmean     = (mean(Energies.^2)-mean(Energies)^2)/T^2*L^2
    # Cbin      = Binor2nd(Energies, Nbin, T, L)/T
    EnergyBins=Binor(Energies, Nbin)

    anim = @animate for j in eachindex(Lattices)
        heatmap(matrixcolor(Lattices[j]), aspect_ratio = 1, size = (400,400), colormap = :coolwarm, legend = false, framestyle=:box, title = "T=$T   $(Int(j*N/400))")    
    end
    gif(anim, "C:/Users/yanis/Desktop/Cours/These/Spin_$N.gif", fps=20)
    
    return mean(Energies), std(EnergyBins), mean(Mag), std(Magbin), χmean, Errorpropagation(Magbin,std(Magbin))/T*L^2 #=std(χbin)=#, Cmean, Errorpropagation(EnergyBins,std(EnergyBins))/T^2*L^2#=std(Cbin)=#, vor/(N-start+1), corr/(N-start+1), acceptance/(L^2*N)
end

function MH(Lattice, N, T, start, twopi, Save, AllLattices, InitalT, Skip=10)     # Sampler
    L=Int(sqrt(length(Lattice)))
    σ=.2f0+T^1.5f0+T^.3f0/5
    Nmeasurement = Int((N-start))
    acceptance::Int64=0
    Energies= zeros(Float32, Nmeasurement+1)
    Energies[1] = Energy(Lattice, PBC)
    Mag  = zeros(Float32, Nmeasurement)
    vor  = 0f0
    corr = zeros(length(RowCorr(Lattice, PBC, L)))
    
    for i=1:(start-1)                   # Calculate the energy at each flip
        for j=1:L
            for k=1:L
                Lattice[j,k], a, ΔE = SingleFlip(Lattice,j,k,T,σ,twopi, PBC, L)
                acceptance+=a
                Energies[1]+=ΔE
            end
        end
    end
    for i=1:Nmeasurement                       # 2nd for loop to optimize the code, to not push useless values in vectors
        E=Energies[i]
        for j=1:L
            for k=1:L
                Lattice[j,k], a, ΔE = SingleFlip(Lattice,j,k,T,σ,twopi, PBC, L)
                acceptance+=a
                E+=ΔE
            end
        end
        if mod(i, Skip) ==0                            # to not measure every lattice sweeps
            Energies[i+1] = E
            Mag[i] = abs(mean(ℯ^(im*x) for x in Lattice))
            vor += Vortex(Lattice, twopi, PBC)
            corr += RowCorr(Lattice, PBC, L)
        end
    end
    popfirst!(Energies)     # as there still was the energy of the first configuration
    Energies /= !PBC*L*(L-1) + PBC*L^2
    accept = acceptance/(L^2*N)
    if Save==true; @save "Data/$(AllLattices)_$(InitalT)_$N/$(L)_$T.jld2" Energies Mag vor corr accept;end
    return mean(Energies)#mean(Energies), std(EnergyBins), mean(Mag), std(Magbin), χmean#= mean(Susc)=#, #=Errorpropagation(Magbin,std(Magbin))/T*L^2 std(Binor(Susc, Nbin))=# std(χbin), Cmean, Errorpropagation(EnergyBins,std(EnergyBins))/T^2*L^2#=std(Cbin)=#, vor/(N-start+1), corr/(N-start+1), acceptance/L^2/N
end

# function MHcluster(Lattice,N,T,Nbin,start,twopi)     # Sampler
#     L=Int(sqrt(length(Lattice)))
#     σ=.2f0+T^1.5f0+T^.3f0/5
#     acceptance::Int64=0
#     Energies= zeros(Float32, N-start+2)
#     Energies[1] = Energy(Lattice, PBC)
#     Mag  = zeros(Float32, N-start+1)
#     vor  = 0f0
#     corr = zeros(length(RowCorr(Lattice, PBC, L)))
#     M = Int(start/500)
    
#     for m=1:M                   # Calculate the energy at each flip
#         for i=1:499             # as the 500th one will be calculate after the this loop
#             for j=1:L
#                 for k=1:L
#                     Lattice[j,k], a, ΔE=SingleFlip(Lattice,j,k,T,σ,twopi, PBC)
#                     acceptance+=a
#                     Energies[1]+=ΔE
#                 end
#             end
#         end

#         # macroscopic change

#     end

#     M = Int((N-start)/500)
#     for m=1:M                       # 2nd for loop to optimize the code, to not push useless values in vectors
#         for i=1:499
#             E=Energies[i+500*(M-1)]
#             for j=1:L
#                 for k=1:L
#                     Lattice[j,k], a, ΔE=SingleFlip(Lattice,j,k,T,σ,twopi, PBC)
#                     acceptance+=a
#                     E+=ΔE
#                 end
#             end
#             Energies[i+500*(M-1)+1] = E
#             # magx[i+500*(M-1)] = mean([cos(Lattice[j]) for j=1:L^2])
#             Mag[i+500*(M-1)] = abs(mean(ℯ^(im*x) for x in Lattice))
#             vor += Vortex(Lattice,twopi,PBC)
#             corr += RowCorr(Lattice, PBC, L)    
#         end
#         # macroscopic change
        
#         Lattice, a, ΔE=MultipleFlips(Lattice,T,σ,twopi, PBC)
#         acceptance2+=a
#         Energies[500*M+1] = E+ΔE
#         # magx[500*M] = mean([cos(Lattice[j]) for j=1:L^2])
#         Mag[500*M] = abs(mean(ℯ^(im*x) for x in Lattice))
#         vor += Vortex(Lattice,twopi,PBC)
#         corr += RowCorr(Lattice, PBC, L) 
#     end
#     popfirst!(Energies)
#     Energies /= !PBC*L*(L-1) + PBC*L^2
#     Magbin    = Binor(Mag, Nbin)
#     # χxxmean   = mean(magx.^2)-mean(magx)^2/T*L^2
#     # χxxbin    = Binor2nd(magx, Nbin, T, L)
#     χmean     = (mean(Mag.^2)-mean(Mag)^2)/T*L^2
#     # χbin    = Binor2nd(Mag, Nbin, T, L)
#     Cmean     = (mean(Energies.^2)-mean(Energies)^2)/T^2*L^2
#     # Cbin      = Binor2nd(Energies, Nbin, T, L)/T
#     EnergyBins= Binor(Energies, Nbin)

#     return mean(Energies), std(EnergyBins), mean(Mag), std(Magbin), χmean, Errorpropagation(Magbin,std(Magbin))/T*L^2 #=std(χbin)=#, Cmean, Errorpropagation(EnergyBins,std(EnergyBins))/T^2*L^2#=std(Cbin)=#, vor/(N-start+1), corr/(N-start+1), acceptance/L^2/(N-M)
# end

function matrixcolor(A)     # change the phase by a color
    N=Int(sqrt(length(A)))
    B=zeros(RGB{Float32},N,N)
    for i=1:N^2
        B[i]= RGB{Float32}(.5+cos(A[i])/2, .0, .5+sin(A[i])/2,)
    end
    return B
end

function RowCorr(Lattice, PBC, L)   # Return the average row Correlation of the matrix (in a Correlation vector : first neigbour, second neigbour...)
    corr=[]
    if PBC==false
        for i=1:(L-1)
            a=0
            for j=1:L
                for k=1:(L-i)
                    a+=cos(Lattice[j,k]-Lattice[j,k+i])
                end
            end
            push!(corr, a/(L-i)/L)
        end
    else
        for i=1:round(Int, L/2-1)
            a=0
            for j=1:L
                for k=1:L
                    a+=cos(Lattice[j,k]-Lattice[j,mod(k+i-1,L)+1])
                end
            end
            push!(corr, a)
        end
    end
    return corr/L^2
end

function Vortex(A, twopi, PBC) # A is the lattice, Vortex return a matrix with -1,0 or 1 depending on if there is an anti-vortex, nothing or a vertex   
    N=Int(sqrt(length(A)))
    vortex=0
    n = PBC ? N : N-1
    for i=1:n
        for j=1:n
            v = mod(A[i,j]-A[i,mod(j,N)+1]+twopi/2,twopi) + mod(A[i,mod(j,N)+1]-A[mod(i,N)+1,mod(j,N)+1]+twopi/2,twopi) + mod(A[mod(i,N)+1,mod(j,N)+1]-A[mod(i,N)+1,j]+twopi/2,twopi) + mod(A[mod(i,N)+1,j]-A[i,j]+twopi/2,twopi)-2*twopi
            if round(Int,v) == 3
                vortex += 1
            end
        end
    end
    return vortex
end

function Binor(arr, Nbin) # take them Nbin by Nbin and average each the bins.
    Bins=[]
    Nperbin = round(Int, length(arr)/Nbin)
    for i=1:Nbin
        push!(Bins, mean(arr[((i-1)*Nperbin)+1:i*Nperbin]))
    end
    return Bins
end

function Binor2nd(arr, Nbin, T, L) # USELESS FOR MH (as I use error propag)
    Bins=[]
    Nperbin = round(Int, length(arr)/Nbin)
    for i=1:Nbin
        push!(Bins, (mean(arr[(i-1)*Nperbin+1:i*Nperbin].^2) - mean(arr[(i-1)*Nperbin+1:i*Nperbin])^2)/T*L^2)
    end
    return Bins
end

function SaveAll(t0)  # save data in different files in the same folder
    ΔT = round(T[2]-T[1], digits=4)
    CSV.write("Data/$(L)_$(T[1])-$(ΔT)_$N/mag.csv", DataFrame(Mmean, :auto))
    CSV.write("Data/$(L)_$(T[1])-$(ΔT)_$N/magerror.csv", DataFrame(MΔ, :auto))
    CSV.write("Data/$(L)_$(T[1])-$(ΔT)_$N/e.csv", DataFrame(Emean, :auto))
    CSV.write("Data/$(L)_$(T[1])-$(ΔT)_$N/eerror.csv", DataFrame(EΔ, :auto))
    CSV.write("Data/$(L)_$(T[1])-$(ΔT)_$N/susc.csv", DataFrame(χmean, :auto))
    CSV.write("Data/$(L)_$(T[1])-$(ΔT)_$N/suscerror.csv", DataFrame(χΔ, :auto))
    CSV.write("Data/$(L)_$(T[1])-$(ΔT)_$N/c.csv", DataFrame(Cmean, :auto))
    CSV.write("Data/$(L)_$(T[1])-$(ΔT)_$N/cerror.csv", DataFrame(CΔ, :auto))
    CSV.write("Data/$(L)_$(T[1])-$(ΔT)_$N/vor.csv", DataFrame(vor, :auto))
    CSV.write("Data/$(L)_$(T[1])-$(ΔT)_$N/corr.csv", DataFrame(Corrmean, :auto))
    CSV.write("Data/$(L)_$(T[1])-$(ΔT)_$N/accept.csv", DataFrame(Acceptance, :auto))
    t = round(Int, time()-t0)
    open("Data/$(L)_$(T[1])-$(ΔT)_$N/elapsed_time.txt", "w") do file
        write(file, "$t\n")
    end
    println("  --  $t  --  Saved  --  ")
end

function basicplot(L, T, vect, ytitle="", error=[],  save=false)   # read data from a file and plot it
    p=plot()
    if typeof(vect) == Matrix{Any}
        for t=1:Int(length(T)/8):length(T)
            plot!(collect(1:Int(L[end]/2-1)), vect[end,:][t], label="$(T[t])",legend=:topright)
        end
        title!(p, "Correlation (L=$(L[end])) as a function of distance")
        xlabel!(p, "Distance (site)")
    else
        if error != [] 
            for l in eachindex(L)
                plot!(T, vect[l, :], yerr = Vector(error[l,:]), markerstrokecolor=:auto, label="$(L[l])")#, ylims=(0,2))
            end
        else
            for l in eachindex(L)
                plot!(T, vect[l, :], label="$(L[l])")
            end
        end
        title!(p, ytitle*" as a function of T")
        xlabel!(p, "Temperature")
    end
    ylabel!(p, ytitle)
    if save ==true
        savefig("Plot/"*ytitle*".pdf")
    end
    display(p)
end

# function getdata(L, T, N) # read and save the data in a dictionnary
#     ΔT = round(T[2]-T[1];digits=5)
#     keys = [:"e", :"eerror", :"mag", :"magerror", :"c", :"cerror", :"susc", :"suscerror", :"vor", :"corr"]
#     dict= Dict(k => Tables.matrix(CSV.read("Data/$(L)_$(T[1])-$(ΔT)_$N/"*string(k)*".csv", DataFrame)) for k in keys)
#     dict[:"corr"] = map(s -> parse.(Float64, split(strip(s, ['[', ']']), ",")), dict[:"corr"])
#     return dict
# end

function Errorpropagation(v, Δ) # propagation of error for c and χ from the error on E and M
    return sqrt(sum((v .- mean(v)).^2)/length(v))*2*Δ
end

function interp(L, T, thing) # interpolate and find the max
    Max = []
    xinterp = collect(.1:.001:2)
    p=plot()
    for l in eachindex(L)
        yinterp=Spline1D(T, thing[l,:], k=3)(xinterp)
        push!(Max, [findmax(yinterp[301:end])[1], xinterp[findmax(yinterp[301:end])[2]]])
    end
    return Max
end

function linear(x, p) # linear function
    return p[1]*x .+ p[2]   
end

function crit(L, T, y, graph="") # calculate α & γ
    ymax = interp(L, T, y)
    fit =  curve_fit(linear, log.(L), log.([v[1] for v in ymax]), if (y=="c"); [.1,-.7]; elseif (y=="susc"); [2.,-4.]; else [.9,-.9] end)
    m,p = coef(fit); σ = stderror(fit)
    if graph!=""
        a=plot()
        plot!(a, log.(L), log.([v[1] for v in ymax]), seriestype=:scatter)
        plot!(a, log.(L), m*log.(L) .+p)
        xlabel!(a, "Ln(Lattice length)")
        ylabel!(a, "Ln(max("*graph*"))")
        display(a)
    end
    return m, σ[1]
end

function critlength(T, data, t, ln=false) # Calculate critical exp or correlation length
    closestT = argmin(abs.(T .- t)); println(T[closestT])
    data=data[Int(length(Corrmean)/length(T))*closestT]
    negativ = findfirst(x -> x < 1e-5, data)
    if typeof(negativ) == Nothing; l = length(data); else; l=negativ-1; end
    println(l)
    x=collect(1:l)
    fit =  curve_fit(linear, log.(x)*ln + !ln*x, log.(data[1:l]), [-.09,9.])
    m,p = coef(fit); σ = stderror(fit)

    a=plot()
    plot!(a, log.(x)*ln + !ln*x, log.(data[1:l]), seriestype=:scatter)
    plot!(a, log.(x)*ln + !ln*x, m*(log.(x)*ln + !ln*x) .+p)
    ylabel!(a, "Ln(Correlation)")
    if ln == true; xlabel!(a, "Ln(Distance)"); ξ=σξ=0
    else; xlabel!(a, "Distance"); ξ =-1/m; σξ = σ[1]/m^2;end
    title!("Correlation vs distance at T = $(T[closestT])")
    display(a)
    return m*ln + !ln*ξ, σ[1]*ln+!ln*σξ
end

# α : specific heat
# β : zero field mag
# γ : zero field isothermal suscxx
# δ : Critical isothermal
# ν : corr length

# heatmap(matrixcolor(Lattices[end]), aspect_ratio = 1, size = (375,400), colormap = :coolwarm,legend = false, framestyle=:box)
# using DelimitedFiles; writedlm(stdout, vort)  # print matrix in a nice form


#=
function Gauss(x,p)
    return p[3]*ℯ.^(-((x.-p[1])/p[2]).^2)
end

function fitmax(x, y, p=[1.1, .3, .9, .5])
    return curve_fit(Gauss, x, y, p).param[1]
end
=#


#=
acceptsatu=[[36000000, 36300000,38400000, 39200000, 40000000, 40500000, 43200000, 57600000, 67600000, 90000000, 102400000], [.466, .462, .437, .428, .419, .414, .388, .291, .246, .186, .164]]
plot(acceptsatu[1], acceptsatu[2], label="Data")
plot!(acceptsatu[1],16720000 ./acceptsatu[1], label="Fit: cte * 1/x")
xlabel!("Number of site multiplied by number of lattice sweep")
ylabel!("Saturation of the computed accpetance rate")
=#