function Energy(Lattice, L, PBC, Δz)    # Calculate the Energy of a configuration
    E=0
    n = PBC ? L : L-1
    for i=1:n
        for j=1:n
            E += -Correlation(Lattice[:,i,j],Lattice[:,i,mod(j,L)+1]) - Correlation(Lattice[:,i,j],Lattice[:,mod(i,L)+1,j]) + Δz*cos(Lattice[2,i,j])^2
        end
    end
    return E
end

function Correlation(spin1, spin2)
    return sin(spin1[2])*sin(spin2[2])*cos(spin1[1]-spin2[1]) + cos(spin1[2])*cos(spin2[2])  
end

function Initial_lattice(L, pi32)
    lattice = Array{Float32}(undef, 2, L, L)
    lattice[1, :, :] .= 2*pi32*rand(TaskLocalRNG(), Float32, L, L)
    lattice[2, :, :] .= acos.(1 .- 2*rand(TaskLocalRNG(), Float32, L, L))
    return lattice
end

function SingleFlip(Lattice, i, j, T, L, σ, pi32, PBC, Δz, p)      # Propose a new configuration and either accept or reject it
    newspin, Which = NewSpin(Lattice[:,i,j], σ, pi32, p)
    Δ = 0f0
    if PBC==false
        if i!=1
            Δ += Correlation(newspin,Lattice[:,i-1,j]) - Correlation(Lattice[:,i,j],Lattice[:,i-1,j]) # This is the right sign as the bounding energy is negative (and is calculate positively)
        end
        if i != L
            Δ += Correlation(newspin,Lattice[:,i+1,j]) - Correlation(Lattice[:,i,j],Lattice[:,i+1,j])
        end
        if j != 1
            Δ += Correlation(newspin,Lattice[:,i,j-1]) - Correlation(Lattice[:,i,j],Lattice[:,i,j-1])
        end
        if j != L
            Δ += Correlation(newspin,Lattice[:,i,j+1]) - Correlation(Lattice[:,i,j],Lattice[:,i,j+1])
        end
    else
        Δ += Correlation(newspin,Lattice[:,mod(i-2,L)+1,j]) - Correlation(Lattice[:,i,j],Lattice[:,mod(i-2,L)+1,j])   +   Correlation(newspin,Lattice[:,mod(i,L)+1,j]) - Correlation(Lattice[:,i,j],Lattice[:,mod(i,L)+1,j])   +   Correlation(newspin,Lattice[:,i,mod(j-2,L)+1]) - Correlation(Lattice[:,i,j],Lattice[:,i,mod(j-2,L)+1])   +   Correlation(newspin,Lattice[:,i,mod(j,L)+1]) - Correlation(Lattice[:,i,j],Lattice[:,i,mod(j,L)+1])
    end
    Δ += Δz * (cos(Lattice[2,i,j])^2 - cos(newspin[2])^2)
    if exp(Δ/T) > rand(TaskLocalRNG())
        Lattice[:,i,j] .= [mod(newspin[1], 2*pi32), newspin[2]]
        acceptance=true
        E=-Δ
    else
        acceptance=false
        E=0
    end
    return Lattice[:,i,j], E, acceptance, Which
end

function NewSpin(spin, σ, pi32, p)
    Which = rand(TaskLocalRNG()) < p
    if Which == true
        spin[2] = pi32 - spin[2]
    else
        alpha   = rand(TaskLocalRNG(), Float32)*2*pi32
        beta    = min(abs(rand(TaskLocalRNG(), Normal(0f0, σ))), pi32)
        spin[1] = atan(sin(spin[1])*(cos(spin[2])*sin(beta)*cos(alpha) + sin(spin[2])*cos(beta)) + cos(spin[1])*sin(beta)*sin(alpha),   cos(spin[1])*(cos(spin[2])*sin(beta)*cos(alpha) + sin(spin[2])*cos(beta)) - sin(spin[1])*sin(beta)*sin(alpha)) # alpha c'est le phi'
        spin[2] = acos(max(-1,min(1,cos(spin[2])*cos(beta) - sin(spin[2])*sin(beta)*cos(alpha))))
    end
    return spin, Which+1
end

function MultipleIsingFlips(Lattice, L, T, E0, pi32, Δz, PBC)      # Propose a new configuration and either accept or reject it
    R = min(rand(TaskLocalRNG(), 1:6), round(Int,L/4-1))
    k, l = rand(TaskLocalRNG(), 1:L, 2)
    newLattice = Lattice
    for i=k-R:k+R
        for j=l-R:l+R
            newLattice[2, mod(i-1,L)+1, mod(j-1,L)+1] = pi32 - Lattice[2, mod(i-1,L)+1, mod(j-1,L)+1]
        end
    end
    E1 = Energy(newLattice, L, PBC, Δz)
    Δ = (E1 - E0)/(4*R^2) # check the difference energy per site
    # if exp(-Δ/T) > rand(TaskLocalRNG())
    acceptance=true
    E = E1
    # else
    #     acceptance=false
    #     E = E0
    #     newLattice = Lattice
    # end
    return newLattice, acceptance, E
end

function MHvideo(Lattice, L, N, T, Nbin, start, pi32, PBC, Δz, p)     # Sampler    
    σ=.2f0+T^1.5f0+T^.3f0/5
    acceptance=[0,0,0]                  # around, Ising, multiflips
    Try = [0,0,0]
    Energies= zeros(Float32, N-start+2)
    Energies[1] = Energy(Lattice, L, PBC, Δz)
    Mag  = zeros(Float32, N-start+1)
    Energies[1] = Energy(Lattice, L, PBC, Δz)
    corr = zeros(Float32, length(RowCorr(Lattice, L, PBC)))
    Lattices=[]
    
    for i=1:(start-1)                   # Calculate the energy at each flip
        IsingFlip = mean(cos.(Lattice[2,:,:])^2) > .75       # require that in average [theta < pi/6 or theta > 5 pi/6] to have the possibility (with proba p) to Ising-flip
        if mod(i,400)==0# && IsingFlip==true
            Lattice, b, Energies[1] = MultipleIsingFlips(Lattice, L, T, Energies[1], pi32, Δz, PBC)
            acceptance[3] += b
            Try[3] += 1
        else
            for j=1:L
                for k=1:L
                    Lattice[:,j,k], ΔE, a, Which = SingleFlip(Lattice,j,k,T,L,σ,pi32, PBC, Δz, p*IsingFlip)
                    acceptance[Which]+=a
                    Try[Which]+=1
                    Energies[1]+=ΔE
                end
            end
        end
        if mod(i,Int(N/200))==0
            push!(Lattices, copy(Lattice))
        end
    end
    for i=start:N                       # 2nd for loop to optimize the code, to not push useless values in vectors
        E=Energies[i-start+1]
        IsingFlip = mean(cos.(Lattice[2,:,:]).^2) > .75       # require that in average [theta < pi/6 or theta > 5 pi/6] to have the possibility (with proba p) to Ising-flip
        for j=1:L
            for k=1:L
                Lattice[:,j,k], ΔE, a, Which = SingleFlip(Lattice,j,k,T,L,σ,pi32, PBC, Δz, p*IsingFlip)
                acceptance[Which]+=a
                Try[Which]+=1
                E+=ΔE
            end
        end
        Energies[i-start+2] = E
        Mag[i-start+1] = sqrt(mean(cos(xy)*sin(z) for xy in Lattice[1,:,:], z in Lattice[2,:,:])^2 + mean(sin(xy)*sin(z) for xy in Lattice[1,:,:], z in Lattice[2,:,:])^2 + mean(cos(z) for z in Lattice[2,:,:])^2)
        # vor += Vortex(Lattice,pi32, PBC)
        corr += RowCorr(Lattice, L, PBC)
        if mod(i,Int(N/200))==0
            push!(Lattices, copy(Lattice))
        end
    end
    popfirst!(Energies)
    Energies /= !PBC*L*(L-1) + PBC*L^2
    Magbin = Binor(Mag, Nbin)
    # χmean     = (mean(Mag.^2)-mean(Mag)^2)/T*L^2
    # Cmean     = (mean(Energies.^2)-mean(Energies)^2)/T^2*L^2
    EnergyBins=Binor(Energies, Nbin)
    println(typeof(Lattices))

    anim = @animate for j in eachindex(Lattices)
        heatmap(matrixcolor(Lattices[j], L), aspect_ratio = 1, size = (400,400), colormap = :coolwarm, legend = false, framestyle=:box, title = "T=$T  Δ=$Δz  $(Int(j*N/200))")    
    end
    # anim = @animate for j in eachindex(Lattices)
    #     p1 = heatmap(matrixcolor(Lattices[j], L),
    #         aspect_ratio = 1, size = (288,288),
    #         colormap = :coolwarm, legend = false,
    #         framestyle=:box, title = "matrixcolor")
    #     p2 = heatmap(matrixcolorXY(Lattices[j], L),
    #         aspect_ratio = 1, size = (288,288),
    #         colormap = :coolwarm, legend = false,
    #         framestyle=:box, title = "matrixcolorXY")
    #     p3 = heatmap(matrixcolorZ(Lattices[j], L),
    #         aspect_ratio = 1, size = (288,288),
    #         colormap = :coolwarm, legend = false,
    #         framestyle=:box, title = "matrixcolorZ")
    #     plot(p1, p2, p3, layout = (1,3), size = (864,288), title = "T=$T   $j")
    # end

    gif(anim, "C:/Users/yanis/Desktop/Cours/These/Heisenberg/Spin_$N.gif", fps=10)

    return mean(Energies), std(EnergyBins), mean(Mag), std(Magbin)#=, χmean, Errorpropagation(Magbin,std(Magbin))/T*L^2 #=std(χbin)=#, Cmean, Errorpropagation(EnergyBins,std(EnergyBins))/T^2*L^2, #=vor/(N-start+1),=# corr/(N-start+1)=#, acceptance ./ Try
end

function MH(Lattice, L, N, T, start, pi32, AllLattices, Δz, p, Save, Skip=10)     # Sampler
    σ = .2f0 + T^1.5f0 + T^.3f0/5
    Nmeasurement = Int((N-start)/Skip)
    acceptance = [0,0,0]              # around, Ising, multiflips
    Try = [0,0,0]
    Energies = zeros(Float32, Nmeasurement+1)
    Energies[1] = Energy(Lattice, L, PBC, Δz)
    Mag  = zeros(Float32, Nmeasurement)
    # vor  = 0f0
    corr = zeros(length(RowCorr(Lattice, L, PBC)))
    
    for i=1:(start-1)                   # Calculate the energy at each flip
        IsingFlip = mean(cos.(Lattice[2,:,:])^2) > .75       # require that in average [theta < pi/6 or theta > 5 pi/6] to have the possibility (with proba p) to Ising-flip
        if mod(i,2000)==0 && IsingFlip==true
            Lattice, b, Energies[1] = MultipleIsingFlips(Lattice, L, T, Energies[1], pi32, Δz, PBC)
            acceptance[3] += b
            Try[3] += 1
        else
            for j=1:L
                for k=1:L
                    Lattice[:,j,k], ΔE, a, Which = SingleFlip(Lattice,j,k,T,L,σ,pi32, PBC, Δz, p*IsingFlip)
                    acceptance[Which]+=a
                    Try[Which]+=1
                    Energies[1]+=ΔE
                end
            end
        end
    end
    for i=1:Nmeasurement                       # 2nd for loop to optimize the code, to not push useless values in vectors
        E=Energies[i]
        IsingFlip = mean(cos.(Lattice[2,:,:])^2) > .75       # require that in average [theta < pi/6 or theta > 5 pi/6] to have the possibility (with proba p) to Ising-flip
        for j=1:L
            for k=1:L
                Lattice[:,j,k], ΔE, a, Which = SingleFlip(Lattice,j,k,T,L,σ,pi32, PBC, Δz, p*IsingFlip)
                acceptance[Which]+=a
                Try[Which]+=1
                E+=ΔE
            end
        end
        if mod(i, Skip) == 0                            # to not measure every lattice sweeps
            Energies[i+1] = E
            Mag[i] = sqrt(mean(cos(xy)*sin(z) for xy in Lattice[1,:,:], z in Lattice[2,:,:])^2 + mean(sin(xy)*sin(z) for xy in Lattice[1,:,:], z in Lattice[2,:,:])^2 + mean(cos(z) for z in Lattice[2,:,:])^2)
            # vor += Vortex(Lattice, pi32, PBC) #    ?????    A passer en 3D   !!!    ????? ça a du sens ?
            corr += RowCorr(Lattice, L, PBC)
        end
    end
    popfirst!(Energies)     # as there still was the energy of the first configuration
    Energies /= !PBC*L*(L-1) + PBC*L^2
    accept = acceptance ./ Try
    if Save==true; @save "Data/$(AllLattices)_$N/$(L)_$(T)_$(Δz).jld2" Energies Mag #=vor=# corr accept
    end
    return mean(Energies)
end

function matrixcolor(Lattice, L)     # change the phase by a color
    B=zeros(RGB{Float32},L,L)
    for i=1:L
        for j=1:L
            B[i,j]= RGB{Float32}(.5+cos(Lattice[1,i,j])*sin(Lattice[2,i,j])/2, .5+cos(Lattice[2,i,j])/2, .5+sin(Lattice[1,i,j])*sin(Lattice[2,i,j])/2)
        end
    end
    return B
end

function matrixcolorXY(Lattice, L)     # change the phase by a color
    B=zeros(RGB{Float32},L,L)
    for i=1:L
        for j=1:L
            B[i,j]= RGB{Float32}(.5+cos(Lattice[1,i,j])/2, 0, .5+sin(Lattice[1,i,j])/2)
        end
    end
    return B
end

function matrixcolorZ(Lattice, L)     # change the phase by a color
    B=zeros(RGB{Float32},L,L)
    for i=1:L
        for j=1:L
            B[i,j]= RGB{Float32}(0, .5+cos(Lattice[2,i,j])/2, 0)
        end
    end
    return B
end

function RowCorr(Lattice, L, PBC)   # Return the average row Correlation of the matrix (in a Correlation vector : first neigbour, second neigbour...)
    corr=[]
    if PBC==false
        for i=1:(L-1)
            a=0
            for j=1:L
                for k=1:(L-i)
                    a+=Correlation(Lattice[:,j,k]-Lattice[:,j,k+i])
                end
            end
            push!(corr, a/(L-i)/L)
        end
    else
        for i=1:round(Int, L/2-1)
            a=0
            for j=1:L
                for k=1:L
                    a+=cos(Lattice[1,j,k]-Lattice[1,j,mod(k+i-1,L)+1]) * cos(Lattice[2,j,k]-Lattice[2,j,mod(k+i-1,L)+1])
                end
            end
            push!(corr, a)
        end
    end
    return corr/L^2
end

function Vortex(A, pi32, PBC) # A is the lattice, Vortex return a matrix with -1,0 or 1 depending on if there is an anti-vortex, nothing or a vertex   
    N=Int(sqrt(length(A)))
    vortex=zeros(Int, N, N)
    n = PBC ? N : N-1
    for i=1:n
        for j=1:n
            v=mod(A[i,j]-A[i,mod(j,N)+1]+pi32,2*pi32) + mod(A[i,mod(j,N)+1]-A[mod(i,N)+1,mod(j,N)+1]+pi32,2*pi32) + mod(A[mod(i,N)+1,mod(j,N)+1]-A[mod(i,N)+1,j]+pi32,2*pi32) + mod(A[mod(i,N)+1,j]-A[i,j]+pi32,2*pi32)-4*pi32
            vortex[i,j] = round(Int,v/pi32/2)
        end
    end
    return #=vortex, =#Float32(count(x -> x==1, vortex)/N^2)#, count(x -> x==-1, vortex)/N^2
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

function basicplot1(L, T, vect, ytitle="", error=[],  save=false)   # read data from a file and plot it
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

function basicplotΔ(T, Δ, vect, ytitle="", error=[],  save=false)   # read data from a file and plot it
    p=plot()
    if error != [] 
        for z in eachindex(Δ)
            plot!(T, vect[end, :, z], yerr = Vector(error[end, :, z]), markerstrokecolor=:auto, label="$(Δ[z])")#, ylims=(0,2))
        end
    else
        for z in eachindex(Δ)
            plot!(T, vect[end, :, z], label="$(Δ[z])")
        end
    end
    title!(p, ytitle*" as a function of T")
    xlabel!(p, "Temperature")
    ylabel!(p, ytitle)
    if save ==true
        savefig("Plot/"*ytitle*".pdf")
    end
    display(p)
end

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