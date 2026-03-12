function Get_T_d(N,L)  # for the first lattice size (supposing they all have the same T and d), gives which T have each d
    files = glob("$(L[1])*.jld2", ("Data/$(L)_$N"))
    for f in files
        name = splitext(last(splitpath(f)))[1]
        parts = split(name, "_")
        if length(parts) != 3
            @warn "Bad filename format" f parts
        end
    end

    parsed = map(files) do f
        name = splitext(last(splitpath(f)))[1]   # robust path handling
        nums = parse.(Float32, split(name, "_"))
        (; n2 = nums[2], n3 = nums[3])
    end

    # all_L = [p.n1 for p in parsed]  # Vector of all L
    all_d = unique(sort([p.n3 for p in parsed]))    # Vector of all d
    all_T = unique(sort([p.n2 for p in parsed]))    # Vector of all T
    T_for_d = Dict{Float32, Vector{Float32}}() # Dictionary mapping each #3 → list of #2
    for p in parsed
        push!(get!(T_for_d, p.n3, Float32[]), p.n2)
    end
    T_for_d = Dict(k => sort(v) for (k,v) in T_for_d)
    return all_T, all_d, T_for_d
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

function basicplot1(L, T_for_d, dvect, d, vect, ytitle="", error=[], save=false)   # read data from a file and plot it
    i = findfirst(==(d), dvect)
    endT = length(vect[1,:,i])
    while endT >= 1 && vect[1,endT,i]==0
        endT -= 1
    end
    p=plot()
    if typeof(vect) == Matrix{Any}
        for t=1:max(round(Int, length(T_for_d[d])/8),1):length(T_for_d[d])
            plot!(collect(1:Int(L[end]/2-1)), vect[end,:][t], label="$(T[t])",legend=:topright) # yet not update to work
        end
        title!(p, "Correlation (L=$(L[end])) as a function of distance")
        xlabel!(p, "Distance (site)")
    else
        if error != []
            if ytitle=="C" || ytitle=="Susc"; lim=(max(0, minimum(vect[:,:,i])),maximum(vect[:,3:end,i])+maximum(error[:,3:end,i]))
            else; lim = :auto
            end
            for l in eachindex(L)
                plot!(T_for_d[d], vect[l, 1:endT,i], yerr = Vector(error[l,1:endT,i]), markerstrokecolor=:auto, label="$(L[l])", ylims = lim)
            end
        else
            for l in eachindex(L)
                plot!(T_for_d[d], vect[l, 1:endT,i], label="$(L[l])", ylims=(max(0, minimum(vect[:,1:endT,i])),maximum(vect[:,1:endT,i])))
            end
        end
        title!(p, ytitle*" as a function of T with d = $d")
        xlabel!(p, "Temperature")
    end
    ylabel!(p, ytitle)
    if save ==true
        savefig("Plot/"*ytitle*".pdf")
    end
    println("size vect: ", size(vect))
    println("size vect[end,:,i]: ", size(vect[end,:,i]))
    println("endT: ", endT)
    display(p)
end

function basicplotd(T_for_d, d, vect, ytitle="", error=[], i=length(vect[:,1,1]), save=false)   # read data from a file and plot it
    p=plot()
    pal = cgrad([RGB(.4,.6,1), RGB(0,0,.5), RGB(0,0,0), RGB(.6,0,0), RGB(1,.6,.6)], [0., .4999, .5, .5001, 1.], categorical = false)
    if error != [] 
        if ytitle=="C"#= || ytitle=="Susc"=#; lim=(max(0, minimum(vect[i,:,:])),maximum(vect[i,3:end,:])+maximum(error[i,3:end,:]))
        else; lim = :auto
        end
        for z in eachindex(d)
            rescaleE=0          # To rescale the Energy (only when y = Energy)
            endT = length(vect[i,:,z])
            while endT >= 1 && vect[i,endT,z]==0
                endT -= 1
            end
            if ytitle=="E" && d[z]<0; rescaleE = d[z]; end
            x = PlotColord(d[z], d[1], d[end])
            plot!(T_for_d[d[z]], vect[i, 1:endT, z].-rescaleE, yerr = Vector(error[i, 1:endT, z]), label="$(d[z])", ylims=lim, seriescolor = pal[x], linecolor = pal[x], markercolor = pal[x], markerstrokecolor = pal[x], ecolor = pal[x])
        end
    else
        for z in eachindex(d)
            endT = length(vect[i,:,z])
            while endT >= 1 && vect[i,endT,z]==0
                endT -= 1
            end
            plot!(T_for_d[d[z]], vect[i, 1:endT, z], label="$(d[z])", seriescolor = pal[z])
        end
    end
    title!(p, ytitle*" as a function of T")
    xlabel!(p, "Temperature")
    ylabel!(p, ytitle)
    if save == true
        savefig("Plot/"*ytitle*".pdf")
    end
    display(p)
end

function PlotColord(x, dmin, dmax)
    if x < 0
        return .5 - .5* sqrt(x / dmin)
    elseif x > 0
        return .5 + .5 * sqrt(x / dmax)
    else
        return .5
    end
end

function Errorpropagation(v, Δ) # propagation of error for c and χ from the error on E and M
    return sqrt(sum((v .- mean(v)).^2)/length(v))*2*Δ
end

function interpmax(L, T, y) # interpolate and find the max
    Max = []
    xinterp = collect(.1:.001:2)
    p=plot()
    for l in eachindex(L)
        yinterp=Spline1D(T, y[l,:], k=3)(xinterp)
        push!(Max, (findmax(yinterp[301:end])[1], xinterp[findmax(yinterp[301:end])[2]]))
    end
    return Max
end

function linear(x, p) # linear function
    return p[1]*x + p[2]   
end

function crit(L, T, y, graph="") # calculate α & γ
    ymax = interpmax(L, T, y)
    println(ymax)
    fit =  curve_fit(linear, log.(L), log.(first.(ymax)), if (y=="c"); [.1,-.7]; elseif (y=="susc"); [2.,-4.]; else [.9,-.9] end)
    println(typeof(fit), "\n", fit)
    m,p = coef(fit)
    σ = stderror(fit)
    if graph!=""
        a=plot()
        plot!(a, log.(L), log.(first.(ymax)), seriestype=:scatter)
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

function Correlation(spin1::Vector, spin2::Vector)  # Correlation between two spins
    return sin(spin1[2])*sin(spin2[2])*cos(spin1[1]-spin2[1]) + cos(spin1[2])*cos(spin2[2])  
end

function MeanCorrTime1(Lattices)
    L = length(Lattices[1][:,1,1])
    n = length(Lattices)
    Corr = zeros(n-1)
    for k=1:n-1
        for l=k+1:n
            C = 0
            for i=1:L
                for j=1:L
                    C += Correlation(Lattices[k][:,i,j],Lattices[l][:,i,j])
                end
            end
            Corr[l-k] += abs(C)/(n-l+k)
        end
    end
    return Corr*L^(-2)
end

function MeanCorrTime(Lattices)
    L = length(Lattices[1][:,1,1])
    n = length(Lattices)
    Corr = zeros(n-1)
    for k=2:n
        C = 0
        for i=1:L
            for j=1:L
                C += Correlation(Lattices[k][:,i,j],Lattices[1][:,i,j])
            end
        end
        Corr[k-1] = C
    end
    return Corr*L^(-2)
end

function CorrPlot(AllLattices, T_for_d, dvect, d, L, l,save=false)
    p=plot()
    pal = cgrad([:lightblue, :green, :red])
    for i in eachindex(T_for_d[dvect[d]])
        MCT = MeanCorrTime(AllLattices[l,i,d])
        x = T_for_d[dvect[d]][i]/T_for_d[dvect[d]][end]
        plot!(collect(1:length(MCT))*10, MCT, label="$(T_for_d[dvect[d]][i])", seriescolor = pal[x], linecolor = pal[x], markercolor = pal[x], markerstrokecolor = pal[x], ecolor = pal[x])
    end
    title!(p, "Correlation time, L=$(L[l]), d=$(dvect[d])")
    if save ==true
        savefig("Plot/CorrelationTime$(dvect[d]).pdf")
    end
    display(p)
end

# α : specific heat                     Ising : 0       XY : NOP Essential singularity
# β : zero field mag                    Ising : 1/8     XY : NOP No magnetization (To have a nonzero 𝑀, the correlation function must approach a constant at large distance. But in the 2D XY model: For 𝑇>𝑇BKT: correlations decay exponentially. For 𝑇<𝑇BKT: correlations decay as a power law)
# γ : zero field isothermal suscxx      Ising : 7/4     XY : NOP Essential divergence
# δ : Critical isothermal               Ising : 15      XY : 15
# ν : corr length                       Ising : 1       XY : NOP 𝜉 diverges exponentially

# continuous symmetries cannot be spontaneously broken in 2D