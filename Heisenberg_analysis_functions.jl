function Get_T_d(N::Int64, L::Vector{Int64}, End::String="")  # for the first lattice size (supposing they all have the same T and d), gives which T have each d
    files = glob("$(L[1])*.jld2", ("Data/$(L)_$(N)$End"))
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

    all_d = unique(sort([p.n3 for p in parsed]))    # Vector of all d
    all_T = unique(sort([p.n2 for p in parsed]))    # Vector of all T
    T_for_d = Dict{Float32, Vector{Float32}}() # Dictionary mapping each #3 -> list of #2
    for p in parsed
        push!(get!(T_for_d, p.n3, Float32[]), p.n2)
    end
    T_for_d = Dict(k => sort(v) for (k,v) in T_for_d)
    return all_T, all_d, T_for_d
end

function Binor(arr::Vector{Float32}, Nbin::Int64, Nperbin::Int64) # take them Nbin by Nbin and average each the bins.
    Bins=[]
    for i=1:Nbin
        push!(Bins, mean(arr[((i-1)*Nperbin)+1:i*Nperbin]))
    end
    return Bins
end

function Binor2nd(arr::Vector{Float32}, Nbin::Int64, Nperbin::Int64, T::Float32, L::Int64) # errorbar for c and Χ
    Bins=[]
    for i=1:Nbin
        push!(Bins, (mean(arr[(i-1)*Nperbin+1:i*Nperbin].^2) - mean(arr[(i-1)*Nperbin+1:i*Nperbin])^2)/T*L^2)
    end
    return Bins
end

function Errorpropagation(v::Vector{Any}, Δ::Float32) # propagation of error for c and χ from the error on E and M
    return sqrt(sum((v .- mean(v)).^2)/length(v))*2*Δ
end

function plotL(L::Vector{Int64}, T_for_d::Dict{Float32, Vector{Float32}}, dvect::Vector{Float32}, x::Array{Float32, 3}, d::Int64, ytitle::String="", error=[], save::Bool=false)   # read data from a file and plot it
    endT = length(x[1,:,d])
    while endT >= 1 && x[1,endT,d]==0
        endT -= 1
    end
    p=plot()
    if typeof(x) == Matrix{Any}     # For Correlation length
        for t=1:max(round(Int, length(T_for_d[dvect[d]])/8),1):length(T_for_d[dvect[d]])
            plot!(collect(1:Int(L[end]/2-1)), x[end,:][t], label="$(T[t])",legend=:topright) # yet not update to work
        end
        title!(p, "Correlation (L=$(L[end])) as a function of distance")
        xlabel!(p, "Distance (site)")
    else
        if error != []  # if there are error bars
            if ytitle=="C" || ytitle=="Susc"; lim=(max(0, minimum(x[:,:,d])),maximum(x[:,3:end,d])+maximum(error[:,3:end,d]))
            else; lim = :auto
            end
            for l in eachindex(L)
                plot!(T_for_d[dvect[d]], x[l, 1:endT,d], yerr = Vector(error[l,1:endT,d]), markerstrokecolor=:auto, label="$(L[l])", ylims = lim)
            end
        else            # if there are no error bars
            for l in eachindex(L)
                plot!(T_for_d[dvect[d]], x[l, 1:endT,d], label="$(L[l])", ylims=(max(0, minimum(x[:,1:endT,i])),maximum(x[:,1:endT,i])))
            end
        end
        title!(p, ytitle*" as a function of T with d = $(dvect[d])")
        xlabel!(p, "Temperature")
    end
    ylabel!(p, ytitle)
    if save ==true
        savefig("Plot/"*ytitle*".pdf")
    end
    # println("size vect: ", size(x))
    # println("endT: ", endT)
    display(p)
end

function plotd(L::Vector{Int64}, T_for_d::Dict{Float32, Vector{Float32}}, d::Vector{Float32}, x::Array{Float32, 3}, ytitle::String="", error=[], i::Int64=length(x[:,1,1]), save::Bool=false)   # read data from a file and plot it
    p=plot()
    pal = cgrad([RGB(.4,.6,1), RGB(0,0,.5), RGB(0,0,0), RGB(.6,0,0), RGB(1,.6,.6)], [0., .4999, .5, .5001, 1.], categorical = false)
    if error != [] 
        if ytitle=="C"; lim=(max(0, minimum(x[i,:,:])),maximum(x[i,3:end,:])+maximum(error[i,3:end,:]))
        else; lim = :auto
        end
        for z in eachindex(d)
            rescaleE=0          # To rescale the Energy (only when y = Energy)
            endT = length(x[i,:,z])
            while endT >= 1 && x[i,endT,z]==0
                endT -= 1
            end
            if ytitle=="E" && d[z]<0; rescaleE = d[z]; end
            n = PlotColord(d[z], d[1], d[end])
            plot!(T_for_d[d[z]], x[i, 1:endT, z].-rescaleE, yerr = Vector(error[i, 1:endT, z]), label="$(d[z])", ylims=lim, seriescolor = pal[n], linecolor = pal[n], markercolor = pal[n], markerstrokecolor = pal[n], ecolor = pal[n])
        end
    else
        if ytitle=="C"; lim=(max(0, minimum(x[i,:,:])),maximum(x[i,3:end,:]))
        end
        for z in eachindex(d)
            endT = length(x[i,:,z])
            while endT >= 1 && x[i,endT,z]==0
                endT -= 1
            end
            rescaleE=0
            if ytitle=="E" && d[z]<0; rescaleE = d[z]; end
            n = PlotColord(d[z], d[1], d[end])
            plot!(T_for_d[d[z]], x[i, 1:endT, z].-rescaleE, label="$(d[z])", seriescolor = pal[n])
        end
    end
    title!(p, ytitle*" as a function of T for $(L[i])x$(L[i]) lattices")
    xlabel!(p, "Temperature")
    ylabel!(p, ytitle)
    if save; savefig("Plot/"*ytitle*".pdf"); end
    display(p)
end

function PlotColord(x::Float32, dmin::Float32, dmax::Float32)
    if x < 0
        return .5 - .5* sqrt(x / dmin)
    elseif x > 0
        return .5 + .5 * sqrt(x / dmax)
    else
        return .5
    end
end

function interpmax(l::Int64, T::Vector{Float32}, y::Array{Float32,2}) # interpolate and find the Tmax (for all Lattice size)
    xinterp = collect(.5:.001:2)
    yinterp=Spline1D(T, y[l,:], k=3)(xinterp)
    return (findmax(yinterp)[1], xinterp[findmax(yinterp)[2]])
end

function linear(x::Vector, p::Vector{Float64}) # linear function
    return p[1]*x .+ p[2]
end

function crit(L::Vector{Int64}, T::Vector{Float32}, y::Array{Float32,2}, title::String="") # calculate α & γ, plot if title is a non-empty string
    ymax = Vector{Tuple}(undef, length(L))
    for l in eachindex(L)
        ymax[l] = interpmax(l, T, y)
    end
    println(ymax)
    fit =  curve_fit(linear, log.(L), log.(first.(ymax)), if (y=="c"); [.1,-.7]; elseif (y=="susc"); [2.,-4.]; else [.9,-.9] end)
    println(#=typeof(fit), "\n",=# fit)
    m, p = coef(fit)
    σm, σp = stderror(fit)
    if title != ""
        a = plot()
        plot!(a, log.(L), log.(first.(ymax)), seriestype=:scatter)
        plot!(a, log.(L), m*log.(L) .+p)
        xlabel!(a, "Ln(Lattice length)")
        ylabel!(a, "Ln(max("*title*"))")
        display(a)
    end
    return m, σm
end

function critlength(T::Vector{Float32}, data::Vector, t::Real, PLOT::Bool=true, ln::Bool=false) # Calculate critical exp or correlation length (by default)
    closestT = argmin(abs.(T .- t)); println(T[closestT])
    data=data[closestT]
    negativ = findfirst(x -> x < 4e-4, data)
    if typeof(negativ) == Nothing; n = length(data); else; n=negativ-1; end
    println(n)
    x=collect(1:n)
    fit =  curve_fit(linear, log.(x)*ln + !ln*x, log.(data[1:n]), [-.09,9.])
    m,p = coef(fit); σ = stderror(fit)

    a=plot()
    plot!(a, log.(x)*ln + !ln*x, log.(data[1:n]), seriestype=:scatter, label="data")
    plot!(a, log.(x)*ln + !ln*x, m*(log.(x)*ln + !ln*x) .+p, label="fit")
    ylabel!(a, "Ln(Correlation)")
    if ln == true; xlabel!(a, "Ln(Distance)"); ξ=σξ=0
    else; xlabel!(a, "Distance"); ξ =-1/m; σξ = σ[1]/m^2
    end
    title!("Correlation vs distance at T = $(T[closestT])")
    if PLOT==true;    display(a);   end
    return m*ln + !ln*ξ, σ[1]*ln+!ln*σξ
end

function CritLength(corr::Vector, L::Int64)
    model(x, p) = CorrelationFunction(L, x, p)
    xdata = Int.(collect(1:L/2-1))
    fit = curve_fit(model, xdata, corr, [3.1, 3.1])    # [..,..] random guess
    p=plot(xdata, corr)
    plot!(xdata, CorrelationFunction(L, xdata, coef(fit)))
    display(p)
    return coef(fit), stderror(fit)
end

function CorrelationFunction(L::Int64, x::Vector{Int64}, p::Vector{Float64})
    return exp.(2π*p[2].*G(x)).*(exp.(-x/p[1])+exp.(-(L.+x)/p[1]))
end

function G(x::Vector{Int64}, Nk::Int64=100)
    y=zeros(length(x))
    for k=1:Nk
        y += (cos.(k*x).-1) / (1-cos(k))
    end
    return y/Nk
end

function Correlation(spin1::Vector, spin2::Vector)  # Correlation between two spins
    return sin(spin1[2])*sin(spin2[2])*cos(spin1[1]-spin2[1]) + cos(spin1[2])*cos(spin2[2])  
end

function MeanCorrTime(Lattices::Vector{Array{Float32, 3}})
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

function CorrTimePlot(AllLattices::Array{Vector{Array{Float32, 3}},3}, T_for_d, dvect::Vector{Float32}, d::Int64, L::Vector{Int64}, l::Int64, save::Bool=false)
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

function Plot_Max_C_Χ(d::Vector{Float32}, T_for_d::Dict{Float32, Vector{Float32}}, x::Array{Float32, 3}, L::Vector{Int64}, title::String="C/χ", save::Bool=false)
    p=plot()
    xmax, Tmax = Array{Float64}(undef, length(d), length(L)), Array{Float64}(undef, length(d), length(L))
    for l in eachindex(L)
        for z in eachindex(d)
            xmax[z,l], Tmax[z,l] = interpmax(l, T_for_d[d[z]], x[:,:,z])
        end
        plot!(d, Tmax[:,l], seriestype=:scatter, label="$(L[l])")
    end
    title!("Tpic of $title vs d")
    ylabel!("Temperature of the maximum of $title")
    xlabel!("⬅ Ising                              Value of d                              XY ➡")
    if save; savefig("Plot/Tpic_$title.pdf"); end
    display(p)
    p=plot()
    for l in eachindex(L)
        plot!(d, xmax[:,l], seriestype=:scatter, label="$(L[l])")
    end
    title!("$title max vs d")
    ylabel!("Mximum of $title")
    xlabel!("⬅ Ising                              Value of d                              XY ➡")
    if save; savefig("Plot/$(title)_max.pdf"); end
    display(p)
    return Tmax, xmax
end

function Plot_Max_ξ(d::Vector{Float32}, T_for_d::Dict{Float32, Vector{Float32}}, Corr::Array{Vector{Float32}, 2}, ln::Bool=false)
    Tmax = zeros(length(d))
    for i in eachindex(d)
        critic = zeros(length(T_for_d[d[i]]))
        for j in eachindex(T_for_d[d[i]])
            critic[j], a = critlength(T_for_d[d[i]], Corr[:,i], T_for_d[d[i]][j], true, ln)
        end
        println(critic)
        Tmax[i] = round(T_for_d[d[i]][argmax(critic)];digits=3)
    end
    p=plot(d, Tmax)
    xlabel!("⬅ Ising                              Value of d                              XY ➡")
    ylabel!("Temperature of the maximum of ξ")
    display(p)
    return Tmax
end
# α : specific heat                     Ising : 0       XY : NOP Essential singularity
# β : zero field mag                    Ising : 1/8     XY : NOP No magnetization (To have a nonzero 𝑀, the correlation function must approach a constant at large distance. But in the 2D XY model: For 𝑇>𝑇BKT: correlations decay exponentially. For 𝑇<𝑇BKT: correlations decay as a power law)
# γ : zero field isothermal suscxx      Ising : 7/4     XY : NOP Essential divergence
# δ : Critical isothermal               Ising : 15      XY : 15
# ν : corr length                       Ising : 1       XY : NOP 𝜉 diverges exponentially

# continuous symmetries cannot be spontaneously broken in 2D