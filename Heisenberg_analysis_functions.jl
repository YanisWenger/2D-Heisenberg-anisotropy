function Get_T_DL(N::Int64, L::Vector{Int64}, End::String="")  # for each lattice size gives which T have each d
    T_for_DL = Dict{Tuple{Int64, Float32}, Vector{Float32}}()
    all_D = []
    all_T = []
    for l in L
        files = glob("$l*.jld2", ("Data/$(L)_$(N)$End"))
        for f in files
            name = splitext(last(splitpath(f)))[1]
            parts = split(name, "_")
            if length(parts) != 3
                @warn "Bad filename format" f parts
            end
        end
        parsed = map(files) do f
            name = splitext(last(splitpath(f)))[1]
            nums = parse.(Float32, split(name, "_"))
            (; n2 = nums[2], n3 = nums[3])  # 2: T, 3: D
        end
        all_D = unique(vcat(all_D, sort([p.n3 for p in parsed])))    # Vector of all d
        all_T = unique(vcat(all_T, sort([p.n2 for p in parsed])))    # Vector of all T
        for p in parsed
            push!(get!(T_for_DL, (l, p.n3), Float32[]), p.n2)
        end
    end
    T_for_DL = Dict(k => sort(v) for (k,v) in T_for_DL)
    return Float32.(all_T), Float32.(all_D), T_for_DL
end

function Binor(arr::Vector{Float32}, Nbin::Int64, Nperbin::Int64) # take measured quantities (E, M) Nbin by Nbin, average each the bins and get the error bars
    Bins=[]
    for i=1:Nbin
        push!(Bins, mean(arr[((i-1)*Nperbin)+1:i*Nperbin]))
    end
    return Bins
end

function Binor2nd(arr::Vector{Float32}, Nbin::Int64, Nperbin::Int64, T::Float32, L::Int64) # error bars for c and Χ
    Bins=[]
    for i=1:Nbin
        push!(Bins, (mean(arr[(i-1)*Nperbin+1:i*Nperbin].^2) - mean(arr[(i-1)*Nperbin+1:i*Nperbin])^2)/T*L^2)
    end
    return Bins
end

function Errorpropagation(v::Vector{Any}, Δ::Float32) # propagation of error for c and χ from the error on E and M
    return sqrt(sum((v .- mean(v)).^2)/length(v))*2*Δ
end

function plotL(L::Vector{Int64}, T_for_DL::Dict{Tuple{Int64, Float32}, Vector{Float32}}, D::Vector{Float32}, x::Array{Float32, 3}, d::Int64, ytitle::String="", error=[], save::Bool=false)   # read data from a file and plot it vs L
    xinterp = collect(T[1]:.002:T[end])
    p=Plots.plot()
    hex_codes = ["#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#999999", "#F781BF", "#A65628", "#FFFF33"]
    pal = [parse(Colorant, hex_codes[i]) for i in 1:length(L)]
    if typeof(x) == Matrix{Any}     # For Correlation length
        for t=1:max(round(Int, length(T_for_DL[(L[1], D[d])])/8),1):length(T_for_DL[(L[1], D[d])])
            Plots.plot!(collect(1:Int(L[end]/2-1)), x[end,:][t], label="$(T[t])",legend=:topright) # yet not update to work
        end
        Plots.title!(p, "Correlation (L=$(L[end])) as a function of distance")
        Plots.xlabel!(p, "Distance (site)")
    else
        # if ytitle=="C"; lim=(max(0, minimum(x[:,:,d])),maximum(x[:,3:end,d]))
        # else
            lim = :auto
        # end
        for l in eachindex(L)
            T = T_for_DL[(L[l], D[d])]
            endT = length(T)
            yinterp = Spline1D(T, x[l, 1:endT, d], k=3)(xinterp)
            Plots.plot!(xinterp, yinterp, label="", ylims=lim, seriescolor = pal[l], linecolor = pal[l], markercolor = pal[l], markerstrokecolor = pal[l], ecolor = pal[l])
            if error != []  # if there are error bars
                Plots.scatter!(T, x[l, 1:endT,d], yerr = Vector(error[l,1:endT,d]), label="$(L[l])", ylims = lim, seriescolor = pal[l], linecolor = pal[l], markercolor = pal[l], markerstrokecolor = pal[l], ecolor = pal[l])
            else        # if there are no error bars
                Plots.scatter!(T, x[l, 1:endT,d], label="$(L[l])", ylims=lim, seriescolor = pal[l])
            end
        end
        Plots.title!(p, ytitle*" as a function of T with D = $(D[d])")
        Plots.xlabel!(p, "Temperature")
    end
    Plots.ylabel!(p, ytitle)
    if save ==true
        Plots.savefig("Plot/"*ytitle*".pdf")
    end
    # println("size vect: ", size(x))
    # println("endT: ", endT)
    display(p)
end

function plotD(pal::PlotUtils.ContinuousColorGradient, L::Vector{Int64}, T_for_DL::Dict{Tuple{Int64, Float32}, Vector{Float32}}, D::Vector{Float32}, x::Array{Float32, 3}, ytitle::String="", error=[], l::Int64=length(x[:,1,1]), save::Bool=false)   # read data from a file and plot it vs D
    p2 = Plots.plot(load("Colorbar.png"), framestyle = :none, axis = nothing)
    p1 = Plots.plot()
    # if ytitle=="C"; lim=(max(0, minimum(x[l,:,:])),maximum(x[l,3:end,:]))
    #=else;=# lim = :auto
    # end
    for d in eachindex(D)
        rescaleE=0          # To rescale the Energy (only when y = Energy)
        if ytitle=="E" && D[d]<0; rescaleE = D[d]; end
        n = PlotcolorD(D[d], D[1], D[end])
        T = T_for_DL[(L[l], D[d])]
        endT=length(T)
        xinterp = collect(T[1]:.002:T[end])
        yinterp = Spline1D(T, x[l, 1:endT, d], k=3)(xinterp)
        Plots.plot!(xinterp, yinterp, label=""#=$(D[d])=#, ylims=lim, seriescolor = pal[n], linecolor = pal[n], markercolor = pal[n], markerstrokecolor = pal[n], ecolor = pal[n])
        if error != [] 
            Plots.plot!(T, x[l, 1:endT, d].-rescaleE, yerr = Vector(error[l, 1:endT, d]), label=""#=$(D[d])=#, ylims=lim, seriescolor = pal[n], linecolor = pal[n], markercolor = pal[n], markerstrokecolor = pal[n], ecolor = pal[n], seriestype=:scatter, ms=3)
        else
            Plots.plot!(T, x[l, 1:endT, d].-rescaleE, label=""#=$(D[d])=#, seriescolor = pal[n], seriestype=:scatter, ms=3)
        end
    end
    Plots.title!(p1, ytitle*" as a function of T for $(L[l])x$(L[l]) lattices")
    Plots.xlabel!(p1, "Temperature")
    Plots.ylabel!(p1, ytitle)
    p=Plots.plot(p1, p2, layout= @layout [a{.88w} b{.12w}])
    if save; Plots.savefig("Plot/"*ytitle*".pdf"); end
    display(p)
end

function PlotcolorD(x::Float32, dmin::Float32, dmax::Float32)   # associate a D value with a color
    if x < 0
        return .5 - .5* (x / dmin)^.3
    elseif x > 0
        return .5 + .5 * (x / dmax)^.3
    else
        return .5
    end
end

function interpmax(T::Vector{Float32}, y::Vector{Float32}) # interpolate and find the Tmax (for all Lattice size)
    xinterp = collect(.5:.001:2)
    yinterp=Spline1D(T, y, k=3)(xinterp)
    return (findmax(yinterp)[1], xinterp[findmax(yinterp)[2]])
end

function interpmax_with_error(T::Vector{Float32}, y::Vector{Float32}, y_err::Vector{Float32}; n_bootstrap=1000) # interpolate and find the Tmax (for all Lattice size) then bootstrap resampling to get error bars
    xinterp = collect(T[1]:.002:T[end])
    max_values = Float32[]
    max_positions = Float32[]
    
    for i in 1:n_bootstrap
        y_perturbed = y .+ randn(length(y)) .* y_err # Perturb y values according to their uncertainties
        try # Fit spline and find max
            yinterp = Spline1D(T, y_perturbed, k=3)(xinterp)
            idx = argmax(yinterp)
            push!(max_values, yinterp[idx])
            push!(max_positions, xinterp[idx])
        catch
            continue # Skip failed fits (rare)
        end
    end    
    return (median(max_values), median(max_positions), std(max_values), std(max_positions))
end

function crit(L::Vector{Int64}, T::Vector{Float32}, y::Array{Float32,2}, title::String="") # calculate α & γ, plot if title is a non-empty string
    ymax = Vector{Tuple}(undef, length(L))
    for l in eachindex(L)
        ymax[l] = interpmax(T, y[l,:])
    end
    println(ymax)
    fit =  curve_fit(linear_fit, log.(L), log.(first.(ymax)), if (y=="c"); [.1,-.7]; elseif (y=="susc"); [2.,-4.]; else [.9,-.9] end)
    println(#=typeof(fit), "\n",=# fit)
    m, p = coef(fit)
    σm, σp = stderror(fit)
    if title != ""
        a = Plots.plot()
        Plots.plot!(a, log.(L), log.(first.(ymax)), seriestype=:scatter)
        Plots.plot!(a, log.(L), m*log.(L) .+p)
        Plots.xlabel!(a, "Ln(Lattice length)")
        Plots.ylabel!(a, "Ln(max("*title*"))")
        display(a)
    end
    return m, σm
end

function critlength(T::Vector{Float32}, data::Vector, t::Real, PLOT::Bool=true, ln::Bool=false) # Calculate critical exp or correlation length (by default)
    closestT = argmin(abs.(T .- t))     # Compute the temperature closest to the requested one
    data=data[closestT] 
    negativ = findfirst(x -> x < 4e-4, data)    # don't take in account when the correlation is below than 4e-4 (0.07)
    if typeof(negativ) == Nothing; n = length(data)
    else; n=negativ-1; end
    println(n)
    x=collect(1:n)
    fit =  curve_fit(linear_fit, log.(x)*ln + !ln*x, log.(data[1:n]), [-.09,9.])
    m,p = coef(fit); σ = stderror(fit)

    a=Plots.plot()
    Plots.plot!(a, log.(x)*ln + !ln*x, log.(data[1:n]), seriestype=:scatter, label="data")
    Plots.plot!(a, log.(x)*ln + !ln*x, m*(log.(x)*ln + !ln*x) .+p, label="fit")
    Plots.ylabel!(a, "Ln(Correlation)")
    if ln == true; Plots.xlabel!(a, "Ln(Distance)"); ξ=σξ=0
    else; Plots.xlabel!(a, "Distance"); ξ =-1/m; σξ = σ[1]/m^2
    end
    Plots.title!("Correlation vs distance at T = $(T[closestT])")
    if PLOT==true;    display(a);   end
    return m*ln + !ln*ξ, σ[1]*ln+!ln*σξ
end

function CritLength(corr::Vector, L::Int64)
    model(x, p) = CorrelationFunction(L, x, p)
    xdata = Int.(collect(1:L/2-1))
    fit = curve_fit(model, xdata, corr, [3.1, 3.1])    # [..,..] random guess
    p=Plots.plot(xdata, corr, label="data")
    Plots.plot!(xdata, CorrelationFunction(L, xdata, coef(fit)), label="fit")
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

function Correlation(spin1::Vector, spin2::Vector)  # Correlation between two spins (i.e. dot product)
    return sin(spin1[2])*sin(spin2[2])*cos(spin1[1]-spin2[1]) + cos(spin1[2])*cos(spin2[2])  
end

function MeanCorrTime(Lattices::Vector{Array{Float32, 3}})  # calculate the mean correlation time
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

function CorrTimePlot(AllLattices::Array{Vector{Array{Float32, 3}},3}, T_for_DL::Dict{Tuple{Int64, Float32}, Vector{Float32}}, D::Vector{Float32}, d::Int64, L::Vector{Int64}, l::Int64, save::Bool=false)
    p=Plots.plot()
    pal = cgrad([:lightblue, :green, :red])
    for i in eachindex(T_for_DL[(L[l], D[d])])
        MCT = MeanCorrTime(AllLattices[l,i,d])
        x = T_for_DL[(L[l], D[d])][i]/T_for_DL[(L[l], D[d])][end]
        Plots.plot!(collect(1:length(MCT))*10, MCT, label="$(T_for_DL[(L[l], D[d])][i])", seriescolor = pal[x], linecolor = pal[x], markercolor = pal[x], markerstrokecolor = pal[x], ecolor = pal[x])
    end
    Plots.title!(p, "Correlation time, L=$(L[l]), D=$(D[d])")
    if save ==true
        Plots.savefig("Plot/CorrelationTime$(D[d]).pdf")
    end
    display(p)
end

function Plot_Max_C_Χ(pal::PlotUtils.ContinuousColorGradient, D::Vector{Float32}, T_for_DL::Dict{Tuple{Int64, Float32}, Vector{Float32}}, Y::Array{Float32, 3}, Yerr::Array{Float32, 3}, L::Vector{Int64}, title::String="C/χ", save::Bool=false)
    p=Plots.plot() # xmax & T of xmax vs D
    xmax, Tmax, xmaxerr, Tmaxerr = Array{Float64}(undef, length(D), length(L)), Array{Float64}(undef, length(D), length(L)), Array{Float64}(undef, length(D), length(L)), Array{Float64}(undef, length(D), length(L))
    for l in eachindex(L)
        for d in eachindex(D)
            xmax[d,l], Tmax[d,l], xmaxerr[d,l], Tmaxerr[d,l] = interpmax_with_error(T_for_DL[(L[l], D[d])], filter(!iszero, Y[l,:,d]), filter(!iszero, Yerr[l,:,d]))
        end
        Plots.plot!(D, Tmax[:,l], seriestype=:scatter, label="$(L[l])", yerr=Tmaxerr[:,l])
    end
    Plots.title!("Tpic of $title vs D")
    Plots.ylabel!("Temperature of the maximum of $title")
    Plots.xlabel!("⬅ Ising                              Value of D                              XY ➡")
    if save; Plots.savefig("Plot/Tpic_$title.pdf"); end
    display(p)
    p=Plots.plot()
    for l in eachindex(L)
        Plots.plot!(D, xmax[:,l], seriestype=:scatter, label="$(L[l])", yerr=xmaxerr[:,l])
    end
    Plots.title!("$title max vs D")
    Plots.ylabel!("Maximum of $title")
    Plots.xlabel!("⬅ Ising                              Value of D                              XY ➡")
    if save; Plots.savefig("Plot/$(title)_max.pdf"); end
    display(p)
    
    fit_Tc = []
    fit_ln = []
    fit_power = []
    x = 8:0.1:70
    p2 = Plots.plot(load("Colorbar.png"), framestyle = :none, axis = nothing)
    n=zeros(length(D))
    for d in eachindex(D)
        n[d] = PlotcolorD(D[d], D[1], D[end])
    end
    p1=Plots.plot() # xmax & T of xmax vs L
    for d in eachindex(D)
        for l in eachindex(L)
            xmax[d,l], Tmax[d,l], xmaxerr[d,l], Tmaxerr[d,l] = interpmax_with_error(T_for_DL[(L[l], D[d])], filter(!iszero, Y[l,:,d]), filter(!iszero, Yerr[l,:,d]))
        end
        if title=="Susc"
            fitTc =  curve_fit(power_fit_plus, L, Tmax[d,:], [-1.1, .9, .9])
            push!(fit_Tc, (round.(coef(fitTc);digits=3), round.(stderror(fitTc);digits=3)))
            Plots.plot!(x, coef(fitTc)[3]*x.^coef(fitTc)[1].+coef(fitTc)[2], seriescolor=pal[n[d]], label="")
        end
        Plots.plot!(L, Tmax[d,:], seriestype=:scatter, seriescolor=pal[n[d]], label="", yerr=Tmaxerr[d,:])
    end
    Plots.title!("Tpic of $title vs L")
    Plots.ylabel!("Temperature of the maximum of $title")
    Plots.xlabel!("Lattice length")

    p=Plots.plot(p1, p2, layout= @layout [a{.88w} b{.12w}])
    if save; Plots.savefig("Plot/TpicvsL_$title.pdf"); end
    display(p)
    p1=Plots.plot()
    for d in eachindex(D)
        fitln =  curve_fit(ln_fit, L, xmax[d,:], [1.1, 1.8])
        fitpower = curve_fit(power_fit, L, xmax[d,:], [.8, 1.1])
        push!(fit_ln, (round(coef(fitln)[1];digits=3), round(stderror(fitln)[1];digits=3)))
        push!(fit_power, (round.(coef(fitpower);digits=3), round.(stderror(fitpower);digits=3)))

        n = PlotcolorD(D[d], D[1], D[end])
        Plots.plot!(L, xmax[d,:], seriestype=:scatter, label="", z=n, seriescolor=pal[n], yerr=xmaxerr[d,:], xaxis=:log)
        if title=="C"
            Plots.plot!(x, coef(fitln)[1]*log.(x) .+coef(fitln)[2], label="", z=n, seriescolor=pal[n], xaxis=:log)
        else
            Plots.plot!(x, coef(fitpower)[2]*x.^coef(fitpower)[1], label="", z=n, seriescolor=pal[n])
        end
    end
    Plots.plot!(p1, colorbar=true, colorbar_title="D Value")
    Plots.title!("$(title)_max vs L")
    Plots.ylabel!("Maximum of $title")
    Plots.xlabel!("Lattice length")

    p=Plots.plot(p1, p2, layout= @layout [a{.88w} b{.12w}])
    if save; Plots.savefig("Plot/$(title)_max_vsL.pdf"); end
    display(p)
    return Tmax, xmax, fit_ln, fit_power, fit_Tc
end

function Plot_Max_ξ(L::Vector{Int64}, D::Vector{Float32}, T_for_DL::Dict{Tuple{Int64, Float32}, Vector{Float32}}, Corr::Array{Vector{Float32}, 3}, l::Int64, ln::Bool=false, only_one_d::Int64=0)
    Corr = Corr[l,:,:]
    if only_one_d == 0
        Tmax = zeros(length(D))
        for d in eachindex(D)
            critic = zeros(length(T_for_DL[(L[l], D[d])]))
            for j in eachindex(T_for_DL[(L[l], D[d])])
                critic[j], a = critlength(T_for_DL[(L[l], D[d])], Corr[:,d], T_for_DL[(L[l], D[d])][j], true, ln)
            end
            println(critic)
            Tmax[d] = round(T_for_DL[(L[l], D[d])][argmax(critic)];digits=3)
        end
        p=Plots.plot(D, Tmax)
        Plots.xlabel!("⬅ Ising                              Value of D                              XY ➡")
        Plots.ylabel!("Temperature of the maximum of ξ")
        display(p)
    else    # To be able to use this function for one D value at a time
        d = only_one_d
        critic = zeros(length(T_for_DL[(L[l], D[d])]))
        for j in eachindex(T_for_DL[(L[l], D[d])])
            critic[j], a = critlength(T_for_DL[(L[l], D[d])], Corr[:,d], T_for_DL[(L[l], D[d])][j], true, ln)
        end
    end
    return Tmax
end

function linear_fit(x::Vector, p::Vector{Float64}) # linear function
    return p[1]*x .+ p[2]
end

function ln_fit(x::Vector{Int64}, p::Vector{Float64})       # a * ln(x) + b
    return p[1]*log.(x) .+ p[2]
end

function power_fit(L::Vector{Int64}, p::Vector{Float64})    # a * L ^ b
    return p[2]*L.^p[1]
end

function power_fit_plus(L::Vector{Int64}, p::Vector{Float64})  # a * L ^ b + c
    return p[3]*L.^p[1].+p[2]
end

function colorbar(pal::PlotUtils.ContinuousColorGradient, D)
    fig = Figure(size = (62, 330))          # Make the wanted colorbar requires other packages. The colorbar is then saved and combined with the plot
    Colorbar(fig[1, 1]; colormap = pal, colorrange = (D[1], D[end]), ticks = ([D[1], D[1]/2, 0, D[end]/2, D[end]], string.(round.([D[1], D[1]*2^(-10/3), 0, D[end]*2^(-10/3), D[end]]; digits=2))), label = "← Ising                    D value                    XY →")
    Makie.save("Colorbar.png", fig; px_per_unit = 8)
end

# α : specific heat                     Ising : 0       XY : NOP Essential singularity
# β : zero field mag                    Ising : 1/8     XY : NOP No magnetization (To have a nonzero 𝑀, the correlation function must approach a constant at large distance. But in the 2D XY model: For 𝑇>𝑇BKT: correlations decay exponentially. For 𝑇<𝑇BKT: correlations decay as a power law)
# γ : zero field isothermal suscxx      Ising : 7/4     XY : NOP Essential divergence
# δ : Critical isothermal               Ising : 15      XY : 15
# ν : corr length                       Ising : 1       XY : NOP 𝜉 diverges exponentially

# continuous symmetries cannot be spontaneously broken in 2D






#    To analyse single D simulation, have to be suppressed at some point

function Plot_Max_C_Χ_singleD(D::Vector{Float32}, T_for_DL::Dict{Tuple{Int64, Float32}, Vector{Float32}}, Y::Array{Float32, 3}, Yerr::Array{Float32, 3}, L::Vector{Int64}, title::String="C/χ", save::Bool=false)
    p=Plots.plot() # xmax & T of xmax vs D
    xmax, Tmax, xmaxerr, Tmaxerr = Array{Float64}(undef, length(D), length(L)), Array{Float64}(undef, length(D), length(L)), Array{Float64}(undef, length(D), length(L)), Array{Float64}(undef, length(D), length(L))
    for l in eachindex(L)
        xmax[1,l], Tmax[1,l], xmaxerr[1,l], Tmaxerr[1,l] = interpmax_with_error(T_for_DL[(L[l], D[1])], Y[l,:,1], filter(!iszero, Yerr[l,:,1]))
        Plots.plot!(D, Tmax[:,l], seriestype=:scatter, label="$(L[l])", yerr=Tmaxerr[:,l])
    end
    Plots.title!("Tpic of $title vs D")
    Plots.ylabel!("Temperature of the maximum of $title")
    Plots.xlabel!("⬅ Ising                              Value of D                              XY ➡")
    if save; Plots.savefig("Plot/Tpic_$title.pdf"); end
    display(p)
    p=Plots.plot()
    for l in eachindex(L)
        Plots.plot!(D, xmax[:,l], seriestype=:scatter, label="$(L[l])", yerr=xmaxerr[:,l])
    end
    Plots.title!("$title max vs D")
    Plots.ylabel!("Maximum of $title")
    Plots.xlabel!("⬅ Ising                              Value of D                              XY ➡")
    if save; Plots.savefig("Plot/$(title)_max.pdf"); end
    display(p)
    
    fit_Tc = []
    fit_ln = []
    fit_power = []
    x = L[1]:0.1:L[end]
    n=zeros(length(D))
    p1=Plots.plot() # xmax & T of xmax vs L
    for d in eachindex(D)
        for l in eachindex(L)
            xmax[d,l], Tmax[d,l], xmaxerr[d,l], Tmaxerr[d,l] = interpmax_with_error(T_for_DL[(L[l], D[d])], Y[l,:,d], filter(!iszero, Yerr[l,:,d]))
        end
        if title!=""
            fitTc =  curve_fit(power_fit_plus, L, Tmax[d,:], [-1.1, .9, .9])
            push!(fit_Tc, (round.(coef(fitTc);digits=3), round.(stderror(fitTc);digits=3)))
            Plots.plot!(x, coef(fitTc)[3]*x.^coef(fitTc)[1].+coef(fitTc)[2], label="")
            println("power: $(coef(fitTc)[1])\tTc: $(coef(fitTc)[2])")
        end
        Plots.plot!(L, Tmax[d,:], seriestype=:scatter, label="", yerr=Tmaxerr[d,:])
    end
    Plots.title!("Tpic of $title vs L")
    Plots.ylabel!("Temperature of the maximum of $title")
    Plots.xlabel!("Lattice length")
    if save; Plots.savefig("Plot/TpicvsL_$title.pdf"); end
    display(p1)
    p1=Plots.plot()
    for d in eachindex(D)
        fitln =  curve_fit(ln_fit, L, xmax[d,:], [1.1, 1.8])
        fitpower = curve_fit(power_fit, L, xmax[d,:], [.8, 1.1])
        push!(fit_ln, (round(coef(fitln)[1];digits=3), round(stderror(fitln)[1];digits=3)))
        push!(fit_power, (round.(coef(fitpower);digits=3), round.(stderror(fitpower);digits=3)))

        Plots.plot!(L, xmax[d,:], seriestype=:scatter, label="", z=n, yerr=xmaxerr[d,:])
        if title=="C"
            Plots.plot!(x, coef(fitln)[1]*log.(x) .+coef(fitln)[2], label="", z=n)
        else
            Plots.plot!(x, coef(fitpower)[2]*x.^coef(fitpower)[1], label="", z=n)
        end
    end
    Plots.plot!(p1, colorbar=true, colorbar_title="D Value")
    Plots.title!("$(title)_max vs L")
    Plots.ylabel!("Maximum of $title")
    Plots.xlabel!("Lattice length")

    if save; Plots.savefig("Plot/$(title)_max_vsL.pdf"); end
    display(p1)
    return Tmax, xmax, fit_ln, fit_power, fit_Tc
end