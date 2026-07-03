function Get_T_LD(N::Int64, L::Vector{Int64}, End::String="")  # for each lattice size gives which T have each d
    T_for_LD = Dict{Tuple{Int64, Float32}, Vector{Float32}}()
    D_for_L = Dict{Int64, Vector{Float32}}()
    all_T = []
    for l in L
        files = glob("$(l)_*.jld2", ("Data/$(L)_$(N)$End"))
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
        all_T = unique(vcat(all_T, sort([p.n2 for p in parsed])))    # Vector of all T
        for p in parsed
            push!(get!(D_for_L, l, Float32[]), p.n3)
        end
        for p in parsed
            push!(get!(T_for_LD, (l, p.n3), Float32[]), p.n2)
        end
    end
    D_for_L = Dict(k => unique(sort(v)) for (k,v) in D_for_L)
    T_for_LD = Dict(k => sort(v) for (k,v) in T_for_LD)
    return Float32.(all_T), D_for_L, T_for_LD
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

# function Errorpropagation(v::Vector{Any}, Δ::Float32) # propagation of error for c from the error on E (binor2nd is used instead)
#     return sqrt(sum((v .- mean(v)).^2)/length(v))*2*Δ
# end

function plotL(L::Vector{Int64}, T_for_LD::Dict{Tuple{Int64, Float32}, Vector{Float32}}, D_for_L::Dict{Int64, Vector{Float32}}, x::Array{Float32, 3}, Dvalue::Real, ytitle::String="", error=[], save::Bool=false)   # read data from a file and plot it vs L
    p=Plots.plot()
    Dvalue = Float32(Dvalue)
    hex_codes = ["#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#999999", "#F781BF", "#A65628", "#4790FF","#BBBB33"]
    pal = [parse(Colorant, hex_codes[i]) for i in 1:length(L)]
    # if ytitle=="C"; lim=(max(0, minimum(x[:,:,d])),maximum(x[:,3:end,d]))
    # else
        lim = :auto
    # end
    for l in eachindex(L)
        if Dvalue ∈ D_for_L[L[l]]
            d = findfirst(==(Dvalue), D_for_L[L[l]])
            T = T_for_LD[(L[l], Dvalue)]
            endT = length(T)
            if error != []  # if there are error bars
                Plots.plot!(T, x[l, 1:endT,d], yerr = Vector(error[l,1:endT,d]), label="$(L[l])", ylims = lim, seriescolor = pal[l], linecolor = pal[l], markercolor = pal[l], markerstrokecolor = pal[l], ecolor = pal[l])
            else        # if there are no error bars
                Plots.plot!(T, x[l, 1:endT,d], label="$(L[l])", ylims=lim, seriescolor = pal[l])
            end
        else
            println("For L=$(L[l]), valid D are: ", D_for_L[L[l]])
        end
    end
    Plots.title!(p, ytitle*" as a function of T with D = $(Dvalue)")
    Plots.xlabel!(p, "Temperature")
    Plots.ylabel!(p, ytitle)
    if save ==true
        Plots.savefig("Plot/"*ytitle*"_d=$(Dvalue).pdf")
    end
    display(p)
end

function plotD(pal::PlotUtils.ContinuousColorGradient, L::Vector{Int64}, T_for_LD::Dict{Tuple{Int64, Float32}, Vector{Float32}}, D_for_L::Dict{Int64, Vector{Float32}}, x::Array{Float32, 3}, ytitle::String="", error=[], l::Int64=length(x[:,1,1]), save::Bool=false)   # read data from a file and plot it vs D
    p2 = Plots.plot(load("Colorbar.png"), framestyle = :none, axis = nothing)
    p1 = Plots.plot()
    # if ytitle=="C"; lim=(max(0, minimum(x[l,:,:])),maximum(x[l,3:end,:]))
    #=else;=# lim = :auto
    # end
    D = D_for_L[L[l]]
    for d in eachindex(D)
        rescaleE=0          # To rescale the Energy (only when y = Energy)
        if ytitle=="E" && D[d]<0; rescaleE = D[d]; end
        n = PlotcolorD(D[d], D[1], D[end])
        T = T_for_LD[(L[l], D[d])]
        endT=length(T)
        if error != [] 
            Plots.plot!(T, x[l, 1:endT, d].-rescaleE, yerr = Vector(error[l, 1:endT, d]), label=""#=$(D[d])=#, ylims=lim, seriescolor = pal[n], linecolor = pal[n], markercolor = pal[n], markerstrokecolor = pal[n], ecolor = pal[n])#, seriestype=:scatter, ms=3)
        else
            Plots.plot!(T, x[l, 1:endT, d].-rescaleE, label=""#=$(D[d])=#, seriescolor = pal[n])#, seriestype=:scatter, ms=3)
        end
    end
    Plots.title!(p1, ytitle*" as a function of T for $(L[l])x$(L[l]) lattices")
    Plots.xlabel!(p1, "Temperature")
    Plots.ylabel!(p1, ytitle)
    p=Plots.plot(p1, p2, layout= @layout [a{.88w} b{.12w}])
    if save; Plots.savefig("Plot/"*ytitle*".pdf"); end
    display(p)
end

function plotLCorrelation(L::Vector{Int64}, T_for_LD::Dict{Tuple{Int64, Float32}, Vector{Float32}}, D_for_L::Dict{Int64, Vector{Float32}}, Corr::Array{Vector{Float32}, 3}, Dvalue::Real, t::Real)
    p=Plots.plot()
    Dvalue = Float32(Dvalue)
    for l in eachindex(L)
        T=T_for_LD[L[l], Dvalue]
        if Dvalue ∈ D_for_L[L[l]]
            d = findfirst(==(Dvalue), D_for_L[L[l]])
            closestT = argmin(abs.(T .- t))     # Compute the temperature closest to the requested one
            Plots.plot!(collect(1:Int(L[l]/2-1)), Corr[l, closestT, d], label="L=$(L[l]), T=$(T[closestT])", legend=:topright)
        else
            println("For L=$(L[l]), valid D are: ", D_for_L[L[l]])
        end
    end
    Plots.title!(p, "Correlation  D=$(Dvalue), T close to $t")
    Plots.xlabel!(p, "Distance (site)")
    Plots.ylabel!(p, "average Correlation")
end

function plotDCorrelation(L::Vector{Int64}, T_for_LD::Dict{Tuple{Int64, Float32}, Vector{Float32}}, D_for_L::Dict{Int64, Vector{Float32}}, Corr::Array{Vector{Float32}, 3}, l::Int64, t::Real)
    p=Plots.plot()
    for d in 1:3:length(D_for_L[L[l]])
        T=T_for_LD[L[l],D_for_L[L[l]][d]]
        closestT = argmin(abs.(T .- t))     # Compute the temperature closest to the requested one
        Plots.plot!(collect(1:Int(L[l]/2-1)), Corr[l, closestT, d], label="D=$(D_for_L[L[l]][d]), T=$(T[closestT])", legend=:topright)
    end
    Plots.title!(p, "Correlation  L=$(L[l]), T close to $t")
    Plots.xlabel!(p, "Distance (site)")
    Plots.ylabel!(p, "average Correlation")
end

function plotTCorrelation(L::Vector{Int64}, T_for_LD::Dict{Tuple{Int64, Float32}, Vector{Float32}}, D_for_L::Dict{Int64, Vector{Float32}}, Corr::Array{Vector{Float32}, 3}, Dvalue::Real, l::Int64)
    p=Plots.plot()
    if Dvalue ∈ D_for_L[L[l]]
        d = findfirst(==(Dvalue), D_for_L[L[l]])
        T = T_for_LD[(L[l], Dvalue)]
        nT = length(T)
        for t=1:max(round(Int, nT/8),1):nT
            Plots.plot!(collect(1:Int(L[l]/2-1)), Corr[l, t, d], label="$(T[t])",legend=:topright)
        end
        Plots.title!(p, "Correlation  L=$(L[l]), D=$(Dvalue)")
        Plots.xlabel!(p, "Distance (site)")
        Plots.ylabel!(p, "average Correlation")
    else
        println("For L=$(L[l]), valid D are: ", D_for_L[L[l]])
    end
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

# function interpmax(T::Vector{Float32}, y::Vector{Float32}) # interpolate and find the Tmax (for all Lattice size)
#     xinterp = collect(.5:.001:2)
#     yinterp=Spline1D(T, y, k=3)(xinterp)
#     return (findmax(yinterp)[1], xinterp[findmax(yinterp)[2]])
# end

# function interpmax_with_error(T::Vector{Float32}, y::Vector{Float32}, y_err::Vector{Float32}; n_bootstrap=1000) # interpolate and find the Tmax (for all Lattice size) then bootstrap resampling to get error bars
#     xinterp = collect(T[1]:.002:T[end])
#     max_values = Float32[]
#     max_positions = Float32[]
    
#     for i in 1:n_bootstrap
#         y_perturbed = y .+ randn(length(y)) .* y_err # Perturb y values according to their uncertainties
#         try # Fit spline and find max
#             yinterp = Spline1D(T, y_perturbed, k=3)(xinterp)
#             idx = argmax(yinterp)
#             push!(max_values, yinterp[idx])
#             push!(max_positions, xinterp[idx])
#         catch
#             continue # Skip failed fits (rare)
#         end
#     end    
#     return (median(max_values), median(max_positions), std(max_values), std(max_positions))
# end

function Gaussmax(T::Vector{Float32}, y::Vector{Float32})
    ymax = maximum(y)
    i=1
    while y[i] < .9*ymax
        i+=1
    end
    i_min = copy(i)
    i=length(T)
    while y[i] < .9*ymax
        i-=1
    end
    i_max = i
    if T[i_max] - T[i_min] < .02
        println("The pic is well defined, no need to fit a gaussian")
        return 0
    else
        fit = curve_fit(Gauss_fit, T[i_min:i_max], y[i_min:i_max], [0.8, 0.2, 0.5, 0.9])
        x = collect(T[i_min]:0.001f0:T[i_max])
        p=Plots.plot(T, y[1:length(T)])
        println(coef(fit))
        Plots.plot!(x, Gauss_fit(x, coef(fit)))
        display(p)
        return coef(fit)[1]
    end
end

# function crit(L::Vector{Int64}, T::Vector{Float32}, y::Array{Float32,2}, title::String="") # calculate α & γ, plot if title is a non-empty string
#     ymax = Vector{Tuple}(undef, length(L))
#     for l in eachindex(L)
#         ymax[l] = interpmax(T, y[l,:])
#     end
#     println(ymax)
#     fit =  curve_fit(linear_fit, log.(L), log.(first.(ymax)), if (y=="c"); [.1,-.7]; elseif (y=="susc"); [2.,-4.]; else [.9,-.9] end)
#     println(#=typeof(fit), "\n",=# fit)
#     m, p = coef(fit)
#     σm, σp = stderror(fit)
#     if title != ""
#         a = Plots.plot()
#         Plots.plot!(a, log.(L), log.(first.(ymax)), seriestype=:scatter)
#         Plots.plot!(a, log.(L), m*log.(L) .+p)
#         Plots.xlabel!(a, "Ln(Lattice length)")
#         Plots.ylabel!(a, "Ln(max("*title*"))")
#         display(a)
#     end
#     return m, σm
# end

function critlength(data::Array{Vector{Float32},3}, T_for_LD::Dict{Tuple{Int64, Float32}, Vector{Float32}}, L::Vector{Int64}, D_for_L::Dict{Int64, Vector{Float32}}, l::Int64, Dvalue::Real, t::Real, PLOT::Bool=true, ln::Bool=false) # Calculate critical exp or correlation length (by default)
    Dvalue = Float32(Dvalue)
    if Dvalue ∈ D_for_L[L[l]]
        d = findfirst(==(Dvalue), D_for_L[L[l]])
        closestT = argmin(abs.(T_for_LD[L[l],Dvalue] .- t))     # Compute the temperature closest to the requested one
        data=data[l,closestT,d] 
        negativ = findfirst(x -> x < 4e-4, data)    # don't take in account when the correlation is below than 4e-4 (0.07)
        if typeof(negativ) == Nothing; n = length(data)
        else; n=negativ-1; end
        x=collect(1:n)
        fit =  curve_fit(linear_fit, log.(x)*ln + !ln*x, log.(data[1:n]), [-.09,9.])
        m,p = coef(fit); σ = stderror(fit)

        if ln; lab="power law fit"; else; lab="exponential fit";end
        a=Plots.plot()
        Plots.plot!(a, log.(x)*ln + !ln*x, log.(data[1:n]), seriestype=:scatter, label="data at T=$(T_for_LD[L[l],Dvalue][closestT])")
        Plots.plot!(a, log.(x)*ln + !ln*x, m*(log.(x)*ln + !ln*x) .+p, label=lab)
        Plots.ylabel!(a, "Ln(Correlation)")
        if ln == true; Plots.xlabel!(a, "Ln(Distance)"); ξ=σξ=0
        else; Plots.xlabel!(a, "Distance"); ξ =-1/m; σξ = σ[1]/m^2
        end
        if PLOT==true;  title!(a, "Correlation at T=$(T_for_LD[L[l],Dvalue][closestT]), L=$(L[l]), D=$(Dvalue)");   display(a);   end
        return m*ln + !ln*ξ, σ[1]*ln+!ln*σξ, a
    else
        println(Dvalue, L[l])
        println("For L=$(L[l]), valid D are: ", D_for_L[L[l]])
    end
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

function CorrTimePlot(AllLattices::Array{Vector{Array{Float32, 3}},3}, T_for_LD::Dict{Tuple{Int64, Float32}, Vector{Float32}}, D::Vector{Float32}, d::Int64, L::Vector{Int64}, l::Int64, save::Bool=false)
    p=Plots.plot()
    pal = cgrad([:lightblue, :green, :red])
    for i in eachindex(T_for_LD[(L[l], D[d])])
        MCT = MeanCorrTime(AllLattices[l,i,d])
        x = T_for_LD[(L[l], D[d])][i]/T_for_LD[(L[l], D[d])][end]
        Plots.plot!(collect(1:length(MCT))*10, MCT, label="$(T_for_LD[(L[l], D[d])][i])", seriescolor = pal[x], linecolor = pal[x], markercolor = pal[x], markerstrokecolor = pal[x], ecolor = pal[x])
    end
    Plots.title!(p, "Correlation time, L=$(L[l]), D=$(D[d])")
    if save ==true
        Plots.savefig("Plot/CorrelationTime$(D[d]).pdf")
    end
    display(p)
end

function Plot_Max_C_Χ(pal::PlotUtils.ContinuousColorGradient, D_for_L::Dict{Int64, Vector{Float32}}, T_for_LD::Dict{Tuple{Int64, Float32}, Vector{Float32}}, Y::Array{Float32, 3}, Yerr::Array{Float32, 3}, L::Vector{Int64}, title::String="C/χ", save::Bool=false)
    p=Plots.plot() # xmax & T of xmax vs D
    xmax, Tmax, xmaxerr, Tmaxerr = Dict{Tuple{Int64, Float32}, Float32}(), Dict{Tuple{Int64, Float32}, Float32}(), Dict{Tuple{Int64, Float32}, Float32}(), Dict{Tuple{Int64, Float32}, Float32}()
    for l in eachindex(L)
        ll = L[l]
        for d in eachindex(D_for_L[L[l]])
            dd = D_for_L[ll][d]
            xmax[ll,dd], ArgMax = maximum(Y[l,:,d]), argmax(Y[l,:,d])
            Tmax[ll,dd], xmaxerr[ll,dd], Tmaxerr[ll,dd] = T_for_LD[ll, dd][ArgMax], Yerr[l,ArgMax,d], 0
        end
        Plots.plot!(D_for_L[L[l]], [Tmax[ll,i] for i in D_for_L[ll]], seriestype=:scatter, label="$(L[l])", yerr=[Tmaxerr[ll,i] for i in D_for_L[ll]])
    end
    Plots.title!("Tpic of $title vs D")
    Plots.ylabel!("Temperature of the maximum of $title")
    Plots.xlabel!("⬅ Ising                              Value of D                              XY ➡")
    if save; Plots.savefig("Plot/Tpic_$title.pdf"); end
    display(p)
    p=Plots.plot()
    for l in eachindex(L)
        Plots.plot!(D_for_L[L[l]], [xmax[L[l],i] for i in D_for_L[L[l]]], seriestype=:scatter, label="$(L[l])", yerr=[xmaxerr[L[l],i] for i in D_for_L[L[l]]])
    end
    Plots.title!("$title max vs D")
    Plots.ylabel!("Maximum of $title")
    Plots.xlabel!("⬅ Ising                              Value of D                              XY ➡")
    if save; Plots.savefig("Plot/$(title)_max.pdf"); end
    display(p)
    
    allD=Vector{Float32}()
    for ll in L
        allD = sort(unique(vcat(allD, D_for_L[ll])))
    end
    fit_Tc = []
    fit_ln = Dict(); fit_ln_err = Dict()
    fit_power = []
    x = L[1]:0.1:L[end]
    p2 = Plots.plot(load("Colorbar.png"), framestyle = :none, axis = nothing)
    n=Dict()
    for d in allD
        n[d] = PlotcolorD(d, allD[1], allD[end])
    end
    p1=Plots.plot() # xmax & T of xmax vs L
    for d in eachindex(allD)
        dd = allD[d]
        L_for_D = [k for (k, v) in pairs(D_for_L) if allD[d] ∈ v]
        for l in [findfirst(==(L_for_D[i]), L) for i in eachindex(L_for_D)]
            ll = L[l]
            xmax[ll,dd], ArgMax = maximum(Y[l,:,d]), argmax(Y[l,:,d])
            Tmax[ll,dd], xmaxerr[ll,dd], Tmaxerr[ll,dd] = T_for_LD[ll, dd][ArgMax], Yerr[l,ArgMax,d], 0
        end
        vect_Tmax = [Tmax[i,dd] for i in L_for_D]
        vect_Tmaxerr = [Tmaxerr[i,dd] for i in L_for_D]
        if title=="Susc"
            try
                fitTc = curve_fit(power_fit_plus, L, vect_Tmax, [-1.1, .9, .9])
                push!(fit_Tc, (round.(coef(fitTc);digits=3), round.(stderror(fitTc);digits=3)))
                Plots.plot!(x, coef(fitTc)[3]*x.^coef(fitTc)[1].+coef(fitTc)[2], seriescolor=pal[n[dd]], label="")
            catch
                println("can't fit D = $dd")
            end
        end
        Plots.plot!(L, vect_Tmax, seriestype=:scatter, seriescolor=pal[n[dd]], label="", yerr=vect_Tmaxerr)
    end
    Plots.title!("Tpic of $title vs L")
    Plots.ylabel!("Temperature of the maximum of $title")
    Plots.xlabel!("Lattice length")

    p=Plots.plot(p1, p2, layout= @layout [a{.88w} b{.12w}])
    if save; Plots.savefig("Plot/TpicvsL_$title.pdf"); end
    display(p)
    p1=Plots.plot()
    for d = eachindex(allD)
        dd = allD[d]
        L_for_D = [k for (k, v) in pairs(D_for_L) if allD[d] ∈ v]
        vect_xmax = [xmax[i,dd] for i in L_for_D]
        vect_xmaxerr = [xmaxerr[i,dd] for i in L_for_D]
        if title=="C"
            fitln =  curve_fit(fit_C2, L, vect_xmax, [1.1, 0.7, 0.1])
            fit_ln[dd] = round.(coef(fitln);digits=3)
            fit_ln_err[dd] = round.(stderror(fitln);digits=3)
        end
        fitpower = curve_fit(power_fit, L, vect_xmax, [.8, 1.1])
        push!(fit_power, (round.(coef(fitpower);digits=3), round.(stderror(fitpower);digits=3)))

        Plots.plot!(L, vect_xmax, seriestype=:scatter, label="", #=z=n,=# seriescolor=pal[n[dd]], yerr=vect_xmaxerr, xaxis=:log, xticks=(L, string.(L)))
        if title=="C"
            Plots.plot!(x, coef(fitln)[1]*log.(x) .+ coef(fitln)[2] .+ coef(fitln)[3]./x, label="", #=z=n,=# seriescolor=pal[n[dd]], xaxis=:log, xticks=(L, string.(L)))
        else
            Plots.plot!(x, coef(fitpower)[2]*x.^coef(fitpower)[1], label="", #=z=n,=# seriescolor=pal[n[dd]])
        end
    end
    Plots.plot!(p1, colorbar=true, colorbar_title="D Value")
    Plots.title!("$(title)_max vs L")
    Plots.ylabel!("Maximum of $title")
    Plots.xlabel!("Lattice length")

    p=Plots.plot(p1, p2, layout= @layout [a{.88w} b{.12w}])
    if save; Plots.savefig("Plot/$(title)_max_vsL.pdf"); end
    display(p)
    return Tmax, xmax, fit_ln, fit_ln_err, fit_power, fit_Tc
end

function Plot_Max_ξ(L::Vector{Int64}, D::Vector{Float32}, T_for_LD::Dict{Tuple{Int64, Float32}, Vector{Float32}}, Corr::Array{Vector{Float32}, 3}, l::Int64, ln::Bool=false, only_one_d::Int64=0)
    if only_one_d == 0
        Tmax = zeros(length(D))
        for d in eachindex(D)
            critic = zeros(length(T_for_LD[(L[l], D[d])]))
            for j in eachindex(T_for_LD[(L[l], D[d])])
                critic[j], a = critlength(Corr, T_for_LD, L, D_for_L, l, d, T_for_LD[(L[l], D[d])][j], false, ln)
            end
            # println(critic)
            Tmax[d] = round(T_for_LD[(L[l], D[d])][argmax(critic)];digits=3)
        end
        p=Plots.plot(D, Tmax)
        Plots.xlabel!("⬅ Ising                              Value of D                              XY ➡")
        Plots.ylabel!("Temperature of the maximum of ξ")
        display(p)
    else    # To be able to use this function for one D value at a time
        d = only_one_d
        critic = zeros(length(T_for_LD[(L[l], D[d])]))
        for j in eachindex(T_for_LD[(L[l], D[d])])
            critic[j], a = critlength(Corr, T_for_LD, L, D_for_L, l, d, T_for_LD[(L[l], D[d])][j], true, ln)
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

function Gauss_fit(x::Vector{Float32}, p::Vector{Float64})
    return p[3]*exp.(-((x .-p[1])/p[2]).^2) .+p[4]
end

function fit_C(x::Vector{Int64}, p::Vector{Float64})       # a * ln(L/L0)(1+L'/L)
    return p[1]*log.(x/p[2]).*(1 .+ p[3]./x)
end

function fit_C2(x::Vector{Int64}, p::Vector{Float64})       # a * ln(L/L0)(1+L'/L)
    return p[1]*log.(x) .+p[2] .+ p[3]./x
end

function colorbar(pal::PlotUtils.ContinuousColorGradient, D)
    fig = Figure(size = (62, 330))          # Make the wanted colorbar requires other packages. The colorbar is then saved and combined with the plot
    Colorbar(fig[1, 1]; colormap = pal, colorrange = (D[1], D[end]), ticks = ([D[1], D[1]/2, 0, D[end]/2, D[end]], string.(round.([D[1], D[1]*2^(-10/3), 0, D[end]*2^(-10/3), D[end]]; digits=2))), label = "← Ising                    D value                    XY →")
    Makie.save("Colorbar.png", fig; px_per_unit = 8)
end

function scaling_plot_C_Ising(L::Vector{Int64}, T_for_LD::Dict{Tuple{Int64, Float32}, Vector{Float32}}, D_for_L::Dict{Int64, Vector{Float32}}, C::Array{Float32, 3}, Dvalue::Real, Tpic::Dict{Tuple{Int64, Float32}, Float32}, error::Array{Float32, 3}=[], save::Bool=false)
    Dvalue = Float32(Dvalue)
    hex_codes = ["#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#999999", "#F781BF", "#A65628", "#FFFF33"]
    pal = [parse(Colorant, hex_codes[i]) for i in 1:length(L)]
    p=Plots.plot()
    for l in eachindex(L)
        ll = L[l]
        if haskey(T_for_LD, (ll, Dvalue))
            T = T_for_LD[ll, Dvalue]
            endT = length(T)
            d = findfirst(==(Dvalue), D_for_L[ll])
            Plots.plot!((T.-Tpic[ll, Dvalue])*ll, C[l, 1:endT,d]/log(ll), yerr = error[l,1:endT,d]/log(ll), label=ll, seriescolor = pal[l], linecolor = pal[l], markercolor = pal[l], markerstrokecolor = pal[l], ecolor = pal[l])
        end
    end
    Plots.title!(p, "C, with D = $Dvalue")
    Plots.xlabel!(p, "(T-Tc(L))*L")
    Plots.ylabel!(p, "C/ln(L)")
    if save ==true
        Plots.savefig("Plot/Scaling_C_Ising.pdf")
    end
    display(p)
end

function scaling_plot_χ_Ising(L::Vector{Int64}, T_for_LD::Dict{Tuple{Int64, Float32}, Vector{Float32}}, D_for_L::Dict{Int64, Vector{Float32}}, χ::Array{Float32, 3}, Dvalue::Real, Tpic::Dict{Tuple{Int64, Float32}, Float32}, error::Array{Float32, 3}=[], save::Bool=false)
    Dvalue = Float32(Dvalue)
    hex_codes = ["#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#999999", "#F781BF", "#A65628", "#FFFF33"]
    pal = [parse(Colorant, hex_codes[i]) for i in 1:length(L)]
    p=Plots.plot()
    for l in eachindex(L)
        ll = L[l]
        if haskey(T_for_LD, (ll, Dvalue))
            T = T_for_LD[ll, Dvalue]
            endT = length(T)
            d = findfirst(==(Dvalue), D_for_L[ll])
            Plots.plot!((T.-Tpic[ll, Dvalue])*ll, χ[l, 1:endT,d]*ll.^(-7/4), yerr = error[l,1:endT,d]*ll^(-7/4), label=ll, seriescolor = pal[l], linecolor = pal[l], markercolor = pal[l], markerstrokecolor = pal[l], ecolor = pal[l])
        end
    end
    Plots.title!(p, "χ, with D = $Dvalue")
    Plots.xlabel!(p, "(T-Tc(L))*L")
    Plots.ylabel!(p, "χ*L^(-7/4))")
    if save ==true
        Plots.savefig("Plot/Scaling_C_Ising.pdf")
    end
    display(p)
end

# α : specific heat                     Ising : 0       XY : NOP Essential singularity
# β : zero field mag                    Ising : 1/8     XY : NOP No magnetization (To have a nonzero 𝑀, the correlation function must approach a constant at large distance. But in the 2D XY model: For 𝑇>𝑇BKT: correlations decay exponentially. For 𝑇<𝑇BKT: correlations decay as a power law)
# γ : zero field isothermal suscxx      Ising : 7/4     XY : NOP Essential divergence
# δ : Critical isothermal               Ising : 15      XY : 15
# ν : corr length                       Ising : 1       XY : NOP 𝜉 diverges exponentially

# continuous symmetries cannot be spontaneously broken in 2D