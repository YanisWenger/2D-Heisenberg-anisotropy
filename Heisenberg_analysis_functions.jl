function Get_T_LD(L::Vector{Int64}, FolderName::String="")  # for each lattice size gives which T have each d
    T_for_LD = Dict{Tuple{Int64, Float32}, Vector{Float32}}()
    D_for_L = Dict{Int64, Vector{Float32}}()
    all_T = []
    for l in L
        files = glob("$(l)_*.jld2", ("Data/$FolderName"))
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
    p=Plots.plot(framestyle = :box)
    Dvalue = Float32(Dvalue)
    pal = cgrad([RGB(.3,1,.3), RGB(.3,.3, 1), RGB(1,.3,.3)], length(L), categorical = true)
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
    # Plots.title!(p, ytitle*" as a function of T with D = $(Dvalue)")
    Plots.xlabel!(p, "Temperature")
    Plots.ylabel!(p, ytitle)
    if save ==true
        Plots.savefig("Plot/"*ytitle*"_d=$(Dvalue).pdf")
    end
    display(p)
end

function plotD(pal::PlotUtils.ContinuousColorGradient, L::Vector{Int64}, T_for_LD::Dict{Tuple{Int64, Float32}, Vector{Float32}}, D_for_L::Dict{Int64, Vector{Float32}}, x::Array{Float32, 3}, ytitle::String="", error=[], l::Int64=length(x[:,1,1]), save::Bool=false)   # read data from a file and plot it vs D
    p2 = Plots.plot(load("Colorbar.png"), framestyle = :none, axis = nothing)
    p1 = Plots.plot(framestyle = :box)
    ll=L[l]; 
    D = D_for_L[ll];
    # if ytitle=="C"; lim=(max(0, minimum(x[l,:,:])),maximum(x[l,3:end,:]))
    #=else;=# lim = :auto
    # end
    for d in eachindex(D)
        dd = D[d]
        rescaleE=0          # To rescale the Energy (only when y = Energy)
        if ytitle=="E" && dd<0; rescaleE = dd; end
        n = PlotcolorD(dd, D[1], D[end])
        T = T_for_LD[(ll, dd)]
        endT=length(T)
        if error != [] 
            Plots.plot!(T, x[l, 1:endT, d].-rescaleE, yerr = Vector(error[l, 1:endT, d]), label=""#=$(dd)=#, ylims=lim, seriescolor = pal[n], linecolor = pal[n], markercolor = pal[n], markerstrokecolor = pal[n], ecolor = pal[n])#, seriestype=:scatter, ms=3)
        else
            Plots.plot!(T, x[l, 1:endT, d].-rescaleE, label=""#=$(dd)=#, seriescolor = pal[n])#, seriestype=:scatter, ms=3)
        end
    end
    # Plots.title!(p1, ytitle*" as a function of T for $(ll)x$(ll) lattices")
    Plots.xlabel!(p1, "Temperature")
    Plots.ylabel!(p1, ytitle)
    p=Plots.plot(p1, p2, layout= @layout [a{.84w} b{.16w}])
    if save; Plots.savefig("Plot/"*ytitle*".pdf"); end
    display(p)
end

function plotLCorrelation(L::Vector{Int64}, T_for_LD::Dict{Tuple{Int64, Float32}, Vector{Float32}}, D_for_L::Dict{Int64, Vector{Float32}}, Corr::Array{Vector{Float32}, 3}, Dvalue::Real, t::Real)
    p=Plots.plot(framestyle = :box)
    Dvalue = Float32(Dvalue)
    for l in eachindex(L)
        ll=L[l]
        if (ll, Dvalue) ∈ keys(T_for_LD)
            T=T_for_LD[ll, Dvalue]
            if Dvalue ∈ D_for_L[ll]
                d = findfirst(==(Dvalue), D_for_L[ll])
                closestT = argmin(abs.(T .- t))     # Compute the temperature closest to the requested one
                Plots.plot!(collect(1:Int(ll/2-1)), Corr[l, closestT, d], label="L=$(ll), T=$(T[closestT])", legend=:topright)
            else
                println("For L=$(ll), valid D are: ", D_for_L[ll])
            end
        end
    end
    Plots.title!(p, "Correlation  D=$(Dvalue), T close to $t")
    Plots.xlabel!(p, "Distance (site)")
    Plots.ylabel!(p, "average Correlation")
end

function plotDCorrelation(pal::PlotUtils.ContinuousColorGradient, L::Vector{Int64}, T_for_LD::Dict{Tuple{Int64, Float32}, Vector{Float32}}, D_for_L::Dict{Int64, Vector{Float32}}, Corr::Array{Vector{Float32}, 3}, ll::Int64, t::Real)
    p2 = Plots.plot(load("Colorbar.png"), framestyle = :none, axis = nothing)
    p1 = Plots.plot(legend=:topright, framestyle = :box)
    l = findfirst(x -> x==ll, L)
    D = D_for_L[ll]
    for d in 1:3:length(D_for_L[ll])
        dd = D[d]
        T=T_for_LD[ll,dd]
        closestT = argmin(abs.(T .- t))     # Compute the temperature closest to the requested one
        Plots.plot!(collect(1:Int(ll/2-1)), Corr[l, closestT, d], label="D=$(dd), T=$(T[closestT])", seriescolor = pal[PlotcolorD(D[d], D[1], D[end])])
    end
    Plots.title!(p1, "Correlation  L=$ll, T close to $t")
    Plots.xlabel!(p1, "Distance (site)")
    Plots.ylabel!(p1, "average Correlation")
    display(Plots.plot(p1, p2, layout= @layout [a{.84w} b{.16w}]))
end

function plotTCorrelation(L::Vector{Int64}, T_for_LD::Dict{Tuple{Int64, Float32}, Vector{Float32}}, D_for_L::Dict{Int64, Vector{Float32}}, Corr::Array{Vector{Float32}, 3}, Dvalue::Real, ll::Int64)
    Dvalue = Float32(Dvalue)
    l = findfirst(x -> x==ll, L)
    p=Plots.plot(framestyle = :box)
    if Dvalue ∈ D_for_L[ll]
        d = findfirst(==(Dvalue), D_for_L[ll])
        T = T_for_LD[(ll, Dvalue)]
        nT = length(T)
        for t=1:max(round(Int, nT/8),1):nT
            Plots.plot!(collect(1:Int(ll/2-1)), Corr[l, t, d], label="$(T[t])",legend=:topright)
        end
        Plots.title!(p, "Correlation  L=$(L[l]), D=$(Dvalue)")
        Plots.xlabel!(p, "Distance (site)")
        Plots.ylabel!(p, "average Correlation")
    else
        println("For L=$(L[l]), valid D are: ", D_for_L[ll])
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

# function Gaussmax(T::Vector{Float32}, y::Vector{Float32})
#     ymax = maximum(y)
#     i=1
#     while y[i] < .9*ymax
#         i+=1
#     end
#     i_min = copy(i)
#     i=length(T)
#     while y[i] < .9*ymax
#         i-=1
#     end
#     i_max = i
#     if T[i_max] - T[i_min] < .02
#         println("The pic is well defined, no need to fit a gaussian")
#         return 0
#     else
#         fit = curve_fit(Gauss_fit, T[i_min:i_max], y[i_min:i_max], [0.8, 0.2, 0.5, 0.9])
#         x = collect(T[i_min]:0.001f0:T[i_max])
#         p=Plots.plot(T, y[1:length(T)], framestyle = :box)
#         println(coef(fit))
#         Plots.plot!(x, Gauss_fit(x, coef(fit)))
#         display(p)
#         return coef(fit)[1]
#     end
# end

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

function critlength(data::Array{Vector{Float32},3}, T_for_LD::Dict{Tuple{Int64, Float32}, Vector{Float32}}, L::Vector{Int64}, D_for_L::Dict{Int64, Vector{Float32}}, ll::Int64, Dvalue::Real, t::Real, pow::Bool, PLOT::Bool, y="Connected correlation function") # Calculate critical exp or correlation length (by default)
    Dvalue = Float32(Dvalue)
    if Dvalue ∈ D_for_L[ll]
        d = findfirst(==(Dvalue), D_for_L[ll])
        closestT = argmin(abs.(T_for_LD[ll,Dvalue] .- t))     # Compute the temperature closest to the requested one
        data=data[findfirst(x -> x==ll, L),closestT,d] 
        stop = findfirst(x -> x < .002, data)    # doesn't take in account when the correlation is below than 0.002 to avoid going in the too noisy part
        if typeof(stop) == Nothing; n = length(data)
        else; n=stop-1; end
        if y == "Correlation function";  n=min(n, round(Int,ll/4));    end
        xx = [2^j for j in 0:Int(trunc(log2(n)))]
        if y == "Correlation function"; yy=collect(trunc(data[n]+.1; digits=1):.1:trunc(data[1]+.1;digits=1))
        else;   yy = [2.0^j/1000 for j in Int(trunc(log2(data[n]*1000))+1):2:Int(trunc(log2(data[1]*1000)))]
        end
        if n == 1; println("no: D = ", Dvalue, "  L = ", ll); return 0,0
        else
            x=collect(1:n)
            fit = curve_fit(linear_fit, log.(x)*pow + !pow*x, log.(data[1:n]), [-.09,9.])
            m, p = coef(fit)
            σ = stderror(fit)
            if pow;  ξ=σξ=0
            else;    ξ=-1/m;    σξ = σ[1]/m^2
            end
            if PLOT
                a=Plots.plot(yticks=(log.(yy), string.(yy)), xaxis="Distance", yaxis=y)# yticks=collect(max(data[1], exp(p)) :1: min(data[n], n^m*exp(p))), string.(1))
                if pow;  Plots.plot!(xticks=(log.(xx), string.(xx))); end
                if pow;  lab="\$\\propto d^{"*nice_result(m, σ[1])*"}\$"; else; lab="\$\\propto e^{-d/"*nice_result(-m, σ[1])*"}\$";    end
                Plots.plot!(a, log.(x)*pow + !pow*x, log.(data[1:n]), seriestype=:scatter, markersize=3, label="T=$(T_for_LD[ll,Dvalue][closestT]), D=$Dvalue")
                Plots.plot!(a, log.(x)*pow + !pow*x, m*(log.(x)*pow + !pow*x) .+p, label=lab)
                # title!(a, "Correlation at T=$(T_for_LD[ll,Dvalue][closestT]), L=$(ll), D=$(Dvalue)")
                display(a)
            end
            return m*pow + !pow*ξ, σ[1]*pow+!pow*σξ
        end
    else
        println(Dvalue, ll)
        println("For L=$(ll), valid D are: ", D_for_L[ll])
    end
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
    p=Plots.plot(framestyle = :box)
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
    pal_L = cgrad([RGB(.3,1,.3), RGB(.3,.3, 1), RGB(1,.3,.3)], length(L), categorical = true)
    p=Plots.plot(framestyle = :box) # xmax & T of xmax vs D
    xmax, Tmax, xmaxerr, Tmaxerr = Dict{Tuple{Int64, Float32}, Float32}(), Dict{Tuple{Int64, Float32}, Float32}(), Dict{Tuple{Int64, Float32}, Float32}(), Dict{Tuple{Int64, Float32}, Float32}()
    for l in eachindex(L)
        ll = L[l]
        for d in eachindex(D_for_L[ll])
            dd = D_for_L[ll][d]
            xmax[ll,dd], ArgMax = maximum(Y[l,:,d]), argmax(Y[l,:,d])
            Tmax[ll,dd], xmaxerr[ll,dd], Tmaxerr[ll,dd] = T_for_LD[ll, dd][ArgMax], Yerr[l,ArgMax,d], 0
        end
        Plots.plot!(D_for_L[ll], [Tmax[ll,i] for i in D_for_L[ll]], seriestype=:scatter, label=ll, yerr=[Tmaxerr[ll,i] for i in D_for_L[ll]], seriescolor=pal_L[l])
    end
    Plots.title!("Tpic of $title vs D")
    Plots.ylabel!("Temperature of the maximum of $title")
    Plots.xlabel!("⬅ Ising                Value of D                XY ➡")
    if save; Plots.savefig("Plot/Tpic_$title.pdf"); end
    display(p)
    p=Plots.plot(framestyle = :box)
    for l in eachindex(L)
        Plots.plot!(D_for_L[L[l]], [xmax[L[l],i] for i in D_for_L[L[l]]], seriestype=:scatter, label="$(L[l])", yerr=[xmaxerr[L[l],i] for i in D_for_L[L[l]]], seriescolor=pal_L[l])
    end
    Plots.title!("$title max vs D")
    Plots.ylabel!("Maximum of $title")
    Plots.xlabel!("⬅ Ising                Value of D                XY ➡")
    if save; Plots.savefig("Plot/$(title)_max.pdf"); end
    display(p)
    
    allD=Vector{Float32}()
    for ll in L
        allD = sort(unique(vcat(allD, D_for_L[ll])))
    end
    fit_Tc = Dict{Float32, Tuple{Vector{Float64}, Vector{Float64}}}()
    fit_power_or_ln = Dict{Float32, Tuple{Vector{Float64}, Vector{Float64}}}()
    p2 = Plots.plot(load("Colorbar.png"), framestyle = :none, axis = nothing)
    n=Dict()
    for d in allD
        n[d] = PlotcolorD(d, allD[1], allD[end])
    end
    p1=Plots.plot(framestyle = :box) # xmax & T of xmax vs L
    for d in eachindex(allD)
        dd = allD[d]
        L_for_D = sort([k for (k, v) in pairs(D_for_L) if allD[d] ∈ v])
        x = L_for_D[1]:0.1:L_for_D[end]
        vect_Tmax = [Tmax[i,dd] for i in L_for_D]
        vect_Tmaxerr = [Tmaxerr[i,dd] for i in L_for_D]
        try
            fitTc = curve_fit(power_fit_plus, L_for_D, vect_Tmax, [-1.1, .9, .9])
            fit_Tc[dd] = (coef(fitTc), stderror(fitTc))
            Plots.plot!(x, coef(fitTc)[3]*x.^coef(fitTc)[1].+coef(fitTc)[2], seriescolor=pal[n[dd]], label="")
        catch
            println("can't fit D = $dd")
        end
        Plots.plot!(L_for_D, vect_Tmax, seriestype=:scatter, seriescolor=pal[n[dd]], label="", yerr=vect_Tmaxerr)
    end
    Plots.title!("Tpic of $title vs L")
    Plots.ylabel!("Temperature of the maximum of $title")
    Plots.xlabel!("Lattice length")
    p=Plots.plot(p1, p2, layout= @layout [a{.84w} b{.16w}])
    if save; Plots.savefig("Plot/TpicvsL_$title.pdf"); end
    display(p)

    p1=Plots.plot(framestyle = :box)
    for d = eachindex(allD)
        dd = allD[d]
        L_for_D = sort([k for (k, v) in pairs(D_for_L) if allD[d] ∈ v])
        x = L_for_D[1]:0.1:L_for_D[end]
        vect_xmax = [xmax[i,dd] for i in L_for_D]
        vect_xmaxerr = [xmaxerr[i,dd] for i in L_for_D]
        if title=="Heat capacity"
            fitpowerorln =  curve_fit(fit_C2, L_for_D, vect_xmax, [1.1, 0.7, 0.1])
        else
            fitpowerorln = curve_fit(power_fit, L_for_D, vect_xmax, [.8, 1.1])
        end
        fit_power_or_ln[dd] = (coef(fitpowerorln), stderror(fitpowerorln))

        Plots.plot!(L_for_D, vect_xmax, seriestype=:scatter, label="", #=z=n,=# seriescolor=pal[n[dd]], yerr=vect_xmaxerr, xaxis=:log, xticks=(L, string.(L)))
        if title=="Heat capacity"
            Plots.plot!(x, coef(fitpowerorln)[1]*log.(x) .+ coef(fitpowerorln)[2] .+ coef(fitpowerorln)[3]./x, label="", #=z=n,=# seriescolor=pal[n[dd]], xaxis=:log, xticks=(L, string.(L)))
        else
            Plots.plot!(x, coef(fitpowerorln)[2]*x.^coef(fitpowerorln)[1], label="", #=z=n,=# seriescolor=pal[n[dd]])
        end
    end
    Plots.plot!(p1, colorbar=true, colorbar_title="D Value")
    Plots.title!("$title max vs L")
    Plots.ylabel!("Maximum of $title")
    Plots.xlabel!("Lattice length")

    p=Plots.plot(p1, p2, layout= @layout [a{.84w} b{.16w}])
    if save; Plots.savefig("Plot/$(title)_max_vsL.pdf"); end
    display(p)
    return Tmax, xmax, fit_power_or_ln, fit_Tc
end

function Fit_ξ(L::Vector{Int64}, D_for_L::Dict{Int64, Vector{Float32}}, T_for_LD::Dict{Tuple{Int64, Float32}, Vector{Float32}}, Corr::Array{Vector{Float32}, 3}, only_one_l::Int64=123456789, only_one_d::Float32=123456789f0)
    pal = cgrad([RGB(.3,1,.3), RGB(.3,.3, 1), RGB(1,.3,.3)], length(L), categorical = true)
    critic_exp = Dict{Tuple{Int64, Float32, Float32}, Tuple{Float64, Float64}}()
    critic_pow = Dict{Tuple{Int64, Float32, Float32}, Tuple{Float64, Float64}}()
    if only_one_d == 123456789f0
        Tmax = Dict{Tuple{Int64, Float32}, Tuple{Float64, Tuple{Float64, Float64}}}()
        N = maximum([length(d) for d in values(D_for_L)])
        Tmax_array = zeros(length(L), N); Tmax_err_lower = zeros(length(L), N); Tmax_err_upper = zeros(length(L), N)
        for l in eachindex(L)
            ll=L[l]
            for d in eachindex(D_for_L[ll])
                dd = D_for_L[ll][d]
                Exp_vect = []; Exp_err_vect = []
                Pow_vect = []; Pow_err_vect = []
                for tt in T_for_LD[(ll, dd)]
                    Exp, Exp_err = critlength(Corr, T_for_LD, L, D_for_L, ll, dd, tt, false, false)
                    critic_exp[ll,tt,dd] = (Exp, Exp_err)
                    push!(Exp_vect, Exp)
                    push!(Exp_err_vect, Exp_err)
                    Pow, Pow_err = critlength(Corr, T_for_LD, L, D_for_L, ll, dd, tt, true, false)
                    critic_pow[ll,tt,dd] = (Pow, Pow_err)
                    push!(Pow_vect, Pow)
                    push!(Pow_err_vect, Pow_err)
                end
                argMAX = argmax(Exp_vect)
                MAX = round(T_for_LD[(ll, dd)][argMAX];digits=3)
                Terr = (MAX-T_for_LD[ll,dd][findfirst(x -> x>Exp_vect[argMAX]-Exp_err_vect[argMAX], Exp_vect)], T_for_LD[ll,dd][findlast(x -> x>Exp_vect[argMAX]-Exp_err_vect[argMAX], Exp_vect)]-MAX)
                Tmax[ll,dd] = (MAX, Terr)
                Tmax_array[l,d] = MAX; Tmax_err_lower[l,d] = Terr[1]; Tmax_err_upper[l,d] = Terr[2]
            end
        end
        p=Plots.scatter(xaxis="⬅ Ising                Value of D                XY ➡", yaxis = "Temperature of the maximum of ξ")
        for l in eachindex(L); Plots.scatter!(D_for_L[L[l]], Tmax_array[l,1:length(D_for_L[L[l]])], yerr=(Tmax_err_lower[l,1:length(D_for_L[L[l]])], Tmax_err_upper[l,1:length(D_for_L[L[l]])]), seriescolor = pal[l], label="\$\\xi\$"); end
        display(p)
        return critic_exp, critic_pow, Tmax
    else    # To be able to use this function for one D value at a time
        dd = only_one_d
        ll = only_one_l
        Exp_vect = []; Exp_err_vect = []
        Pow_vect = []; Pow_err_vect = []
        for tt in T_for_LD[(ll, dd)]
            Exp, Exp_err = critlength(Corr, T_for_LD, L, D_for_L, ll, dd, tt, false, false)
            critic_exp[ll,tt,dd] = (Exp, Exp_err)
            push!(Exp_vect, Exp)
            push!(Exp_err_vect, Exp_err)
            Pow, Pow_err = critlength(Corr, T_for_LD, L, D_for_L, ll, dd, tt, true, false)
            critic_pow[ll,tt,dd] = (Pow, Pow_err)
            push!(Pow_vect, Pow)
            push!(Pow_err_vect, Pow_err)
        end
        display(Plots.plot(T_for_LD[ll,dd], Exp_vect, yerr=Exp_err_vect, xaxis="Temperature", yaxis="Correlation length", label="$dd $ll", linecolor =:purple, markercolor =:purple))#, ylim=(0,80)))
        display(Plots.plot(T_for_LD[ll,dd], Pow_vect, yerr=Pow_err_vect, xaxis="Temperature", yaxis="Pow", label="$dd $ll", linecolor =:purple, markercolor =:purple))
        return critic_exp, critic_pow
    end
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

function exp_plus(x::Vector{Int64}, p::Vector{Float64})       # a * ln(L/L0)(1+L'/L)
    return exp.(-x/p[1]) .+p[2]
end

function colorbar(pal::PlotUtils.ContinuousColorGradient, D_for_L::Dict{Int64, Vector{Float32}})    #   Make the wanted colorbar requires other packages. The colorbar is then saved and combined with the plot
    set_theme!(fonts = (; regular="CMU Serif"), fontsize = 20)
    fig = Figure(size = (85, 350))
    D = sort(unique(vcat(values(D_for_L)...)))
    Colorbar(fig[1, 1]; colormap = pal, colorrange = (D[1], D[end]), ticks = ([D[1], D[1]/2, 0, D[end]/2, D[end]], string.(round.([D[1], D[1]*2^(-10/3), 0, D[end]*2^(-10/3), D[end]]; digits=2))), label = "← Ising           D value           XY →")
    Makie.save("Colorbar.png", fig; px_per_unit = 8)
end

function scaling_plot_C_Ising(L::Vector{Int64}, T_for_LD::Dict{Tuple{Int64, Float32}, Vector{Float32}}, D_for_L::Dict{Int64, Vector{Float32}}, C::Array{Float32, 3}, Dvalue::Real, Tpic::Dict{Tuple{Int64, Float32}, Float32}, error::Array{Float32, 3}, FIT::Dict{Float32, Tuple{Vector{Float64}, Vector{Float64}}})
    Dvalue = Float32(Dvalue)
    pal = cgrad([RGB(.3,1,.3), RGB(.3,.3, 1), RGB(1,.3,.3)], length(L), categorical = true)
    p=Plots.plot()
    Lprime = FIT[Float32(Dvalue)][1][3]
    L0 = exp(-FIT[Float32(Dvalue)][1][2] / FIT[Float32(Dvalue)][1][1])
    for l in eachindex(L)
        ll = L[l]
        if haskey(T_for_LD, (ll, Dvalue))
            T = T_for_LD[ll, Dvalue]
            endT = length(T)
            d = findfirst(==(Dvalue), D_for_L[ll])
            Plots.plot!((T.-Tpic[ll, Dvalue])*ll, (C[l, 1:endT,d] .-Lprime/ll)/log(ll/L0), yerr = error[l,1:endT,d]/log(ll), label=ll, seriescolor = pal[l], linecolor = pal[l], markercolor = pal[l], markerstrokecolor = pal[l], ecolor = pal[l])
        end
    end
    # Plots.title!(p, "C, with D = $Dvalue")
    Plots.xlabel!(p, "(T-Tc(L))\$\\,\\cdot\$L")
    Plots.ylabel!(p, "\$ (C-\\frac{L'}{L})/ln(\\frac{L}{L_0})\$")
    display(p)
end

function scaling_plot_χ_Ising(L::Vector{Int64}, T_for_LD::Dict{Tuple{Int64, Float32}, Vector{Float32}}, D_for_L::Dict{Int64, Vector{Float32}}, χ::Array{Float32, 3}, Dvalue::Real, Tpic::Dict{Tuple{Int64, Float32}, Float32}, error::Array{Float32, 3}, fitpower::Dict{Float32, Tuple{Vector{Float64}, Vector{Float64}}})
    Dvalue = Float32(Dvalue)
    pal = cgrad([RGB(.3,1,.3), RGB(.3,.3, 1), RGB(1,.3,.3)], length(L), categorical = true)
    p=Plots.plot(framestyle = :box)
    pow  = fitpower[Float32(Dvalue)][1][1]
    Δpow = fitpower[Float32(Dvalue)][2][1]
    for l in eachindex(L)
        ll = L[l]
        if haskey(T_for_LD, (ll, Dvalue))
            T = T_for_LD[ll, Dvalue]
            endT = length(T)
            d = findfirst(==(Dvalue), D_for_L[ll])
            Plots.plot!((T.-Tpic[ll, Dvalue])*ll, χ[l, 1:endT,d]*ll.^(-pow), yerr = error[l,1:endT,d]*ll^(-7/4), label=ll, seriescolor = pal[l], linecolor = pal[l], markercolor = pal[l], markerstrokecolor = pal[l], ecolor = pal[l])
        end
    end
    # Plots.title!(p, "\$\\chi\$, with D = $Dvalue")
    Plots.xlabel!(p, raw"$(T-Tc(L))\cdotL$")
    Plots.ylabel!(p, raw"$\chi \cdot L^{-"*nice_result(pow, Δpow)*"}\$")
    display(p)
end

function scaling_plot_C_XY(L::Vector{Int64}, T_for_LD::Dict{Tuple{Int64, Float32}, Vector{Float32}}, D_for_L::Dict{Int64, Vector{Float32}}, C::Array{Float32, 3}, Dvalue::Real, Tpic::Dict{Tuple{Int64, Float32}, Float32}, error::Array{Float32, 3})
    Dvalue = Float32(Dvalue)
    pal = cgrad([RGB(.3,1,.3), RGB(.3,.3, 1), RGB(1,.3,.3)], length(L), categorical = true)
    p=Plots.plot(framestyle = :box)
    for l in eachindex(L)
        ll = L[l]
        if haskey(T_for_LD, (ll, Dvalue))
            T = T_for_LD[ll, Dvalue]
            endT = length(T)
            d = findfirst(==(Dvalue), D_for_L[ll])
            # Plots.plot!((T.-Tpic[ll, Dvalue])*ll, χ[l, 1:endT,d]*ll.^(-7/4), yerr = error[l,1:endT,d]*ll^(-7/4), label=ll, seriescolor = pal[l], linecolor = pal[l], markercolor = pal[l], markerstrokecolor = pal[l], ecolor = pal[l])
            Plots.plot!(T.-Tpic[ll, Dvalue], C[l, 1:endT,d], yerr = error[l,1:endT,d]*ll^(-7/4), label=ll, seriescolor = pal[l], linecolor = pal[l], markercolor = pal[l], markerstrokecolor = pal[l], ecolor = pal[l])
        end
    end
    # Plots.title!(p, "C, with D = $Dvalue")
    Plots.xlabel!(p, "\$T-Tc(L)\$")
    Plots.ylabel!(p, "Heat capacity")
    display(p)
end

function scaling_plot_χ_XY(L::Vector{Int64}, T_for_LD::Dict{Tuple{Int64, Float32}, Vector{Float32}}, D_for_L::Dict{Int64, Vector{Float32}}, χ::Array{Float32, 3}, Dvalue::Real, Tpic::Dict{Tuple{Int64, Float32}, Float32}, error::Array{Float32, 3}, fitpower::Dict{Float32, Tuple{Vector{Float64}, Vector{Float64}}})
    Dvalue = Float32(Dvalue)
    pal = cgrad([RGB(.3,1,.3), RGB(.3,.3, 1), RGB(1,.3,.3)], length(L), categorical = true)
    p=Plots.plot(framestyle = :box)
    pow  = fitpower[Float32(Dvalue)][1][1]
    Δpow = fitpower[Float32(Dvalue)][2][1]
    for l in eachindex(L)
        ll = L[l]
        if haskey(T_for_LD, (ll, Dvalue))
            T = T_for_LD[ll, Dvalue]
            endT = length(T)
            d = findfirst(==(Dvalue), D_for_L[ll])
            x = T.-Tpic[ll, Dvalue]
            Plots.plot!(sign.(x).*sqrt.(x.*sign.(x)).^3*ll, χ[l, 1:endT,d]*ll.^(-pow), yerr = error[l,1:endT,d]*ll^(-7/4), label=ll, seriescolor = pal[l], linecolor = pal[l], markercolor = pal[l], markerstrokecolor = pal[l], ecolor = pal[l], xlim = (-2.5,4.5))
        end
    end
    Plots.xlabel!(p, "\$L \\cdot (T-Tc(L))^{3/2}\$")
    Plots.ylabel!(p, "\$\\chi\\cdot L ^{-"*nice_result(pow, Δpow)*"}\$")
    display(p)
end

function LogLogPlotAndFit(x::Vector, y::Vector, guess::Vector, xax::String, yax::String, lab::String)
    P=Plots.plot(framestyle = :box)
    xlog = log.(x)
    ylog = log.(y)
    # yy = round.(y;digits=1)
    # yy = Int.([2^j for j in trunc(log2(y[1])):1:trunc(log2(y[end]))])
    yy = [1, 3, 10, 30, 100, 300]
    fit  = curve_fit(linear_fit, xlog, ylog, guess)
    m, p = coef(fit)
    Δm, Δp = stderror(fit)
    Plots.plot!(xlog, ylog, label=" "*lab, seriestype=:scatter, xticks=(xlog, string.(x)), yticks=(log.(yy), string.(yy)))
    Plots.plot!(xlog, xlog*m .+p, label = "\$\\propto L^{"*nice_result(m, Δm)*"}\$")
    Plots.xlabel!(xax)
    Plots.ylabel!(yax)
    display(P)
    return m, p, Δm, Δp
end

function xlogPlotAndFit(x::Vector, y::Vector, guess::Vector)
    P=Plots.plot()
    xlog = log.(x)
    fit  = curve_fit(linear_fit, xlog, y, guess)
    m, p, = coef(fit)
    Δm = stderror(fit)[1]
    Plots.plot!(xlog, y, label=" \$C_{max}\$", seriestype=:scatter, xticks=(xlog, string.(x)))
    Plots.plot!(xlog, xlog*m .+p, label = "\$\\propto\\ln(L)\$")
    Plots.xlabel!("Lattice length L")
    Plots.ylabel!("Maximum of heat capacity")
    display(P)
    return m, p, Δm, stderror(fit)[2]
end

function first_nonzero_digit(x::Float64)  # give the first non-zero digit of a number and its position
    if x >= 10 || x <= 0
        error("first_nonzero_decimal requires an input must be between 0 and 1, given : $x")
    end
    position = -floor(Int, log10(x))
    shifted = x * (10 ^ position)
    digit = floor(Int, shifted)
    if digit == 0           # Rare case due to precision, try next position
        return first_nonzero_decimal(x * 10)
    elseif digit == 10      # Precision issue where it rounded up to 10
        digit = 1
        position += 1
    end
    return digit, position
end

function nice_result(x::Float64, err::Float64)  #   Input: result and error margin      Output: result with proper number of decimal according to the error margin
    pow = floor(Int, log10(abs(x)))
    if pow > 0
        error("x=$x ... Cannont handle number smaller than -10 nor bigger than 10")
    end
    er, digit = first_nonzero_digit(err)
    if er == 9; digit -= 1; er = 1
    else; er +=1
    end
    if -pow > digit
        return "0."*"0"^digit*"($er)"
    end
    x = round(x, digits=digit)
    s = string(x)
    if '.' ∉ s;     s *= ".";   end
    parts = split(s, '.')
    current_digits = length(parts[2])
    if current_digits < digit
        s *= "0"^ (digit - current_digits)
    end
    return s*"($er)"
end

function Convergence_check(L::Vector{Int64}, FolderName::String="", ll::Int64=0)
    T, D_for_L, T_for_LD = Get_T_LD(L, FolderName)   # Obtain Temperature and anisotropic term in Data/$(L)_$N i.e. Data/[8, 12, 20]_1000000
    conv=[]
    if ll == 0  # if the lattice size is not given, check for all lattice size
        for ll in L
            for d in eachindex(D_for_L[ll])
                dd = D_for_L[ll][d]
                for t in eachindex(T_for_LD[(ll, dd)])
                    tt = T_for_LD[(ll, dd)][t]
                    Data = load("Data/$FolderName/$(ll)_$(tt)_$(dd).jld2")
                    x = rhat(Chains(Data["Energies"], ["my param"])).nt[2][1]   # Gelman-Rubin convergnce test
                    push!(conv, (x, ll, tt, dd))
                end
            end
        end
    else        # else check for the given latice size
        for d in eachindex(D_for_L[ll])
            dd = D_for_L[ll][d]
            for t in eachindex(T_for_LD[(ll, dd)])
                tt = T_for_LD[(ll, dd)][t]
                Data = load("Data/$FolderName/$(ll)_$(tt)_$(dd).jld2")
                x = rhat(Chains(Data["Energies"], ["my param"])).nt[2][1]   # Gelman-Rubin convergnce test
                push!(conv, (x, tt, dd))
            end
        end
    end
    return conv
end

# α : specific heat                     Ising : 0       XY : NOP Essential singularity
# β : zero field mag                    Ising : 1/8     XY : NOP No magnetization (To have a nonzero 𝑀, the correlation function must approach a constant at large distance. But in the 2D XY model: For 𝑇>𝑇BKT: correlations decay exponentially. For 𝑇<𝑇BKT: correlations decay as a power law)
# γ : zero field isothermal suscxx      Ising : 7/4     XY : NOP Essential divergence
# δ : Critical isothermal               Ising : 15      XY : 15
# ν : corr length                       Ising : 1       XY : NOP 𝜉 diverges exponentially

# continuous symmetries cannot be spontaneously broken in 2D