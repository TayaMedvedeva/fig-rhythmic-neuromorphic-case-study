module Analysis

using ..Model   # <- THIS makes Model visible inside Analysis
using ..Utils
using LinearAlgebra, Statistics
export G_func, N_func, find_best_beta_w, HB_amplitude, func_th_beta
export amplitude_vs_beta

# Frequency response G(jw) for second-order system with gain K
function G_func(w, K, ωn, ξ)
    return K / (ωn^2 - w^2 + im * 2 * ξ * ωn * w)
end

# describing function N (if you need it)
function N_func(A, w, β, k)
    a1 = 2 * k / π * sin(w * β)
    b1 = 2 * k / π * (1 - cos(w * β))
    return (a1 + im * b1) / A
end

"""
    find_best_beta_w(β_range, w_range; K, ωn, ξ, k, θ_ref)
Search grid for (w, β) minimizing the error metric.
"""
function find_best_beta_w(β_range, w_range; K=15.0, ωn=8.0, ξ=0.1, k=1.0, θ_ref=0.5)
    errors = zeros(length(w_range), length(β_range))
    for i in eachindex(w_range)
        w = w_range[i]
        G = G_func(w, K, ωn, ξ)
        for j in eachindex(β_range)
            β = β_range[j]
            err1 = abs((π - w*β)/2 + angle(G))
            err2 = abs(4*k/π * sin(w*β/2) * abs(G) - θ_ref)
            errors[i, j] = err1 + err2
        end
    end

    # convert linear index -> (i,j) robustly using CartesianIndices
    li = argmin(errors)                    # linear index of minimum
    ci = CartesianIndices(errors)[li]      # CartesianIndex(i, j)
    i, j = Tuple(ci)

    w_best = w_range[i]
    β_best = β_range[j]
    return w_best, β_best, errors
end

"""
    HB_amplitude(w, β; K, ωn, ξ, k)
Analytic amplitude: A_hat = 4*k/pi * sin(w*β/2) * |G(jw)|
"""
function HB_amplitude(w, β; K=15.0, ωn=8.0, ξ=0.1, k=1.0)
    G = G_func(w, K, ωn, ξ)
    return 4*k/π * sin(w*β/2) * abs(G)
end

function amplitude_vs_beta(β_arr; ξ=0.1, ωn=8.0, k=1.0, K=15.0,
                           θ_ref=0.5, c=0.2, γ=0.0, dt=1e-3, Tfinal=50.0)

    A_arr_sim, w_arr_sim = Float64[], Float64[]
    A_arr_anal, w_arr_anal_HB = Float64[], Float64[]
    A_arr_anal_1 = Float64[]

    ωd = ωn * sqrt(1 - ξ^2)

    for β in β_arr
        sol, _ = Model.run_simulation_const(; dt=dt, Tfinal=Tfinal,
                                      ξ=ξ, ωn=ωn, c=c, γ=γ,
                                      β_ref=β, k=k, K=K, θ_ref=θ_ref)
        Amax = maximum(sol[1, floor(Int, length(sol[1,:])/3):end])
        push!(A_arr_sim, Amax)

        idxs = findall(diff(sign.(sol[1, floor(Int, length(sol[1,:])/3):end])) .!= 0)
        t_crossings = sol.t[idxs]
        mean_dt = mean(diff(t_crossings))
        push!(w_arr_sim, 2π/mean_dt)

        errors = [abs((π - w*β)/2 + angle(G_func(w,K,ωn,ξ))) for w in range(0.1, stop=10, length=500)]
        w_best = range(0.1, stop=10, length=500)[argmin(errors)]
        push!(w_arr_anal_HB, w_best)

        Aanal = 4*k/π * sin(w_best*β/2) * abs(G_func(w_best,K,ωn,ξ))
        push!(A_arr_anal, Aanal)

        Aanal1 = 4*k/π * sin(ωd*β/2) * abs(G_func(ωd,K,ωn,ξ))
        push!(A_arr_anal_1, Aanal1)
    end

    return β_arr, A_arr_sim, w_arr_sim, A_arr_anal, w_arr_anal_HB, A_arr_anal_1
end

function func_th_beta(β; K=15.0, ωn=8.0, ωd=0, ξ=0.1, k=1.0) 
    absG = K / sqrt((ωn^2 - ωd^2)^2 + (2*ξ*ωn*ωd)^2)
    return 4*k / pi * sin(ωd*β/2) * absG
end

end # module
