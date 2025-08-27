module Model

using DifferentialEquations
using ..Utils: heaviside, pulse, delta
export second_order_sys_const!, second_order_sys!, build_callbacks, run_simulation, run_simulation_const

# The core ODE. p is a tuple:
# p = (ξ, ωn, ti_arr_ref, tk_arr_ref, θk_arr_ref, Ti_arr_ref, c, γ, Tc, dt, k, K, θ_ref)
function second_order_sys_const!(du,u,p,t)
    ξ, ωn, ti_arr, tk_arr, θk_arr, Ti_arr, c, γ, Tc, dt, k, K, θ_ref = p

    x1 = u[1]
    x2 = u[2]
    ti = u[3]
    Ti = u[4]
    sc2 = u[5]
    T   = u[6]

    # control input: pulse of width Tc starting at ti and scaled by sc2
    u_ctrl = k * sc2 * pulse(t, ti, ti + Tc) 

    du[1] = x2
    du[2] = -2*ξ*ωn * x2 - ωn^2 * sin(x1) + K*u_ctrl
    du[3] = 0.0
    du[4] = 0.0
    du[5] = 0.0

    # accumulate instantaneous sensory events using delta comb at saved tk_arr times
    Σ = 0.0
    for tk in tk_arr
        Σ += sign(abs(x1) - θ_ref) * delta(t, tk, dt)
    end

    du[6] = -c*T - γ*Σ
end

function second_order_sys!(du,u,p,t)
    ξ, ωn, ti_arr, tk_arr, θk_arr, Ti_arr, c, γ, Tc, dt, k, K, θ_ref = p

    x1 = u[1]
    x2 = u[2]
    ti = u[3]
    Ti = u[4]
    sc2 = u[5]
    T   = u[6]

    # control input: pulse of width Tc starting at ti and scaled by sc2
    u_ctrl = k * sc2 * pulse(t, ti, ti + T) 

    du[1] = x2
    du[2] = -2*ξ*ωn * x2 - ωn^2 * sin(x1) + K*u_ctrl
    du[3] = 0.0
    du[4] = 0.0
    du[5] = 0.0

    # accumulate instantaneous sensory events using delta comb at saved tk_arr times
    Σ = 0.0
    for tk in tk_arr
        Σ += sign(abs(x1) - θ_ref) * delta(t, tk, dt)
    end

    du[6] = -c*T - γ*Σ
end

# Callbacks similar to your notebook: sensory and control events.
# These are simple examples that push time / state into the external arrays (held in p).
function condition_sensory_event(u,t,integrator)
    # sensory event triggers whenever x2 (velocity) crosses anything - here just return u[2]
    return u[2]
end
function affect_sensory_event!(integrator)
    # record current T into Ti and store time into provided arrays
    integrator.u[4] = integrator.u[6]
    push!(integrator.p[4], integrator.t)   # tk_arr (third element of p is ti_arr; 4 is tk_arr)
    push!(integrator.p[5], integrator.u[1]) # θk_arr store x1
end

function condition_control_event(u,t,integrator)
    # control event trigger on x1 (position) crossing
    return u[1]
end
function affect_control_event!(integrator)
    integrator.u[3] = integrator.t
    integrator.u[4] = integrator.u[6]
    integrator.u[5] = sign(integrator.u[2])
    push!(integrator.p[3], integrator.t)   # ti_arr
    push!(integrator.p[6], integrator.u[6]) # Ti_arr
end

# returns a CallbackSet ready to be given to ODEProblem
function build_callbacks()
    cb_s = ContinuousCallback(condition_sensory_event, affect_sensory_event!)
    cb_c = ContinuousCallback(condition_control_event, affect_control_event!)
    return CallbackSet(cb_s, cb_c)
end

"""
run_simulation(; kwargs...)

Runs one simulation and returns the solution object.
Keyword arguments:
  dt, Tfinal, ξ, ωn, c, γ, β_ref, k, K, θ_ref
"""
function run_simulation_const(; dt=1e-3, Tfinal=50.0,
                          ξ=0.1, ωn=8.0, c=0.2, γ=0.0, β_ref=0.1,
                          k=1.0, K=15.0, θ_ref=0.5)

    # arrays used by callbacks (mutable containers)
    ti_arr = []
    tk_arr = []
    θk_arr = []
    Ti_arr = []

    p = (ξ, ωn, ti_arr, tk_arr, θk_arr, Ti_arr, c, γ, β_ref, dt, k, K, θ_ref)

    x0 = [1.1; 0.0; 0.0; 0.0; 0.0; 0.0]
    tspan = (0.0, Tfinal)
    prob = ODEProblem(second_order_sys_const!, x0, tspan, p, callback = build_callbacks())
    sol = solve(prob, Euler(), dt=dt, saveat=0.01)
    return sol, (ti_arr=ti_arr, tk_arr=tk_arr, θk_arr=θk_arr, Ti_arr=Ti_arr)
end

function run_simulation(; dt=1e-3, Tfinal=50.0,
                          ξ=0.1, ωn=8.0, c=0.2, γ=0.0, β_ref=0.1,
                          k=1.0, K=15.0, θ_ref=0.5)

    # arrays used by callbacks (mutable containers)
    ti_arr = []
    tk_arr = []
    θk_arr = []
    Ti_arr = []

    p = (ξ, ωn, ti_arr, tk_arr, θk_arr, Ti_arr, c, γ, β_ref, dt, k, K, θ_ref)

    x0 = [1.1; 0.0; 0.0; 0.0; 0.0; 0.0]
    tspan = (0.0, Tfinal)
    prob = ODEProblem(second_order_sys!, x0, tspan, p, callback = build_callbacks())
    sol = solve(prob, Euler(), dt=dt, saveat=0.01)
    return sol, (ti_arr=ti_arr, tk_arr=tk_arr, θk_arr=θk_arr, Ti_arr=Ti_arr)
end

end # module
