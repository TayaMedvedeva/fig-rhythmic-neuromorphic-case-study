module Utils

export heaviside, pulse, delta, range_float

# heaviside / pulse helpers
heaviside(t) = t < 0 ? 0.0 : 1.0
pulse(t, ti, tf) = heaviside(t - ti) - heaviside(t - tf)

# approximate delta: returns 1/dt when t in [a-dt/2, a+dt/2)
function delta(t, a, dt)
    return (a - dt/2 <= t < a + dt/2) ? 1.0/dt : 0.0
end

# small helper: create an array range similar to your code
range_float(a; stop=1.0, length=100) = range(a, stop=stop, length=length)

end # module
