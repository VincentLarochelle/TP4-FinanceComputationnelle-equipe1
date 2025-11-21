using Plots

function american_put_binomial_boundary(S0, E, T, r, σ, M)
    dt, p = T/M, 0.5
    u = exp(σ*√dt + (r-0.5σ^2)*dt)
    d = exp(-σ*√dt + (r-0.5σ^2)*dt)
    
    dpowers = [d^(M-i) for i in 0:M]
    upowers = [u^i for i in 0:M]
    S_values = S0 .* dpowers .* upowers
    V = max.(E .- S_values, 0)
    
    boundary = zeros(M+1)
    boundary[M+1] = E
    
    for i = M-1:-1:0
        S_current = S0 .* [d^(M-i+j) * u^j for j in 0:i]
        V_cont = exp(-r*dt) * (p*V[2:i+2] + (1-p)*V[1:i+1])
        V_ex = max.(E .- S_current, 0)
        V = max.(V_cont, V_ex)
        
        boundary_found = false
        for j = 1:length(S_current)
            if V_ex[j] ≥ V_cont[j] && !boundary_found
                boundary[i+1] = S_current[j]
                boundary_found = true
            end
        end
        !boundary_found && (boundary[i+1] = E)
    end
    
    return boundary, [i*dt for i in 0:M], V[1]
end

# Parameters
E, r, σ, T, M = 10.0, 0.06, 0.3, 1.0, 1000
boundary, times, option_value = american_put_binomial_boundary(9.0, E, T, r, σ, M)

# Output
println("Option value: ", option_value)
println("\nExercise boundary at selected times (M = $M):")
for i in [1, 125, 250, 375, 500, 625, 750, 875, 1001]
    i = min(i, length(times))
    println("t = $(round(times[i], digits=3)): S* = $(round(boundary[i], digits=4))")
end

# Plot
plot(times, boundary, fillrange=0, fillalpha=0.4, fillcolor=:red, linewidth=0, label="")
plot!(times, boundary, linewidth=3, color=:blue, label="Optimal Exercise Boundary S*(t)")
hline!([E], linestyle=:dash, color=:black, linewidth=2, label="Exercise Price E = $E")
plot!(title="American Put Exercise Boundary", xlabel="Time t", ylabel="Asset Price S",
      legend=:topright, xlims=(0,T), ylims=(0,12), grid=true, gridwidth=1, gridalpha=0.3)