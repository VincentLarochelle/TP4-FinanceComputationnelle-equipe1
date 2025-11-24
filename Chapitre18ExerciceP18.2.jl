using Plots, Random, Statistics, Distributions

# Calcul EXACT du put européen
function exact_european_put(S0, K, r, σ, T)
    d1 = (log(S0/K) + (r + σ^2/2)*T) / (σ*sqrt(T))
    d2 = d1 - σ*sqrt(T)
    K*exp(-r*T)*cdf(Normal(), -d2) - S0*cdf(Normal(), -d1)
end

# Calcul EXACT du put américain par arbre binomial
function exact_american_put(S0, K, r, σ, T, N=100)
    dt = T/N
    u = exp(σ*sqrt(dt))
    d = 1/u
    p = (exp(r*dt) - d)/(u - d)
    
    # Initialisation des valeurs à l'échéance
    V = zeros(N+1)
    for i in 0:N
        S = S0 * u^i * d^(N-i)
        V[i+1] = max(K - S, 0)
    end
    
    # Rétropopagation
    for step in (N-1):-1:0
        for i in 0:step
            S = S0 * u^i * d^(step-i)
            continuation = exp(-r*dt) * (p*V[i+2] + (1-p)*V[i+1])
            exercise = max(K - S, 0)
            V[i+1] = max(continuation, exercise)
        end
    end
    return V[1]
end

# Test des stratégies
function figure_18_5()
    S0, K, r, σ, T = 40, 39, 0.05, 0.3, 1
    
    # Calcul des valeurs de référence
    euro_value = exact_european_put(S0, K, r, σ, T)
    amer_value = exact_american_put(S0, K, r, σ, T)
    
    println("VALEURS DE RÉFÉRENCE:")
    println("Put européen: $(round(euro_value, digits=4))")
    println("Put américain: $(round(amer_value, digits=4))")
    println()
    
    # Test des stratégies
    results = []
    for α in 0:0.5:10
        vals = Float64[]
        for _ in 1:3000
            S, t = S0, 0.0
            dt = 0.01
            exercised = false
            
            while t < T && !exercised
                S *= exp((r - σ^2/2)*dt + σ*sqrt(dt)*randn())
                t += dt
                payoff = max(K - S, 0)
                if payoff > α
                    push!(vals, payoff * exp(-r*t))
                    exercised = true
                end
            end
            
            if !exercised
                push!(vals, max(K - S, 0) * exp(-r*T))
            end
        end
        value = mean(vals)
        push!(results, (α, value))
        println("α = $α: $(round(value, digits=4))")
    end
    
    # Graphique
    αs = [r[1] for r in results]
    vals = [r[2] for r in results]
    
    plt = plot(αs, vals, 
               marker=:circle, 
               linewidth=2, 
               label="Stratégie α", 
               title="Stratégie d'exercice - Put Américain",
               xlabel="α",
               ylabel="Valeur du Put",
               legend=:bottomright,
               grid=true)
    
    # Lignes de référence 
    hline!([amer_value], 
           linestyle=:dash, 
           color=:red, 
           linewidth=2,
           label="Put Américain ($(round(amer_value, digits=3)))")
    
    hline!([euro_value], 
           linestyle=:dash, 
           color=:blue, 
           linewidth=2,
           label="Put Européen ($(round(euro_value, digits=3)))")
    
    # Point optimal
    opt_idx = argmax(vals)
    opt_α = αs[opt_idx]
    opt_val = vals[opt_idx]
    
    scatter!([opt_α], [opt_val], 
             color=:green,
             markersize=6,
             label="Optimum α=$(opt_α)")
    
    return plt, euro_value, amer_value, opt_α, opt_val
end

# Exécution
plt, euro, amer, opt_α, opt_val = figure_18_5()
display(plt)

println("\nAvec paramètres S0=40, K=39, r=0.05, σ=0.3, T=1")
