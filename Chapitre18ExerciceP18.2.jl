using Plots

# Arbre binomial pour référence
function binomial_am_put(S0=40, K=39, r=0.05, σ=0.3, T=1, N=100)
    dt = T/N
    u, d = exp(σ*√dt), exp(-σ*√dt)
    p = (exp(r*dt) - d)/(u - d)
    
    # Calcul récursif
    function node_value(i, j, S)
        i == N && return max(K-S, 0)
        cont = exp(-r*dt) * (p*node_value(i+1, j+1, S*u) + (1-p)*node_value(i+1, j, S*d))
        return max(cont, max(K-S, 0))
    end
    return node_value(0, 0, S0)
end

# Test de stratégies d'exercice
function test_strategies()
    S0, K, r, σ, T = 40, 39, 0.05, 0.3, 1
    ref = binomial_am_put()
    println("Valeur référence (arbre binomial): $(round(ref, digits=4))")
    
    strategies = [
        ("Seuil 0.2K", (S,t) -> max(K-S,0) >= 0.2*K),
        ("Seuil 0.3K", (S,t) -> max(K-S,0) >= 0.3*K), 
        ("Niveau 0.8K", (S,t) -> S <= 0.8*K),
        ("Niveau 0.9K", (S,t) -> S <= 0.9*K),
        ("Temps 0.8", (S,t) -> t >= 0.8*T)
    ]
    
    results = []
    for (name, strat) in strategies
        values = Float64[]
        for _ in 1:10000
            S, t, exercised = S0, 0.0, false
            while t < T && !exercised
                dt = 0.01
                S *= exp((r-σ^2/2)*dt + σ*√dt*randn())
                t += dt
                if strat(S, t)
                    push!(values, max(K-S,0)*exp(-r*t))
                    exercised = true
                end
            end
            !exercised && push!(values, max(K-S,0)*exp(-r*T))
        end
        value = mean(values)
        push!(results, (name, value))
        println("$name: $(round(value, digits=4))")
    end
    
    # Graphique
    names = [r[1] for r in results]
    vals = [r[2] for r in results]
    bar(names, vals, legend=false, title="Stratégies d'exercice - Put Américain")
    hline!([ref], linestyle=:dash, color=:red, label="Valeur référence")
    ylabel!("Valeur de l'option")
end

test_strategies()