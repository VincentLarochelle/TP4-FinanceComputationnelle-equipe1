using Plots, Random, Statistics

# Fonction potentiel de Lennard-Jones
Φ(r) = 4000 * ((0.1/r)^12 - (0.1/r)^6)

# Fonction énergie totale
function energie(x)
    E = 0.0
    n = length(x)
    for i in 1:n
        for j in (i+1):n
            r_ij = abs(x[i] - x[j])
            E += Φ(r_ij)
        end
    end
    return E
end

# Fonction réflexion aux bornes
function reflexion(x, a, b)
    if x < a
        return 2*a - x
    elseif x > b
        return 2*b - x
    else
        return x
    end
end

# Algorithme de recuit simulé
function recuit_simule(T, n_particules=4, n_iter=100000)
    # Initialisation aléatoire
    x = rand(n_particules) * 2 .- 1
    energies = Float64[]
    E_courant = energie(x)
    
    for iter in 1:n_iter
        # Choisir une particule au hasard
        i = rand(1:n_particules)
        
        # Appliquer perturbation gaussienne avec réflexion
        x_nouveau = copy(x)
        x_nouveau[i] = reflexion(x_nouveau[i] + randn() * 0.2, -1.0, 1.0)
        
        # Calculer nouvelle énergie
        E_nouveau = energie(x_nouveau)
        ΔE = E_nouveau - E_courant
        
        # Critère d'acceptation de Metropolis
        if ΔE < 0 || rand() < exp(-ΔE/T)
            x = x_nouveau
            E_courant = E_nouveau
        end
        
        # Enregistrer l'énergie périodiquement
        if iter % 100 == 0
            push!(energies, E_courant)
        end
    end
    
    return energies, x
end

# Simulations pour les deux températures
println("=== Simulation de recuit simulé pour 4 particules ===")
println("Fonction d'énergie: Φ(r) = 4000[(0.1/r)¹² - (0.1/r)⁶]")

Random.seed!(1234)  # Pour reproductibilité

println("\n--- Simulation à T = 0.04 ---")
energies_T1, positions_finales_T1 = recuit_simule(0.04)
println("Énergie moyenne: ", round(mean(energies_T1), digits=2))
println("Énergie minimale: ", round(minimum(energies_T1), digits=2))
println("Positions finales: ", round.(positions_finales_T1, digits=3))

println("\n--- Simulation à T = 0.08 ---")
energies_T2, positions_finales_T2 = recuit_simule(0.08)
println("Énergie moyenne: ", round(mean(energies_T2), digits=2))
println("Énergie minimale: ", round(minimum(energies_T2), digits=2))
println("Positions finales: ", round.(positions_finales_T2, digits=3))

# Tracé des histogrammes
p1 = histogram(energies_T1, bins=40, alpha=0.7, color=:blue, label="T = 0.04",
               title="Distribution des énergies - T = 0.04", 
               xlabel="Énergie E(x)", ylabel="Fréquence",
               xlims=(-8000, 2000))

p2 = histogram(energies_T2, bins=40, alpha=0.7, color=:red, label="T = 0.08",
               title="Distribution des énergies - T = 0.08",
               xlabel="Énergie E(x)", ylabel="Fréquence",
               xlims=(-8000, 2000))

# Affichage côte à côte
plot(p1, p2, layout=(1,2), size=(800, 600))
