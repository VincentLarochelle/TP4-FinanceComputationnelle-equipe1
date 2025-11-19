using Random, Statistics, Plots, Printf

# Fonction d'énergie de Lennard-Jones
function Φ(r)
    return 4000 * ((0.1/r)^12 - (0.1/r)^6)
end

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

# Fonction de réflexion aux bornes
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
function recuit_simule(T, n_particules=4, n_iter=100000, bounds=(-1.0, 1.0))
    # Initialisation aléatoire des positions
    x = rand(n_particules) * 2 .- 1  # positions entre -1 et 1
    energies = Float64[]
    positions_hist = Vector{Float64}[]
    
    a, b = bounds
    σ = 0.2  # écart-type pour les perturbations
    
    E_courant = energie(x)
    push!(energies, E_courant)
    push!(positions_hist, copy(x))
    
    acceptations = 0
    
    for iter in 1:n_iter
        # Sauvegarder l'état courant
        x_nouveau = copy(x)
        
        # Choisir une particule au hasard
        i = rand(1:n_particules)
        
        # Appliquer une perturbation gaussienne
        perturbation = randn() * σ
        x_nouveau[i] += perturbation
        
        # Appliquer la réflexion aux bornes
        x_nouveau[i] = reflexion(x_nouveau[i], a, b)
        
        # Calculer la nouvelle énergie
        E_nouveau = energie(x_nouveau)
        
        # Critère d'acceptation de Metropolis
        ΔE = E_nouveau - E_courant
        
        if ΔE < 0 || rand() < exp(-ΔE/T)
            x = x_nouveau
            E_courant = E_nouveau
            acceptations += 1
        end
        
        # Enregistrer périodiquement
        if iter % 100 == 0
            push!(energies, E_courant)
            push!(positions_hist, copy(x))
        end
    end
    
    taux_acceptation = acceptations / n_iter
    println("Taux d'acceptation pour T = $T: $(taux_acceptation)")
    
    return energies, positions_hist, x
end

# Simulations pour les deux températures
Random.seed!(1234)  # Pour la reproductibilité

println("Simulation pour T = 0.04")
energies_T1, positions_T1, etat_final_T1 = recuit_simule(0.04)

println("\nSimulation pour T = 0.08")
energies_T2, positions_T2, etat_final_T2 = recuit_simule(0.08)

# Analyse des résultats
function analyser_resultats(energies, T, positions_hist)
    println("\n=== Analyse pour T = $T ===")
    println("Énergie moyenne: $(mean(energies))")
    println("Énergie minimale: $(minimum(energies))")
    println("Énergie maximale: $(maximum(energies))")
    println("Écart-type: $(std(energies))")
    
    # État final
    etat_final = positions_hist[end]
    println("Positions finales: $(round.(etat_final, digits=3))")
    
    # Calcul des distances entre paires
    distances = Float64[]
    for i in 1:4
        for j in (i+1):4
            push!(distances, abs(etat_final[i] - etat_final[j]))
        end
    end
    println("Distances moyennes entre particules: $(round(mean(distances), digits=3))")
    
    return energies
end

# Analyses
energies_T1_clean = analyser_resultats(energies_T1, 0.04, positions_T1)
energies_T2_clean = analyser_resultats(energies_T2, 0.08, positions_T2)

# Tracé des histogrammes
p1 = histogram(energies_T1_clean, bins=50, alpha=0.7, label="T = 0.04", 
               xlabel="Énergie", ylabel="Fréquence", title="Histogramme des énergies - T = 0.04")

p2 = histogram(energies_T2_clean, bins=50, alpha=0.7, label="T = 0.08",
               xlabel="Énergie", ylabel="Fréquence", title="Histogramme des énergies - T = 0.08")

plot(p1, p2, layout=(2,1), size=(800, 600))
savefig("histogrammes_energies.png")

# Tracé de l'évolution des énergies
p3 = plot(energies_T1, label="T = 0.04", xlabel="Itération (×100)", ylabel="Énergie", title="Évolution des énergies")
plot!(energies_T2, label="T = 0.08")
savefig("evolution_energies.png")

# Affichage des positions finales
println("\n=== Positions finales ===")
println("T = 0.04: $(round.(etat_final_T1, digits=3))")
println("T = 0.08: $(round.(etat_final_T2, digits=3))")

# Calcul de l'énergie potentielle pour différentes distances
println("\n=== Analyse de la fonction Φ(r) ===")
distances_test = 0.15:0.05:0.5
for r in distances_test
    println("r = $(r): Φ(r) = $(round(Φ(r), digits=2))")
end