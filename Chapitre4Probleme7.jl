using Random, LinearAlgebra

function permanent_approx(A)
    n = size(A, 1)
    total = 0.0
    for _ in 1:500  
        perm_val = 1.0
        for i in 1:n
            available = findall(A[i,:] .> 0)
            isempty(available) && (perm_val = 0.0; break)
            perm_val *= A[i, rand(available)]
        end
        total += perm_val
    end
    return total / 500
end

# Génère une matrice 8×8 avec 20 uns
function random_matrix_8_20()
    A = zeros(Int, 8, 8)
    A[randperm(64)[1:20]] .= 1
    return A
end

# Propose un voisin en échangeant un 1 et un 0
function neighbor_matrix(A)
    new_A = copy(A)
    ones_pos = findall(==(1), A)
    zeros_pos = findall(==(0), A)
    if !isempty(ones_pos) && !isempty(zeros_pos)
        new_A[rand(ones_pos)] = 0
        new_A[rand(zeros_pos)] = 1
    end
    return new_A
end

# Recuit simulé principal
function simulated_annealing_8_20(;iterations=50000, T_start=10.0, T_end=0.1)
    current = random_matrix_8_20()
    current_perm = permanent_approx(current)
    best_perm, best_matrix = current_perm, copy(current)
    
    for i in 1:iterations
        T = T_start * (T_end/T_start)^(i/iterations)
        
        candidate = neighbor_matrix(current)
        candidate_perm = permanent_approx(candidate)
        
        ΔE = candidate_perm - current_perm
        if ΔE > 0 || rand() < exp(ΔE/T)
            current, current_perm = candidate, candidate_perm
            if candidate_perm > best_perm
                best_perm, best_matrix = candidate_perm, copy(candidate)
            end
        end
        
        # Affichage pour T≈2
        if abs(T - 2.0) < 0.2 && i % 100 == 0
            println("T=$((round(T, digits=2))) - Permanent: $(round(candidate_perm, digits=4))")
        end
    end
    
    println("\nRésultat final:")
    println("Permanent maximal: $(round(best_perm, digits=4))")
    println("Matrice optimale (format 8×8):")
    display(best_matrix)
    println("Vérification: $(sum(best_matrix)) uns sur 20 requis")
    
    return best_matrix, best_perm
end

# Exécution
best_matrix, best_permanent = simulated_annealing_8_20()