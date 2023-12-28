include("TSP_Etoile.jl")
include("TSP_Etoile_NC.jl")
include("MetroMeta.jl")

medians = 6

# Dans les expériences ici on ne s'interesse qu'au temps de resolution

# input = "../Instances_TSP/berlin52.tsp"
# cluster = "../Clustering/berlin52.txt"
input = "../Instances_TSP/pr107.tsp"
cluster = "../Clustering/6pr107.txt"


# Utiliser @elapsed pour recuperer seulement les temps et @time si on a besoin de la solution
# Attention pour l'heuristique il faut utiliser le bon fichier resultant du clustering
println("\nDebut methode heuristique\n")
# resHeur = @time @CPUtime solve(input,cluster)
resHeur = @elapsed @CPUelapsed solve(input,cluster)
println("Fin heuristique")

println("\nDebut methode PLNE compact\n")
#resCompact = @time @CPUtime solveE(input,medians)
resCompact = @elapsed @CPUelapsed solveE(input,medians)
println("Fin PLNE compact")
# sleep(10)

println("\nDebut methode PLNE non compact")
#resNonC = @time @CPUtime solveE_NC(input,medians)
resNonC = @elapsed @CPUelapsed solveE_NC(input,medians)
println("Fin PLNE non compact")

println("Les temps sont ")
println("Heur (sans le temps du clustering) : ",resHeur)
println("Compact : ",resCompact)
println("Non Compact : ",resHeur)

