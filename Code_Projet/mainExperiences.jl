include("TSP_Etoile.jl")
include("TSP_Etoile_NC.jl")
include("MetroMeta.jl")

medians = 6

input = "../Instances_TSP/berlin52.tsp"
cluster = "../Clustering/berlin52.txt"

# Utiliser @elapsed pour recuperer seulement les temps et @time si on a besoin de la solution
# Attention pour l'heuristique il faut utiliser le bon fichier resultant du clustering
println("\nDebut methode heuristique\n")
resHeur = solve(input,cluster)
println("Fin heuristique")

println("\nDebut methode PLNE compact\n")
nuageComp, tspComp, liensComp, coutComp, meanComp, ratioMarcheComp  = solveE(input,medians)
println("Fin PLNE compact")

println("\nDebut methode PLNE non compact")
nuageNC, tspNC, liensNC, coutNC, meanNC, ratioMarcheNC  = solveE_NC(input,medians)
println("Fin PLNE non compact")

# resHeur est a une liste de solutions non dominees avec les couts, moyennes et ratio
println("Les resultats :")
println("Heuristique\n")
println(resHeur)

println("\nCompact\n")
println("Le nuage de points ",nuageComp)
println("\n le Tour ",tspComp)
println("\n Les stations en 1 et les affectations",liensComp)

println("\nLes critères : ")
println("Le coût de la solution (critere 1): ",round(coutComp;digits=2))
println("Le temps moyen de la solution (critere 2): ",round(meanComp;digits=2))
println("Les ratios de marche et metro (critère 3): Marche - ",round(ratioMarcheComp;digits=3), " et Metro - ",round(1-ratioMarcheComp;digits=3))
println()

println("\nNon Compact\n")
println("Le nuage de points ",nuageNC)
println("\n le Tour ",tspNC)
println("\n Les stations en 1 et les affectations",liensNC)

println("\nLes critères : ")
println("Le coût de la solution (critere 1): ",round(coutNC;digits=2))
println("Le temps moyen de la solution (critere 2): ",round(meanNC;digits=2))
println("Les ratios de marche et metro (critère 3): Marche - ",round(ratioMarcheNC;digits=3), " et Metro - ",round(1-ratioMarcheNC;digits=3))
println()

