include("TSP_Etoile.jl")
include("TSP_Etoile_NC.jl")
include("MetroMeta.jl")

medians = 6

input = "../Instances_TSP/burma14.tsp"
cluster = "../Clustering/6burma14.txt"
# input = "../Instances_TSP/berlin52.tsp"
# cluster = "../Clustering/6berlin52.txt"
# input = "../Instances_TSP/eil101.tsp"
# cluster = "../Clustering/20eil101.txt"
# input = "../Instances_TSP/a280.tsp"
# cluster = "../Clustering/6a280.txt"
# input = "../Instances_TSP/att532.tsp"
# cluster = "../Clustering/6att532.txt"
# cluster2 = "../Clustering/6att532_2.txt"
# cluster3 = "../Clustering/6att532_3.txt"

nuage = Read_undirected_TSP(input)

#Utiliser @elapsed pour recuperer seulement les temps et @time si on a besoin de la solution
# Pour visualiser le résultat -> décommenter les lignes WritePdf_visualization_solution_projet et changer le nom du dernier paramètre !

println("\nDebut methode PLNE compact\n")
nuageComp, tspComp, liensComp, coutComp, meanComp, ratioMarcheComp  = @time @CPUtime solveE(input,medians)
println("Fin PLNE compact")

println("\nDebut methode PLNE non compact")
nuageNC, tspNC, liensNC, coutNC, meanNC, ratioMarcheNC  = @time @CPUtime solveE_NC(input,medians)
println("Fin PLNE non compact")
# sleep(10)

# # Attention pour l'heuristique il faut utiliser le bon fichier resultant du clustering
# # Note : Les nuages de points sont les memes pour les 3
println("\nDebut methode heuristique\n")
critH, tourH, liensH = @time @CPUtime solve(input,cluster)
# critH, tourH, liensH = @time @CPUtime solve(input,cluster2)
# critH, tourH, liensH = @time @CPUtime solve(input,cluster3)


# resHeur est a une liste de solutions non dominees avec les couts, moyennes et ratio
println("Les resultats :")

println("\nCompact\n")
println("Le nuage de points ",nuageComp)
println("\n le Tour ",tspComp)
println("\n Les stations en 1 et les affectations",liensComp)
#La solution compact donne la meme solution que la version Non compacte
# WritePdf_visualization_solution_projet(nuageComp,tspComp,liensComp,"6berlin52solutionC")

println("\nLes critères : ")
println("Le coût de la solution (critere 1): ",round(coutComp;digits=2))
println("Le temps moyen de la solution (critere 2): ",round(meanComp;digits=2))
println("Les ratios de marche et metro (critère 3): Marche - ",round(ratioMarcheComp;digits=3), " et Metro - ",round(1-ratioMarcheComp;digits=3))
println()


println("\nNon Compact\n")
# println("Le nuage de points ",nuageNC)
println("\n le Tour ",tspNC)
println("\n Les stations en 1 et les affectations",liensNC)
# WritePdf_visualization_solution_projet(nuageNC,tspNC,liensNC,"6eil101solutionNC")

println("\nLes critères : ")
println("Le coût de la solution (critere 1): ",round(coutNC;digits=2))
println("Le temps moyen de la solution (critere 2): ",round(meanNC;digits=2))
println("Les ratios de marche et metro (critère 3): Marche - ",round(ratioMarcheNC;digits=3), " et Metro - ",round(1-ratioMarcheNC;digits=3))
println()

println("\nHeuristique\n")
println("\n le Tour ",tourH)
println("\n Les stations en 1 et les affectations",liensH)

println("\nLes critères : ")
println("Le coût de la solution (critere 1): ",round(critH[1];digits=2))
println("Le temps moyen de la solution (critere 2): ",round(critH[2];digits=2))
println("Les ratios de marche et metro (critère 3): Marche - ",round(critH[3];digits=3), " et Metro - ",round(1-critH[3];digits=3))
println()

# WritePdf_visualization_solution_projet(nuage,tourH,liensH,"20eil101solutionHeur")
# println("Fin heuristique")