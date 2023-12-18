using JuMP
using CPLEX
using CPUTime
using Graphs
using GraphPlot
using Statistics

struct Solutions
	c    # critere 1,2,3
	tours # les tours
  liens # les liens
end

include("TSP_IO.jl")
include("TSP_NC.jl")


function compare(allc,alltours,allLiens, solc,soltour,solLien)
  # La liste des elements a retirer
  currC = allc
  currtours = alltours
  currLien = allLiens
  toremove = Int64[]
  for i in 1:size(currC,1)
    # Si c'est = 1 pour chaque composante alors allc[i] est dominé par solc
    curr = currC[i]
    bincomp = zeros(Int64,size(curr,1)-1)
    for j in 1:size(curr,1)-1
      # solc[j] est mieux 
      if solc[j] < curr[j]
        bincomp[j] = 1
      end
    end

    # Si c'est = 3 alors solc domine allc[i]
    if sum(bincomp) == 3
      # ON devra retirer l'element a l'indice i
      push!(toremove,i)
    end

    # Un element domine au contraire solc
    if sum(bincomp) == 0
      return allc,alltours,allLiens,false
    end
  end

  # On retire tous les elements domines
  if size(toremove) != 0
    deleteat!(currC, toremove)
    deleteat!(currtours, toremove)
    deleteat!(currLien, toremove)
  end

  # On ajoute la nouvelle solution
  push!(currC,solc)
  push!(currtours,soltour)
  push!(currLien,solLien)
    
  return currC,currtours,currLien,true
end

function solve(filename,clusterfile)
  I = Read_undirected_TSP(filename)
  Liens = read_Clustering(clusterfile)
  Stations = map(x->x[1],Liens)

  #filename_inst = replace(filename, ".tsp" => "_inst")
  #WritePdf_visualization_TSP(I, "filename_inst")
  #WritePdf_visualization_TSP(I, "test")

  c = calcul_dist(I)

  @time @CPUtime TSP_Sol =BandC_TSP(I,Stations, c) # on le résout
  # On doit changer par les stations
  RealSol = Int64[]
  for e in TSP_Sol
    push!(RealSol,Stations[e])
  end

  # val_STAR=Compute_value_TSP(I, S_STAR)
  # println("Solution TSP :S=",RealSol)
  # println("Solution des villes (par clustering) : L=",Liens)
  # println("Les stations : Stations =",Stations)
  alld = calcul_dist(I)

  # PARTIE CALCUL DES CRITERES
  # Les liens entre stations, l'anneau et les distances
  cost = calc_cost(Liens,RealSol, alld)
  # println("Le coût de la solution (critere 1): ",round(cost;digits=2))
  # # println("Valeur: ",val_STAR)
  mean, ratioMarche, ratioMetro = calc_meanTime(Liens,RealSol,alld,I)
  # println("Le temps moyen de la solution (critere 2): ",round(mean;digits=2))
  # println("Les ratios de marche et metro (critère 3): Marche - ",round(ratioMarche;digits=3), " et Metro - ",round(ratioMetro;digits=3))
  # println()

  #filename_STAR = replace(filename, ".tsp" => "_HEURI")

  # META-HEURISTIQUE VOISINAGE 1-1 (Pareto Local Search sans reiterer plusieurs fois car ce serait TROP long)
  # Stations = les stations
  
  allc = []
  push!(allc,[round(cost;digits=2),round(mean;digits=2),round(ratioMarche;digits=3),round(ratioMetro;digits=3)])
  tours = []
  push!(tours, RealSol)
  aLiens = []
  push!(aLiens,Liens)

  # Pour chaque cluster
  println(" \nPHASE iteration unique du PLS\n")
  # println(" Initialement ",allc, " ", tours," ",aLiens)
  # println(Stations)
  # println(Liens)
  
  for i in 1:size(Liens,1)
    # On essaye de changer le centre du cluster i
    for j in 2:size(Liens[i],1)
      copyL = copy(Liens)
      copyL[i] = copy(Liens[i])
      newSta = copy(Stations)
      # On change le centre la station est toujours en 1
      newSta[i] =copyL[i][j]
      # On place le nouveau centre au debut pour facilite la reutilisation des fonctions de calcul de criteres
      copyL[i][j] = copyL[i][1]
      copyL[i][1] = newSta[i]


      @time @CPUtime TSP_Sol =BandC_TSP(I,newSta, c) # on le résout avec le nouveau centre
      
      RealSol = Int64[]
      for e in TSP_Sol
        push!(RealSol,newSta[e])
      end

      # PARTIE CALCUL DES CRITERES
      # Les liens entre stations, l'anneau et les distances
      cost = calc_cost(copyL,RealSol, alld)
      # println("Valeur: ",val_STAR)

      mean, ratioMarche, ratioMetro = calc_meanTime(copyL,RealSol,alld,I)


      # Vecteur contenant les nouveaux criteres
      tmpSol = [round(cost;digits=2),round(mean;digits=2),round(ratioMarche;digits=3),round(ratioMetro;digits=3)]

      # Il faut maintenant faire la comparaison pour verifier qu'il n'y a pas de solutions dominées au sens de pareto
      allc, tours, aLiens, bool = compare(allc,tours,aLiens,tmpSol,TSP_Sol,copyL)

    end
  end

  finalSol = Solutions(allc,tours,aLiens)

  #WritePdf_visualization_solution_projet(I,RealSol,Liens,filename_STAR)
  return finalSol
end

# file1 = "../Instances_TSP/berlin52.tsp"
# file2 = "../Clustering/berlin52.txt"
# res = solve(file1,file2)

