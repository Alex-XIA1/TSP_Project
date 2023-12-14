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

# Calcul du cout de l'anneau + des etoiles
function calc_cost(liens,ring,d)
    somme = 0
    
    # Somme pour l'anneau
    for i in 1:size(ring,1)
      # si i == i+1, on a atteint la station de depart
      if ring[1] == ring[i]
        break
      end
      source = ring[i]
      dest = ring[i+1]
      somme = somme + d[source,dest]
    end
  
    # Somme pour les zones d'habitation
    for i in 1:size(liens,1)
      tmp = liens[i]
      # La presentation fait qu'on a toujours la station en premier
      median = tmp[1]
      for j in 2:size(tmp,1)
        # Facteur 10 pour marcher jusqu'a la station
        somme = somme + 10*d[tmp[j],median]
      end
    end
    
    return somme
  end
  
function pair_DistRing(ring, d, i, j)
  res = 0
  # indice dans l'anneau
  source = i
  dest = j 
  while ring[source] != ring[dest]
    res = res + d[ring[source], ring[source+1]]
    source = mod(source, size(ring,1)-1) + 1
    #println("SOurce ",source, " Fin ",ring[dest])
  end
  return res
end

# Calcul du second critere et 3e critere
function calc_meanTime(liens, ring, d, G)
  # temps de marche
  res = Float64[]
  # le temps passé en métro
  resmetro = Float64[]
  resratio = Float64[]

  for i in 1:G.nb_points
    for j in i+1:G.nb_points
      if i != j
        # Le temps si on va directement
        tmp = d[i,j]
        tmp = tmp*10
        # si le chemin n'est pas sur la ligne : On peut optimiser le code
        if i in ring && j in ring
          indStart = findall(x->x==i,ring)[1]
          indFin = findall(x->x==j,ring)[1]
          road1 = pair_DistRing(ring,d,indStart,indFin)
          road2 = pair_DistRing(ring,d,indStart,indFin)
          if road1 > road2
            if road2 > tmp
              push!(res,tmp)
            else
              push!(resmetro,road2)
            end
          else
            if road1 > tmp
              push!(res,tmp)
            else
              push!(resmetro,road1)
            end
          end
          continue
        end

        lsource = -1
        ldest = -1
        # calcul a partir de prise du chemin par metro
        for k in liens
          # On cherche a quel station le premier point est dans
          if i in k
            lsource = k
            break
          end
        end
        # on cherche le second point
        for l in liens
          if j in l
            ldest = l
            break
          end
        end

        # On separe la distance metro et marche pour faire plus facilement le ratio
        tmp2 = 0
        # Distance de source vers la station
        tmp2 = tmp2 + 10*d[i, lsource[1]]
        # Distance de la stations de fin vers la destination
        tmp2 = tmp2 + 10*d[j, ldest[1]]
        # Distance entre les 2 stations lsource -> ldest
        indStart = findall(x->x==lsource[1],ring)[1]
        indFin = findall(x->x==ldest[1],ring)[1]

        # Distance entre les 2 stations de metro
        tmp3=0
        # On doit choisir le sens qui minimise la distance
        road1 = pair_DistRing(ring,d,indStart,indFin)
        road2 = pair_DistRing(ring,d,indFin,indStart)
        if road1 > road2
          tmp3 = road2
        else
          tmp3 = road1
        end

        if tmp2 + tmp3 > tmp
          push!(res,tmp)
          # c'est mieux de marcher
          push!(resratio,1)
        else
          # C'est mieux de marcher et de prendre le metro
          push!(res,tmp2+tmp3)
          push!(resratio,tmp2/(tmp2+tmp3))
        end
      end
    end
  end

  ratioMarche = mean(resratio)
  ratioMetro = 1 - ratioMarche

  return mean(res), ratioMarche, ratioMetro

end


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

file1 = "../Instances_TSP/burma14.tsp"
file2 = "../Clustering/burma14.txt"
res = solve(file1,file2)

