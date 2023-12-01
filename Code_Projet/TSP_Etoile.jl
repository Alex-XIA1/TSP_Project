using JuMP
using CPLEX
using CPUTime
using Graphs
using GraphPlot

include("TSP_IO.jl")
include("TSP_heuristique.jl")
include("TSP_meta_heuristique.jl")
include("TSP_PLNE_MTZ.jl")
include("TSP_BandC.jl")


# PLNE compacte pour l'anneau etoile, p = nombre de stations
function PLNE_compact_star(G,p)

    # formera les dij dans la formulation PLNE
    d = calcul_dist(G)

    println("Creation du PLNE")

    # Optimizer m, not Adam
    m = Model(CPLEX.Optimizer)


    # Variables
    # Les variables xij
    @variable(m, x[1:G.nb_points, 1:G.nb_points], Bin)
    # Les variables yij
    @variable(m, y[1:G.nb_points, 1:G.nb_points], Bin)
    # Les variables zij, i in V, j in V \ {1,i}
    @variable(m, z[1:G.nb_points, 1:G.nb_points])

    # Fonction Objective
    # somme dij*xij + somme di,j*yi,j
    # ij appartient à E et i,j appartient à V*V
    @objective(m, Min, sum(sum(d[i, j] * x[i, j] for j = 1:G.nb_points) for i = 1:G.nb_points)  +  sum((sum(d[i, j] * y[i, j]) for j = 1:G.nb_points) for i = 1:G.nb_points ) )

    # Contraintes
    # p (medians) stations au maximum (1)
    @constraint(m, (sum(y[i, i] for i in 1:G.nb_points))== p)

    # un point i est affecte a un seul median (2)
    for i in 1:G.nb_points
      @constraint(m, (sum(y[i, j] for j in 1:G.nb_points))== 1)
    end

    # contrainte (3) : si j n'est pas une station, on affecte pas i a j
    for i in 1:G.nb_points
      for j in 1:G.nb_points
        # on veut i dans V \ {j}
        if i != j
          @constraint(m, y[i, j] <= y[j, j])
        end
      end
    end

    # contrainte (4) : une station a exactement une arête entrante et sortante (edges constraint)
    for i in 1:G.nb_points
      @constraint(m, (sum(x[j, i] for j in 1:i-1)) + (sum(x[i, j] for j in i+1:G.nb_points))== 2 * y[i,i])
    end

    # contrainte de flot 1 (5)
    @constraint(m, (sum(z[1, j] for j in 2:G.nb_points) == p-1))

    # contrainte de flot 2 (6)
    for i in 2:G.nb_points
      @constraint(m, ( sum(z[j, i] for j in 1:i-1) + sum(z[j, i] for j in i+1:G.nb_points) ) == ( sum(z[i, j] for j in 2:i-1) + sum(z[i, j] for j in i+1:G.nb_points) + y[i, i]))
    end

    # contrainte de flot 3 (7)
    for i in 1:G.nb_points
      for j in 2:G.nb_points
        # on veut i dans V \ {j}
        if i != j
          @constraint(m, z[i, j] + z[j, i] <= (p-1)*x[i, j])
        end
      end
    end

    # contrainte liee a la variable zij
    for i in 1:G.nb_points
      for j in 2:G.nb_points
        # on veut i dans V \ {j}
        if i != j
          # works ?
          @constraint(m, 0 <= z[i, j] <= p-1 )
        end
      end
    end

    # contrainte y11 = 1
    @constraint(m, y[1, 1] == 1)
    # contrainte y1j = 0
    for j in 2:G.nb_points
      @constraint(m, y[1, j] == 0)
    end

  print(m)
	println()
	
	println("Résolution du PLNE par le solveur")
	optimize!(m)
   	println("Fin de la résolution du PLNE par le solveur")
   	
	#println(solution_summary(m, verbose=true))

	status = termination_status(m)
  #println(status)

	# un petit affichage sympathique, A MODIFIER, On prend les yii median, ainsi que les points associes a chaque median ensuite les liens entre stations
  # quelques notes : On peut creer des vecteurs en julia et mettre des vecteurs dans une liste utile pour generaliser
	if status == JuMP.MathOptInterface.OPTIMAL
		println("Valeur optimale = ", objective_value(m))
		println("Solution primale optimale :")

    # Recuperation des stations
    Stations = Int64[]
    for i in 1:G.nb_points
      if (value(y[i, i]) > 0.999)
        push!(Stations, i)
      end
    end

    # Recuperation des points lies a chaque station
    Liens = []
    # On parcours chaque station
    for j in 1:size(Stations)[1]
      # initialisation avec la station 
      Tmp = Int64[]
      push!(Tmp, value(Stations[j]))
      for i in 1:G.nb_points
        # j un median et i un point affecte, i != j, on utilise directement les stations pour aller plus vite
        if (value(y[i, value(Stations[j])]) > 0.999 && i!=value(Stations[j]))
          push!(Tmp,i)
        end
      end
      # On garde les points lies a une station
      push!(Liens,Tmp)
    end

    # println("LES STATIONS")
    # println(Stations)
    # println("Suite")
    #1 - 19 - 21 - 20 - 1
    #12 - 24 - 27 - 12
    #On a un TSP symetrique donc on peut avoir x1j = 1 pour 2 j differents
    # for i in 1:G.nb_points
    #   for j in 2:G.nb_points
    #     if (value(x[i, j]) >0.999)
    #       println(value(x[i, j]))
    #       println("OK")
    #       println("Start ",i, " End ",j)
          
    #     end
    #   end
    # end

    # Recuperation du TSP pour les stations
		S = Int64[]
        # y[1, 1] = 1 donc Stations[1] contient ce point
        i=1
        j=2
        while (value(x[Stations[i], Stations[j]]) < 0.999) 
     	  j=j+1
        
        end
        # Numero de Station
        push!(S,Stations[i])
        push!(S,Stations[j])
        # On cherche depuis la nouvelle Station
        i=j
        # On ne revient pas en arriere (exemple : on veut eviter 1 - 19 - 1)
        tmp = 1
        while (Stations[i]!=Stations[1])
          j=1
          # Les conditions : j==i, j=tmp, si la station precedente = station qu'on verifie, s'il n y a pas d'arete entre les 2 stations
          while  ( j==i || j==tmp || (value(x[Stations[i],Stations[j]]) <0.999 && value(x[Stations[j],Stations[i]]) < 0.999 ) ) 
            j=j+1
          end
          #println("WHEREEEE ", j," and ",tmp, " and ", Stations[j])
          push!(S,Stations[j])
          # On ne revient pas en arriere, avant le changement de i, on met un taboo sur le retour vers i
          tmp = i
          i=j
		    end
		 println("Temps de résolution :", solve_time(m))
		 return S, Stations, Liens
	else
		 println("Problème lors de la résolution")
	end

end

function solve(filename)
  I = Read_undirected_TSP(filename)
	
	#filename_inst = replace(filename, ".tsp" => "_inst")
    #WritePdf_visualization_TSP(I, "filename_inst")
    #WritePdf_visualization_TSP(I, "test")
  

  @time @CPUtime S_STAR, Stations, Liens=PLNE_compact_star(I, 5) # on le résout
 
	# val_STAR=Compute_value_TSP(I, S_STAR)
	println("Solution Etoile :S=",S_STAR)
  println("Solution des villes : L=",Liens)
  println("Les stations : Stations =",Stations)
	# println("Valeur: ",val_STAR)
	 println()
	
	 filename_STAR = replace(filename, ".tsp" => "_STAR")

	 WritePdf_visualization_solution_projet(I,S_STAR,Liens,filename_STAR)
end

input = "./Instances_TSP/eil51.tsp"
solve(input)