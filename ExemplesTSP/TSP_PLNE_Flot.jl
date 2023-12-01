using JuMP
using CPLEX
using CPUTime

include("TSP_IO.jl")


# fonction principale pour notre PLNE MTZ, graphe orienté ici
function PLNE_TSP_MTZ(G)

	c = calcul_dist(G)
	
	println("Création du PLNE")

	m = Model(CPLEX.Optimizer)

	@variable(m, x[1:G.nb_points, 1:G.nb_points], Bin)
	@variable(m, u[1:G.nb_points], Int)
	@variable(m, z[1:G.nb_points, 1:G.nb_points])

	@objective(m, Min, sum((sum(x[i, j] * c[i, j]) for j = 1:G.nb_points) for i = 1:G.nb_points ) )

	# contraite pour dire "une arête entrante dans un sommet"
	for i in 1:G.nb_points
		@constraint(m, (sum(x[i, j] for j in 1:G.nb_points))== 1)
	end
	
	# contraite pour dire "une arête sortante dans un sommet"
	for j in 1:G.nb_points
		@constraint(m, (sum(x[i, j] for i in 1:G.nb_points))== 1)
	end
	
	# pas d'arêtes qui commencent et finissent dans le même point !
	for i in 1:G.nb_points
		@constraint(m, x[i, i] == 0)
	end
	
	# il faut une certaine symétrie...
	for i in 1:G.nb_points
		for j in 1:G.nb_points
			@constraint(m, x[i, j] + x[j, i] <= 1 )
		end
	end


	# contrainte de flot 1 (5)
    @constraint(m, sum(z[1, j] for j in 2:G.nb_points) == G.nb_points-1)

    # contrainte de flot 2 (6)
    for i in 2:G.nb_points
      @constraint(m, ( sum(z[j, i] for j in 1:i-1) + sum(z[j, i] for j in i+1:G.nb_points) ) == ( sum(z[i, j] for j in 2:i-1) + sum(z[i, j] for j in i+1:G.nb_points) + 1))
    end

    # contrainte de flot 3 (7)
    for i in 1:G.nb_points
      for j in 2:G.nb_points
        # on veut i dans V \ {j}
        if i != j
          @constraint(m, z[i, j] + z[j, i] <= (G.nb_points-1)*(x[i, j]+x[j,i]))
        end
      end
    end


	# contrainte liee a la variable zij
    for i in 1:G.nb_points
		for j in 2:G.nb_points
		  # on veut i dans V \ {j}
		  if i != j
			# works ?
			@constraint(m, 0 <= z[i, j] )
		  end
	    end
	end
		

	print(m)
	println()
	
	println("Résolution du PLNE par le solveur")
	optimize!(m)
   	println("Fin de la résolution du PLNE par le solveur")
   	
	#println(solution_summary(m, verbose=true))

	status = termination_status(m)

	# un petit affichage sympathique
	if status == JuMP.MathOptInterface.OPTIMAL
		println("Valeur optimale = ", objective_value(m))
		println("Solution primale optimale :")
		
		S = Int64[]
        i=1
        j=2
        while (value(x[i, j]) < 0.999) 
     	  j=j+1
        end
        push!(S,1)
        push!(S,j)
        i=j
        while (i!=1)
          j=1
          while  ( j==i || value(x[i,j]) < 0.999 ) 
             j=j+1
          end
          push!(S,j)
          i=j
		end
		 println("Temps de résolution :", solve_time(m))
		 return S
	else
		 println("Problème lors de la résolution")
	end

end


function solve(filename)
	I = Read_undirected_TSP(filename)
	  
	  #filename_inst = replace(filename, ".tsp" => "_inst")
	  #WritePdf_visualization_TSP(I, "filename_inst")
	  #WritePdf_visualization_TSP(I, "test")
	
  
	@time @CPUtime S_STAR=PLNE_TSP_MTZ(I) # on le résout
   
	#   val_STAR=Compute_value_TSP(I, S_STAR)
	#   println("Solution Etoile :S=",S_STAR)
	#   println("Valeur: ",val_STAR)
	#   println()
	  
	  filename_STAR = replace(filename, ".tsp" => "_STAR")
  
	  WritePdf_visualization_solution_ordre(I,S_STAR,filename_STAR)
  end
  
  input = "./Instances_TSP/berlin52.tsp"
  solve(input)