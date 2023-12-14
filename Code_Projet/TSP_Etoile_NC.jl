using JuMP
using CPLEX
using CPUTime
using Graphs
using GraphPlot
using Statistics

include("TSP_IO.jl")


function edge_value(x,i,j)
  if (i<j)
     return value(x[i,j])
  else
     return value(x[j,i])
  end
end

# Si x est entier et respecte les contraintes de degre sum_{v=1}^n x[u,v)]=1 pour tout u
# Renvoie la liste des sommets formant un cycle passant par u
function find_cycle_in_integer_x(x, u)
     S = Int64[]
     # le sommet de base du cycle
     #push!(S,u)
     i=u
     prev=-1
     while true
        j=1
        while  ( j==i || j==prev || edge_value(x,i,j) < 0.999 ) 
          j=j+1
        end
        push!(S,j)
        prev=i
        i=j
     #println(i)
        (i!=u) || break   # si i==u alors fin de boucle
     end
     return S
end

# PLNE compacte pour l'anneau etoile, p = nombre de stations
function PLNE_compact_star(G,p)
    # Setting some stat variables
    nbViolatedMengerCut_fromIntegerSep = 0
    nbViolatedMengerCut_fromFractionalSep = 0

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

    # Fonction Objective
    # somme dij*xij + somme di,j*yi,j
    # ij appartient à E et i,j appartient à V*V
    @objective(m, Min, sum(sum(d[i, j] * x[i, j] for j = 1:G.nb_points) for i = 1:G.nb_points)  +  sum((sum(d[i, j] * y[i, j]) for j = 1:G.nb_points) for i = 1:G.nb_points ) )

    # Contraintes
    # p (medians) stations au maximum (1)
    @constraint(m, (sum(y[i, i] for i in 1:G.nb_points))== p)

    # un point i est affecte a un seul median (2)
    for i in 2:G.nb_points
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
      @constraint(m, (sum(x[j, i] for j in 1:i-1)) + (sum(x[i, j] for j in i+1:G.nb_points)) == 2 * y[i,i])
    end


    # contrainte y11 = 1
    @constraint(m, y[1, 1] == 1)
    # contrainte y1j = 0
    for j in 2:G.nb_points
      @constraint(m, y[1, j] == 0)
    end

  # Contrainte pour retirer les tours : Coupe de Menger
  G_sep=complete_digraph(G.nb_points)

  #################
  # our function lazySep_ViolatedMengerCut : SEPARATION ENTIERE
  function lazySep_ViolatedMengerCut(cb_data)
      # cb_data is the CPLEX value of our variables for the separation algorithm
      # In the case of a LazyCst, the value is integer, but sometimes, it is more 0.99999 than 1

      # 1) Reperage des points median
      tmpSta = Int64[]
      for i in 1:G.nb_points
        if (callback_value(cb_data, y[i,i]) >0.999)
          push!(tmpSta,i)
        end
      end

      # par rapport a la solution actuelle, on cherche tous les liens xij
      xsep =zeros(Float64,size(tmpSta)[1], size(tmpSta)[1]);
      for i in 1:size(tmpSta)[1]
          for j in 1:i-1
            # Attention on doit mettre tmpSta[j] < tmpSta[i]
            if tmpSta[j] < tmpSta[i]
              xsep[i,j]=callback_value(cb_data, x[tmpSta[j],tmpSta[i]])
            else
              xsep[i,j]=callback_value(cb_data, x[tmpSta[i],tmpSta[j]])
            end
          end

          for j in i+1:size(tmpSta)[1]
            if tmpSta[j] < tmpSta[i]
              xsep[i,j]=callback_value(cb_data, x[tmpSta[j],tmpSta[i]])
            else
              xsep[i,j]=callback_value(cb_data, x[tmpSta[i],tmpSta[j]])
            end
          end
      end
      # Dans notre cas il faut tout d'abord recuperer les stations !!!! sinon le code ne marche pas

      #println("All the stations are ", tmpSta)

      # for i in 1:G.nb_points
      #   print(xsep[i]," ")
      # end
      # println()
      
#        violated, W = ViolatedMengerCut_IntegerSeparation(G,xsep)
      
      # random (on ne sait pas pourquoi mais c'est dans l'exemple donnee dans les fichiers)
      start=rand(1:size(tmpSta)[1])

      # println("Depart ",start)
      # for i in 1:size(tmpSta)[1]
      #   println(xsep[i,:])
      # end
      # println("OK ",size(tmpSta))

      # On cherche s'il y a un cycle
      W =find_cycle_in_integer_x(xsep, start)
      choices = Int64[]
      
      # Si on trouve un cycle avec 1, on refait (1 ne doit pas etre dans W)
      if 1 in W
        for e in tmpSta
          if e in W
            continue
          push!(choices,e)
          end
        end
        # Si le vecteur choix est nul alors il n'y a pas de sous-tours
        if size(choices,1) != 0
          start=rand(1:size(choices)[1])
          W =find_cycle_in_integer_x(xsep, start)
        end
      end


      # Formulation 1 : (8) et (3)
      # On ajoute les contraintes si il y a une sous-tour. Violation de la contrainte (8)
      if size(W,1)!=size(tmpSta)[1]    # size(W) renvoie sinon (taille,)
        # println("size ",size(W,1), " And ", size(tmpSta,1))
    
          #println(W)
          
          for l in 1:size(tmpSta)[1]
            if (l in W)
              con = @build_constraint(sum(x[i,j] for i ∈ W for j ∈ i+1:G.nb_points if j ∉ W) 
                                      + sum(x[j,i] for i ∈ W for j ∈ 1:i-1 if j ∉ W)  
                                      >= 2 * y[l,l])
              
            #println(con)
              
              MOI.submit(m, MOI.LazyConstraint(cb_data), con) 
            end
          end
          nbViolatedMengerCut_fromIntegerSep=nbViolatedMengerCut_fromIntegerSep+1
          
      end

      # Formulation 2 Optimale (10) dans le sujet
      # On ajoute les contraintes s'il y a une sous-tour (Violation de la contrainte (10))
      # NE FONCTIONNE PAS POUR L'INSTANT
      # if size(W,1)!=size(tmpSta,1)   # size(W) renvoie sinon (taille,)
      #   # println("size ",size(W,1), " And ", size(tmpSta,1))
    
      #     #println(W)
          
      #     for l in tmpSta
      #       if (l in W)
      #         con = @build_constraint(sum(x[i,j] for i ∈ W for j ∈ i+1:G.nb_points if j ∉ W) 
      #                                 + sum(x[j,i] for i ∈ W for j ∈ 1:i-1 if j ∉ W)  
      #                                 >= 2 * (sum(y[l,m] for m in 1:G.nb_points if m ∈ W)))
              
      #       #println(con)
              
      #         MOI.submit(m, MOI.LazyConstraint(cb_data), con) 
      #       end
      #     end
      #     nbViolatedMengerCut_fromIntegerSep=nbViolatedMengerCut_fromIntegerSep+1
          
      # end
      
  end
#
#################




#################
# our function userSep_ViolatedMengerCut
# CBDATA correspond à la solution courante
  function userSep_ViolatedMengerCut(cb_data)
      # cb_data is the CPLEX value of our variables for the separation algorithm
      # In the case of a usercut, the value is fractional or integer (and can be -0.001)
      tmpSta = Int64[]
      for i in 1:G.nb_points
        if (callback_value(cb_data, y[i,i]) >0.1)
          push!(tmpSta,i)
        end
      end

      # Get the x value from cb_data 
      xsep =zeros(Float64,size(tmpSta)[1], size(tmpSta)[1]);
      for i in 1:size(tmpSta)[1]
          for j in 1:i-1
            # Attention on doit mettre tmpSta[j] < tmpSta[i]
            if tmpSta[j] < tmpSta[i]
              xsep[i,j]=callback_value(cb_data, x[tmpSta[j],tmpSta[i]])
            else
              xsep[i,j]=callback_value(cb_data, x[tmpSta[i],tmpSta[j]])
            end
          end

          for j in i+1:size(tmpSta)[1]
            if tmpSta[j] < tmpSta[i]
              xsep[i,j]=callback_value(cb_data, x[tmpSta[j],tmpSta[i]])
            else
              xsep[i,j]=callback_value(cb_data, x[tmpSta[i],tmpSta[j]])
            end
          end
      end

      # On a besoin du digraph complet
      G_sep2=complete_digraph(size(tmpSta)[1])
      
      # mincut pour faire les sous-tours
      Part,valuecut=mincut(G_sep2,xsep)  # Part is a vector indicating 1 and 2 for each node to be in partition 1 or 2
      
      # 1 n'appartient pas a S, pour obtenir une des sous-tours qui nous interessent
      W=Int64[]
      for k in 1:2
        b = true
        for i in 1:size(tmpSta)[1]
          if (Part[i]==k)
            if (tmpSta[i] == 1)
              W=Int64[]
              b = false
            break
            end
          push!(W,tmpSta[i])
          end
        end
        if (b)
        break
        end
      end
      
      if (valuecut<2.0)
    #     println(W)
          
        # i in W
        # Version 1 : contraintes (8) et (3)
        # On peut utiliser G.nb_points parce que les points non utilisees seront = 0
        for l in 1:G.nb_points
          if (l in W)
            con = @build_constraint(sum(x[i,j] for i ∈ W for j ∈ i+1:G.nb_points if j ∉ W) 
                                    + sum(x[j,i] for i ∈ W for j ∈ 1:i-1 if j ∉ W)  
                                    >= 2 * y[l,l])
            
            # println(con)
            
            MOI.submit(m, MOI.UserCut(cb_data), con) 
          end
        end

        # Version 2 : contrainte (10) dans le sujet
        # NE FONCTIONNE PAS POUR L'INSTANT
        # for l in tmpSta
        #   if (l in W)
        #     con = @build_constraint(sum(x[i,j] for i ∈ W for j ∈ i+1:G.nb_points if j ∉ W) 
        #                             + sum(x[j,i] for i ∈ W for j ∈ 1:i-1 if j ∉ W)  
        #                             >= 2 * (sum(y[l,m] for m in 1:G.nb_points if m ∈ W)  ))
            
        #     # println(con)
            
        #     MOI.submit(m, MOI.UserCut(cb_data), con) 
        #   end
        # end
        
        nbViolatedMengerCut_fromFractionalSep=nbViolatedMengerCut_fromFractionalSep+1
        
      end
            
  end
#
#################

#################
# our function primalHeuristicTSP
  function primalHeuristicTSP(cb_data)

    tmpSta = Int64[]
    for i in 1:G.nb_points
      if (callback_value(cb_data, y[i,i]) >0)
        push!(tmpSta,i)
      end
    end
    # println(" Le nombre de stations ",tmpSta[:])

    # Get the x value from cb_data 
    xfrac =zeros(Float64,size(tmpSta)[1], size(tmpSta)[1]);
    for i in 1:size(tmpSta)[1]
        for j in 1:i-1
          # Attention on doit mettre tmpSta[j] < tmpSta[i]
          if tmpSta[j] < tmpSta[i]
            xfrac[i,j]=callback_value(cb_data, x[tmpSta[j],tmpSta[i]])
          else
            xfrac[i,j]=callback_value(cb_data, x[tmpSta[i],tmpSta[j]])
          end
        end

        for j in i+1:size(tmpSta)[1]
          if tmpSta[j] < tmpSta[i]
            xfrac[i,j]=callback_value(cb_data, x[tmpSta[j],tmpSta[i]])
          else
            xfrac[i,j]=callback_value(cb_data, x[tmpSta[i],tmpSta[j]])
          end
        end
    end


    # The global idea is to add the edges one after the other
    # in the order of the x_ij values sorted from the highest to the lowest
    # Adding an edge is valid only 
    # if the edge linked two nodes having a degree < 2 in the solution
    # and if the edge does not form a subtour  (i.e. a cycle of size < n nodes)
    # the detection of the creation of a cycle is done by the techniques     
    # called "union-find" structure where each node is associated with the number
    # of the smallest index of a node linked by a path
    # each time an edge is added, this number (call the connected component)
    # must be updated
      
    sol=zeros(Float64,G.nb_points, G.nb_points);
    
    # Toutes les aretes du TSP sur les stations
    L=[]
    for i in 1:size(tmpSta,1)
        for j in i+1:size(tmpSta,1)
          push!(L,(i,j,xfrac[i,j]))
        end
    end
    # Les aretes sont rangees dans l'ordre decroissante (glouton)
    sort!(L,by = x -> x[3])
    
    CC= zeros(Int64,size(tmpSta,1));  #Connected component of node i
    for i in 1:size(tmpSta,1)
      CC[i]=-1
    end

    # 2 noeuds par station
    tour=zeros(Int64,size(tmpSta,1),2)  # the two neighbours of i in a TSP tour, the first is always filled before de second
    for i in 1:size(tmpSta,1)
        tour[i,1]=-1
        tour[i,2]=-1
    end
    
    cpt=0
    while ( (cpt!=size(tmpSta,1)-1) && (size(L)!=0) )
    
      # On sort la solution la plus grande
      (i,j,val)=pop!(L)   

      if ( ( (CC[i]==-1) || (CC[j]==-1) || (CC[i]!=CC[j]) )  && (tour[i,2]==-1) && (tour[j,2]==-1) ) 
      
          cpt=cpt+1 
          
          # aucun noeud sortant alors on ajoute j
          if (tour[i,1]==-1)  # if no edge going out from i in the sol
          tour[i,1]=j        # the first outgoing edge is j
        CC[i]=i;
          else
        # sinon c'est le second noeud
        tour[i,2]=j        # otherwise the second outgoing edge is j
          end

          # meme chose pour j
          if (tour[j,1]==-1)
        tour[j,1]=i
        CC[j]=CC[i]
          else
        tour[j,2]=i
        
        oldi=i
        k=j

        # On cherche un cycle ?
        # si != -1 on alors un entrant et un sortant
        while (tour[k,2]!=-1)  # update to i the CC of all the nodes linked to j
          if (tour[k,2]==oldi) 
              l=tour[k,1]
            else 
                l=tour[k,2]
            end
          # on connecte a l'element connecte de i
          CC[l]=CC[i]
          oldi=k
          k=l
        end
      end
      end
    end
    
    i1=-1          # two nodes haven't their 2nd neighbour encoded at the end of the previous loop
    i2=0
    for i in 1:size(tmpSta,1)
    if tour[i,2]==-1
      if i1==-1
          i1=i
      else 
          i2=i
      end
    end
    end
    tour[i1,2]=i2
    tour[i2,2]=i1
  
    # Il faut bien faire attention, tour utilise des indices entre 1 et size(tmpSta)
    # tmpSta[indice] pour recuperer le vrai point
    value=0
    for i in 1:size(tmpSta,1)
      for j in 2:size(tmpSta,1)
        if tmpSta[j] > tmpSta[i]
          if ((j!=tour[i,1])&&(j!=tour[i,2]))
            sol[tmpSta[i],tmpSta[j]]=0
          else          
            sol[tmpSta[i],tmpSta[j]]=1      
            value=value+dist(G,tmpSta[i],tmpSta[j])
          end
        end
      end
    end
    
    xvec=vcat([m[:x][i, j] for i = 1:G.nb_points for j = i+1:G.nb_points])
    solvec=vcat([sol[i, j] for i = 1:G.nb_points for j = i+1:G.nb_points])

    MOI.submit(m, MOI.HeuristicSolution(cb_data), xvec, solvec)
  
  end
  #
  #################

  #################
  # Setting callback in CPLEX
    # our lazySep_ViolatedAcyclic function sets a LazyConstraintCallback of CPLEX
    # Entiere
    MOI.set(m, MOI.LazyConstraintCallback(), lazySep_ViolatedMengerCut) 
    
    # our userSep_ViolatedAcyclic function sets a LazyConstraintCallback of CPLEX  
    # Version Exacte
    MOI.set(m, MOI.UserCutCallback(), userSep_ViolatedMengerCut)
    
    # # our primal heuristic to "round up" a primal fractional solution
    # Pas vraiment une heuristique
    MOI.set(m, MOI.HeuristicCallback(), primalHeuristicTSP)
  #
  #################

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

    tmpmat = zeros(Int64,G.nb_points, G.nb_points)

    # for i in 1:G.nb_points
    #   for j in 1:G.nb_points
    #       tmpmat[i,j]=value(y[i, j])
    #   end
    # end
    # println("MATRICE YYYYY")
    # for i in 1:G.nb_points
    #   println(tmpmat[i,:])
    # end
    

    # Recuperation des points lies a chaque station
    Liens = []
    # On parcours chaque station
    for j in 1:size(Stations)[1]
      median = Stations[j]
      # initialisation avec la station 
      Tmp = Int64[]
      push!(Tmp, value(Stations[j]))
      for i in 1:G.nb_points
        # j un median et i un point affecté, i != j, on utilise directement les stations pour aller plus vite
        if (value(y[i, value(median)]) > 0.999 && i!=median)
          push!(Tmp,i)
        end
      end
      # On garde les points lies a une station
      push!(Liens,Tmp)
    end

    for i in 1:size(Liens)[1]
      println(Liens[i,:])
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
    S= find_cycle_in_integer_x(x, 1)
    push!(S,first(S))
    #println("TOUVEEEE ",S)
        # y[1, 1] = 1 donc Stations[1] contient ce point
        # i=1
        # j=2
        # while (value(x[Stations[i], Stations[j]]) < 0.999) 
     	  # j=j+1
        
        # end
        # # Numero de Station
        # push!(S,Stations[i])
        # push!(S,Stations[j])
        # # On cherche depuis la nouvelle Station
        # i=j
        # # On ne revient pas en arriere (exemple : on veut eviter 1 - 19 - 1)
        # tmp = 1
        # while (Stations[i]!=Stations[1])
        #   j=1
        #   # Les conditions : j==i, j=tmp, si la station precedente = station qu'on verifie, s'il n y a pas d'arete entre les 2 stations
        #   while  ( j==i || j==tmp || (value(x[Stations[i],Stations[j]]) <0.999 && value(x[Stations[j],Stations[i]]) < 0.999 ) ) 
        #     j=j+1
        #   end
        #   #println("WHEREEEE ", j," and ",tmp, " and ", Stations[j])
        #   push!(S,Stations[j])
        #   # On ne revient pas en arriere, avant le changement de i, on met un taboo sur le retour vers i
        #   tmp = i
        #   i=j
		    # end
		 println("Temps de résolution :", solve_time(m))
		 return S, Stations, Liens
	else
		 println("Problème lors de la résolution")
     println("Status ",status)
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
  alld = calcul_dist(I)

  # PARTIE CALCUL DES CRITERES
  # Les liens entre stations, l'anneau et les distances
  cost = calc_cost(Liens,S_STAR, alld)
  println("Le coût de la solution (critere 1): ",round(cost;digits=2))
	# println("Valeur: ",val_STAR)
  mean, ratioMarche, ratioMetro = calc_meanTime(Liens,S_STAR,alld,I)
  println("Le temps moyen de la solution (critere 2): ",round(mean;digits=2))
  println("Les ratios de marche et metro (critère 3): Marche - ",round(ratioMarche;digits=3), " et Metro - ",round(ratioMetro;digits=3))
  println()

  filename_STAR = replace(filename, ".tsp" => "_STAR")

  WritePdf_visualization_solution_projet(I,S_STAR,Liens,filename_STAR)
end

input = "../Instances_TSP/berlin52.tsp"
solve(input)