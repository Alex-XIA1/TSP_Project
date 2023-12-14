# Code largement inspiré du travail de Florian Belhassen-Dubois
using Plots

# Sructure de données "Nuage de points"
struct NuagePoints
	nb_points    # nombre de points
	X            # tableau des coordonnées x des points
	Y            # tableau des coordonnées y des points    
	max_distance # plus longue distance entre deux points
end

# Lecture d'un fichier de la TSPLIB .tsp (atttention uniquement les instances symétriques données par des coordonnées)
# Renvoie l'instance du TSP c
function Read_undirected_TSP(filename)

       I = NuagePoints(0, [], [],0)

	open(filename) do f
			
		node_coord_section = 0 # repère dans quelle section du fichier .tsp on est en train de lire
		nbp = 0
		X = Array{Float64}(undef, 0)
		Y = Array{Float64}(undef, 0)
			
		# on lit une par une nos ligne du fichier	
		for (i,line) in enumerate(eachline(f))
			
			# on sépare cette ligne en mots
			x = split(line," ") # For each line of the file, splitted using space as separator
             
			# on supprime les mots vides, en effet une espace suivi d'un autre espace renvoie le mot vide
			deleteat!(x, findall(e->e=="", x))

			
			if(node_coord_section == 0)       # If it's not a node section
					
				# si on est dans les coordonnées on va s'arrêter et remplir notre instance sinon il nous reste des labels à lire
				if(x[1] == "NODE_COORD_SECTION")
					node_coord_section = 1
				# si on est dans le label dimension, on le récupère
				elseif(x[1] == "DIMENSION")
					nbp = parse(Int, x[3])
				end
			
			# on est enfin dans nos coordonnées ! On les lit et on remplit notre instance avec
			elseif(node_coord_section == 1 && x[1] != "EOF")
				
				push!(X, parse(Float64, x[2]))
				push!(Y, parse(Float64, x[3]))
				
			else
				
				node_coord_section = 2
				
			end
		end
		
		
		# Calcule la plus longue distance entre deux points
		max_distance=0
	       for i in 1:nbp
	          for j in 1:nbp
	             if (max_distance < ( (X[i] - X[j])^2 + (Y[i] - Y[j])^2 ) )
	                max_distance = (X[i] - X[j])^2 + (Y[i] - Y[j])^2
	             end
	          end
	       end
	       
      		# on construit notre nuage de points
		I = NuagePoints(nbp, X, Y, max_distance)

		
	end
	return I
end



# Visualisation d'une instance comme un nuage de points dans un fichier pdf dont le nom est passé en paramètre
# On peut visualiser les données du PROJET ici
function WritePdf_visualization_TSP(I, filename)

	filename_splitted_in_two_parts = split(filename,".") # split to remove the file extension
	filename_with_pdf_as_extension= filename_splitted_in_two_parts[1]*".pdf"
	# save to pdf
	
	# un simple plot en fonction des coordonnées 
	p = plot(I.X, I.Y, seriestype = :scatter)
	savefig(p, filename_with_pdf_as_extension)

end




# Renvoie la distance euclidienne entre deux points du nuage
function dist(I, i, j) 
	
	return ((I.X[i] - I.X[j])^2 + (I.Y[i] - I.Y[j])^2)^(0.5)
	
end

# Crée une matrice de toutes les distances point à point
function calcul_dist(I)

	c = Array{Float64}(undef, (I.nb_points, I.nb_points))
	
	for i in 1:I.nb_points
	
		for j in 1:I.nb_points
		
			c[i, j] = dist(I, i, j)
		
		end
	
	end
	
	return c
	
end

 
# calcule la somme des coûts de notre arête solution
function Compute_value_TSP(I, S)
	
	res = ((I.X[S[1]] - I.X[S[I.nb_points]])^2 + (I.Y[S[1]] - I.Y[S[I.nb_points]])^2)^(0.5)
	for i = 1:(I.nb_points -1)	
		res = res + ((I.X[S[i]] - I.X[S[i+1]])^2 + (I.Y[S[i]] - I.Y[S[i+1]])^2)^(0.5)
	end
	
	return res

end


# permet de visualiser notre solution (un circuit / cycle) dans un fichier pdf dont le nom est spécifié en paramètres
# La solution est donnée par la liste ordonné des points à visiter commençant par 1
function WritePdf_visualization_solution_ordre(I, S, filename)

	filename_splitted_in_two_parts = split(filename,".") # split to remove the file extension
	filename_with_pdf_as_extension= filename_splitted_in_two_parts[1]*".pdf"
	# save to pdf
	
	tabX= Float16[]
	tabY= Float16[]
	
    for i in S
       push!(tabX, I.X[i])
       push!(tabY, I.Y[i])
    end
	
	# on ajoute le premier point pour le plot, c'est important sinon il manque l'arête entre 1 et n...
	push!(tabX, I.X[1])
	push!(tabY, I.Y[1])
	
	p = plot(I.X, I.Y, seriestype = :scatter,legend = false)
	plot!(p, tabX, tabY,legend = false)
	savefig(p, filename_with_pdf_as_extension)

end

# Visualisation de la solution pour le PROJET MAOA
# entree : I graphe, S solution, Liens (les zones vers stations), filename, nom du fichier en sortie
function WritePdf_visualization_solution_projet(I, S, Liens, filename)

	filename_splitted_in_two_parts = split(filename,".") # split to remove the file extension
	filename_with_pdf_as_extension= filename_splitted_in_two_parts[1]*".pdf"
	# save to pdf
	
	# tableau des liens entre solutions
	tabX= Float16[]
	tabY= Float16[]
	
    for i in S
       push!(tabX, I.X[i])
       push!(tabY, I.Y[i])
    end
	
	# on ajoute le premier point pour le plot, c'est important sinon il manque l'arête entre 1 et n...
	push!(tabX, I.X[1])
	push!(tabY, I.Y[1])
	
	p = plot(I.X, I.Y, seriestype = :scatter,legend = false)
	plot!(p, tabX, tabY,legend = false)

	tab = []
	for i in 1:size(Liens)[1]
		l = Liens[i]
		# La station est mise en premier
		xs, ys = I.X[Liens[i][1]], I.Y[Liens[i][1]]
		for j in 2:size(l)[1]
			# Pour chaque zone en lien avec la station
			x, y = I.X[Liens[i][j]], I.Y[Liens[i][j]]
			# On trace un trait entre la station et la zone
			plot!(p,[xs, x], [ys, y], line=:solid, color=:gray, legend = false)
		end
	end
	savefig(p, filename_with_pdf_as_extension)

end


function read_Clustering(filename)
	Liens = []
	open(filename) do f
		for (i,line) in enumerate(eachline(f))
			
			# on sépare cette ligne en mots
			x = split(line," ") # For each line of the file, splitted using space as separator
			# on supprime les mots vides, en effet une espace suivi d'un autre espace renvoie le mot vide
			deleteat!(x, findall(e->e=="", x))

			# Un cluster avec ses liens
			tmp = Array{Int}(undef, 0)
			for j in 2:size(x,1)
				push!(tmp,parse(Int,x[j]))
			end
			push!(Liens,tmp)
		end
	end
	return Liens
end

# println(read_Clustering("../Clustering/burma14.txt"))


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
	#source = source + 1
	source = mod(source, size(ring,1)-1) + 1
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