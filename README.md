# Modèles et Applications en Ordonnancementet optimisation combinatoire
## Projet 2023-2024
### Autour du tracé d’un métro circulaire

On s’intéresse au problème du tracé d’une ligne de transport public circulaire (cela peut être un
bus, un tramway ou un RER: toute ressemblance avec les métros du grand Paris serait parfaitement
fortuite).


Pour ce projet, il est demandé:
- une visualisation des données et des solutions
- une heuristique et méta-heuristique
- une formulation compacte
- une formulation non-compacte
- de comparer les résultats obtenus pour les 3 méthodes: temps de calcul et durée d’exécution
- de tester les méthodes sur un lot suffisant d’instances de tailles croissantes: déterminant ainsi
les tailles maximales pouvant être résolues versus la qualité des solutions.
- un rendu sous forme d’un mini-rapport décrivant vos idées, vos expérimentations (tableaux,
courbes) avec une analyse critique
- un rendu des archives de votre code
- une mini-soutenance en janvier pour présenter vos idées et résultats (quelques slides issus de
votre rapport pendant 5 à 10 min).


Les méthodes qu'on utilise dans le projet sont :
- Heuristique et méta-heuristique -> Heuristique Kmeans + TSP (À la façon d'un TSP Généralisé), Métaheuristique voisinage 1-1 au sein d'un cluster.
- Formulation compacte -> Problème de l'anneau étoile
- Formulation non-compacte -> Problème de l'anneau étoile pour comparer avec la formulation compacte

Pour pouvoir faire le Kmeans clustering : installer sci-kitlearn et numpy. Puis pour lancer : python (ou python3) Kmeans.py (en changeant si besoin les paramètres du fichier)

Le dossier imgReport contient des affichages de solutions avec par exemple 6berlin52solutionNC.png qui signifie 6 stations, berlin52, solution PLNE non compact. Si on numéro n'est pas présent au début, ce sera par défaut 6 stations.

Le dossier Clustering contient des clusters déjà fait (pour les expériences)
Le dossier Instances_TSP contient les instances fournies de TSP
Le dossier Code_Projet contient le code du projet :
- TSP_IO.jl contient des codes pour le calcul des critères ainsi que de lecture de fichier (certaines fonctions données dans le TME). 
- TSP_Etoile.jl correspond à la formulation compact de l'anneau étoile. Le fichier TSP_Etoile_NC.jl correspond à la formulation non compact de l'anneau étoile. 
- mainExperiences.jl correspond au fichier avec les initialisations et expériences
- MetroMeta.jl correspond au code de la metaheuristique
- Kmeans.py correspond au code du clustering Kmeans
- TSP_NC.jl correspond au code fourni dans le TME du TSP.

Donc pour lancer une expérience il suffit de lancer le fichier mainExperiences.jl avec Julia en précisant le nom des instances dans le fichier (cluster pour la métaheuristique et l'instance pour tous).


