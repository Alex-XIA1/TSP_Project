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


