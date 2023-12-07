from sklearn.cluster import KMeans
import numpy as np


nCluster = 5
fname = "../Instances_TSP/berlin52.tsp"
points = []

f = open(fname, "r")
x = ""
while x != "NODE_COORD_SECTION\n":
    x = f.readline()
    x = x.split(" ")
    x = x[0]
    print(x)
    


print("OK ", f.readline())