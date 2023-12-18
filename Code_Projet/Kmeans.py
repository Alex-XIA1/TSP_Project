from sklearn.cluster import KMeans
import numpy as np
from sklearn.metrics import pairwise_distances_argmin_min
from scipy.cluster.vq import vq
import matplotlib.pyplot as plt

def read(): return list(map(float, readline().split()))

def kmeans(file_input, file_output, p):
    M = []
    with open(file_input, 'r') as f:
        for line in f:
            try: 
                line = list(map(float, line.split()))
                if line: M.append(line)
            except: continue

    M = np.array(M)
    kmeans = KMeans(n_clusters=p, n_init=1).fit(M[:,1:])
    affect = kmeans.labels_
    centers = np.array(kmeans.cluster_centers_)

    # Les centres sont plus proches des points selon l'index
    closest, distances = vq(centers, M[:,1:])
    reste = [i for i in range(len(M)) if i not in closest]
    test = list(closest) + reste

    affect[closest] = -1
    
    # for i in range(p):
    #     print(i, closest[i], *np.where(affect==i)[0])
    
    with open(file_output, 'w') as f:
        for i in range(p):
            s = ' '.join(map(str,[i+1, closest[i]+1, *(np.where(affect==i)[0]+1)])) + '\n'
            f.write(s)

p = 6

kmeans("../Instances_TSP/pr107.tsp","../Clustering/"+str(p)+"pr107.txt",6)