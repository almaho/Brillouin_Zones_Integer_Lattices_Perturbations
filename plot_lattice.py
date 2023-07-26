import numpy as np

def read_file (point_file, edge_file) :

    file = open(point_file, 'r')
    points = file.read()
    points_str = points.split(" ")
    points = []
    for i in range(int(len(points_str)/3)):
         points.append([float(points_str[3*i+1]), float(points_str[3*i+2]) ])
    file.close() 

    file = open(edge_file, 'r') 
    edges = file.read()
    edges_str = edges.split(" ")
    edges = []
    for i in range(int(len(edges_str)/2)):
         edges.append([int(edges_str[2*i]), int( edges_str[2*i+1] )])
    file.close() 
    return points,edges
def calc_R(order) :
    return (np.sqrt(order/np.pi) + np.sqrt(2)/2) * 10000
    
R= calc_R (57)
points,edges = read_file('points-lattice.txt','edges-lattice.txt')

points = np.array(points)
edges = np.array(edges)

from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection
lc = LineCollection(points[edges])
fig = plt.figure(figsize=(50, 50))

plt.gca().add_collection(lc)
plt.xlim(-R, R)
plt.ylim(-R, R)
plt.plot(points[:,0], points[:,1], 'ro', )
fig.savefig('full_figure.png')
