import math as m
from collections import defaultdict
 def order_calc(dc,last_order):
    
        
        self = dc
        
        out_bound = self.out_bound
        counter = self.counter
        order = self.order

        
        
        

        
        counter_next = 0
        while order != last_order:
            order += 1
            print(order)
            while counter != 0 : 
#                 print(out_bound)
                for s in self.findRegionGivenSegment(out_bound.twin) :
                    
                    s.order = order
                    
                    if  s.twin.order == 0 :
                        if counter_next  != 0 :
                            
                            out_bound_next.next_bound  = s
                        else :
                            start = s
                        
                        out_bound_next = s
                        counter_next += 1
                counter -= 1
                out_bound = out_bound.next_bound
                    
            out_bound_next.next_bound = start 
            
            out_bound = out_bound_next
            counter = counter_next
            counter_next = 0 
            
        self.out_bound = out_bound
        self.counter = counter
        self.order = order
        
# This class represents a directed graph
# using adjacency list representation
class Graph:
 # Constructor
    def __init__(self):
 # default dictionary to store graph
        self.graph = defaultdict(list)
    # function to add an edge to graph
 
    
    def addEdge(self,u,v):
        self.graph[u].append(v)
 # Function to print a BFS of graph
    def BFS(self, s):
 
        # Mark all the vertices as not visited
        visited = [False] * (max(self.graph) + 1)
        level = [0] * (max(self.graph) + 1)
        # Create a queue for BFS
        queue = []
 
        # Mark the source node as
        # visited and enqueue it
        queue.append(s)
        visited[s] = True
        level[s] = 0
 
        while queue:
 
            # Dequeue a vertex from
            # queue and print it
            s = queue.pop(0)
            print (s, end = " ")
 
            # Get all adjacent vertices of the
            # dequeued vertex s. If a adjacent
            # has not been visited, then mark it
            # visited and enqueue it
            for i in self.graph[s]:
                if visited[i] == False:
                    queue.append(i)
                    visited[i] = True
                    level[i] =  level[s] +1
                    
        return level
# Utils
def findHAngle(dx, dy):
  """Determines the angle with respect to the x axis of a segment
  of coordinates dx and dy
  """
  l = m.sqrt(dx*dx + dy*dy)
  if dy > 0:
    return m.acos(dx/l)
  else:
    return 2*m.pi - m.acos(dx/l)



class Vertex:
  def __init__(self, x, y):
    self.x = x
    self.y = y
    self.hedges = []  # list of halfedges whose tail is this vertex

  def __eq__(self, other):
    if isinstance(other, Vertex):
      return self.x == other.x and self.y == other.y
    return NotImplemented

  def sortHedges(self):
    self.hedges.sort(key=lambda a: a.angle, reverse=True)

  def __repr__(self):
    return "({0},{1})".format(self.x, self.y)


class Hedge:
  # v1 -> v2
  def __init__(self, v1, v2):
    self.prev = None
    self.twin = None
    self.next = None
    self.tail = v1
    self.face = None
    self.angle = findHAngle(v2.x-v1.x, v2.y-v1.y)
    self.order = 0
    self.next_bound = None

  def __eq__(self, other):
    return self.tail == other.tail and \
        self.next.tail == other.next.tail

  def __repr__(self):
    if self.next is not None:
      return "({0},{1})->({2},{3})".format(self.tail.x, self.tail.y,
                                           self.next.tail.x,
                                           self.next.tail.y)
    else:
      return "({0},{1})->()".format(self.tail.x, self.tail.y)


class Face:
  def __init__(self):
    self.halfEdge = None
    self.name = None
    self.order = 0
    self.angles = []
    self.area = 0
    self.number_of_edges = 0
    self.neighbours = []
    self.id = 0

class DCEL:
  def __init__(self):
    self.vertices = []
    self.hedges = []
    self.faces = []
    self.out_bound = None
    self.order = 1
    self.counter = 0
    self.number_of_faces = 0

  # Returns vertex object given x and y
  def findVertex(self, x, y):
    for v in self.vertices:
      if v.x == x and v.y == y:
        return v
    return None

  # Returns Halfedge whole vertices are v1 and v2
  # v1 and v2 are tuples
  def findHalfEdge(self, v1, v2):
    for halfEdge in self.hedges:
      nextEdge = halfEdge.next
      if (halfEdge.tail.x == v1[0] and halfEdge.tail.y == v1[1]) and (nextEdge.tail.x == v2[0] and nextEdge.tail.y == v2[1]):
        return halfEdge

    return None

  def findRegionGivenSegment(self, segment):
    segments = []
#     # We need to find the half edge whose vertices
#     # are that of the passed segment
#     v1 = segment[0]
#     v2 = segment[1]
    startEdge = segment

    h = startEdge
    while (not h.next == startEdge):
#       print(h, end="--->")
      segments.append(h)  
      h = h.next
#     print(h, '--->', startEdge)

    segments.append(h)
    return segments

  def build_dcel(self, points, segments):

    #  For each point create a vertex and add it to vertices
    for point in points:
      self.vertices.append(Vertex(point[0], point[1]))

    # For each input segment, create to hedges and assign their
    # tail vertices and twins

    # Structures of segment is [(0, 5), (2, 5)]
    for segment in segments:
      startVertex = segment[0]
      endVertex = segment[1]

      v1 = self.findVertex(startVertex[0], startVertex[1])
      v2 = self.findVertex(endVertex[0], endVertex[1])

      h1 = Hedge(v1, v2)
      h2 = Hedge(v2, v1)

      h1.twin = h2
      h2.twin = h1

      v1.hedges.append(h1)
      v2.hedges.append(h2)

      self.hedges.append(h1)
      self.hedges.append(h2)

    # For each endpoint, sort the half-edges whose
    # tail vertex is that endpoint in clockwise order.

    for vertex in self.vertices:
      vertex.sortHedges()

      noOfHalfEdges = len(vertex.hedges)

      if noOfHalfEdges < 2:
        return Exception("Invalid DCEL. There should be at least two half edges for a vertex")

      # For every pair of half-edges e1, e2 in clockwise order,
      # assign e1->twin->next = e2 and e2->prev = e1->twin.
      for i in range(noOfHalfEdges - 1):
        e1 = vertex.hedges[i]
        e2 = vertex.hedges[i+1]

        e1.twin.next = e2
        e2.prev = e1.twin

      # for the last and first halfedges pair
      e1 = vertex.hedges[noOfHalfEdges - 1]
      e2 = vertex.hedges[0]

      e1.twin.next = e2
      e2.prev = e1.twin

    # For every cycle, allocate and assign a face structure.
    faceCount = 0
    for halfEdge in self.hedges:

      if halfEdge.face == None:
#         print('here')
        faceCount += 1

        f = Face()
        f.name = "f" + str(faceCount)

        f.halfEdge = halfEdge
        halfEdge.face = f

        h = halfEdge
        while (not h.next == halfEdge):
          h.face = f
          h = h.next
        h.face = f

        self.faces.append(f)
        f.id = len(self.faces )-1
        
        
    a =  self.findHalfEdge([0.5,0.5],[-0.5,0.5])
    b =  self.findHalfEdge([-0.5,0.5],[-0.5,-0.5])
    c =  self.findHalfEdge([-0.5,-0.5],[0.5,-0.5])
    d =  self.findHalfEdge([0.5,-0.5],[0.5,0.5])
        
    a.next_bound = b
    b.next_bound = c
    c.next_bound = d
    d.next_bound = a
        
    a.order = 1
    b.order = 1
    c.order = 1
    d.order = 1
    self.counter = 4 
    self.out_bound = a

  





    ####################

def __main__():
  
    file = open('points.txt', 'r') #read the file
    points = file.read()
        
    points_str = points.split(" ")
        
        
        
    points = []
    for i in range(int(len(points_str)/3)):
         points.append((float(points_str[3*i+1]), float(points_str[3*i+2]) ))
        
    file.close() #close file
            
    file = open('edges.txt', 'r') #read file
    edges = file.read()
    edges_str = edges.split(" ")
    edges = []
    for i in range(int(len(edges_str)/2)):
         edges.append((int(edges_str[2*i]), int( edges_str[2*i+1] )))
    file.close() #close file
    
   


   
    
#     myDCEL.findRegionGivenSegment(starting_segment)
     



    segments = []
    for edge in edges :
        if edge[0] < edge[1] :
            segments.append([points[edge[0]],points[edge[1]]])



    myDCEL = DCEL()
    myDCEL.build_dcel(points, segments)



    g=Graph()
    
    for hedge in myDCEL.hedges :
        hedge.face.neighbour = hedge.twin.face
        g.addEdge(hedge.face.id,hedge.twin.face.id)
    
    level = g.BFS( myDCEL.findHalfEdge([0.5,0.5],[-0.5,0.5]).face.id )
    
    for face in myDCEL.faces :
        face.order = level[face.id]




    dc = myDCEL
    
    array = []
    for i in range(300) :
        array.append([])
    


    
    for face in dc.faces :
        
        array[face.order+1].append(face)
    
        
    order = 60
    
    # number of regions for every k :
    
    
    array[0] = []
    number_of_regions = []
    b_zone =[ ]
    for i in range(order) :
        b_zone.append(i)
    
        number_of_regions.append(len(array[i]))
    import pandas as pd 
    dataframe = pd.DataFrame()
    dataframe['b_zone'] = b_zone
    dataframe['number_of_regions'] = number_of_regions
    dataframe['karan'] = 6*dataframe.b_zone - 6

# plt.figure(figsize=(9, 3*1.6))
# plt.plot(dataframe.loc[2:].b_zone, dataframe.loc[2:].number_of_regions)
# plt.plot(dataframe.loc[2:].b_zone, dataframe.loc[2:].karan)

# print(dataframe.number_of_regions)

# plt.plot(dataframe.b_zone, dataframe.number_of_regions.cumsum())

    areas = []
    for i in range(order) :
        for face in array[i] :
            import numpy as np
            x = []
            y = []
            hedges = dc.findRegionGivenSegment(face.halfEdge)
            for he in hedges :
                x.append(he.tail.x)
                y.append(he.tail.y)
            x = np.array(x)
            y = np.array(y)
    
            face.area =  0.5*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1))) 
            areas.append(face.area)

    max_arr = []
    for i in range(order) :
        maximum = 0
        for face in array[i] :
            if face.area > maximum :
                maximum = face.area
        max_arr.append(maximum)
        
    dataframe['maximum_area'] =  max_arr


    
    df = pd.DataFrame()
    df['areas'] = areas




    dataframe['bound'] = 1/ np.sqrt(dataframe.b_zone)


# plt.plot(dataframe.loc[2:].b_zone, dataframe.loc[2:].bound)



    dataframe['bound'] = 1/ np.sqrt(dataframe.b_zone)







