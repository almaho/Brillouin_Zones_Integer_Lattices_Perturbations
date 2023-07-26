import numpy as np
import math
import sympy as sym
import random
import time
import arrangement.geometricTraits as trts 
import networkx as nx

class Node:
    def __init__ (self, selfIdx, point, traitIdx, traits):
        self.attributes = {}
        self.selfIdx = selfIdx   # self-index
        self.point = point
        self._traitIdx = traitIdx
        self._traitTval = []
        self.update_tval(traits)

    def transform_sequence(self, operTypes, operVals, operRefs, traits):


        # transforming the self.poin
        for opIdx, opType in enumerate(operTypes):

            if opType == 'T':# and all(operVals[opIdx]!=(0,0)):
                tx,ty = operVals[opIdx]
                self.point = self.point.translate(tx,ty)

            elif opType == 'R':# and operVals[opIdx]!=0:
                theta = operVals[opIdx]
                ref = operRefs[opIdx]
                self.point = self.point.rotate(theta,ref)

            elif opType == 'S':# and all(operVals[opIdx]!=(1,1)):
                sx,sy = operVals[opIdx]
                ref = operRefs[opIdx]
                self.point = self.point.scale(sx,sy,ref)

        # updating self._traitTval after transforming the self.point
        self.update_tval(traits)

    def update_tval (self, traits):
        '''
        Node class
        '''
        self._traitTval = [traits[cIdx].IPE(self.point) for cIdx in self._traitIdx]
class HalfEdge:
    def __init__ (self,
                  selfIdx, twinIdx,
                  traitIdx, direction):
        '''
        '''
        self.attributes = {}
        self.selfIdx = selfIdx   # self-index (startNodeIdx, endNodeIdx, pathIdx)
        self.twinIdx = twinIdx   # index to the twin half-edge
        self.succIdx = None      # index to the successor half-edge

        # half edge Trait's attributes:
        self.traitIdx = traitIdx    # Index of the trait creating the edge
        self.direction = direction  # defines the direction of t-value (positive: t2>t1, negative: t1>t2)


    def get_tvals(self, traits, nodes):

        (s,e,_) = self.selfIdx # (s,e,k) = self.selfIdx

        sTVal = nodes[s]['obj']._traitTval[nodes[s]['obj']._traitIdx.index(self.traitIdx)]
        eTVal = nodes[e]['obj']._traitTval[nodes[e]['obj']._traitIdx.index(self.traitIdx)]

        # for the case of circles and the half-edge that crosses the theta=0=2pi
        if (self.direction=='positive') and not(sTVal < eTVal):
            eTVal += 2*np.pi
        if (self.direction=='negative') and not(sTVal > eTVal):
            sTVal += 2*np.pi

        return (sTVal, eTVal)

class Arrangement:
    def __init__ (self, traits, point_file, edge_file, flag_file,flag, R):

        self.traits = traits
        self.graph = nx.MultiDiGraph()


        # STAGE A: construct nodes
        tic = time.time()
        self._construct_nodes(point_file, flag_file,flag, R)
        
        #### STAGE B: construct edges
        tic = time.time()
        self._construct_edges()
            
        out = ""
        for halfEdgeIdx in self.graph.edges(keys=True):
            (s,e,k) = (startNodeIdx, endNodeIdx, path) = halfEdgeIdx
#             print ( (s,e,k), ': ', self.graph[s][e][k]['obj'].attributes )
            out = out + str(s) +' ' +  str(e) + ' '
            
        with open(edge_file, 'w') as f:
            f.write(out)
            
    def _construct_nodes(self, point_file,flag_file,flag, R):
        from fractions import Fraction as frac
        
        def equal_f (x1,x2, y1,y2) :
            if (x1.numerator == x2.numerator) and (y1.numerator == y2.numerator) :
                if (x1.denominator == x2.denominator) and (y1.denominator == y2.denominator) :
                    return True
            return False
                
#         flag_dists = []
        print('start')
        intersection = [ [ 0
                                for col in range(len(self.traits)) ]
                              for row in range(len(self.traits)) ]


        print('step 1')

        dicts = {}
        intersectionsFlat = []
        import time
        begin = time.time()
        for row, trait_row in enumerate(self.traits):
        
            print(row)
            for col in range(row):
                
                if intersection[col][row] == 0 :
                    obj1 = trait_row.obj
                    obj2 = self.traits[col].obj
                    if len( sym.intersection(obj1,obj2) ) > 0 : 

                        ip_tmp = sym.intersection(obj1,obj2)[0]
                        
                        intersection[col][row] = ip_tmp 
                        intersection[row][col] = ip_tmp 
                        fx2 = frac(ip_tmp.x)
                        fy2 = frac(ip_tmp.y)
#                         x2 = np.float(ip_tmp.x)
#                         y2 = np.float(ip_tmp.y)
#                         if fx2**2 + fy2**2 > R**2 :
#                             continue
                        intersectionsFlat.append(ip_tmp)
                        
                        dicts[ip_tmp] = [row,col]
                        
                        n= len(self.traits)
                        for ind in range(col+1,n) :
                            
                            if ind == row :
                                continue
                            obj2 = self.traits[ind].obj
                            if len( sym.intersection(obj1,obj2) ) == 0 :
                                continue
                            from fractions import Fraction
                            fx1 = frac(sym.intersection(obj1,obj2)[0].x) 
                            fy1 = frac(sym.intersection(obj1,obj2)[0].y) 

                            
#                             if fx1**2 + fy1**2 > R**2 :
#                                 continue
#                             dist = np.sqrt((x1-x2)**2 + (y1-y2)**2)
#                             if dist < 10*flag :
#                                     flag_dists.append(dist)
#                             if  dist < np.spacing(10**10):
#                                 print(dist)
                            if equal_f(fx1,fx2,fy1,fy2) :
                                print(fx1.numerator)
                                print(fx2.numerator)
                                for line in dicts[ip_tmp] :
                                     
                                    intersection[line][ind] = ip_tmp
                                    intersection[ind][line] = ip_tmp

                                    
                                
                                    
                                dicts[ip_tmp].append(ind)
                                    
        ipsTraitIdx =  []
        for ind in range(len(intersectionsFlat )) :
            ipsTraitIdx.append(dicts[intersectionsFlat[ind]])
        if not (len(intersectionsFlat) == len(ipsTraitIdx)): raise AssertionError()

 
        nodes = [ [ pIdx, {'obj': Node(pIdx, intersectionsFlat[pIdx], ipsTraitIdx[pIdx], self.traits)} ]
                      for pIdx in range(len(intersectionsFlat)) ]
#         print(nodes)
        out = ""
        for i in range(len(intersectionsFlat)) :
            out = out+ str(i) + ' ' + str(np.float(intersectionsFlat[i].x)) + ' ' +str(np.float(intersectionsFlat[i].y)) + ' ' 
            
        with open(point_file, 'w') as f:
            f.write(out)
            
        self.graph.add_nodes_from( nodes )

        if not (len(self.graph.nodes()) == len(intersectionsFlat)): raise AssertionError()

    
    def _construct_edges(self):
        
        ########################################
        # step 1_a: find intersection points of traits
        # indeces of intersection points corresponding to each traits
        traitIpsIdx = [[] for i in range(len(self.traits))]
        traitIpsTVal = [[] for i in range(len(self.traits))]

        for nodeIdx in self.graph.nodes():
            for (tVal,cIdx) in zip(self.graph.nodes[nodeIdx]['obj']._traitTval,
                                   self.graph.nodes[nodeIdx]['obj']._traitIdx) :
                traitIpsIdx[cIdx].append(nodeIdx)
                traitIpsTVal[cIdx].append(tVal)


        for cIdx in range(len(self.traits)):
            # in python2 zip returns a list, in pothon3 it returns a zip() object
            tmp = sorted( list( zip( traitIpsTVal[cIdx], traitIpsIdx[cIdx] )))
            traitIpsIdx[cIdx] = [pIdx for (tVal,pIdx) in tmp]
            traitIpsTVal[cIdx].sort()


        # ########################################
        # step 3: half-edge construction
        for (cIdx,trait) in enumerate(self.traits):

            ipsIdx = traitIpsIdx[cIdx]
            tvals = traitIpsTVal[cIdx]

            # step a:
            # for each trait, create all edges (half-edges) located on it
            if isinstance(trait.obj, ( sym.Line, sym.Segment, sym.Ray) ):

                startIdxList = ipsIdx[:-1]
                startTValList = tvals[:-1]

                endIdxList = ipsIdx[1:]
                endTValList = tvals[1:]

            

            

            # in python2 zip returns a list, in pothon3 it returns a zip() object
            l = list( zip (startIdxList, startTValList, endIdxList, endTValList) )

            # create a half-edge for each pair of start-end point
            for ( sIdx,_, eIdx,eTVal ) in l: # (sIdx,sTVal, eIdx,eTVal) in l

                newPathKey1 = len(self.graph[sIdx][eIdx]) if eIdx in self.graph[sIdx].keys() else 0
                newPathKey2 = len(self.graph[eIdx][sIdx]) if sIdx in self.graph[eIdx].keys() else 0

                # in cases where sIdx==eIdx, twins will share the same key ==0
                # this will happen if there is only one node on a circle
                # also a non-intersecting circles with one dummy node
                # next line will take care of that only
                if sIdx==eIdx: newPathKey2 += 1

                idx1 = (sIdx, eIdx, newPathKey1)
                idx2 = (eIdx, sIdx, newPathKey2)

                # Halfedge(selfIdx, twinIdx, cIdx, side)

                # first half-edge
                direction = 'positive'
                he1 = HalfEdge(idx1, idx2, cIdx, direction)
                e1 = ( sIdx, eIdx, {'obj':he1} )

                # second half-edge
                direction = 'negative'
                he2 = HalfEdge(idx2, idx1, cIdx, direction)
                e2 = ( eIdx, sIdx, {'obj': he2} )

                self.graph.add_edges_from([e1, e2])
            

def load_perturbed_data (length, xs, ys) :
    def calc_line(x,y) :
            return [x/2 , y/2 , x/2 - y/2 , x/2 + y/2]
    Point = sym.Point
    Line = trts.LineModified 
    result = {}
    traits = []
    
    all_lattice_points = []
    traits = []
    Point = sym.Point
    Line = trts.LineModified
    
    data = {'lines':[] }
    points = []
    for i in range(2*length+1) :
        for j in range(2*length+1) :
            if (i-length == 0 ) and (j-length == 0) :
                point = [0,0]
            else :
                
                x = xs[i][j]
                y = ys[i][j]
                point = [((i-length)*10000+x),((j-length)*10000+y)] 
                all_lattice_points.append(point)
                l = calc_line(point[0], point[1])
                traits += [ Line( args=(Point(l[0],l[1]), Point(l[2],l[3]))) ]
                

    return traits


def calc_flag(order) :
    return  ((np.pi / (16*order + 16*np.sqrt(2*order*np.pi ) + 8*np.pi) )**2)*10000
    
def calc_R(order) :
    return (np.sqrt(order/np.pi) + np.sqrt(2)/2)*10000
    
    
if __name__ == '__main__':
    
    length = 3
    order = 10
    
    flag  = calc_flag(order) 
    R = calc_R(order)
#     print(flag)
    
#     print ('weak perturbation')
    
    
#     perturbation = 10
    
    xs = np.zeros((2*length+1, 2*length+1))
    ys = np.zeros((2*length+1, 2*length+1))
#     for i in range (2*length+1) :
#         for j in range(2*length+1) :
#             xs[i][j] = random.randint(-perturbation,perturbation)
#             ys[i][j] = random.randint(-perturbation,perturbation)
    
#     traits10 = load_perturbed_data(length,xs, ys)
#     Arrangement(traits10, 'points-run7-weak.txt', 'edges-run7-weak.txt','flags_run7-weak.txt',flag, R)
    
#     print ('medium perturbation')
    
    
#     perturbation = 100
    
   
#     for i in range (2*length+1) :
#         for j in range(2*length+1) :
#             xs[i][j] *= 10
#             ys[i][j] *= 10
    
#     traits100 = load_perturbed_data(length,xs, ys)
#     Arrangement(traits100, 'points-run7-medium.txt', 'edges-run7-medium.txt','flags_run7-medium.txt',flag, R)
    
#     print ('strong perturbation')
#     for i in range (2*length+1) :
#         for j in range(2*length+1) :
#             xs[i][j] = xs[i][j]*10
#             ys[i][j] = ys[i][j]*10
            
#     traits1000 = load_perturbed_data(length,xs, ys) 
#     Arrangement(traits1000, 'points-run7-strong.txt', 'edges-run7-strong.txt', 'flags_run7_strong.txt',flag, R)
    
    print('lattice')
    for i in range (2*length+1) :
        for j in range(2*length+1) :
            xs[i][j] = 0
            ys[i][j] = 0
    traits = load_perturbed_data(length,xs, ys) 
    Arrangement(traits, 'points-run8-lattice.txt', 'edges-run8-lattice.txt', 'flags_run7_lattice.txt',flag, R)
