# Brillouin_Zones_Integer_Lattices_Perturbations

We study the k-th Brillouin zone of a point in a locally finite set which is the region in which the point
is the k-closest in the set. If the set is a lattice, the k-th Brillouin zones of different points are
translations of each other, which tile space.
Here is the link to our paper: https://arxiv.org/abs/2204.01077

![photo1634583830](https://github.com/almaho/Brillouin_Zones_Integer_Lattices_Perturbations/assets/126069668/d45ae885-4f00-4ccb-b99d-158758921b38)


In our paper "Brillouin Zones of Integer Lattices and Their
Perturbations" We study characteristics of perturbed lattices as well. Here is an example:
![download (7)](https://github.com/almaho/Brillouin_Zones_Integer_Lattices_Perturbations/assets/126069668/b61e85d5-163f-4275-bd14-2d4ae448d747)

Here are some instructions for using the code:

Libraries:
arrangement: https://github.com/saeedghsh/arrangement 
This library is used for creating files containing coordinates for points created in a given line arrangemnet.
Using lattice.py we create the edges and points related to Brillouin Zones. They are created by drawing the perpendicular bisectors of any 2 points of the lattice.

The input to the file is the line arrangement created by the perpendicular bisectors. The size of the input lattice is 19x19, insuring the accuracy of orders at most 39 to 57 for Brillouin Zone of the center point in Perturbed lattices of weak to strong.
The output is two files, one with the coordinates of points and their indices. The other is the indices of any two points creating an edge.

Using dcel.py faces, the number of regions, areas, and angles are explored. We assign order values to each face created in the Brillouin zone.

The "dcel" class is used to create objects saving the features of faces created in Brillouin Zones. The function "calc_order", assigns the order number to each face.
Using this assignment we calculate the minimum and maximum areas, distances and number of regions of any order.
