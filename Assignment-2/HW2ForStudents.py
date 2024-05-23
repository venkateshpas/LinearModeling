import HW2functions_5963540 as HW2f
import numpy as np
import math
# input
E = 70e3 # E-modulus
BC = np.array([[2,5,7,9,10],[0,0,0,0,0]]) #boundary conditions; form: [degrees of freedom], [applied displacements]
F = np.array([0, 0 ,4250, 0, -3000, -5000, 0,0,0,0]) # force vector
Con = np.array([[1,2],[2,3],[3,4],[4,5],[1,5],[1,4],[2,4]]) # connectivity matrix; form: [[element 1], [element 2],etc]
NodePos = np.array([[-500, 0], [-350, 150], [0, 200],[0,0],[0,-100]]) # node position; form: [[element 1],[element 2],etc]
Area = np.array([20, 20, 30, 30, 20, 30, 30]) # vector with the area for each element
# functions here

u, strain, stress, reactions = HW2f.calcUSSR(E,BC,F,Con,NodePos,Area)

#output: do not change!
for el in range(len(Con)):
	print ("the stress in element ", (el+1) , " equals ", "%.2f" % stress[el], "MPa")

for node in range(len(NodePos)):
	print("The displacement of node", node+1 , " in x-direction is ", "%.3f" % u[2*(node+1)-2], "mm")
	print("The displacement of node", node+1 , " in y-direction is ", "%.3f" % u[2*(node+1)-1], "mm")
	print("The reaction force at node", node+1 , " in x-direction is ", "%.0f" % reactions[2*(node+1)-2], "N")
	print("The reaction force at node", node+1 , " in y-direction is ", "%.0f" % reactions[2*(node+1)-1], "N")
