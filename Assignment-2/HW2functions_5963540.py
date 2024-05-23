import numpy as np
import math
from numpy.linalg import inv

def calcStrain(Nelem,length,u, Con,T):
	strain = np.zeros((Nelem,1))
	counter = 0
	for i in Con:
		ulocal = np.zeros((4,1))
		a = (i[0]-1) * 2
		b = (i[0]-1) * 2 + 1
		c = (i[1]-1) * 2
		d = (i[1]-1) * 2 + 1
		ulocal[0][0] = u[a][0]
		ulocal[1][0] = u[b][0]
		ulocal[2][0] = u[c][0]
		ulocal[3][0] = u[d][0]
		ulocal = np.dot(inv(T[counter]), ulocal)
		strain[counter] = (ulocal[2][0] - ulocal[0][0])/length[counter]
		counter += 1
	return strain 
	
def calcStress(strain,E):
	stress = E * strain
	return stress

def calcUSSR(E,BC,F,Con,NodePos,Area):
	Nelem = Area.shape[0]
	F = F[...,None]
	DOF = NodePos.shape[0] * 2
	theta,length = thetaLength(NodePos,Con)
	u = np.zeros((DOF,1))
	Kglobal,T = stiffnessMatrix(E,theta,Area,length,Con,DOF)
	Kglobal_BC, F_BC, u_BC = BoundaryConditionImplementation(BC,Kglobal,F,u)
	u_BC = np.dot(inv(Kglobal_BC), F_BC)
	node_values = BC[0,:]
	disp_values = BC[1,:]
	dof_no = 1
	count1 = 0
	count2 = 0
	while dof_no <= DOF:
		if dof_no in node_values:
			u[dof_no-1] = disp_values[count1]
			count1 +=1
		else:
			u[dof_no-1] = u_BC[count2]
			count2 += 1
		dof_no += 1
	reactions = np.zeros((DOF,1))
	reactions = F - np.dot(Kglobal,u)
	strain = calcStrain(Nelem, length, u, Con,T) # just an example of a function, input is also just to have an input
	stress = calcStress(strain,E) # just an example of a function, input is also just to have an input
	reactions = np.zeros((DOF,1))
	reactions = np.dot(Kglobal,u) - F  # just a placeholder to have a vector of the right size
	return u, strain, stress, reactions

def thetaLength(NodePos,Con):
	theta = []
	length = []
	for element in Con:
		deltaX = NodePos[element[1]-1] [0] - NodePos[element[0]-1][0]
		deltaY = NodePos[element[1]-1] [1] - NodePos[element[0]-1][1]
		length.append(math.sqrt(deltaX**2 + deltaY ** 2))
		theta.append(math.atan2(deltaY,deltaX))
	return theta,length

def stiffnessMatrix(E,theta,Area,length,Con,DOF):
	Nelem = Area.shape[0]
	Kglobal = np.zeros((DOF,DOF))
	Klocal = [np.zeros((4, 4)) for _ in range(Nelem)]
	Krotated = []
	T = [np.zeros((4, 4)) for _ in range(Nelem)]
	j = 0
	while j < Nelem:
		Klocal[j][0][0] = E * Area[j]/length[j]
		Klocal[j][0][2] = -E * Area[j]/length[j]
		Klocal[j][2][0] = -E * Area[j]/length[j]
		Klocal[j][2][2] = E * Area[j]/length[j]
		T[j][0][0] = math.cos((theta[j]))
		T[j][1][1] = math.cos((theta[j]))
		T[j][2][2] = math.cos((theta[j]))
		T[j][3][3] = math.cos((theta[j]))
		T[j][0][1] = -math.sin((theta[j]))
		T[j][2][3] = -math.sin((theta[j]))
		T[j][1][0] = math.sin((theta[j]))
		T[j][3][2] = math.sin((theta[j]))
		j = j + 1
	j = 0
	while j < Nelem:
		Krotated.append(np.dot(np.dot(T[j],Klocal[j]),inv(T[j])))
		j = j + 1
	j =0
	for i in Con:
		x1 = (i[0] - 1) * 2
		x2 = (i[1] - 1) * 2
		y1 = (i[0] - 1) * 2 + 1
		y2 = (i[1] - 1) * 2 + 1
		Kglobal[x1,x1] += Krotated[j][0,0]
		Kglobal[x1,y1] += Krotated[j][0,1]
		Kglobal[x1,x2] += Krotated[j][0,2]
		Kglobal[x1,y2] += Krotated[j][0,3]
		Kglobal[y1,x1] += Krotated[j][1,0]
		Kglobal[y1,y1] += Krotated[j][1,1]
		Kglobal[y1,x2] += Krotated[j][1,2]
		Kglobal[y1,y2] += Krotated[j][1,3]
		Kglobal[x2,x1] += Krotated[j][2,0]
		Kglobal[x2,y1] += Krotated[j][2,1]
		Kglobal[x2,x2] += Krotated[j][2,2]
		Kglobal[x2,y2] += Krotated[j][2,3]
		Kglobal[y2,x1] += Krotated[j][3,0]
		Kglobal[y2,y1] += Krotated[j][3,1]
		Kglobal[y2,x2] += Krotated[j][3,2]
		Kglobal[y2,y2] += Krotated[j][3,3]
		j = j + 1	
	return Kglobal,T


def BoundaryConditionImplementation(BC,Kglobal,F,u):
    node_values = []
    disp_values = []
    k = 0
    for i in BC:
        if k == 0:
            for j in i:
                node_values.append(j)
            k = k + 1
        else:
            for j in i:
                disp_values.append(j)
            k = k + 1
    for i in range(len(node_values)):
        F = F - Kglobal[:, [node_values[i]-1]] * disp_values[i]
    
    for i in range(len(node_values)):
        Kglobal = np.delete(Kglobal,node_values[i]-1,0)
        Kglobal = np.delete(Kglobal,node_values[i]-1,1)
        F = np.delete(F,node_values[i]-1,0)
        u = np.delete(u,node_values[i]-1,0)
        j = i
        while j < len(node_values)-1:
            if node_values[j+1] > node_values[j]:
                node_values[j+1] = node_values[j+1] - 1
            j = j + 1
	    
    return Kglobal, F, u

