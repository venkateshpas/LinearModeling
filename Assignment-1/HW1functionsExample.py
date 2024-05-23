import numpy as np
from numpy.linalg import inv

def length():
	L = float(input('Enter the length of the bar in millimeters'))
	return L

def Nelem():
	Nelem = int(input('Enter the number of elements'))
	return Nelem

def Force(Nelem):
	DOF = (Nelem + 1) * 1
	forceInput = float(input('Enter the applied force on the final node in Newtons'))
	F = np.zeros((DOF,1))
	F[DOF-1][0] = forceInput
	return F

def EModulus():
	E = float(input('Enter the elastic modulus in N/mm^2'))
	return E

def BoundaryCondition():
	Rows = 2  
	Columns = int(input("Mention the number of nodes that have boundary conditions as a number. Example: 4"))
	BC = []
	for i in range(Rows):
		if i == 0:
			single_row = list(map(int, input("Enter the node values as a single row followed by space. Example: 1 2 3 4").split()))
			BC.append(single_row)
		else:
			single_row = list(map(int, input("Enter the displacement of the above node values as a single row. Example: 0 1 2 3").split()))
			BC.append(single_row)
	BC = np.array(BC)
	return BC

def ConnectivityMatrix(Nelem):
	response = input("Do you want to define the connectivity matrix ? Yes or No")
	Rows = Nelem
	Colums = 2
	CM = []
	if response == "Yes":
		for i in range(Rows):
			single_row = list(map(int, input("Enter the row of the connectivity matrix one by one").split()))
			CM.append(single_row)
	else:
		for i in range(Rows):
			CM.append([i+1,i+2])
	CM = np.array(CM)
	return CM

def stiffnessMatrix(L,E,Nelem,CM):
	DOF = (Nelem + 1) * 1
	Kglobal = np.zeros((DOF,DOF))
	W1 = float(input("Enter W1 in millimeters"))
	W2 = float(input("Enter W2 in millimeters"))
	t = float(input("Enter the thickness in millimeters"))
	counter = 1
	Kelements = []
	while counter <= Nelem:
		if Nelem == 1:
			Kelements.append(E * (W1 + W2) * t/(2*L))
		else:
			Kelements.append(E * ( W2 + (W1-W2)* (Nelem - counter)/Nelem ) * t *Nelem /L)
		counter = counter + 1
	j = 0
	for i in CM:
		r = i[0] - 1
		c = i[1] - 1
		Kglobal[r,r] += Kelements[j]
		Kglobal[r,c] += -Kelements[j]
		Kglobal[c,r] += -Kelements[j]
		Kglobal[c,c] += Kelements[j]
		j = j + 1	
	return Kglobal

def BoundaryImplementation(BC,Kglobal, F, u):
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
	
	print(Kglobal)

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
	print(Kglobal)
	return Kglobal, F, u, node_values, disp_values


def calcStrain(Nelem,u,L,CM):
	strain = np.zeros((Nelem,1))
	counter = 0
	while counter < Nelem:
		strain[counter] = (u[counter+1] - u[counter])*Nelem/L
		counter += 1
	return strain 
	
def calcStress(E,strain):
	stress = E*strain
	return stress

def calcDisp(L, Nelem, F, E, BC):
	# this is where the input is tranferred to; the functions defined above are called from line 15 and 16 in this example
	CM = ConnectivityMatrix(Nelem)
	Kglobal = stiffnessMatrix(L,E,Nelem,CM)
	print(Kglobal) #Remove
	u = np.zeros((Nelem+1,1))
	Kglobal_BC, F_BC, u_BC, node_values, disp_values = BoundaryImplementation(BC,Kglobal,F,u)
	u_BC = np.dot(inv(Kglobal_BC), F_BC)
	node_no = 1
	count = 0
	while node_no <= Nelem+1:
		if node_no in node_values:
			u[node_no-1] = disp_values[count]
			count += 1
			node_no += 1
		else:
			u[node_no-1] = u_BC[node_no-1-count]
			node_no += 1

	strain = calcStrain(Nelem,u,L,CM) # just an example of a function, input is also just to have an input
	stress = calcStress(E,strain) # just an example of a function, input is also just to have an input
	return u, strain, stress 

