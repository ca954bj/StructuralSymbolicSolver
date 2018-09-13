execfile('PyStrSymbolic.py')
E1, E2, A1, A2, I1, I2, L1, L2, dh, F1, F2, inf1, k3 = symbols('E1, E2, A1, A2, I1, I2, L1, L2, dh, F1, F2, inf1, k3', positive=True)
M = Symbol('M')

Node1 = node(1, 0, -L2)
Node3 = node(3, 0, L2)
Node2 = node(2, 0, 0)
Node4 = node(4, -L1, 0)
Node5 = node(5, L1, 0)
#Node6 = node(6, 0, 0)
#Node7 = node(7, 0, 0)

ColumnDown = ElasticBeam(1, Node1, Node2, A2, E2, I2)
ColumnUp = ElasticBeam(2, Node2, Node3, A2, E2, I2)
BeamLeft = ElasticBeam(3, Node4, Node2, A1, E1, I1)
BeamRight = ElasticBeam(4, Node2, Node5, A1, E1, I1)
#SpringLeft = twonodelink(5, Node6, Node2, inf1, inf1, k3)
#SpringRight = twonodelink(6, Node2, Node7, inf1, inf1, k3)

Node4.defineload(0, -F2, 0)
Node5.defineload(0, F2, 0)
Node1.defineimposeddisp(0, 0, 0)
Node3.defineimposeddisp(0, -dh, 'free')

'''Node1.defineimposeddisp(0, 0, 0)
Node3.defineimposeddisp(0, -dh, 0)
Node4.defineload(0, -F1, 0)
Node5.defineload(0, F1, 0)'''

Analyze()

for item in node.instances:
	print("Node No.%d" % item.n)
	print("Displacements")
	print(item.postdisp)
	#print("Reactions")
	#print(item.postreaction)

for item in ElasticBeam.instances:
	print("ElasticBeam element No.%d" % item.n)
	print("Local Internal Force")
	print(item.postlocforce)


print("End")

	
		