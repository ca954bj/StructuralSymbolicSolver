execfile('PyStrSymbolic.py')
E, I, A, L1, L2, q = symbols('E, I, A, L1, L2, q', positive=True)
M = Symbol('M')

Node1 = node(1, 0, 0)
Node2 = node(2, L1, 0)
Node3 = node(3, L1 + L2, 0)


Ele1 = ElasticBeam(1, Node1, Node2, A, E, I)
Ele2 = ElasticBeam(2, Node2, Node3, A, E, I)

Node1.defineimposeddisp(0, 0, 0)
Node2.defineimposeddisp(0, 0, 'free')
Node3.defineimposeddisp(0, 0, 'free')

Ele2.defineeleload(0, -q)

Analyze()

dimension = {}
dimension[L1] = 4
dimension[L2] = 6
dimension[q] = 20

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

	
		