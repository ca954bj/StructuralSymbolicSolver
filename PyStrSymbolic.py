##########################

import os, sys
file_path = 'mpmath-1.0.0/mpmath'
sys.path.append(os.path.dirname(file_path))

file_path = 'sympy-1.2/sympy'
sys.path.append(os.path.dirname(file_path))

from sympy import *
import time

# ================================= Select Solver ==============================
solver = "Sympy"
#solver = "Giacpy"
if solver == "Giacpy":
	from giacpy import *

# ================================= Structure Solver ===========================

# ========= Define classes =========

class node:
	instances = []
	Counter = 0

	def __init__(self, n, x, y):
		self.n = n
		if type(x) == int:
			self.x = float(x)
		else:
			self.x = x
		if type(y) == int:
			self.y = float(y)
		else:
			self.y = y
		self.load = [0, 0, 0]
		self.imposeddisp = []
		self.postdisp = []
		self.postreaction = []
		node.instances.append(self)
		self.number = node.Counter
		node.Counter += 1

	def defineload(self, Fx, Fy, Mz):
		if self.load[0] != "unknown":
			self.load[0] += Fx
		if self.load[1] != "unknown":
			self.load[1] += Fy
		if self.load[2] != "unknown":
			self.load[2] += Mz

	def defineimposeddisp(self, Dx, Dy, Rz):
		self.imposeddisp = [Dx, Dy, Rz]
		if Dx != "free":
			self.load[0] = "unknown"
		if Dy != "free":
			self.load[1] = "unknown"
		if Rz != "free":
			self.load[2] = "unknown"

	def imposeddisp(self):
		return self.imposeddisp

	def setpostdisp(self, Ux, Uy, Rz):
		self.postdisp = [Ux, Uy, Rz]

	def setpostreaction(self, Fx, Fy, Mz):
		self.postreaction = [Fx, Fy, Mz]

class ElasticBeam:
	instances = []
	Counter = 0
	def __init__(self, n, node1, node2, A, E, I, StartRelease=False, EndRelease=False):
		self.n = n
		self.node1 = node1
		self.node2 = node2
		self.A = A
		self.E = E
		self.I = I
		self.eleload = []
		self.release = [StartRelease, EndRelease]
		self.postlocforce = []
		self.postglbforce = []

		VE = [self.node2.x - self.node1.x, self.node2.y - self.node1.y]
		self.lx = VE[0]
		self.ly = VE[1]
		self.l = sqrt((VE[0]**2+VE[1]**2))
		
		# ================ Setup Local Stiffness Matrix =======
		
		k1 = E*A/self.l
		k2 = 12*E*I/(self.l**3)
		k3 = 6*E*I/(self.l**2)
		k4 = 4*E*I/self.l
		k5 = k4/2
		
		if self.release == [False, False]:
			self.Kl = Matrix([[k1, 0, 0, -k1, 0, 0], [0, k2, k3, 0, -k2, k3], [0, k3, k4, 0, -k3, k5], [-k1, 0, 0, k1, 0, 0], [0, -k2, -k3, 0, k2, -k3], [0, k3, k5, 0, -k3, k4]])
		elif self.release == [True, False]:
			self.Kl = Matrix([[k1, 0, 0, -k1, 0, 0], [0, k2/4, 0, 0, -k2/4, k3/2], [0, 0, 0, 0, 0, 0], [-k1, 0, 0, k1, 0, 0], [0, -k2/4, 0, 0, k2/2, -k3/2], [0, k3/2, 0, 0, -k3/2, 3*k4/4]])
		elif self.release == [False, True]:
			self.Kl = Matrix([[k1, 0, 0, -k1, 0, 0], [0, k2/4, k3/2, 0, -k2/4, 0], [0, k3/2, 3*k4/4, 0, -k3/2, 0], [-k1, 0, 0, k1, 0, 0], [0, -k2/4, -k3/2, 0, k2/4, 0], [0, 0, 0, 0, 0, 0]])
		elif self.release == [True, True]:
			self.Kl = Matrix([[k1, 0, 0, -k1, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [-k1, 0, 0, k1, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0]])


		# ================ Setup Transformation Matrix ========
		
		cosa = VE[0]/self.l 
		sina = VE[1]/self.l 
		T11 = Matrix([[cosa, sina, 0], [-sina, cosa, 0], [0, 0, 1]])
		T12 = Matrix([[0, 0, 0], [0, 0, 0], [0, 0, 0]])
		temp1 = T11.row_join(T12)
		temp2 = T12.row_join(T11)
		self.transfermatrix = temp1.col_join(temp2)

		# ================ Transfer from local stiffness matrix to global stiffness matrix ===========

		self.Kg = Transpose(self.transfermatrix)*self.Kl*self.transfermatrix

		# ================ Add to element set =================
		ElasticBeam.instances.append(self)
		self.number = ElasticBeam.Counter
		ElasticBeam.Counter += 1

	def defineeleload(self, wx, wy):
		# in global axis
		self.eleload = [wx, wy]

	def setpostlocforce(self, Fx1, Fy1, Mz1, Fx2, Fy2, Mz2):
		self.postlocforce = [Fx1, Fy1, Mz1, Fx2, Fy2, Mz2]

	def setpostglbforce(self, Fx1, Fy1, Mz1, Fx2, Fy2, Mz2):
		self.postglbforce = [Fx1, Fy1, Mz1, Fx2, Fy2, Mz2]

class twonodelink:
	instances = []
	counter = 0
	def __init__(self, n, node1, node2, k1, k2, k3, VE=0):
		self.n = n
		self.node1 = node1
		self.node2 = node2
		self.k1 = k1
		self.k2 = k2
		self.k3 = k3
		self.VE = VE
		self.l = 0
		self.Kl = Matrix([[k1, 0, 0, -k1, 0, 0], [0, k2, 0, 0, -k2, 0], [0, 0, k3, 0, 0, -k3], [-k1, 0, 0, k1, 0, 0], [0, -k2, 0, 0, k2, 0], [0, 0, -k3, 0, 0, k3]])

		twonodelink.instances.append(self)
		self.number = twonodelink.counter
		twonodelink.counter += 1

		if VE == 0:
			self.Kg = self.Kl
			self.transfermatrix = eye(6)
		else:
			# ============== setup transpose matrix =======================
			VEl = sqrt(VE[0]*VE[0] + VE[1]*VE[1])
			cosa = VE[0] / VEl
			sina = VE[1] / VEl
			T11 = Matrix([[cosa, sina, 0], [-sina, cosa, 0], [0, 0, 1]])
			T12 = Matrix([[0, 0, 0], [0, 0, 0], [0, 0, 0]])
			temp1 = T11.row_join(T12)
			temp2 = T12.row_join(T11)
			self.transfermatrix = temp1.col_join(temp2)
			self.Kg = Transpose(self.transfermatrix) * self.Kl * self.transfermatrix

class inverse:
	def __init__(self, A):
		size = A.shape[0]
		if size == 1:
			A[0] = simplify(S(1)/A[0])
		else:
			h = ones(1, size)
			for ki in range(size, 0, -1):
				print("Inversing: %d out of %d" % (ki, size))
				p = simplify(expand(A[0, 0]))
				for ii in range(2, size + 1):
					q = simplify(expand(A[ii - 1, 0]))
					if ki < ii:
						h[ii - 1] = simplify(expand(q / p))
					else:
						h[ii - 1] = simplify(expand(-q / p))
					for j in range(2, ii + 1):
						step1 = cancel(q * h[j - 1])
						step11, step12 = fraction(step1)
						A_1, A_2 = fraction(A[ii - 1, j - 1])
						A[ii - 2, j - 2] = expand(step11*A_2 + A_1*step12)/(A_2*step12)
					A[ii - 1, ii - 1] = simplify(expand(1 / p))
				for ii in range(2, size + 1):
					A[size - 1, ii - 2] = h[ii - 1]
			for ii in range(0, size):
				for j in range(0, ii):
					A[j, ii] = A[ii, j]
		print("Inversed")

class doubleinverse:
	def __init__(self, A):
		size = A.shape[0]
		if size <= 5:
			inverse(A)
		else:
			size1 = size/2
			size2 = size - size1
			MA = A[0: size1, 0: size1]
			MD = A[size1:, size1:]
			MB = A[0: size1, size1:]
			MC = A[size1:, 0:size1]
			MAI = MA.inv()
			print(MAI)
			print("DoubleInverse: 1 out of 4")
			MDI = MD.inv()
			print(MDI)
			print("DoubleInverse: 2 out of 4")
			M1 = simplify((MA-MB*MDI*MC))
			M1I = M1.inv()
			print(M1I)
			print("DoubleInverse: 3 out of 4")
			M2 = simplify(MD-MC*MAI*MB)
			print(M2)
			M2I = M2.inv(method='LU', try_block_diag=True)
			print(M2I)
			print("DoubleInverse: 4 out of 4")
			A[0: size1, 0: size1] = M1I
			A[0: size1, size1:] = -MAI*MB*M2I
			A[size1:, 0:size1] = -MDI*MC*M1I
			A[size1:, size1:] = M2I

class detectzeros:
	def __init__(self, A):
		size = A.shape[0]
		self.scope = []
		for i in range(1, size):
			if A[i:size, 0:i] == zeros(size - i, i):
				iRecord = i
				self.scope = iRecord
				break

class Analyze:
	def __init__(self):
		inf1 = Symbol('inf1')
		# =================== Convert ELement Load to Nodal Load =========================================
		for item in ElasticBeam.instances:
			if len(item.eleload) == 2:
				# Input shoud be in global coordinate system
				nodalmoment = item.eleload[1]*item.l*item.lx/12 + item.eleload[0]*item.l*item.ly/12
				nodalloady = item.eleload[1]*item.l/2
				nodalloadx = item.eleload[0]*item.l/2
				item.node1.defineload(nodalloadx, nodalloady, nodalmoment)
				item.node2.defineload(nodalloadx, nodalloady, -nodalmoment)

		# =================== Determine the size of the total global matrix first ========================
		size = len(node.instances)*3
		self.Kstiff = zeros(size)
		self.U = zeros(size, 1)
		self.F = zeros(size, 1)
		self.UnknownSet = symbols('a0:%d' % size)

		# =================== Begin to assemble ==========================================================
		AllElement = ElasticBeam.instances + twonodelink.instances
		for item in AllElement:
			Node1Start = item.node1.number*3
			Node2Start = item.node2.number*3

			self.Kstiff[Node1Start:Node1Start+3, Node1Start:Node1Start+3] = simplify(self.Kstiff[Node1Start:Node1Start+3, Node1Start:Node1Start+3] + item.Kg[0:3, 0:3])
			self.Kstiff[Node1Start:Node1Start+3, Node2Start:Node2Start+3] = simplify(self.Kstiff[Node1Start:Node1Start+3, Node2Start:Node2Start+3] + item.Kg[0:3, 3:6])
			self.Kstiff[Node2Start:Node2Start+3, Node1Start:Node1Start+3] = simplify(self.Kstiff[Node2Start:Node2Start+3, Node1Start:Node1Start+3] + item.Kg[3:6, 0:3])
			self.Kstiff[Node2Start:Node2Start+3, Node2Start:Node2Start+3] = simplify(self.Kstiff[Node2Start:Node2Start+3, Node2Start:Node2Start+3] + item.Kg[3:6, 3:6])

		# ===================== Construct Equations =======================================================
		Vars = symbols('x0:%d'%size)
		self.U1 = zeros(size, 1)
		for item in node.instances:
			if len(item.imposeddisp) > 0:
				for i in range(0, 3):
					if item.imposeddisp[i] == 'free':
						self.U1[item.number*3 + i] = Vars[item.number*3 + i]
					else:
						self.U1[item.number * 3 + i] =  item.imposeddisp[i]
			else:
				self.U1[item.number * 3] = Vars[item.number * 3]
				self.U1[item.number * 3 + 1] = Vars[item.number * 3 + 1]
				self.U1[item.number * 3 + 2] = Vars[item.number * 3 + 2]

		self.F1 = zeros(size, 1)
		for item in node.instances:
			for i in range(0, 3):
				if item.load[i] == "unknown":
					self.F1[item.number*3 + i] = Vars[item.number*3 + i]
				else:
					self.F1[item.number * 3 + i] = item.load[i]

		eqsm = self.Kstiff*self.U1 - self.F1
		eqs = []
		for i in range(0, size):
			eqs.append(eqsm[i])

		for i, obj in enumerate(eqs):
			jd = map(lambda x: x in obj.free_symbols, Vars)
			jdt = filter(lambda x: x==True, jd)
			print(jd, len(jdt))

		starttime = time.time()
		sol = solve(eqs, Vars)

		endtime = time.time()
		print("Solving time %f seconds" % (endtime - starttime))
		
		USolved = self.U1.subs(sol)
		FSolved = self.F1.subs(sol)

		# ================== Case 1: Only Apply Displacements to the structure ======================================
		'''starttime = time.time()
		self.U1 = zeros(size, 1)
		for item in node.instances:
			for i in range(0, 3):
				if len(item.imposeddisp) > 0:
					if item.imposeddisp[i] != 'free' and item.imposeddisp[i] != 0:
						self.U1[item.number*3 + i] = item.imposeddisp[i]

		self.F1 = self.Kstiff*self.U1

		# ================= Case 2: Only Apply Force to the structure ==============================================
		self.F2 = zeros(size, 1)

		for item in node.instances:
			for i in range(0, 3):
				if item.load[i] != "unknown" and item.load[i] != 0:
					self.F2[item.number*3 + i] = item.load[i]
				if len(item.imposeddisp) > 0:
					if item.imposeddisp[i] == 0:
						NodeStart = item.number*3 + i
						for ii in range(0, size):
							if ii == NodeStart:
								self.Kstiff[NodeStart, ii] = 1
							else:
								self.Kstiff[NodeStart, ii] = 0
								self.Kstiff[ii, NodeStart] = 0

		MatrixD = self.Kstiff
		sco = detectzeros(MatrixD).scope
		ljsco = 0
		while sco != []:
			ljsco += sco
			MatrixA = MatrixD[0: sco, 0: sco]
			print("AShape", MatrixA.shape)
			doubleinverse(MatrixA)
			self.Kstiff[ljsco - sco: ljsco, ljsco - sco: ljsco] = MatrixA

			MatrixD = MatrixD[sco:, sco:]
			print("DShape", MatrixD.shape)
			sco = detectzeros(MatrixD).scope
		for i in range(0, MatrixD.shape[0]):
			print(MatrixD[i, 0:])
		doubleinverse(MatrixD)
		self.Kstiff[ljsco:, ljsco:] = MatrixD

		#doubleinverse(self.Kstiff)
		self.U2 = self.Kstiff*self.F2
		endtime = time.time()
		print("Solving time %f seconds" % (endtime - starttime))

		USolved = self.U1 + self.U2
		FSolved = self.F1 + self.F2
		# USolved is displacements, FSolved is not reaction force !!!'''

		# ============================ Post Processing ===================================


		for item in node.instances:
			NodeStart = item.number * 3
			item.setpostdisp(USolved[NodeStart], USolved[NodeStart + 1], USolved[NodeStart + 2])

		for item in ElasticBeam.instances:
			Node1Start = item.node1.number * 3
			Node2Start = item.node2.number * 3
			temp = zeros(6, 1)
			temp[0] = USolved[Node1Start]
			temp[1] = USolved[Node1Start + 1]
			temp[2] = USolved[Node1Start + 2]
			temp[3] = USolved[Node2Start]
			temp[4] = USolved[Node2Start + 1]
			temp[5] = USolved[Node2Start + 2]
			GlbEleFocTemp = simplify(item.Kg * temp)
			GlbEleFoc = zeros(6, 1)
			for i in range(0, 6):
				GlbEleFoc[i] = GlbEleFocTemp[i]
			if len(item.eleload) == 2:
				nodalmoment = item.eleload[1] * item.l * item.lx / 12 + item.eleload[0] * item.l * item.ly / 12
				nodalloady = item.eleload[1] * item.l / 2
				nodalloadx = item.eleload[0] * item.l / 2
				GlbEleFoc[0] = GlbEleFoc[0] - nodalloadx
				GlbEleFoc[3] = GlbEleFoc[3] - nodalloadx
				GlbEleFoc[1] = GlbEleFoc[1] - nodalloady
				GlbEleFoc[4] = GlbEleFoc[4] - nodalloady
				GlbEleFoc[2] = GlbEleFoc[2] - nodalmoment
				GlbEleFoc[5] = GlbEleFoc[5] + nodalmoment

			for i in range(0, 6):
				if inf1 in GlbEleFoc[i].free_symbols:
					GlbEleFoc[i] = limit(GlbEleFoc[i], inf1, oo)
			item.setpostglbforce(GlbEleFoc[0], GlbEleFoc[1], GlbEleFoc[2], GlbEleFoc[3], GlbEleFoc[4], GlbEleFoc[5])
			LocEleFoc = simplify(item.transfermatrix.inv()*GlbEleFoc)
			for i in range(0, 6):
				if inf1 in LocEleFoc[i].free_symbols:
					LocEleFoc[i] = limit(LocEleFoc[i], inf1, oo)
			item.setpostlocforce(LocEleFoc[0], LocEleFoc[1], LocEleFoc[2], LocEleFoc[3], LocEleFoc[4], LocEleFoc[5])
