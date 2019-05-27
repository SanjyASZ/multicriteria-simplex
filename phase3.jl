__precompile__()

using LinearAlgebra
using JuMP
# using GLPKMathProgInterface

#PHASE 3 ARTHUR
# Problème auxiliaire qui, si il est non borné nous donne l'efficiency d'une variable entrante xj
# il résoud le PL min{ e'v : Ry - r^jd + Iv = 0,y,d,v <= 0 } où R est la matrice des coef des var hors base des fonctions objectifs
function pj(solverSelected, C::Vector{Float64}, A::Array{Float64,2})
	m = Model(solver = solverSelected)
	nb = length(C)

	nb_obj=size(A,1)

	@variable(m,x[1:nb] >= 0)

	@objective(m,Min,sum(x[i]*C[i] for i in 1:nb))
	@constraint(m,ctr[ictr = 1:nb_obj], sum(x[iv]*A[ictr,iv] for iv in 1:nb) <= 0)
	return m
end

# fonction qui test si la base B est dans la liste de bases L
function isin(L,B)
	isinlist = false
	for i in 1:length(L)
		isinbase = true
		for j in 1:length(B)
			isinbase = isinbase && (B[j] in L[i])
		end
		isinlist = isinlist || isinbase
	end
	return isinlist
end

function p3(C,A,b,B1)
	# construction de la liste de bases et la liste de solutions
	LB = [B1]
	Lx = []

	cb = 1 # indice de la base courante
	# tant qu'il nous reste des bases dans la liste L à explorer
	while cb <= length(LB)
		B = LB[cb] # base courante
		N = [n for n in 1:(nbvar+nbneq) if !(n in B)] # N : variables hors base

		println("B = ",B)

		AB = A[:,B] # matrice A sous la base B

		# décomposition LU
		L,U = lu(AB)# p est le vecteur des permutations que fais julia pour la decomp LU

		# Calcul de R, matrice des coef des var en bases des fonctions objectifs
		# le calcul est fait depuis les donnés du probleme avec la forme révisée + LU
		R = Matrix{Float64}(undef,nbobj,length(N))
		for k in 1:nbobj
			t = U'\C[k,B]
			y = L'\t
			R[k,:] = (C[k,N])' - y'*A[:,N]
		end
		#R = roundm(R);printm(R)
		#construction de la solution xb de la base B
		y = L\b
		xb = U\y
		push!(Lx,xb)

		# énumération des variables hors bases pour tester si elles sont efficientes
		for j in 1:length(N)
			#print("    x",N[j],":")

			# Modélisation et résolution du PL auxiliaire pour savoir si une var est intéressante
			rj = R[:,j]

			A3 = hcat(R,-rj)

			# la fonction objectif est la somme des colonnes de R-rj
			C3 = [sum(A3[:,o]) for o in 1:size(A3,2)]

			m = pj(ClpSolver(),C3,A3)
			status = solve(m, suppress_warnings = true)

			if status != :Optimal
				#println("Inefficient variable")
			else
				# calcul de la colone de la variable entrante xj dans le tableau simplex révisé
				y = L\A[:,N[j]]
				d = U\y

				# Identification de la variable sortante
				ratios = map(x -> if x<0 Inf else x end,(xb./d))
				s,si = findmin(ratios)

				if s == Inf || round(d[si],digits=12) == 0
					if s == Inf
						#println("Rayon infini")
					else
						#println("Pivot impossible")
					end
				else
					#print("x",B[si]," ==> ")
					# construcion de la nouvelle base Bj
					Bj = copy(B)
					Bj[si] = N[j]

					if isin(LB,Bj)
						#println("Base déja vue")
					else
						# ajout de la nouvelle base et de la nouvelle solution dans les listes
						push!(LB,Bj)
						#println("Nouvelle base ")
					end
				end
			end
		end
		cb = cb + 1
	end
	return LB,Lx
end
