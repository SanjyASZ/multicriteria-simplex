__precompile__()

using LinearAlgebra
using JuMP
#using CPLEX
using GLPKMathProgInterface
using Clp

#choix solver
#solverSelected = GLPKSolverLP()
GLPK=false
#solverSelected = CplexSolver()
#solverSelected = GurobiSolver()
solverSelected= ClpSolver()
#GLPK=false

include("fonctions.jl")
include("data.jl")
include("phase1.jl")
include("phase2.jl")
include("phase3.jl")

nbctr = size(constraint)[1]
nbvar = size(constraint)[2]
nbobj = length(objectiv)

function completion(eq)
	global nbneq = nbctr - length(findall(!isodd,eq))
	I = zeros(Int,nbctr,nbneq)
	global i = 1
	for j in 1:nbctr
		if eq[j] != 0
			I[j,i] = eq[j]
			global i = i + 1
		end
	end
    i=1
	return I,zeros(Int,nbobj,nbneq)
end

#constraint=verif_combi_lin(constraint,equ_const,b)
#PHASE 1 --------------------------------------------------
res= phase1(objectiv,constraint,b,equ_const,solverSelected)
println("\n \n")
x0=res[2]

#PHASE 2 -------------------------------------------------
res2=phase2(objectiv,constraint,b,x0,solverSelected,equ_const)
LABASE=enBase(res2)
println(LABASE)

#PHASE 3 ------------------------------------------------------
	println("-P3-")
	II,OO = completion(equ_const)
	A3 = hcat(A,II)
	C3 = hcat(C,OO)
	L,Lx = @time p3(C3,A3,b,LABASE)

	println("")
	if test
		#Ls,Lxs = @time p3test(C,A,b,B1)
		Ls,Lxs = @time p3sorted(C,A,b,B1)

		println("")
		println("")

		same = true
		for i in 1:length(Ls)
			same = same && isin(L,Ls[i])
		end
		for i in 1:length(L)
			same = same && isin(Ls,L[i])
		end
		println(" sames : ",same)
		println("Bases test : ",Ls)
		println("Sol associées : ",Lxs,"\n")
	end
	println("Bases : ",L)
	println("Sol associées : ",Lx)
