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


#PHASE 3 Arthur à optimiser il faudrait enlever les globals ----------------------------------------------------
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
	global nbvar = size(A,2) # nb de variables
	global nbctr = size(A,1) # nb de contraintes
	global nbobj = size(C,1) # nb d' objectifs
#---------------------------------------------------------------------------------


#PHASE 1 --------------------------------------------------
A,b=verif_combi_lin(A,equ_const,b)
res= phase1(C,A,b,equ_const,solverSelected)
println("\n \n")
x0=res[2]

#PHASE 2 -------------------------------------------------
A,equ_cont,b,C,M = formalise(A,C,MM,equ_const,b)
res2=phase2(C,A,b,x0,solverSelected,equ_const)
LABASE=enBase(res2)
println(LABASE)

#PHASE 3 ------------------------------------------------------
println("-P3-")
II,OO = completion(equ_const)
A3 = hcat(A,II)
C3 = hcat(C,OO)
L,Lx = @time p3(C3,A3,b,LABASE)

println("Bases : ",L)
println("Sol associées : ",Lx)
