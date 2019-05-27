__precompile__()

using LinearAlgebra
using JuMP
#using CPLEX
using GLPKMathProgInterface
using Clp

GLPK=false

#choix solver
# solverSelected = GLPKSolverLP()
#GLPK=false
#solverSelected = CplexSolver()
#solverSelected = GurobiSolver()
solverSelected= ClpSolver()
#GLPK=false

include("fonctions.jl")
include("data/molp_10_779_10174_entropy.jl")
include("phase1.jl")
include("phase2.jl")
include("phase3.jl")
include("parser.jl")

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
#----------------------------------------------------------

C=float(C)
A=float(A)
b=float(b)

#PHASE 1 --------------------------------------------------
# A,b=verif_combi_lin(A,eq,b)
res= phase1(C,A,b,eq,solverSelected)
println("\n \n")
x0=res[2]

#PHASE 2 -------------------------------------------------
# A,equ_cont,b,C,M = formalise(A,C,MM,eq,b)
res2=phase2(C,A,b,x0,solverSelected,eq)
LABASE=enBase(res2)
println("LaBase = ",LABASE)

#PHASE 3 -------------------------------------------------
println("-P3-")
Id,O = completion(eq)
A3 = hcat(A,Id)
C3 = hcat(C,O)
L,Lx = @time p3(C3,A3,b,LABASE)

println("Bases : ",L)
println("Sol associées : ",Lx)
