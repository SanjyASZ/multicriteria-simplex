__precompile__()

using LinearAlgebra
using JuMP
# using GLPKMathProgInterface

#Mise en place du LP phase1
function setLP1(Cout::Array{Float64,1},A::Matrix,b::Array{Float64,1},solverSelected, equ_const::Array{Int64})
    # facilitation creation LP
    lp = Model(solver=solverSelected)
    #lp = vModel(solver=solverSelected) #VOPTSOLVER :)

    #nb colon
    nb_var=size(A,2)
    #nb ligne
    nb_ligne=size(A,1)

    #définition variable
    @variable(lp, x[1:nb_var] >= 0 )

    #définition objectif
    @objective(lp, Min, 0)

    #@constraint(lp , cte[i=1:nb_ligne], sum( A[i,j] * x[j] for j=1:nb_var ) <= b[i] )
    #définition contraintes
    for i in 1:nb_ligne
        if equ_const[i]==1
            @constraint(lp , sum( A[i,j] * x[j] for j=1:nb_var ) <= b[i] )
        elseif equ_const[i]==2
            @constraint(lp , sum( A[i,j] * x[j] for j=1:nb_var ) >= b[i] )
        else
            @constraint(lp , sum( A[i,j] * x[j] for j=1:nb_var ) == b[i] )
        end
    end

    #return le lp et les valeurs de x
    return lp, x
end

#PHASE 1 à donner sous forme standard
function phase1(objectiv::Array{Array{Float64,1},1}, constraint::Array{Float64,2}, b::Array{Float64,1}, equ_const::Array{Int64},solverSelected)
    #verification des combinaisons linéaire

    #nombre de contraintes
    nb_contr=size(constraint)[1]
    #transformation des objectiv rajout des var d'écart
#    objectiv=transfo_objectiv(objectiv,nb_contr)
    #nombre de variables
    nb_var=length(objectiv[1])

    #constru cout objectif nul car on veut la phase 1 du simplexe 1ere sol admissible
    zeroo=zeros(Float64,nb_var)

    #matrice des contraintes correspond à A
    A=constraint

    #mise en place LP Phase1
    lp, lp_x = setLP1(zeroo,A,b,solverSelected,equ_const)
    println("\n \n Resolution: Phase 1... Le MOLP est il faisable?"); @time solve(lp)

    println("\n the LP correspond to : ", lp)
    println(" x values : ", lp_x)

    #affichage et résolution du lp

    #println("Resolution: Phase 1..."); @time solve(lp , method = :lex)
    println("\n Resolution finie: ")

    #solution
    x=getvalue(lp_x)
    z= dot(x,make_eT(objectiv))
    WOW = getrawsolver(lp)

    #affichage
    if (isnan(z))
        println("\n Your linear problem is infeasible \n there isn't any basic feasible solution")
        var_ecart = -(constraint*x-b)
    else
        var_ecart = -(constraint*x-b)
        println("There is a feasible set of solution")
        println( "\n var de decision + var d'écart \n",
        x ,  var_ecart)
    end

    #result
    return(z,vcat(x,var_ecart))
end
