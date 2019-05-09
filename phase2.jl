__precompile__()

using LinearAlgebra
using JuMP
# using GLPKMathProgInterface
using MathProgBase

#choix solver

function setLP2_1(objectiv::Array{Array{Float64,1},1},constraint::Matrix,b::Array{Float64,1},x0::Array{Float64},solverSelected)
    lp= Model(solver=solverSelected)

    #nbligne(contr) et Col(var)
    nb_var=size(constraint,2)
    nb_contr=size(constraint,1)
    nb_obj=length(objectiv)

    A=make_A(constraint)
    C=make_C(objectiv,nb_contr)

    #définition variable
    @variable(lp, u[1:nb_contr])
    @variable(lp, w[1:nb_obj])

    println("DOT DOT ",dot(u,b)) #las

    #définition objectif
    @objective(lp, Min, dot(u,b) + dot(w,C*x0) )

    println( "u ",u ) #las
    println( "A ",A ) #las
    println( "C ",C ) #las

    #définition contraintes
    @constraint(lp, cte[j=1:nb_var+nb_obj], sum(A[i,j] * u[i] for i=1:nb_contr) + sum(C[i,j] * w[i] for i=1:nb_obj) >= 0 )
    @constraint(lp, cte[i=1:nb_obj], w[i] >= 1 )

    #return
    return lp,u,w

end

function setLP2_2(poids::Array{Float64,1},A::Matrix,b::Array{Float64,1},solverSelected,equ_const::Array{Int64,1})
    # facilitation creation LP
    lp = Model(solver=solverSelected )
    #lp = vModel(solver=solverSelected) #VOPTSOLVER :)

    #nb colon
    nb_var=size(A,2)
    #nb ligne
    nb_ligne=size(A,1)

    #définition variable

    libre = Array{Int,1}()
    pas_libre = Array{Int,1}()

    for i in 1:nb_var
        if equ_const[i]== 3
            push!(libre,i)
        else
            push!(pas_libre,i)
        end
    end

    println("Libre = ", libre, "\n pas_libre= ", pas_libre) #LAS

    #A
    if length(pas_libre)!=0
        @variable(lp, x[ 1:length(pas_libre) ] >= 0)
    end
    if length(libre)!=0
        @variable(lp, x[ (length(pas_libre)+1):(length(pas_libre) + length(libre)) ])
    end

    #définition objectif
    #@objective(lp, Min, 0)
    @objective(lp,Min,dot(poids, x)) #LAS

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

function phase2(objectiv::Array{Array{Float64,1},1}, constraint::Array{Float64,2}, b::Array{Float64,1}, x0::Array{Float64} , solverSelected,equ_const::Array{Int64})
    nb_contr=length(constraint)
    nb_obj=length(objectiv)

    lp, lp_u, lp_w = setLP2_1(objectiv,constraint,b,x0,solverSelected)

    println("\n \n Resolution: Phase 2-1 existenz efficience solution? "); @time solve(lp)

    #solution
    w = getvalue(lp_w)

    println("the LP correspond to : ", lp)
    println(" u,w  values : ", lp_u,lp_w)

    #affichage
    if (isnan(w[1]))
        println("Your MOLP doesn't have any efficient solution")
        basic_or_not=0
    else

        println("There is a feasible set of efficient solution")
        println(" la valeur de w = ",w)

        A=constraint
        poids = scalair(objectiv,w)
        lp, x = setLP2_2(poids,A,b,solverSelected,equ_const)
        println("\n \n Resolution: Phase 2-2 première solution de base efficiente"); @time solve(lp)

        sol_eff = getvalue(x)

        println("the LP correspond to : ", lp)
        println(" u,w  values : ", lp_u,lp_w)

        #affichage
        if (isnan(sol_eff[1]))
            println("Your MOLP doesn't have any efficient solution")
        else
            var_ecart = -(constraint*sol_eff-b)
            println("There is a Basic Feasible efficient solution")
            println(" la valeur de sol = ",sol_eff,var_ecart)
            var_total=vcat(sol_eff,var_ecart)
            enBase=Array{Int64}(undef,0)
            for i in 1:length(var_total)
                if var_total[i]!=0
                    push!(enBase,i)
                end
            end
        end

        if GLPK
            println("INTERAL MODEL 2 GLPKK !! \n",internalmodel(lp))
            #GLPK.get_bhead(lp, 2)
            basic_or_not = Base.getindex(lp,x)
            println("BBBB=", basic_or_not)
        else
            println("INTERAL MODEL 2 !! \n",internalmodel(lp))
            cbasis, rbasis = MathProgBase.getbasis(internalmodel(lp))
            basic_or_not = vcat(cbasis,rbasis)
        end
    end

    return basic_or_not
end
