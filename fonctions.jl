__precompile__()

#ensemble de fonctions

using LinearAlgebra

## fonctions de prétraitement pour phase1 phase2

## fonctions utiles à la phase 1 problème admissible
#----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------#

#fonction qui verifie l'existence de combinaison linéaire
function verif_combi_lin(A::Matrix,equ_const::Array{Int64},b::Array{Float64})
    #initialisation matrice resultat
    res=Array{Array{Float64,1},1}(undef,0)
    #indice des contraintes d'égalité possiblement redondantes (vide au début)
    ind=Array{Int64}(undef,0)
    egal=false
    #boucle pour rajout ligne contr égal == (parcous equ_const indique contr ==)
    for i in 1:length(equ_const)
        #si == 3 alors on a une contrainte d'égalité
        if equ_const[i]==3
            #creation ligne coef A + membre b
            ligne=push!(A[i,:],b[i])
            #rajoute ligne dans res
            push!(res,ligne)
            #rajout le num contr d'= concernée
            push!(ind,i)
            #il existe une contr =
            egal=true

            #suppr dans A et equ_const
            A = A[1:size(A,1) .!= i,: ]
            deleteat!(equ_const,i)
        end
    end
    #s'il existe contr d'égalité
    if egal
        #transfo du res en matrix
        res=Arr_to_Mat(res)
        #calcul du rang de res
        res=ordonne(res)
        rang = rank(res)
    else
        #pas de contraintes d'égalité rang nul
        rang = 0
    end

    #si la diff>0 alors il existe contr combinaison linéaire
    for i in 1:size(res)[1] - rang
        #suppr contrainte redondante (supprime les i dernières lignes)
        res = res[1:size(res,1) .!= size(res,1),: ]
    end
    #retour de la matrice A sans les contraintes redondantes avec le vecteur b second membre

    ##############################################""
    #A COMPLETER
    return A,res
end

function echelonnage(M::Matrix)
    M=float(M)
    i=1;j=1
    while i <= size(M,1) && j <= size(M,2)
        if M[i,j]!=0
            M=pivote(M,i,j)
        else
            while M[i,j]==0
                j=j+1
            end
        end
        i=i+1
        j=j+1
    end
    return M
end

function pivote(M::Matrix,ligne_p::Int64,col_p::Int64)
    nb_ligne=size(M,1)
    #mets des 0 sauf sur le pivot
    for i in (ligne_p+1):nb_ligne
        #Ligne i = ligne i - a(ind i, ind col_p)/a(ind ligne_p, ind col_p) * ligne_p
        M[i,:]=M[i,:]- M[i,col_p]/M[ligne_p,col_p]*M[ligne_p,:]
    end
    return M
end

function ordonne(M::Matrix)
    res=Array{Array{Float64,1},1}(undef,0)
    trouve=false
    pivot=0
    c=1
    l=1
    while c <= size(M,2) && pivot < min(size(M,1),size(M,2))
        while l  <= size(M,1) && c <= size(M,2)
            if M[l,c] != 0
                push!(res,M[l,:])
                #suppr contrainte redondante (supprime la ligne l)
                M = M[1:size(M,1) .!= l,: ]
                c=c+1
                l=1
                pivot=pivot+1
            else
                l=l+1
            end
        end
        c=c+1
    end
    for i in 1:size(M,1)
        push!(res,M[i,:])
    end
    res=Arr_to_Mat(res)
    return res
end


##########################################################################################

#fonction objectif avec poids = 1 retourne ::Array{Array{Float64,1},1}
function make_eT( objectiv::Array{Array{Float64,1},1} )
    #nombre de var
    nb_var=length(objectiv[1])
    #nombre de fonctions objectiv
    nb_obj=length(objectiv)
    #result
    res=zeros(nb_var)

    #calcul somme des fonctions objectiv
    for j in 1:nb_var
        for i in 1:nb_obj
            res[j]=res[j]+objectiv[i][j]
        end
    end
    return res
end

#fonction qui transforme un array d'array en matrix
function Arr_to_Mat(tab2D::Array{Array{Float64,1},1})
    res=hcat(tab2D...)
    res=transpose(res)
    res=convert(Matrix,res)
    return res
end

#fonction qui transforme une matrix en array d'array
function Mat_to_Arr(mat::Array{Float64,2})
    #res = array de array vide
    res=Array{Array{Float64,1},1}(undef,0)
    #pour le nombre de lignes dans la matrice
    for i in 1:size(mat)[1]
        #on ajoute dans le res (array) la ligne matrix transform en array
        push!(res,mat[i,:])
    end
    return res
end

#----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------#




## fonctions utiles à la phase 2 recherche 1ere base efficiente
#----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------#

function normalise(constraint::Array{Float64,2},objectiv::Array{Array{Float64,1},1},MinMax::Array{Int64},equ_const::Array{Int64})
    constraint=canonise(constraint,equ_const,b)
    res=convertObj(objectiv,MinMax)
    MinMax=res[2]
    objectiv=res[1]
    return objectiv,constraint,MinMax
end

#fonction qui convertit les maximisation en minimisation retourne Array{Array{Float64,1},1}
function convertObj(objectiv::Array{Array{Float64,1},1},MinMax::Array{Int64})
    for i in 1:length(MinMax)
        if MinMax[i] == 1
            println("- - ", -objectiv[i])
            objectiv[i] = -objectiv[i]
            println("- - ", -objectiv[i])
            MinMax[i]=0
        end
    end
    return objectiv,MinMax
end

#fonction d'attribution de poids aux objectifs retourne ::Array{Array{Float64,1},1}
function scalair(objectiv::Array{Array{Float64,1},1},w::Array{Float64} )
    #nombre de var
    nb_var=length(objectiv[1])
    #nombre de fonctions objectiv
    nb_obj=length(objectiv)
    #result
    res=zeros(nb_var)

    #calcul somme des fonctions objectiv
    for j in 1:nb_var
        for i in 1:nb_obj
            res[j]=res[j]+ w[i]*objectiv[i][j]
        end
    end

    return res
end

#fonction qui change les égalités en inégalités pour construire le dual
function canonise( constraint::Matrix, equ_const::Array{Int64,1} , b::Array{Int64,1}  )
    for i in 1:length(equ_const)
        if equ_const[i]==3
            equ_const[i]=1
            v= constraint[i,:]
            v=-v
            constraint = vcat(constraint,v')
            println(constraint) #LAS
            push!(equ_const,2)
        end
        if equ_const[i]==2
            equ_const=1
            constraint[i,:]=-constraint[i,:]
            b[i]=-b[i]
        end
    end
    return constraint,equ_const,b
end

#construction matrice contr partie gauche A du dual
function make_A( constraint::Array{Float64,2} )
    #creation matrice Identité
    Id = Matrix{Float64}(I,size(constraint,1),size(constraint,1))
    #rajout des var d'écart dans A
    return hcat(constraint,Id)
end

#construction matrice contr partie droite C du dual
function make_C( objectiv::Array{Array{Float64,1},1} , nb_contr::Int64 )
    #creation matrice nb_objectiv
    C = Arr_to_Mat(objectiv)
    println(C) #LAS
    #rajout dans C var_ecart (nb_contr) tout à 0
    Var_add = zeros(Float64,length(objectiv),nb_contr)
    return hcat(C,Var_add)
end

#fonction qui vérifie les variables en base
function enBase(vec)
    res=Array{Int64,1}(undef,0)
    for i in 1:length(vec)
        if string((vec[i]))=="Basic"
            push!(res,i)
        end
    end
    return res
end
