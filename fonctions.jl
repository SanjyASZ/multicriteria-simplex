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
    res_A=Array{Array{Float64,1},1}(undef,0)
    #indice des contraintes d'égalité possiblement redondantes (vide au début)
    ind=Array{Int64}(undef,0)
    egal=false
    #boucle pour rajout ligne contr égal == (parcous equ_const indique contr ==)
    for i in 1:length(equ_const)
        #si equ_const[i]== 3 alors on a une contrainte d'égalité
        if equ_const[i]==3
            #creation ligne coef A + membre b
            ligne=push!(A[i,:],b[i])
            #rajoute ligne dans res
            push!(res_A,ligne)
            #rajout le num contr d'= concernée
            push!(ind,i)
            #il existe une contr =
            egal=true
        end
    end
        #suppr dans A et equ_const et b
        sort!(ind, rev=true)
        for i in ind
            A = A[1:size(res_A,1) .!= i,: ]
        end
        deleteat!(equ_const,ind)
        deleteat!(b,ind)

    #s'il existe contr d'égalité
    if egal
        #transfo du res en matrix
        res_A=Arr_to_Mat(res_A)
        #ordonne sur la famille génératrice
        res_A=ordonne(res_A)
        #calcul du rang de res
        rang = rank(res_A)
    else
        #pas de contraintes d'égalité rang nul
        rang = 0
    end

    #si la diff>0 alors il existe contr combinaison linéaire
    for i in 1:size(res_A)[1] - rang
        #suppr contrainte redondante (supprime la dernière ligne)
        res_A = res_A[1:size(res,1) .!= size(res,1),: ]
    end

    if rang !=0
    #la dernière colonne correspond au second membre
        res_b=res_A[:,size(A)[2]]
        #sppr la dernière colonne dans res_A
        res_A = res_A[:, 1:size(res,1) .!= size(res,1) ]
        #rajoute des contraintes non redondantes parmi les contraintes initiales
        res_A = vcat(A,res_A)
        #rajoute des seconds membres b
        res_b = vcat(b,res_b)
    else
         res_A = A
         res_b = b
    end
    #retour de la matrice A sans les contraintes redondantes avec le vecteur b second membre
    return res_A,res_b
end

#fonction d'échelonnage pas utile dans notre cas mais pourrait servir à d'autres
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

#fonction qui fait un pivot au sens de Gauss pas utilisé non plus
function pivote(M::Matrix,ligne_p::Int64,col_p::Int64)
    nb_ligne=size(M,1)
    #mets des 0 sauf sur le pivot
    for i in (ligne_p+1):nb_ligne
        #Ligne i = ligne i - a(ind i, ind col_p)/a(ind ligne_p, ind col_p) * ligne_p
        M[i,:]=M[i,:]- M[i,col_p]/M[ligne_p,col_p]*M[ligne_p,:]
    end
    return M
end

#fonction qui ordonne la famille libre à récup
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
function make_eT( objectiv::Array{Float64,2} )
    #nombre de var
    nb_var=size(objectiv)[2]
    #nombre de fonctions objectiv
    nb_obj=size(objectiv)[1]
    #result
    res=zeros(nb_var)

    #calcul somme des fonctions objectiv
    for j in 1:nb_var
        for i in 1:nb_obj
            res[j]=res[j]+objectiv[i,j]
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



## fonctions pré-traitement utiles à la phase 2
#----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------#

function formalise(constraint::Array{Float64,2},objectiv::Array{Float64,2},MinMax::Array{Int64},equ_const::Array{Int64}, b::Array{Float64,1})
    constraint,equ_const,b = inferiorise(constraint,equ_const,b)
    res=convertObj(objectiv,MinMax)
    #par soucis de clarté
    MinMax=res[2]
    objectiv=res[1]
    #retourne la forme final du problème en minimisation et inférieure ou égale
    return constraint,equ_const,b,objectiv,MinMax
end

#fonction qui convertit les maximisation en minimisation retourne Array{Array{Float64,1},1}
function convertObj(objectiv::Array{Float64,2},MinMax::Array{Int64})
    for i in 1:length(MinMax)
        if MinMax[i] == 1
            println("- - ", -objectiv[i,:])
            objectiv[i,:] = -objectiv[i,:]
            MinMax[i]=0
        end
    end
    return objectiv,MinMax
end

#fonction d'attribution de poids aux objectifs retourne ::Array{Array{Float64,1},1}
function scalair(objectiv::Array{Float64,2},w::Array{Float64} )
    #nombre de var
    nb_var=size(objectiv)[2]
    #nombre de fonctions objectiv
    nb_obj=size(objectiv)[1]
    #result
    res=zeros(nb_var)
    #calcul somme des fonctions objectiv
    for j in 1:nb_var
        for i in 1:nb_obj
            res[j]=res[j]+ w[i]*objectiv[i,j]
        end
    end
    return res
end

# dans equ_const "3 signifie = | 2 signifie >= et 1 signifie <="

#fonction qui change les supérieure en inégalités inférieures
function inferiorise( constraint::Array{Float64,2} , equ_const::Array{Int64,1} , b::Array{Float64,1}  )
    #parcourt les contraintes
    for i in 1:length(equ_const)
        #si on a une contrainte supérieure ou égale on la transforme en inf
        if equ_const[i]==2
            #simplement pour respecter la notation dans equ_const
            equ_const=1
            #multiplie par -1 la contrainte
            constraint[i,:]=-constraint[i,:]
            #multiplie par -1 second membre
            b[i]=-b[i]
        end
    end
    return constraint,equ_const,b
end

#construction matrice contr partie gauche A
function make_A( constraint::Array{Float64,2},equ_const )
    #creation matrice Identité
    Id = Matrix{Float64}(I,size(constraint,1),size(constraint,1))
    #attention au var d'écart des égalités sont nulles
    for i in 1:size(constraint)[1]
        if equ_const[i]==3
            Id[i,i]=0
        end
    end
    #rajout des var d'écart dans A
    return hcat(constraint,Id)
end

#construction matrice contr partie droite C du dual
function make_C( objectiv::Array{Float64,2} , nb_contr::Int64 )
    #creation matrice nb_objectiv
    println(C) #LAS
    #rajout dans C var_ecart (nb_contr) tout à 0
    Var_add = zeros(size(objectiv)[1],nb_contr)
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
