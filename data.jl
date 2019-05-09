__precompile__()

# test and example

# dans equ_const "3 signifie = | 2 signifie >= et 1 signifie <="

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------

#exemple Mathias
MM=[0,0,0]
objectiv=[ [-1.0, -2.0, 0.0] , [ -1.0,0.0,2.0], [ 1.0,0.0,-1.0] ]
C= Arr_to_Mat(objectiv)
#matrice contraintes
constraint=[ 1.0 1.0 0.0 ; 0.0 1.0 0.0; 1.0 -1.0 1.0]
equ_const=[1,1,1]
#membre de droite
b=[1.0,2.0,4.0]

# # # exemple dégénéréscence cours
# equ_const=[1,1,1]
# MM=[1]
# objectiv=[ [10.0, -57.0, -9.0, -24.0 ] ]
# #matrice contraintes
# constraint=[ 0.5 -5.5 -2.5 9.0 ; 0.5 1.5 -0.5 1.0 ; 1.0 0.0 0.0 0.0]
# #membre de droite
# b=[0.0,0.0,1.0]

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------

# #ex STEUER481
# MM=[0,0,0,0]
# equ_const=[1,1,1,1,1,1]
# C =
# [
#     -6.0 -5.0 -4.0 5.0 -4.0 -2.0 ;
#     -5.0 0.0 0.0 -3.0 -2.0 2.0 ;
#     -6.0 3.0 1.0 5.0 -1.0 -6.0 ;
#     1.0 -3.0 3.0 -1.0 -1.0 2.0
# ]
#
# objectiv=Mat_to_Arr(C)
# println(" objectiv = ",objectiv)
#
# constraint =
# [
#     4.0 5.0 0.0 4.0 5.0 6.0;
#     2.0 4.0 -1.0 0.0 0.0 0.0;
#     4.0 3.0 0.0 5.0 6.0 3.0;
#     0.0 1.0 -1.0 2.0 1.0 0.0;
#     0.0 -1.0 -1.0 0.0 1.0 0.0;
#     0.0 0.0 1.0 0.0 0.0 0.0
# ]
#
# println("constraint = ",constraint)
#
# b = [9.0,8.0,9.0,8.0,10.0,6.0]

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------

#ex STEUER482
# C =
# [
# -3.0 -3.0 1.0 0.0 -2.0 1.0 1.0 0.0 0.0 1.0;
# 2.0 1.0 -2.0 0.0 1.0 0.0 0.0 0.0 2.0 1.0;
# 0.0 -2.0 -1.0 -1.0 2.0 0.0 -2.0 -3.0 0.0 1.0;
# -3.0 -3.0 1.0 0.0 0.0 0.0 -1.0 1.0 -1.0 -3.0;
# 2.0 0.0 0.0 0.0 -2.0 -3.0 -1.0 -3.0 1.0 1.0
# ]
#
# M=[0,0,0,0,0]
# equ_const=[1,1,1,1,1,1,1,1,1,1,1,1,1,1]
#
# objectiv=Mat_to_Arr(C)
#
# A =
# [
# 0 0 0 0 0 0 0 0 5 0;
# 5 0 0 0 4 0 0 0 0 4;
# 0 5 0 -1 0 0 0 0 0 0;
# 5 0 3 0 0 0 0 0 0 0;
# 0 5 3 0 3 0 0 5 4 0;
# 3 0 4 0 2 0 -1 0 0 0;
# 4 0 0 0 0 0 0 0 0 0;
# 0 4 0 0 1 5 0 1 0 0;
# 3 0 0 0 0 0 -1 4 0 0;
# 0 0 0 0 0 0 2 0 0 0;
# 0 0 0 0 5 0 0 0 5 0;
# 0 -1 0 -1 5 4 0 1 0 0;
# 2 5 0 4 0 3 0 -1 0 0;
# 0 4 0 5 0 0 0 1 0 0
# ]
# constraint=float(constraint)
#
# b = [10,10,10,5,10,6,9,7,10,5,5,6,6,8]
# b=float(b)

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
