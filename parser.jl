
function parsemop(fname)
	f=open(fname)
	for i in 1:5
		readline(f)
	end
	nobj = parse.(Int, split(readline(f))[3] )
	nctr = parse.(Int, split(readline(f))[3] )
	nvar = parse.(Int, split(readline(f))[3] )

	for i in 1:5
		readline(f)
	end

	C = zeros(Int,nobj,nvar)
	A = zeros(Int,nctr,nvar)
	b = zeros(Int,nctr)
	eq = ones(Int,nctr)

	line = split(readline(f))
	while line[1] != "ENDATA"
		if length(line)>1 && line[1][1] == 'C'
			j = parse(Int,line[1][2:end])
			i = parse(Int,line[2][2:end])
			if line[2][1] == 'R'
				A[i,j] = parse(Int,line[3])
			else
				C[i,j] = parse(Int,line[3])
			end
		elseif length(line)>1 && line[1][1] == 'R'
			b[parse(Int,line[2][2:end])] = parse.(Int, line[3])
		end
       		line = split(readline(f))
	end
	close(f)
	return C, A, b, eq
end

 
