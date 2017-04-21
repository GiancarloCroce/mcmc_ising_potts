# input - output functions



########## INPUT #########
function offset(i::Int64,
		N::Int64)
	return Int(2+i*N-i*(i-1)/2)
end

function import_h_J_ising(filename::AbstractString)
	x=readdlm(filename)
	N=Int(x[1])
	J=zeros(N,N)
	for i=1:N
		J[i,i]=x[i+1]
	end
	for i=1:N
		off=offset(i,N)
		c=0
		for j=(i+1):N
			J[i,j]=x[off+c]
			J[j,i]=x[off+c]
			c+=1
		end
	end
	return J
end

function import_h_J_potts(filename::AbstractString)

	h=zeros(1,1)
	J=zeros(1,1,1,1)
	N=0
	q=0

	first_line=true
	open(filename) do f
		while !eof(f)
			if(first_line==true)
				line=readline(f)
				line_split=split(line)
				N=parse(Int64,line_split[1])
				q=parse(Int64,line_split[2])
				#println(N,q)
				h=zeros(N,q)
				J=zeros(N,N,q,q)
				first_line=false
			end

		        line = readline(f)
			line_split=split(line)
			#println(line_split)
			ind=map(x->parse(Int64,x),line_split[2:end-1])
			val=parse(Float64,line_split[end])
			if(line_split[1]=="h")
				h[ind[1],ind[2]]=val
			end
			if(line_split[1]=="J")
				J[ind[1],ind[2],ind[3],ind[4]]=val
				J[ind[2],ind[1],ind[4],ind[3]]=val
			end
		end
	end
	return N,q,h,J
end

########## OUTPUT
# stampa magn e corr su file 
function print_magn(magn::Vector)
	f_magn=open("out.ms", "w+") 
	for a=1:length(magn)
		@printf(f_magn,"%f %d \n", magn[a],a)
	end
	close(f_magn)
end


function print_corr(corr::Matrix)
	N=size(corr)[1]
	f_corr=open("out.ccs","w+")
	for a=1:N
		for b=(a+1):N
			@printf(f_corr, "%f %d %d \n", corr[a,b],a,b)
		end
	end
	close(f_corr)
end







