


include("io.jl") #dove ci sono le funzioni di input e di output 

function compute_energy_potts(h::Array{Float64,2},
		 	      J::Array{Float64,4},
		              S::Vector{Int64})

	N=size(h)[1]
	q=size(h)[2]
	E=0.0
	for i=1:N
		E-=h[i,S[i]]
		for j=(i+1):N 
			E-=J[i,j,S[i],S[j]]
		end
	end
return E
end


######### MONTECARLO-Metropolis ######

@compat function mcmc_potts(h::Array{Float64,2},
			    J::Array{Float64,4},
			    beta::Union{Int64,Float64}=1.0,
			    N_term::Union{Float64,Int64}=10e5,
			    N_sweep::Union{Float64,Int64}=10e4,
			    dec::Union{Float64,Int64}=10e2)

	N=size(J)[1]
	q=size(h)[2]

	#initial (random) configuration
	S=rand(1:q,N)
	
	E=compute_energy_potts(h,J,S) # compute initial energy

	#### equilibration at givet T
	E_old=E
	E_stationary=false
	for i=1:N_term
		deltaE=0
		k=rand(1:N) #take a random spin
		#energy before flip
		s1=S[k]
		e1=-h[k,s1]
		for j=1:N
			if(j!=k)
				e1-=J[k,j,s1,S[j]]
			end
		end
		#flip spin
		list_q=collect(1:q)
		deleteat!(list_q, findin(list_q, s1))
		s2=rand(list_q)
		e2=-h[k,s2]
		for j=1:N
			if(j!=k)
				e2-=J[k,j,s2,S[j]]
			end
		end

		deltaE=e2-e1
		#println("E2 ",E1," E2 ",E2)	
		gamma=exp(-beta*deltaE)
		if(rand()<gamma)
			E+=deltaE
			S[k]=s2 #accetta new conf
		end
	end

	#### sweep
	N_data=div(N_sweep,dec) #number of uncorrelated data, div is truncated division
	sampled_conf=Array{Int64}(Int64(N_data),Int64(N))
	cont=1
	for i=1:N_sweep
		deltaE=0
		k=rand(1:N) #take a random spin
		#energy before flip
		s1=S[k]
		e1=-h[k,s1]
		for j=1:N
			if(j!=k)
				e1-=J[k,j,s1,S[j]]
			end
		end
		#flip spin
		list_q=collect(1:q)
		deleteat!(list_q, findin(list_q, s1))
		s2=rand(list_q)
		e2=-h[k,s2]
		for j=1:N
			if(j!=k)
				e2-=J[k,j,s2,S[j]]
			end
		end

		deltaE=e2-e1
		#println("E2 ",E1," E2 ",E2)	
		gamma=exp(-beta*deltaE)
		if(rand()<gamma)
			E+=deltaE
			S[k]=s2 #accetta new conf
		end
		if(i%dec==0)
			sampled_conf[cont,:]=S
			cont+=1
		end
	end

return sampled_conf
end

###### DO MONTECARLO ######
function mc_potts(filename::AbstractString,
		  beta::Union{Float64,Int64}=1.0,
		  N_term::Union{Float64,Int64}=10e5,
		  N_sweep::Union{Float64,Int64}=10e4,
		  dec::Union{Float64,Int64}=10e2)

	@printf("%s\n", "Montecarlo-metropolis for (generalized) potts model")
	
	#import couplings and fields
	@printf("%s", "Reading couplings and fields..")
	N,q,h,J=import_h_J_potts(filename)
	@printf("%s \n", "done")

	@printf("%s %d\n", "Num spin=", N)
	@printf("%s %d\n", "Num of colors=", q)
	@printf("%s %f\n", "beta=", beta)
	@printf("%s %e %s", "Thermalisation steps=", N_term, "\\ ")
	@printf("%s %e %s %e\n", "Sweep steps=", N_sweep, " \\ Decorrelation steps=", dec)
	@printf("%s %e \n", "Num of samples=", N_sweep/dec)

	S=mcmc_potts(h,J,beta,N_term,N_sweep,dec)

return S
end

