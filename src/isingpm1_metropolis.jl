
#Montecarlo code for ising spin s_i = \pm 1
#Hamiltonian is H= - \sum J_ij s_i s_j  - \sum h_i s_i

function compute_energy(J::Array{Float64,2},
	                S::Array{Int64,1})

	N=size(J)[1]
	E=0.0
	for i=1:N
		E-=J[i,i]*S[i]
		for j=(i+1):N 
			E-=J[i,j]*S[i]*S[j]
		end
	end
return E
end


######### MONTECARLO-Metropolis ######

@compat function do_mcmc_ising_pm1(J::Array{Float64,2},
				beta::Union{Int64,Float64}=1.0,
				N_term::Union{Float64,Int64}=10e5,
			        N_sweep::Union{Float64,Int64}=10e4,
			        dec::Union{Float64,Int64}=10e2)

	N=size(J)[1]

	#Initial configuration
	S=rand([-1,1],N)
	E=compute_energy(J,S) # compute initial energy

	#### thermalisation
	for i=1:N_term
		deltaE=0
		k=rand(1:N) #take a random spin

		#energy before flip
		E1=-J[k,k]*S[k]	
		for j=1:N
			if(j!=k)
				E1-=J[j,k]*S[k]*S[j]
			end
		end

		#flip spin and compute energy
		S[k]=-S[k]
		E2=-J[k,k]*S[k]	
		for j=1:N
			if(j!=k)
				E2-=J[j,k]*S[k]*S[j]
			end
		end

		deltaE=E2-E1

		gamma=exp(-beta*deltaE)
		if(rand()<gamma)
			E+=deltaE
		else
			S[k]=-S[k] #riflippa indietro
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
		E1=-J[k,k]*S[k]	
		for j=1:N
			if(j!=k)
				E1-=J[j,k]*S[k]*S[j]
			end
		end

		#flip spin and compute energy
		S[k]=-S[k]
		E2=-J[k,k]*S[k]	
		for j=1:N
			if(j!=k)
				E2-=J[j,k]*S[k]*S[j]
			end
		end

		deltaE=E2-E1

		gamma=exp(-beta*deltaE)
		if(rand()<gamma)
			E+=deltaE
		else
			S[k]=-S[k] #riflippa indietro
		end
		if(i%dec==0)   #to decorralate
			sampled_conf[cont,:]=S
			cont+=1
		end
	end

return sampled_conf 
end



###### DO MONTECARLO ###############
function mc_ising_pm1(filename::AbstractString,
			     beta::Union{Float64,Int64}=1.0,
			     N_term::Union{Float64,Int64}=10e5,
			     N_sweep::Union{Float64,Int64}=10e4,
			     dec::Union{Float64,Int64}=10e2)


	@printf("%s\n", "Montecarlo-metropolis for ising spin +1,-1")

	#import couplings and fields
	@printf("%s", "Reading couplings and fields..")
	J=import_h_J_ising(filename)
	@printf("%s \n", "done")
	
	N=size(J)[1]

	@printf("%s %d\n", "Num spin=", N)
	@printf("%s %f\n", "beta=", beta)
	@printf("%s %e %s", "Thermalisation steps=", N_term, " \\ ")
	@printf("%s %e %s %e\n", "Sweep steps=", N_sweep, " \\ Decorrelation steps=", dec)
	@printf("%s %e \n", "Num of samples=", N_sweep/dec)

	S=do_mcmc_ising_pm1(J,beta,N_term,N_sweep,dec)

	return S
end


