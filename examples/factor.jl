using ToQ

model = ToQ.Model("factor_model", "laptop", "c4-sw_sample", "workingdir", "c4")

numbits = 3

@defparam model factor
@defparam model ancillary

@defvar model n1[1:numbits]
@defvar model n2[1:numbits]
@defvar model s[1:numbits, 1:numbits]

N = 35

for i = 1:numbits
	for j = 1:numbits
		#make s[i, j] = n1[j] * n2[i]
		@addterm model ancillary * n1[j] * n2[i]
		@addterm model ancillary * -2 * n1[j] * s[i, j]
		@addterm model ancillary * -2 * n2[i] * s[i, j]
		@addterm model ancillary * 3 * s[i, j]
		#make n1 * n2 = N
		for i2 = 1:numbits
			for j2 = 1:numbits
				@addterm model factor * 2 ^ (i + j + i2 + j2) * s[i, j] * s[i2, j2]
			end
		end
		@addterm model factor * -2 * 2 ^ (i + j) * N * s[i, j]
	end
end


#solve the system
numreads = 10 ^ 4
@time ToQ.solve!(model; factor=100, ancillary=10, param_chain=1, numreads=numreads, doembed=true)

#load the solutions
i = 1
numsolutions = ToQ.getnumsolutions(model)
successcount = 0
for i = 1:numsolutions
	@loadsolution model energy occurrencesi validi i
	if (n1.value == [1., 0., 1.] && n2.value == [1., 1., 1.]) || (n2.value == [1., 0., 1.] && n1.value == [1., 1., 1.])
		successcount += 1
	end
	#=
	@show energy
	@show n1.value
	@show n2.value
	=#
end
@show successcount / numreads

#=
#print the solutions
validcount = 0
for i = 1:length(energies)
	isvalid = true
	for j = 1:length(provinces)
		for k = 1:j - 1
			if provinces[k] in neighbors[provinces[j]] && norm(solutions[i][j, :] - solutions[i][k, :]) == 0
				@show k, provinces[k]
				@show j, provinces[j]
				isvalid = false
			end
		end
	end
	if isvalid
		validcount += 1
	end
	#=
	println("Solution #$i (valid = $isvalid)")
	println("Energy: $(energies[i])")
	println("Occurrences: $(occurrences[i])")
	println("Solution:\n$(solutions[i])\n")
	=#
end
@show validcount / length(solutions)
=#
