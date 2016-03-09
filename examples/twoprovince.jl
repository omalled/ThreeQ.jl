using ToQ

model = ToQ.Model("modelname", "laptop", "c4-sw_sample", "workingdir", "c4")

@defparam model color
@defparam model neighbor

@defvar model province_rgb[1:2, 1:3]

#add color penalties
for i = 1:2
	for j = 1:3
		@addterm model -1 * color * province_rgb[i, j]
		for k = 1:j - 1
			@addterm model 2 * color * province_rgb[i, j] * province_rgb[i, k]
		end
	end
end

#add neighbor penalties
for i = 1:3
	@addterm model neighbor * province_rgb[1, i] * province_rgb[2, i]
end

#solve the system
ToQ.solve!(model; color=1, neighbor=1, param_chain=2, numreads=100)

#load the solutions
i = 1
solutions = Array{Float64, 2}[]
energies = Float64[]
occurrences = Int[]
numsolutions = ToQ.getnumsolutions(model)
for i = 1:numsolutions
	@loadsolution model energy occurrencesi valid i
	push!(solutions, copy(province_rgb.value))
	push!(energies, energy)
	push!(occurrences, occurrencesi)
end
for i = 1:length(energies)
	isvalid = norm(sum(solutions[i], 2) - [1., 1.]) == 0 && sum(solutions[i]) == 2
	println("Solution #$i (valid = $isvalid)")
	println("Energy: $(energies[i])")
	println("Occurrences: $(occurrences[i])")
	println("Solution:\n$(solutions[i])\n")
end
