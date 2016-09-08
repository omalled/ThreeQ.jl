using ToQ

model = ToQ.Model("canada_model", "laptop", "c4-sw_sample", "workingdir", "c4")

@defparam model color
@defparam model neighbor

provinces = ["BC", "YK", "NW", "AB", "SK", "NV", "MT", "ON", "QB", "NB", "NS", "PE", "NL"]
#province2index = Dict(zip(provinces, 1:length(provinces)))
neighbors = Dict()
neighbors["BC"] = ["YK", "NW", "AB"]
neighbors["YK"] = ["BC", "NW"]
neighbors["NW"] = ["YK", "BC", "AB", "SK", "NV"]
neighbors["AB"] = ["BC", "NW", "SK"]
neighbors["SK"] = ["AB", "NW", "MT"]
neighbors["NV"] = ["NW", "MT"]
neighbors["MT"] = ["NV", "SK", "ON"]
neighbors["ON"] = ["MT", "QB"]
neighbors["QB"] = ["ON", "NB", "NL"]
neighbors["NB"] = ["QB", "NS"]
neighbors["NS"] = ["NB"]
neighbors["PE"] = []
neighbors["NL"] = ["QB"]

@defvar model province_rgb[1:length(provinces), 1:3]

#add color penalties
for i = 1:length(provinces)
	for j = 1:3
		@addterm model -1 * color * province_rgb[i, j]
		for k = 1:j - 1
			@addterm model 2 * color * province_rgb[i, j] * province_rgb[i, k]
		end
	end
end

#add neighbor penalties
for j = 1:length(provinces)
	for k = 1:j - 1
		if provinces[k] in neighbors[provinces[j]]
			for i = 1:3
				@addterm model neighbor * province_rgb[j, i] * province_rgb[k, i]
			end
		end
	end
end

#solve the system
solver = DWQMI.defaultsolver
ToQ.solvesapi!(model; solver=solver, color=1, neighbor=5, param_chain=4, num_reads=100, auto_scale=false)

#load the solutions
i = 1
solutions = Array{Float64, 2}[]
energies = Float64[]
occurrences = Int[]
valid = Bool[]
numsolutions = ToQ.getnumsolutions(model)
for i = 1:numsolutions
	@loadsolution model energy occurrencesi validi i
	push!(solutions, copy(province_rgb.value))
	push!(energies, energy)
	push!(occurrences, occurrencesi)
	push!(valid, validi)
end

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
