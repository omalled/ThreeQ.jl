using ToQ

model = ToQ.Model("modelname", "laptop", "c4-sw_sample", "workingdir")

@defparam model color
@defparam model neighbor

@defvar model onevar
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

ToQ.solve(model; color=1, neighbor=1)
