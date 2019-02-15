"""
Compute the Battacharyya distance
"""
function bhatdist(p1s::AbstractDict, p2s::AbstractDict)
	if length(p1s) < length(p2s)
		bhatdist(p2s, p1s)
	end
	bhatcoeff = 0.0
	for k in keys(p2s)
		if haskey(p1s, k)
			bhatcoeff += sqrt(p1s[k] * p2s[k])
		end
	end
	return -log(bhatcoeff)
end

function enumerateprobabilities(numbits, loglikelihood)
	llhoods = Array{Float64}(2^numbits)
	xs = map(collect, reshape(collect(Base.product(fill(0:1, numbits)...)), 2^numbits))
	#for (i, x) in enumerate(xs)
	@sync @Distributed.distributed for i = 1:length(xs)
		llhoods[i] = loglikelihood(xs[i])
	end
	llhoods -= maximum(llhoods)
	probs = exp.(llhoods)
	probs /= sum(probs)
	return xs, probs
end

function sample(x, burnin, numsamples, loglikelihood, thinning=1; probbitflip::Float64=0.25)
	chain = Array{eltype(x)}(length(x), div(numsamples, thinning))
	llhoods = Array{Float64}(div(numsamples, thinning))
	if burnin > 0
		x = sample(x, 0, burnin, loglikelihood, burnin; probbitflip=probbitflip)[1][:, end]
	end
	proposedx = copy(x)
	flipbit = Array{Float64}(length(x))
	lastllhood = loglikelihood(x)
	for i = 1:numsamples
		copy!(proposedx, x)
		rand!(flipbit)
		for j = 1:length(x)
			if flipbit[j] < probbitflip
				proposedx[j] = 1 - proposedx[j]
			end
		end
		thisllhood = loglikelihood(proposedx)
		if thisllhood - lastllhood > log(rand())
			lastllhood = thisllhood
			copy!(x, proposedx)
		end
		if mod(i, thinning) == 0
			chain[:, div(i, thinning)] = x
			llhoods[div(i, thinning)] = lastllhood
		end
	end
	return chain, llhoods
end
