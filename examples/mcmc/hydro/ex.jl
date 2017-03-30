import ThreeQ
import PyPlot

@everywhere module Samples

import ThreeQ
import DataStructures
import PyPlot
#solve 1d groundwater equation assuming the left and right boundary have 0 head, with a spacing of 1 between nodes
function solver(ks::Array{Float64, 1}, rs::Array{Float64, 1})
	I = Array{Int}(3 * (length(ks) - 1) + 2)
	J = Array{Int}(3 * (length(ks) - 1) + 2)
	V = Array{Float64}(3 * (length(ks) - 1) + 2)
	b = Array(Float64, length(ks) + 1)
	I[1] = 1; J[1] = 1; V[1] = 1.; b[1] = 1.
	I[2] = length(ks) + 1; J[2] = length(ks) + 1; V[2] = 1.; b[2] = 0.
	push!(I, length(ks) + 1); push!(J, length(ks) + 1); push!(V, 1.); b[end] = 0.
	l = 3
	dx = 1 / length(ks)
	for i = 2:length(ks)
		b[i] = -rs[i - 1]
		I[l] = i
		J[l] = i - 1
		V[l] = ks[i - 1] / dx^2
		I[l + 1] = i
		J[l + 1] = i + 1
		V[l + 1] = ks[i] / dx^2
		I[l + 2] = i
		J[l + 2] = i
		V[l + 2] = -(ks[i - 1] + ks[i]) / dx^2
		l += 3
	end
	A = sparse(I, J, V)
	return A \ b
end

function b2f(kb, klow, khigh)
	local ks = Array(Float64, length(kb))
	for i = 1:length(kb)
		if kb[i] == 1
			ks[i] = khigh
		else
			ks[i] = klow
		end
	end
	return ks
end

function solver(kb, rs, klow, khigh)
	ks = b2f(kb, klow, khigh)
	return solver(ks, rs)
end

function sample(h, klow, khigh, rs; modelargs=("uq1d_model", "laptop", "c4-sw_sample", "workingdir", "c4"), num_reads=100, solver=ThreeQ.DWQMI.defaultsolver, kwargs...)
    model = ThreeQ.Model(modelargs...)
    @ThreeQ.defvar model q[1:length(h) - 1]
	Q = zeros(length(h) - 1, length(h) - 1)
    for i = 2:length(h) - 1
		constterm = klow * h[i + 1] - klow * h[i] - klow * h[i] + klow * h[i - 1] + rs[i - 1]
		qicoeff = (khigh - klow) * (h[i + 1] - h[i])
		qim1coeff = (khigh - klow) * (h[i - 1] - h[i])
		Q[i, i] += qicoeff ^ 2 + 2 * constterm * qicoeff
		Q[i - 1, i - 1] += qim1coeff ^ 2 + 2 * constterm * qim1coeff
		Q[i, i - 1] += 2 * qicoeff * qim1coeff
		@ThreeQ.addquadratic model (constterm + qicoeff * q[i] + qim1coeff * q[i - 1]) ^ 2
    end
	hi, J, _ = ThreeQ.DWQMI.qubo2ising(ThreeQ.matrix2Q(Q))
	hi, J = ThreeQ.rescale(hi, J, 2, 1)
	fig, axs = PyPlot.subplots(2, 3, figsize=(16, 9))
	axs[1][:plot](sort(hi), "k.")
	axs[2][:plot](sort(collect(values(J))), "k.")
	axs[3][:plot](h)
	axs[3][:plot](Samples.solver(fill(Float64(klow), length(q.value)), rs))
	axs[5][:plot](h - Samples.solver(fill(Float64(klow), length(q.value)), rs), "k.")
	embeddings = ThreeQ.stringembed(length(q), ThreeQ.DWQMI.getadjacency(solver))
	adjacency = ThreeQ.DWQMI.getadjacency(solver)
	kbs = Array(Int, length(h) - 1, num_reads)
	@show vecnorm(h - Samples.solver(klow + (khigh - klow) * Float64.(truekb), rs))
	function finishsolve_helper(m, embans, emb)
		answer = 0.5 * (ThreeQ.DWQMI.unembedanswer(embans["solutions"], emb)' + 1)#convert from ising to qubo
		i = 1
		for j = 1:size(answer, 2)
			if j <= 3
				@show map(Int, answer[:, j]) == truekb
				@show map(Int, answer[:, j])
				@show truekb
				axs[4][:plot](Samples.solver(klow + (khigh - klow) * answer[:, j], rs))
				axs[5][:plot](h - Samples.solver(klow + (khigh - klow) * answer[:, j], rs))
				axs[6][:plot](h - Samples.solver(klow + (khigh - klow) * rand([0.0, 1.0], size(answer, 1)), rs))
			end
			occurrences = embans["num_occurrences"][j]
			for l = 1:occurrences
				kbs[:, i] = answer[:, j]
				i += 1
			end
		end
	end
    ThreeQ.solvesapi!(Q; auto_scale=true, solver=solver, num_reads=num_reads, embeddings=embeddings, finishsolve_helper=finishsolve_helper, kwargs...)
    #_, _, i2varstring = ThreeQ.solvesapi!(model; auto_scale=true, solver=solver, num_reads=num_reads, embeddings=embeddings, kwargs...)
    #ThreeQ.solvesapi!(model; auto_scale=true, solver=solver, num_reads=num_reads, kwargs...)
	#=
	kbs = Array(Int, length(h) - 1, num_reads)
	i = 1
    for j = 1:length(model.bitsolutions)
		@ThreeQ.loadsolution model energy occurrences valid j
		if j <= 3
			@show occurrences
			@show map(Int, q.value)
			@show ThreeQ.evalqubo!(model; q=q.value)
			@show energy
			axs[4][:plot](Samples.solver(klow + khigh * q.value, rs))
			axs[5][:plot](h - Samples.solver(klow + khigh * q.value, rs))
			axs[6][:plot](h - Samples.solver(klow + khigh * rand([0.0, 1.0], length(q.value)), rs))
		end
		for l = 1:occurrences
			kbs[:, i] = q.value
			i += 1
		end
    end
	=#
	display(fig); println()
	PyPlot.close(fig)
	return kbs
end

function sample(h, klow, khigh, rs, myloglikelihood; burnin=100, numsteps=100, num_reads=100, kwargs...)
	kbs = sample(h, klow, khigh, rs; num_reads=num_reads, kwargs...)
	bigkbs = Array(Int, length(h) - 1, num_reads * numsteps)
	i = 1
	for j = 1:size(kbs, 2)
		bigkbs[:, i:i + numsteps - 1], _ = ThreeQ.sample(kbs[:, j], burnin, numsteps, myloglikelihood)
		i += numsteps
    end
    return bigkbs
end


const klow = 1
const khigh = 2
truekb = 0

function makellhood(numhs)
	srand(0)
	global truekb = rand([0, 1], numhs - 1)
	@show truekb
	#truers = rand(numhs - 2)
	#truers = ones(numhs - 2)
	truers = zeros(numhs - 2)
	trueh = solver(truekb, truers, klow, khigh)
	sigmaobs = 2e-4
	obsh = trueh + sigmaobs * randn(length(trueh))
	function myloglikelihood(kb)
		h = solver(kb, truers, klow, khigh)
		return -sum((h - obsh).^2) / (2 * sigmaobs^2)
	end
	return myloglikelihood, truekb, truers, obsh
end

function exactcompare(numhs, numcomparisons, burnin, numsteps)
	myloglikelihood, truekb, truers, obsh = makellhood(numhs)
	dwsolver = ThreeQ.DWQMI.getdw2xsys4(Main.mytoken)
	#dwsolver = ThreeQ.DWQMI.defaultsolver
	trueprobs = Dict(zip(ThreeQ.enumerateprobabilities(length(truekb), myloglikelihood)...))
	classicals = SharedArray{Float64}(numcomparisons)
	dwaves = SharedArray{Float64}(numcomparisons)
	for j = 1:numcomparisons
		#do the classical sampling
		classicalsamples, _ = ThreeQ.sample(rand([0, 1], length(truekb)), burnin, numsteps, myloglikelihood)
		classicalmcmcprobs = DataStructures.DefaultDict{Array{eltype(classicalsamples), 1}, Float64}(0.0)
		for i = 1:size(classicalsamples, 2)
			classicalmcmcprobs[classicalsamples[:, i]] += 1 / size(classicalsamples, 2)
		end
		classicals[j] = ThreeQ.bhatdist(trueprobs, classicalmcmcprobs)
	end
	#do the d-wave sampling
	alldwavesamples = sample(obsh, klow, khigh, truers, myloglikelihood; num_reads=numcomparisons, solver=dwsolver, burnin=burnin, numsteps=numsteps, token=Main.mytoken, timeout=Inf)
	for j = 1:numcomparisons
		dwavesamples = alldwavesamples[:, 1 + (j - 1) * numsteps:j * numsteps]
		dwavemcmcprobs = DataStructures.DefaultDict{Array{eltype(dwavesamples), 1}, Float64}(0.0)
		for i = 1:size(dwavesamples, 2)
			dwavemcmcprobs[dwavesamples[:, i]] += 1 / size(dwavesamples, 2)
		end
		dwaves[j] = ThreeQ.bhatdist(trueprobs, dwavemcmcprobs)
	end
	return dwaves, classicals
end

function countsamples(numhs, maxsteps, numreads)
	myloglikelihood, truekb, truers, obsh = makellhood(numhs)
	dwsolver = ThreeQ.DWQMI.getdw2xsys4(Main.mytoken)
	kbs = sample(obsh, klow, khigh, truers; num_reads=numreads, solver=dwsolver, token=Main.mytoken, timeout=Inf)
	targetllhood = myloglikelihood(kbs[:, 1])
	samples, llhoods = ThreeQ.sample(rand([0, 1], length(truekb)), 0, maxsteps, myloglikelihood)
	bestsofar = copy(llhoods)
	for i = 2:length(llhoods)
		bestsofar[i] = max(bestsofar[i - 1], llhoods[i])
	end
	return targetllhood, bestsofar
end

end#end Samples Module

#dw, cl = Samples.exactcompare(20, 100, 10^2, 10^2)
#TODO figure out why the d-wave sampling seems to be worse than random as the number of bits gets large...why isn't the energy in 1-1 correspondence with the evalqubo?...they shouldn't be equal because of the rescaling, but they should be proportional
tllh, bestsofar = Samples.countsamples(100, 10^5, 10000)
@show sum(bestsofar .> tllh - 10)
@show bestsofar[1], bestsofar[end]
@show tllh
#bestsofar
