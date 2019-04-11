module DWQMI

import PyCall

dwlocal = PyCall.PyNULL()
dwcore = PyCall.PyNULL()
dwfv = PyCall.PyNULL()
dwutil = PyCall.PyNULL()
dwembed = PyCall.PyNULL()
dwremote = PyCall.PyNULL()
defaultsolver = PyCall.PyNULL()
defaultadjacency = []
defaulturl = "https://127.0.0.1:10443/sapi/"
function __init__()
	global dwlocal
	global dwcore
	global dwfv
	global dwutil
	global dwembed
	global dwremote
	global defaultsolver
	global defaultadjacency
	copy!(dwlocal, PyCall.pyimport("dwave_sapi2.local"))
	copy!(dwcore, PyCall.pyimport("dwave_sapi2.core"))
	copy!(dwfv, PyCall.pyimport("dwave_sapi2.fix_variables"))
	copy!(dwutil, PyCall.pyimport("dwave_sapi2.util"))
	copy!(dwembed, PyCall.pyimport("dwave_sapi2.embedding"))
	copy!(dwremote, PyCall.pyimport("dwave_sapi2.remote"))
	copy!(defaultsolver, dwlocal.local_connection.get_solver("c4-sw_sample"))
	defaultadjacency = collect(dwutil.get_hardware_adjacency(defaultsolver))
end

include("utils.jl")


function connect(token, url=defaulturl)
	return dwremote.RemoteConnection(url, token)
end

function getadjacency(solver)
	return collect(dwutil.get_hardware_adjacency(solver))
end

function getremotesolvers(connection)
	return connection.solver_names()
end

function getremotesolver(connection, name)
	return connection.get_solver(name)
end

function getdw2q(token)
	connection = connect(token)
	return getremotesolver(connection, "DW_2000Q_LANL")
end

function getvfyc(token)
	connection = connect(token)
	return getremotesolver(connection, "DW_2000Q_VFYC_LANL")
end

function findembeddings(Q, adjacency=defaultadjacency; verbose=0, tries=100, timeout=300, kwargs...)
	pythonembedding = dwembed.find_embedding(collect(keys(Q)), adjacency, verbose=verbose, tries=tries, timeout=timeout)
	if typeof(pythonembedding) == Array{Any, 2}
		embedding = convert(Array{Array{Int, 1}}, map(i->vec(pythonembedding[i, :]), 1:size(pythonembedding, 1)))
	else
		embedding = convert(Array{Array{Int, 1}}, pythonembedding)
	end
	return embedding
end

function qubo2ising(Q::Matrix)
	@assert size(Q, 1) == size(Q, 2)
	hJ = zeros(size(Q))
	for i = 1:size(Q, 1)
		for j = 1:size(Q, 2)
			if i == j
				hJ[i, i] += 0.5 * Q[i, i]
			else
				hJ[i, i] += 0.25 * (Q[i, j] + Q[j, i])
				hJ[i, j] +=  0.25 * Q[i, j]
			end
		end
	end
	return hJ
end

function ising2qubo(hJ::Matrix)
	@assert size(hJ, 1) == size(hJ, 2)
	Q = zeros(size(hJ))
	for i = 1:size(Q, 1)
		for j = 1:size(Q, 2)
			if i == j
				Q[i, i] += 2 * hJ[i, i]
			else
				Q[i, i] -= 2 * (hJ[i, j] + hJ[j, i])
				Q[i, j] += 4 * hJ[i, j]
			end
		end
	end
	return Q
end

function qubo2ising(Q::AbstractDict)
	linearterms = Dict{Int, Float64}()
	j = Dict{Tuple{Int, Int}, Float64}()
	energyshift = 0.0
	for kv in Q
		k = kv[1]::Tuple{Int, Int}
		v = kv[2]::Float64
		m, n = min(k[1], k[2]), max(k[1], k[2])
		if m == n
			if haskey(linearterms, m)
				linearterms[m] += 0.5 * v
				energyshift += 0.5 * v
			else
				linearterms[m] = 0.5 * v
				energyshift += 0.5 * v
			end
		else
			j[(m, n)] = 0.25 * v
			energyshift += .25 * v
			if haskey(linearterms, m)
				linearterms[m] += .25 * v
			else
				linearterms[m] = .25 * v
			end
			if haskey(linearterms, n)
				linearterms[n] += .25 * v
			else
				linearterms[n] = .25 * v
			end
		end
	end
	h = zeros(maximum(keys(linearterms)) + 1)
	for kv in linearterms
		k = kv[1]
		v = kv[2]
		h[k + 1] = v
	end
	return h, j, energyshift
end

const savedvgcsandjcs = Dict()

function getverygoodcouplingsandjc(h, j, embedding, fadj)
	if haskey(savedvgcsandjcs, (embedding, fadj))
		return savedvgcsandjcs[(embedding, fadj)]
	end
	physbit2logbit = Dict{Int, Int}()
	for (i, physbits) in enumerate(embedding)
		for physbit in physbits
			physbit2logbit[physbit] = i - 1
		end
	end
	goodcouplings = Dict{Tuple{Int, Int}, Tuple{Int, Int}}()#these are the couplings we want the machine to have...maps from physical couplings to virtual couplings
	for k in keys(j)
		m, n = min(k...), max(k...)
		for b1 in embedding[m + 1]
			for b2 in embedding[n + 1]
				bl, bh = min(b1, b2), max(b1, b2)
				goodcouplings[(bl, bh)] = (m, n)
			end
		end
	end
	verygoodcouplings = Dict{Tuple{Int, Int}, Set{Tuple{Int, Int}}}()#these are the couplings we want the machine to have and which it actually has...maps from virtual couplings to array of physical couplings
	jc = Dict()
	for x in fadj
		m, n = min(x...), max(x...)
		if haskey(physbit2logbit, m) && haskey(physbit2logbit, n)
			if physbit2logbit[m] == physbit2logbit[n]
				jc[(m, n)] = -1.0
			end
		end
		if haskey(goodcouplings, (m, n))
			if haskey(verygoodcouplings, goodcouplings[(m, n)])
				push!(verygoodcouplings[goodcouplings[(m, n)]], (m, n))
			else
				verygoodcouplings[goodcouplings[(m, n)]] = Set{Tuple{Int, Int}}([(m, n)])
			end
		end
	end
	savedvgcsandjcs[(embedding, fadj)] = (verygoodcouplings, jc)
	return verygoodcouplings, jc
end

const savedfadjs = Dict()

function embed_problem(h, j, embedding, adjacency)
	if haskey(savedfadjs, adjacency)
		fadj = savedfadjs[adjacency]
	else
		fadj = convert(Array{Tuple{Int, Int}, 1}, collect(adjacency))
		savedfadjs[adjacency] = fadj
	end
	numbits = maximum(map(x->max(x...), fadj)) + 1
	newh = Any[0.0 for i in 1:numbits]
	newj = Dict()
	for i = 1:length(h)
		newh[embedding[i] .+ 1] .= h[i] / length(embedding[i])
	end
	verygoodcouplings, jc = getverygoodcouplingsandjc(h, j, embedding, fadj)
	for k in keys(verygoodcouplings)
		numcouplings = length(verygoodcouplings[k])
		for x in verygoodcouplings[k]
			newj[x] = j[k] / numcouplings
		end
	end
	return newh, newj, jc, embedding
end

function embedproblem(Q, embeddings, adjacency=defaultadjacency; param_chain=1)
	jparamchain = -param_chain
	h, j, energyshift = qubo2ising(Q)
	newh, newj, jc, newembeddings = embed_problem(h, j, embeddings, adjacency)
	for k in keys(jc)
		if jc[k] == -1
			newj[k] = jparamchain
		end
	end
	return newh, newj, newembeddings, energyshift
end

function unembedanswer(solutions, embeddings; broken_chains="weighted_random", kwargs...)
	if broken_chains != "weighted_random"
		validkws = [:broken_chains, :h, :j]
		solutionsarray = Array{Array{Any, 1}}(size(solutions, 1))
		for i = 1:length(solutionsarray)
			solutionsarray[i] = vec(solutions[i, :])
		end
		return dwembed.unembed_answer(solutionsarray, embeddings; broken_chains=broken_chains, validkwargs(kwargs, validkws)...)
	else
		unembeddedsolutions = Array{Int}(undef, size(solutions, 1), length(embeddings))
		for j = 1:length(embeddings)
			ks = rand(embeddings[j], size(solutions, 1))
			for i = 1:size(solutions, 1)
				unembeddedsolutions[i, j] = solutions[i, ks[i] + 1]
			end
		end
		return unembeddedsolutions
	end
end

function solvequbo(Q, solver=defaultsolver; fixvars=false, auto_scale=true, kwargs...)
	if fixvars
		fvresult = dwfv.fix_variables(Q)
		newQ = fvresult["new_Q"]
		if solver == defaultsolver
			if auto_scale
				error("can't auto_scale with the simulator")
			end
			solverresult = dwcore.solve_qubo(solver, newQ; auto_scale=auto_scale, kwargs...)
		else
			solverresult = dwcore.solve_qubo(solver, newQ; auto_scale=auto_scale, kwargs...)
		end
		numsolutions = size(solverresult["solutions"], 1)
		for i = 1:numsolutions
			solverresult["energies"][i] -= fvresult["offset"]
			for j in keys(fvresult["fixed_variables"])
				solverresult["solutions"][i, j + 1] = fvresult["fixed_variables"][j]
			end
		end
		return solverresult
	else
		if solver == defaultsolver
			if auto_scale
				error("can't auto_scale with the simulator")
			end
			return dwcore.solve_qubo(solver, Q; kwargs...)
		else
			return dwcore.solve_qubo(solver, Q; auto_scale=auto_scale, kwargs...)
		end
	end
end

function asyncsolveising(h, j, solver=defaultsolver; kwargs...)
	if solver == defaultsolver
		validkws = [:annealing_time, :answer_mode, :beta, :chains, :max_answers, :num_reads, :num_spin_reversal_transforms, :postprocess, :programming_thermalization, :readout_thermalization]
	else
		validkws = [:annealing_time, :answer_mode, :auto_scale, :beta, :chains, :max_answers, :num_reads, :num_spin_reversal_transforms, :postprocess, :programming_thermalization, :readout_thermalization]
	end
	p1 = dwcore.async_solve_ising(solver, h, j; validkwargs(kwargs, validkws)...)
	return p1
end

function solveising(h, j, solver=defaultsolver; kwargs...)
	p1 = asyncsolveising(h, j, solver; kwargs...)
	min_done = 1
	timeout = 60
	done = dwcore.await_completion([p1], min_done, timeout)
	if done
		return p1.result()
	else
		error("timed out awaiting solve_ising...or something: done=$done")
	end
end

end
