module DWQMI

import PyCall

@PyCall.pyimport dwave_sapi2.local as dwlocal
@PyCall.pyimport dwave_sapi2.core as dwcore
@PyCall.pyimport dwave_sapi2.fix_variables as dwfv
@PyCall.pyimport dwave_sapi2.util as dwutil
@PyCall.pyimport dwave_sapi2.embedding as dwembed
@PyCall.pyimport dwave_sapi2.remote as dwremote

include("utils.jl")

const defaultsolver = dwlocal.local_connection[:get_solver]("c4-sw_sample")
const defaultadjacency = collect(dwutil.get_hardware_adjacency(defaultsolver))
const defaulturl = "https://dw2x.dwavesys.com/sapi/"

function connect(token, url=defaulturl)
	return dwremote.RemoteConnection(url, token)
end

function getadjacency(solver)
	return dwutil.get_hardware_adjacency(solver)
end

function getremotesolvers(connection)
	return connection[:solver_names]()
end

function getremotesolver(connection, name)
	return connection[:get_solver](name)
end

function getdw2xsys4(token)
	connection = connect(token)
	return getremotesolver(connection, "DW2X_SYS4")
end

function findembeddings(Q, adjacency=defaultadjacency)
	return dwembed.find_embedding(collect(keys(Q)), adjacency, verbose=0, tries=100, timeout=300)
end

function qubo2ising(Q)
	h, j, energyshift = dwutil.qubo_to_ising(Q)
	return h, j, energyshift
end

function embedproblem(Q, embeddings, adjacency=defaultadjacency; param_chain=1)
	jparamchain = -param_chain
	h, j, energyshift = qubo2ising(Q)
	newh, newj, jc, newembeddings = dwembed.embed_problem(h, j, embeddings, adjacency)
	for k in keys(jc)
		if jc[k] == -1
			newj[k] = jparamchain
		end
	end
	return newh, newj, newembeddings, energyshift
end

function unembedanswer(solutions, embeddings; broken_chains="weighted_random", kwargs...)
	validkws = [:broken_chains, :h, :j]
	solutionsarray = Array(Array{Any, 1}, size(solutions, 1))
	for i = 1:length(solutionsarray)
		solutionsarray[i] = vec(solutions[i, :])
	end
	return dwembed.unembed_answer(solutionsarray, embeddings; broken_chains=broken_chains, validkwargs(kwargs, validkws)...)
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

function solveising(h, j, solver=defaultsolver; kwargs...)
	if solver == defaultsolver
		validkws = [:annealing_time, :answer_mode, :beta, :chains, :max_answers, :num_reads, :num_spin_reversal_transforms, :postprocess, :programming_thermalization, :readout_thermalization]
	else
		validkws = [:annealing_time, :answer_mode, :auto_scale, :beta, :chains, :max_answers, :num_reads, :num_spin_reversal_transforms, :postprocess, :programming_thermalization, :readout_thermalization]
	end
	return dwcore.solve_ising(solver, h, j; validkwargs(kwargs, validkws)...)
end

end
