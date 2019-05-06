using Test
import PlotEmbedding
import ThreeQ
function getpostlabel(line)
	return split(line, ":")[2]
end
function removespaces(line)
	return replace(line, " "=>"")
end
function cleanline(line)
	return removespaces(getpostlabel(line))
end
function parse_dw_geometry(filename)
	lines = readlines(filename)
	L = parse(Int, cleanline(lines[1]))
	M = parse(Int, cleanline(lines[2]))
	N = parse(Int, cleanline(lines[3]))
	#figure out which qubits are there
	numqubits = 2 * L * M * N
	qubitlines = 4:4 + div(numqubits, 32) - 1
	qubitindicators = prod(cleanline.(lines[qubitlines]))
	hasqubit = map(x->x == '1', collect(qubitindicators))
	#figure out which couplers are there
	numcouplers = 5 * div(numqubits, 2) - L * M + div(numqubits, 2) - L * N
	couplerlines = maximum(qubitlines) + 1:length(lines)
	couplerindicators = prod(cleanline.(lines[couplerlines]))
	hascoupler = map(x->x == '1', collect(couplerindicators))
	return hasqubit, hascoupler, L, M, N
end
function couplernum2qubits(couplernum, L, M, N)
	if couplernum <= L * L * M * N#if it is an intracell coupler
		celli, cellj = PlotEmbedding.lin2chimera(1 + 2 * L * div(couplernum - 1, L ^ 2), M, N, L)[1:2]
		cellcouplernum = mod(couplernum - 1, L * L)
		qubit1_k = div(cellcouplernum, L) + 1
		qubit2_k = mod(cellcouplernum, L) + 5
		return PlotEmbedding.chimera2lin(celli, cellj, qubit1_k, M, N, L), PlotEmbedding.chimera2lin(celli, cellj, qubit2_k, M, N, L)
	elseif couplernum <= L * L * M * N + L * M * N - L * M#if it is a horizontal intercell coupler
		celli, cellj = PlotEmbedding.lin2chimera(1 + 2 * (couplernum - 1 - L * L * M * N), M, N - 1, L)#N - 1 because there are no couplers associated with the last column of unit cells
		k = mod(couplernum - 1, L) + 1
		return PlotEmbedding.chimera2lin(celli, cellj, L + k, M, N, L), PlotEmbedding.chimera2lin(celli, cellj + 1, L + k, M, N, L)
	elseif couplernum <= L * L * M * N + 2 * L * M * N - L * (M + N)#if it is a vertical intercell coupler
		celli, cellj = PlotEmbedding.lin2chimera(1 + 2 * (couplernum - 1 - L * L * M * N - L * M * N + L * M), M, N, L)
		k = mod(couplernum - 1, L) + 1
		return PlotEmbedding.chimera2lin(celli, cellj, k, M, N, L), PlotEmbedding.chimera2lin(celli + 1, cellj, k, M, N, L)
	else
		error("bad coupler number $couplernum")
	end
end

hasqubit, hascoupler, _ = parse_dw_geometry("geo/geometry-12x12-2")
@test length(hasqubit) == 1152
@test length(hascoupler) == 3360
@test !hasqubit[8]#qubit 8 is missing
@test !hascoupler[4] && !hascoupler[8] && !hascoupler[12] && !hascoupler[16] && !hascoupler[2308]#the couplers associated with qubit 8 are missing
@test couplernum2qubits(4, 4, 12, 12) == (1, 8)
@test couplernum2qubits(8, 4, 12, 12) == (2, 8)
@test couplernum2qubits(12, 4, 12, 12) == (3, 8)
@test couplernum2qubits(16, 4, 12, 12) == (4, 8)
@test couplernum2qubits(17, 4, 12, 12) == (9, 13)
@test couplernum2qubits(18, 4, 12, 12) == (9, 14)
@test couplernum2qubits(4 * 4 * 12 * 12, 4, 12, 12) == (2 * 4 * 12 * 12 - 4, 2 * 4 * 12 * 12)
#test a horizontal copuler
@test couplernum2qubits(2308, 4, 12, 12) == (8, 16)#this test is from the example denny gave
#test a vertical coupler
hasqubit, hascoupler, _ = parse_dw_geometry("geo/geometry-4x4-1")
@test couplernum2qubits(305, 4, 4, 4) == (1, 33)
@test couplernum2qubits(352, 4, 4, 4) == (92, 124)

function getadjacency(filename)
	hasqubit, hascoupler, L, M, N = parse_dw_geometry(filename)
	adj = Tuple{Int, Int}[]
	for i = 1:length(hascoupler)
		v1, v2 = couplernum2qubits(i, L, M, N)
		if hascoupler[i]
			push!(adj, (v1 - 1, v2 - 1))
			push!(adj, (v2 - 1, v1 - 1))
		end
	end
	return adj, M, N
end

adj, _ = getadjacency("geo/geometry-4x4-1")
simulator_adj = Set(ThreeQ.DWQMI.getadjacency(ThreeQ.DWQMI.defaultsolver))
@test length(Set(adj)) == length(simulator_adj)
for i = 1:length(adj)
	@test (adj[i][1], adj[i][2]) in simulator_adj
end

hasqubit, hascoupler, _ = parse_dw_geometry("geo/geometry-16x16-1")
adj, _ = getadjacency("geo/geometry-16x16-1")
@test length(Set(adj)) == 2 * sum(hascoupler)
