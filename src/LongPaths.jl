module LongPaths

import LightGraphs
import LongestPaths
import PlotEmbedding

function findpath(v1, v2, pairs)
	allbits = sort(unique([map(x->x[1], pairs); map(x->x[2], pairs)]))
	bigi2smalli = Dict(zip(allbits, 1:length(allbits)))
	smalli2bigi = Dict(zip(1:length(allbits), allbits))
	g = LightGraphs.DiGraph(length(allbits))
	for pair in pairs
		LightGraphs.add_edge!(g, bigi2smalli[pair[1]], bigi2smalli[pair[2]])
	end
	longest_path = LongestPaths.find_longest_path(g, bigi2smalli[v1], v2 == 0 ? 0 : bigi2smalli[v2])
	return map(x->smalli2bigi[x], longest_path.longest_path)
end

function get2by2group(i1, i2_unreversed, n, m, p, livebits; reverse=false)
	if mod(i1, 2) == 0 && reverse
		i2 = 1 + div(m, 2) - i2_unreversed#on the even rows, go backwards
	else
		i2 = i2_unreversed
	end
	bitsingroup = Int[]
	for j1 = 1:2
		for j2 = 1:2
			for k = 1:p
				lini = PlotEmbedding.chimera2lin(j1 + (i1 - 1) * 2, j2 + (i2 - 1) * 2, k, m, n, div(p, 2)) - 1 
				if lini in livebits
					append!(bitsingroup, lini)
				end
			end
		end
	end
	return bitsingroup
end

function turnright(direction)
	right = [1, 0]
	left = [-1, 0]
	up = [0, 1]
	down = [0, -1]
	if direction == right
		return down
	elseif direction == left
		return up
	elseif direction == up
		return right
	elseif direction == down
		return left
	else
		error("you're lost")
	end
end

function snakeit(adjacency, n, m, p)
	hn = div(n, 2)#half n
	hm = div(m, 2)#half m
	visited = fill(false, hn + 2, hm + 2)
	visited[1, :] .= true
	visited[end, :] .= true
	visited[:, 1] .= true
	visited[:, end] .= true
	visited[2, 2] = true
	stops = [[2, 2]]
	direction = [1, 0]
	while length(stops) < hn * hm
		nextstop = stops[end] + direction
		if !visited[nextstop...]
			push!(stops, nextstop)
			visited[nextstop...] = true
		else
			direction = turnright(direction)
		end
	end
	livebits = Set([map(x->x[1], adjacency); map(x->x[2], adjacency)])
	groups = map(stop->get2by2group(stop[1] - 1, stop[2] - 1, n, m, p, livebits), stops)
	return groups
end

function findconnector(adjacency, group1, group2)
	for v1 in reverse(group1)
		for v2 in group2
			if (v1, v2) in adjacency
				return v1, v2
			end
		end
	end
end

function coarsen(finegroups)
	extra = mod(length(finegroups), 2)
	coarsegroups = Array{Array{Int, 1}}(undef, div(length(finegroups), 2) + extra)
	if extra == 1
		coarsegroups[end] = finegroups[end]
	end
	for i = 1:length(coarsegroups) - extra
		coarsegroups[i] = vcat(finegroups[1 + 2 * (i - 1)], finegroups[2 + 2 * (i - 1)])
	end
	return coarsegroups
end

function findlongpath(adjacency, n, m, p=8; numcoarsens=0)#p=8 means that there are 8 qubits in a unit cell
	groups = snakeit(adjacency, n, m, p)
	for i = 1:numcoarsens
		if length(groups) > 1
			groups = coarsen(groups)
		end
	end
	connectors = Int[groups[1][1]]
	adj = Set(adjacency)
	for i = 1:length(groups) - 1
		append!(connectors, findconnector(adj, groups[i], groups[i + 1]))
	end
	#push!(connectors, minimum(groups[end]))
	push!(connectors, 0)#0 means it can end anywhere
	partialpaths = Any[]
	totallength = 0
	for i = 1:length(groups)
		thispath = findpath(connectors[2 * (i - 1) + 1], connectors[2 * (i - 1) + 2], filter(x->(x[1] in groups[i] && x[2] in groups[i]), adjacency))
		totallength += length(thispath)
		push!(partialpaths, thispath)
	end
	return vcat(partialpaths...)
end

end
