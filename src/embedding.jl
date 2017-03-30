function stringembed(numbits, adjacency)
	livebits = Set()
	for (i, j) in adjacency
		push!(livebits, i)
		push!(livebits, j)
	end
	neighbors = Dict{Int, Array{Int, 1}}()
	for i in livebits
		neighbors[i] = Int[]
	end
	for (i, j) in adjacency
		push!(neighbors[i], j)
		push!(neighbors[j], i)
	end
	return map(i->Int[i], stringembed(numbits, adjacency, livebits, neighbors))
end

function stringembed(numbits, adjacency, availablebits, neighbors)
	for i in availablebits
		try
			pop!(availablebits, i)
			return Int[Int[i]; stringembed(i, numbits - 1, adjacency, availablebits, neighbors)]
		catch
			push!(availablebits, i)
		end
	end
	error("no embedding found")
end

function stringembed(lastbit, numbits, adjacency, availablebits, neighbors)
	for i in neighbors[lastbit]
		if i in availablebits
			try
				pop!(availablebits, i)
				if numbits == 1
					return Int[i]
				else
					return Int[Int[i]; stringembed(i, numbits - 1, adjacency, availablebits, neighbors)]
				end
			catch
				push!(availablebits, i)
			end
		end
	end
	error()
end
