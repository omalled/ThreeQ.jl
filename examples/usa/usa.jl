using ThreeQ

function colormap(regions, neighbors, numcolors; w1=1, w2=1)
	colornames = ["red", "orange", "yellow", "green", "blue", "indigo", "violet"]
	m = ThreeQ.Model("mapcolor", "laptop", "c4-sw_sample", "workingdir", "c4")
	@defvar m colors[1:length(regions), 1:numcolors]
	for i = 1:length(regions)
		for j = 1:numcolors
			@addterm m -w1 * colors[i, j]
			for k = 1:j - 1
				@addterm m 2 * w1 * colors[i, j] * colors[i, k]
			end
		end
	end
	for j = 1:length(regions)
		for k = 1:j - 1
			if regions[k] in neighbors[regions[j]]
				for i = 1:numcolors
					@addterm m w2 * colors[j, i] * colors[k, i]
				end
			end
		end
	end
	ThreeQ.qbsolv!(m; minval=-length(regions) * w1, S=false, showoutput=true)
	result = Dict(zip(regions, map(i->colornames[indmax(vec(colors.value[i, :]))], 1:length(regions))))
	return result
end

function checkanswer(answer, regions, neighbors)
	if length(answer) < length(regions)
		println("Some regions have not a color!")
		return false
	end
	for r1 in regions
		for r2 in neighbors[r1]
			if answer[r1] == answer[r2]
				println("$r1 and $r2 have the same color!")
				return false
			end
		end
	end
	return true
end

function parseusa(filename)
	f = open(filename)
	lines = readlines(f)
	close(f)
	i = 1
	county_neighbors = Dict()
	while i <= length(lines)
		splitline = split(lines[i], "\t")
		thiscounty = splitline[1][2:end-1]
		firstneighbor = splitline[3][2:end-1]
		theseneighbors = Any[firstneighbor]
		i += 1
		while i <= length(lines) && startswith(lines[i], "\t\t")
			splitline = split(lines[i], "\t")
			thisneighbor = splitline[3][2:end-1]
			push!(theseneighbors, thisneighbor)
			i += 1
		end
		county_neighbors[thiscounty] = Any[]
		for n in theseneighbors
			if n != thiscounty
				push!(county_neighbors[thiscounty], n)
			end
		end
	end
	return collect(keys(county_neighbors)), county_neighbors
end

provinces = ["BC", "YK", "NW", "AB", "SK", "NV", "MT", "ON", "QB", "NB", "NS", "PE", "NL"]
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

answer = colormap(provinces, neighbors, 3)
@show checkanswer(answer, provinces, neighbors)

counties, county_neighbors = parseusa("county_adjacency.txt")

nm_counties = filter(x->endswith(x, ", NM"), counties)
nm_neighbors = Dict(zip(nm_counties, map(c->filter(x->endswith(x, ", NM"), county_neighbors[c]), nm_counties)))
nm_answer = colormap(nm_counties, nm_neighbors, 4)
@show checkanswer(nm_answer, nm_counties, nm_neighbors)

@time usa_answer = colormap(counties, county_neighbors, 4)
checkanswer(usa_answer, counties, county_neighbors)
