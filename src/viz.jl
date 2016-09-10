import GraphPlot
import Graphs
using Colors

function adjacencygraph(adjacency)
	legitbits = unique([map(x->x[1], adjacency); map(x->x[2], adjacency)])
	bnum2i = Dict(zip(legitbits, 1:length(legitbits)))
	g = Graphs.simple_graph(length(legitbits), is_directed=false)
	for (i1, i2) in adjacency
		Graphs.add_edge!(g, bnum2i[i1], bnum2i[i2])
	end
	return g
end

function plotadjacency(adjacency)
	g = adjacencygraph(adjacency)
	return GraphPlot.gplot(g)
end

function plotembedding(model, adjacency)
	legitbits = unique([map(x->x[1], adjacency); map(x->x[2], adjacency)])
	pbnum2i = Dict(zip(legitbits, 1:length(legitbits)))
	lbnum2i = Dict(zip(1:length(model.embedding), length(legitbits) + 1:length(model.embedding) + length(legitbits)))
	g = Graphs.simple_graph(length(legitbits) + length(model.embedding), is_directed=false)
	for (i1, i2) in adjacency
		Graphs.add_edge!(g, pbnum2i[i1], pbnum2i[i2])
	end
	embeddings = collect(values(model.embedding))
	for i = 1:length(model.embedding)
		for j = 1:length(embeddings[i])
			Graphs.add_edge!(g, lbnum2i[i], pbnum2i[embeddings[i][j] - 1])
		end
	end
	embeddedbits = vcat(embeddings...) - 1
	threecolors = [colorant"lightseagreen", colorant"black", colorant"orange"]
	plotcolors = [map(x->x in embeddedbits ? threecolors[1] : threecolors[2], legitbits); map(x->threecolors[3], 1:length(model.embedding))]
	return GraphPlot.gplot(g, nodefillc=plotcolors)
end
