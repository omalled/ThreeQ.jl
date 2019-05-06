using ThreeQ
import PyPlot
include("parse_dw_geometry.jl")

#filename = "geo/geometry-4x4-1"
filename = "geo/geometry-16x16-1"
hasqubit, hascoupler, _ = parse_dw_geometry(filename)
adj, m, n = getadjacency(filename)
@time embeddingbits = LongPaths.findlongpath(adj, m, n; numcoarsens=3)
livebits = Set([map(x->x[1] + 1, adj); map(x->x[2] + 1, adj)])

#plot the results
p = 4
c = PlotEmbedding.Chimera(PlotEmbedding.Cross(), 0.8, 0.8, 0.1, 0.1)
dwi2myi = Array{Int}(undef, m * n * p * 2)
livebits = sort(unique([map(x->x[1], adj); map(x->x[2], adj)]))
usededges = Pair{Int, Int}[]
for i = 1:length(embeddingbits) - 1
	push!(usededges, embeddingbits[i] + 1=>embeddingbits[i + 1] + 1)
end

for i = 1:length(livebits)
	dwi2myi[livebits[i] + 1] = i
end
bottom = -0.75
xs, ys = PlotEmbedding.getchimeraxys(m, n, p, c)
fig, ax = PyPlot.subplots(figsize=(12, 12))
PlotEmbedding.plotgraph(fig, ax, xs[livebits .+ 1], ys[livebits .+ 1], (map(x->dwi2myi[x[1]]=>dwi2myi[x[2]], usededges), "#5a9bd4"); edgealpha=0.8, markeredgewidth=0.5)
display(fig)
fig.savefig("longgraph.png")
println()
PyPlot.close(fig)
