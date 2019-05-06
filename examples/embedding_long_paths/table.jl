using ThreeQ
import JLD
import PyPlot
include("parse_dw_geometry.jl")

if !isfile("table.jld")
	coarsenlevels = [0, 1, 2, 3]
	filenames = ["geo/geometry-4x4-1", "geo/geometry-8x8-1", "geo/geometry-8x8-2", "geo/geometry-12x12-1", "geo/geometry-12x12-2", "geo/geometry-12x12-3", "geo/geometry-12x12-4", "geo/geometry-12x12-5", "geo/geometry-12x12-6", "geo/geometry-16x16-1", "geo/geometry-16x16-2", "geo/geometry-16x16-3", "geo/geometry-16x16-4", "geo/geometry-16x16-5", "geo/geometry-16x16-6", "geo/geometry-16x16-7", "geo/geometry-16x16-8"]
	#filenames = ["geo/geometry-4x4-1", "geo/geometry-8x8-1", "geo/geometry-8x8-2", "geo/geometry-12x12-1"]

	timings = Array{Float64}(undef, length(filenames), length(coarsenlevels))
	numused = Array{Float64}(undef, length(filenames), length(coarsenlevels))
	numavail = Array{Float64}(undef, length(filenames), length(coarsenlevels))
	couplersavail = Array{Float64}(undef, length(filenames), length(coarsenlevels))
	for (i, filename) in enumerate(filenames)
		hasqubit, hascoupler, _ = parse_dw_geometry(filename)
		for (j, numcoarsens) in enumerate(coarsenlevels)
			@show filename, numcoarsens
			adj, m, n = getadjacency(filename)
			t = @elapsed embeddingbits = LongPaths.findlongpath(adj, m, n; numcoarsens=numcoarsens)
			timings[i, j] = t
			numused[i, j] = length(embeddingbits)
			numavail[i, j] = sum(hasqubit)
			couplersavail[i, j] = sum(hascoupler)
		end
	end

	@JLD.save "table.jld" filenames coarsenlevels timings numused numavail couplersavail
end


@JLD.load "table.jld" filenames coarsenlevels timings numused numavail couplersavail
#=
println("File name\tCoarsenings\tTime\tFraction qubits used\tQubits used\tQubits available\tCouplers available")
k = 0
for (i, filename) in enumerate(filenames)
	for(j, numcoarsens) in enumerate(coarsenlevels)
		global k
		k += 1
		println("$filename\t$numcoarsens\t$(timings[i, j])\t$(numused[i, j] / numavail[i, j])\t$(numused[i, j])\t$(numavail[i, j])\t$(couplersavail[i, j])")
	end
end
@show k
=#
println("\\begin{longtable}{|l|r|r|r|r|r|r|}")
println("\\hline")
println("Graph & \\begin{tabular}{@{}c@{}}Number of\\\\Coarsenings\\end{tabular} & Time (s) & \\begin{tabular}{@{}c@{}}Fraction of\\\\qubits used\\end{tabular} & \\begin{tabular}{@{}c@{}}Qubits\\\\used\\end{tabular} & \\begin{tabular}{@{}c@{}}Qubits\\\\available\\end{tabular} & \\begin{tabular}{@{}c@{}}Couplers\\\\available\\end{tabular} \\\\")
println("\\hline")
for (i, filename) in enumerate(filenames)
	#if !(filename in ["geo/geometry-4x4-1", "geo/geometry-8x8-1", "geo/geometry-12x12-1", "geo/geometry-16x16-3"])
		for(j, numcoarsens) in enumerate(coarsenlevels)
			println("$(filename[14:end]) & $numcoarsens & $(round(timings[i, j]; digits=1)) & $(round(numused[i, j] / numavail[i, j]; digits=3)) & $(Int(numused[i, j])) & $(Int(numavail[i, j])) & $(Int(couplersavail[i, j])) \\\\")
		end
	#end
end
println("\\hline")
println("\\end{longtable}")
