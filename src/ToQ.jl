module ToQ

export @defparam, @defvar, @addterm, @loadsolution

import Base.getindex
import Base.string

include("types.jl")
include("macros.jl")

function string(p::Param)
	return string(p.name)
end

function getindex(v::Var, args...)
	return string(v.name, "___", join(args, "___"))
end

function getindex{T<:Number}(v::Var{T})
	return string(v.name)
end

function addparam!(model, x)
	push!(model.params, x)
end

function addvar!(model, x)
	push!(model.vars, x)
end

function addterm!(model, term)
	push!(model.terms, term)
end

function assemble(x::Var)
	return "var: $(string(x.name))\n"
end

function assemble{T<:Vector}(x::Var{T})
	strings = ASCIIString[]
	for i = 1:length(x.value)
		push!(strings, "var: $(x[i])")
	end
	return join(strings, "\n") * "\n"
end

function assemble{T<:Matrix}(x::Var{T})
	strings = ASCIIString[]
	for i = 1:size(x.value, 1)
		for j = 1:size(x.value, 2)
			push!(strings, "var: $(x[i, j])")
		end
	end
	return join(strings, "\n") * "\n"
end

function writeqfile(m::Model, filename)
	f = open(filename, "w")
	for x in m.params
		write(f, "param: $(string(x.name))\n")
	end
	for x in m.vars
		write(f, assemble(x))
	end
	for x in m.terms
		write(f, "term: $(join(x.prodstrings, " * "))\n")
	end
	close(f)
end

function writebfile(args, filename)
	f = open(filename, "w")
	for t in args
		write(f, "$(t[1]) = $(t[2])\n")
	end
	close(f)
end

function readsol(filename)
	f = open(filename)
	resultlen = read(f, Int32)
	numsolutions = read(f, Int32)
	solutionlen = div(resultlen, numsolutions)
	solutions = Array(Array{Int32, 1}, numsolutions)
	energies = Array(Float64, numsolutions)
	occurrences = Array(Int32, numsolutions)
	for i = 1:numsolutions
		solutions[i] = read(f, Int32, solutionlen)
	end
	energies = read(f, Float64, numsolutions)
	occurrences = read(f, Int32, numsolutions)
	close(f)
	return solutions, energies, occurrences
end

function parseembedding(filename)
	f = open(filename)
	lines = readlines(f)
	close(f)
	embedding = Dict()
	for line in lines
		if startswith(line, "VAR")
			paramname = split(chomp(line), "\"")[2]
			qsstring = split(chomp(line), ":")[2]
			qstrings = split(qsstring)
			embedding[paramname] = map(x->parse(Int, x[2:end]) + 1, qstrings)
		end
	end
	return embedding
end

function solve!(m::Model; doembed=true, removefiles=false, numreads=10, args...)
	writeqfile(m, m.name * ".q")
	writebfile(args, m.name * ".b")
	bashscriptlines = ASCIIString[]
	push!(bashscriptlines, "set -e")
	push!(bashscriptlines, "dw set connection $(m.connection)")
	push!(bashscriptlines, "dw set solver $(m.solver)")
	push!(bashscriptlines, "dw mkdir -p $(m.workingdir)")
	push!(bashscriptlines, "dw cd $(m.workingdir)")
	if doembed
		push!(bashscriptlines, "dw embed $(m.name).q -o $(m.name).epqmi")
	end
	push!(bashscriptlines, "dw cp $(m.name).epqmi .epqmi")
	push!(bashscriptlines, "dw get embedding >$(m.name).embedding_info")
	push!(bashscriptlines, "dw bind $(m.name).b")
	push!(bashscriptlines, "dw exec num_reads=$numreads $(m.name).qmi")
	push!(bashscriptlines, "rm -f $(m.name).sol_*")
	push!(bashscriptlines, "dw val $(m.name).sol")
	if removefiles
		push!(bashscriptlines, "rm $(m.name).q")
		push!(bashscriptlines, "rm $(m.name).b")
	end
	bashscript = open(m.name * ".bash", "w")
	for line in bashscriptlines
		write(bashscript, string(line, "\n"))
	end
	close(bashscript)
	run(`bash $(m.name).bash`)
	m.embedding = parseembedding("$(m.name).embedding_info")
	m.bitsolutions, m.energies, m.occurrences = readsol(joinpath(ENV["DWAVE_HOME"], string("workspace.", m.workspace), m.workingdir, string(m.name, ".sol")))
	m.valid = Array(Bool, length(m.bitsolutions))
	for j = 1:length(m.valid)
		m.valid[j] = true
		for k in keys(m.embedding)
			firstbit = m.bitsolutions[j][m.embedding[k][1]]
			for i in m.embedding[k]
				thisbit = m.bitsolutions[j][i]
				if firstbit != thisbit || (thisbit != 1 && thisbit != 0)
					m.valid[j] = false
				end
			end
		end
	end
	if removefiles
		run(`rm $(m.name).bash`)
		run(`rm $(m.name).embedding_info`)
	end
end

function getnumsolutions(model)
	return length(model.bitsolutions)
end

end
