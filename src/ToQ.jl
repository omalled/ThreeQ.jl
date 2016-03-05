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

function solve!(m::Model; doembed=true, removefiles=false, numreads=10, args...)
	writeqfile(m, m.name * ".q")
	writebfile(args, m.name * ".b")
	bashscriptlines = ASCIIString[]
	push!(bashscriptlines, "dw set connection $(m.connection)")
	push!(bashscriptlines, "dw set solver $(m.solver)")
	push!(bashscriptlines, "dw mkdir -p $(m.workingdir)")
	push!(bashscriptlines, "dw cd $(m.workingdir)")
	if doembed
		push!(bashscriptlines, "dw embed $(m.name).q -o $(m.name).epqmi")
	end
	push!(bashscriptlines, "dw cp $(m.name).epqmi .epqmi")
	push!(bashscriptlines, "dw bind $(m.name).b")
	push!(bashscriptlines, "dw exec num_reads=$numreads $(m.name).qmi")
	push!(bashscriptlines, "rm -f $(m.name).sol_*")
	push!(bashscriptlines, "dw val $(m.name).sol")
	push!(bashscriptlines, "numsols=\"\$(dw val $(m.name).sol | grep \"DISTINCT SOLUTIONS\" | awk '{print \$4;}')\"")
	push!(bashscriptlines, "for i in `seq 1 \$numsols`; do")
	push!(bashscriptlines, "\tdw val -s \$i $(m.name).sol >$(m.name).sol_\$i")
	push!(bashscriptlines, "done")
	if removefiles
		push!(bashscriptlines, "rm $(m.name).q")
		push!(bashscriptlines, "rm $(m.name).b")
		push!(bashscriptlines, "rm $(m.name).epqmi")
	end
	bashscript = open(m.name * ".bash", "w")
	for line in bashscriptlines
		write(bashscript, string(line, "\n"))
	end
	close(bashscript)
	run(`bash $(m.name).bash`)
	if removefiles
		run(`rm $(m.name).bash`)
	end
end

end
