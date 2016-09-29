module ToQ

export @defparam, @defvar, @addterm, @addquadratic, @loadsolution, DWQMI

import Base.getindex
import Base.string
import Base.length
import Base.size
import Base.==

include("DWQMI.jl")
include("types.jl")
include("macros.jl")
include("viz.jl")

size(v::Var, args...) = size(v.value, args...)
length(v::Var) = length(v.value)

function ==(x::VarRef, y::VarRef)
	return x.v == y.v && x.args == y.args
end

function getindex{T}(v::Var{T}, args...)
	if ndims(T) != length(args)
		error("cannot access var $(v.name) with $(length(args)) dimensions when the dimensions are $(ndims(T))")
	end
	return VarRef(v, args)
end

function string(p::Param)
	return string(p.name)
end

function string(v::Var)
	return string(v.name)
end

function string(vr::VarRef)
	return string(string(vr.v), "___", join(vr.args, "___"))
end

function string(t::LinearTerm)
	return string(join(map(string, [t.realcoeff, t.var]), " * "))
end

function string(t::ParamLinearTerm)
	return string(join(map(string, [t.realcoeff, t.param, t.var]), " * "))
end

function string(t::QuadraticTerm)
	return string(join(map(string, [t.realcoeff, t.var1, t.var2]), " * "))
end

function string(t::ParamQuadraticTerm)
	return string(join(map(string, [t.realcoeff, t.param, t.var1, t.var2]), " * "))
end

function addparam!(model, x)
	push!(model.params, x)
end

function addvar!(model, x)
	push!(model.vars, x)
end

function addterm!(model, term::ConstantTerm)
	#do nothing for constant terms
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
		push!(strings, "var: $(string(x[i]))")
	end
	return join(strings, "\n") * "\n"
end

function assemble{T<:Matrix}(x::Var{T})
	strings = ASCIIString[]
	for i = 1:size(x.value, 1)
		for j = 1:size(x.value, 2)
			push!(strings, "var: $(string(x[i, j]))")
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
		write(f, "term: $(string(x))\n")
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

function varset(t::LinearTerm)
	return Set(map(string, [t.var]))
end

function varset(t::ParamLinearTerm)
	return Set(map(string, [t.param, t.var]))
end

function varset(t::QuadraticTerm)
	return Set(map(string, [t.var1, t.var2]))
end

function varset(t::ParamQuadraticTerm)
	return Set(map(string, [t.param, t.var1, t.var2]))
end

function varset(t::ConstantTerm)
	return Set()
end

function collectterms!(m::Model, paramdict::Associative)
	newterms = Term[]
	for term in m.terms
		if isa(term, ParamLinearTerm)
			push!(newterms, LinearTerm(term.realcoeff * paramdict[term.param.name], term.var))
		elseif isa(term, ParamQuadraticTerm)
			push!(newterms, QuadraticTerm(term.realcoeff * paramdict[term.param.name], term.var1, term.var2))
		else
			push!(newterms, term)
		end
	end
	m.terms = newterms
	collectterms!(m)
end

function collectterms!(m::Model)
	termdict = Dict()
	for term in m.terms
		s = varset(term)
		if haskey(termdict, s)
			push!(termdict[s], term)
		else
			termdict[s] = Term[term]
		end
	end
	newterms = Term[]
	for k in keys(termdict)
		if k != Set()#skip constant terms
			newterm = termdict[k][1]
			for i = 2:length(termdict[k])
				newterm.realcoeff += termdict[k][i].realcoeff
			end
			push!(newterms, newterm)
		end
	end
	m.terms = newterms
end

function model2ab(m::Model)
	diagterms = AbstractLinearTerm[]
	offdiagterms = AbstractQuadraticTerm[]
	varstring2index = Dict{Any, Int}()
	i2varstring = Dict{Int, AbstractString}()
	var2i(var) = varstring2index[string(var)]
	varindex = 1
	for term in m.terms
		if isa(term, AbstractLinearTerm)
			push!(diagterms, term)
			if !haskey(varstring2index, string(term.var))
				varstring2index[string(term.var)] = varindex
				i2varstring[varindex] = string(term.var)
				varindex += 1
			end
		elseif isa(term, AbstractQuadraticTerm)
			push!(offdiagterms, term)
			if !haskey(varstring2index, string(term.var1))
				varstring2index[string(term.var1)] = varindex
				i2varstring[varindex] = string(term.var1)
				varindex += 1
			end
			if !haskey(varstring2index, string(term.var2))
				varstring2index[string(term.var2)] = varindex
				i2varstring[varindex] = string(term.var2)
				varindex += 1
			end
		end
	end
	numvars = length(varstring2index)
	a = zeros(numvars)
	b = zeros(numvars, numvars)
	for term in m.terms
		if isa(term, AbstractLinearTerm)
			a[var2i(term.var)] += term.realcoeff
		elseif isa(term, AbstractQuadraticTerm)
			i = min(var2i(term.var1), var2i(term.var2))
			j = max(var2i(term.var1), var2i(term.var2))
			b[i, j] += term.realcoeff
		end
	end
	return a, b, i2varstring, numvars
end

function writeqbsolvfile(m::Model, filename)
	a, b, i2varstring, numvars = model2ab(m)
	f = open(filename, "w")
	write(f, "p qubo 0 $(numvars) $(sum(a .!= 0)) $(sum(b .!= 0))\n")
	for i = 1:length(a)
		if a[i] != 0
			inn = i - 1
			write(f, "c $(i2varstring[i])\n")
			write(f, "$inn $inn $(a[i])\n")
		end
	end
	for i = 1:size(b, 1)
		for j= 1:i - 1
			if b[j, i] != 0
				inn = i - 1
				jnn = j - 1
				write(f, "c $(i2varstring[j]) * $(i2varstring[i])\n")
				write(f, "$jnn $inn $(b[j, i])\n")
			end
		end
	end
	close(f)
	return i2varstring
end

function model2sapiQ(m::Model)
	a, b, i2varstring, numvars = model2ab(m)
	Q = Dict()
	for i = 1:length(a)
		if a[i] != 0
			inn = i - 1
			Q[(inn, inn)] = a[i]
		end
	end
	for i = 1:size(b, 1)
		for j= 1:i - 1
			if b[j, i] != 0
				inn = i - 1
				jnn = j - 1
				Q[(jnn, inn)] = b[j, i]
			end
		end
	end
	return Q, i2varstring, numvars
end

function rescale(h, j, maxh, maxj)
	j = deepcopy(j)
	hfactor = maxh / maximum(h)
	jfactor = maxj / maximum(map(abs, values(j)))
	goodfactor = min(hfactor, jfactor)
	h *= goodfactor
	for key in keys(j)
		j[key] *= goodfactor
	end
	return h, j
end

function solvesapi!(m::Model, maxh=2, maxj=1; solver=DWQMI.defaultsolver, adjacency=DWQMI.getadjacency(solver), param_chain=1, auto_scale=false, kwargs...)
	paramdict = Dict(kwargs)
	collectterms!(m, paramdict)
	Q, i2varstring, numvars = model2sapiQ(m)
	embeddings = DWQMI.findembeddings(Q, adjacency)
	if length(embeddings) == 0
		error("embedding failed")
	end
	h, j, newembeddings, energyshift = DWQMI.embedproblem(Q, embeddings, adjacency; param_chain=param_chain)
	if auto_scale
		h, j = rescale(h, j, maxh, maxj)
	end
	embeddedanswer = DWQMI.solveising(h, j, solver; auto_scale=false, kwargs...)
	embeddedqubosolutions = map(x->x == -1 ? 0 : x == 1 ? 1 : 3, embeddedanswer["solutions"])
	unembeddedisingsolutions = DWQMI.unembedanswer(embeddedanswer["solutions"], newembeddings; kwargs...)
	m.embedding = Dict(zip(map(i->i2varstring[i], 1:size(unembeddedisingsolutions, 2)), map(i->[newembeddings[i]...] + 1, 1:size(unembeddedisingsolutions, 2))))
	m.bitsolutions = map(i->vec(embeddedqubosolutions[i, :]), 1:size(embeddedqubosolutions, 1))
	m.energies = embeddedanswer["energies"]
	m.occurrences = map(x->convert(Int32, x), embeddedanswer["num_occurrences"])
	fillvalid!(m)
end

function qbsolv!(m::Model; minval=false, S=0, showoutput=false, paramvals...)
	paramdict = Dict(paramvals)
	collectterms!(m, paramdict)
	i2varstring = writeqbsolvfile(m, m.name * ".qbsolvin")
	if minval != false
		targetstring = "-T $minval"
	else
		targetstring = ""
	end
	if is(S, false)
		Sstring = ""
	else
		Sstring = "-S$S"
	end
	qbsolvcommand = `bash -c "dw set connection $(m.connection); dw set solver $(m.solver); qbsolv -i $(m.name * ".qbsolvin") $Sstring $targetstring -v1"`
	output = readlines(qbsolvcommand)
	solutionline = length(output)
	while !contains(output[solutionline - 1], "Number of bits in solution")
		solutionline -= 1
	end
	if showoutput
		for line in output
			print(line)
		end
	end
	bitsolution = map(b->parse(Int, b), collect(output[solutionline][1:end - 1]))
	for i = 1:length(i2varstring)
		fullvarname = i2varstring[i]
		splitvarname = split(fullvarname, "___")
		for j = 1:length(m.vars)
			if splitvarname[1] == string(m.vars[j].name)
				value = bitsolution[i]
				if length(splitvarname) == 1
					m.vars[j].value = value
				else
					indices = map(i->parse(Int, i), splitvarname[2:end])
					m.vars[j].value[indices...] = value
				end
			end
		end
	end
end

function fillvalid!(m::Model)
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
end


function solve!(m::Model; doembed=true, removefiles=false, numreads=10, args...)
	collectterms!(m)
	writeqfile(m, m.name * ".q")
	writebfile(args, m.name * ".b")
	bashscriptlines = ASCIIString[]
	push!(bashscriptlines, "set -e")
	push!(bashscriptlines, "dw set connection $(m.connection)")
	push!(bashscriptlines, "dw set solver $(m.solver)")
	push!(bashscriptlines, "dw mkdir -p $(m.workingdir)")
	push!(bashscriptlines, "dw cd $(m.workingdir)")
	if doembed
		push!(bashscriptlines, "dw embed $(m.name).q -o $(m.name).epqmi 2>$(m.name).embedding_info")
		push!(bashscriptlines, "cat $(m.name).embedding_info")
	end
	push!(bashscriptlines, "dw cp $(m.name).epqmi .epqmi")
	push!(bashscriptlines, "dw get embedding >>$(m.name).embedding_info")
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
	fillvalid!(m)
	if removefiles
		run(`rm $(m.name).bash`)
		run(`rm $(m.name).embedding_info`)
	end
end

function getnumsolutions(model)
	return length(model.bitsolutions)
end

function evalvar(var, vardict::Associative)
	return vardict[var.name]
end

function evalvar(varref::VarRef, vardict::Associative)
	return vardict[varref.v.name][varref.args...]
end

function evalqubo!(m::Model; kwargs...)
	kwdict = Dict(kwargs)
	collectterms!(m, kwdict)
	val = 0.
	for term in m.terms
		if isa(term, LinearTerm)
			val += term.realcoeff * evalvar(term.var, kwdict)
		elseif isa(term, QuadraticTerm)
			val += term.realcoeff * evalvar(term.var1, kwdict) * evalvar(term.var2, kwdict)
		else
			error("innapproprate term type: $(typeof(term))")
		end
	end
	return val
end

end
