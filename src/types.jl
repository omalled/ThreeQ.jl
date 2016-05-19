type Param
	name::Symbol
	value::Float64
end

function Param(s::AbstractString)
	Param(symbol(s), NaN)
end

type Var{T}
	name::Symbol
	value::T
end

function Var(s::AbstractString, value=NaN)
	Var(symbol(s), value)
end

type Term
	prodstrings::Array{ASCIIString, 1}
end

function Term(args...)
	newargs = Any[1.]
	numvars = 0
	for i = 1:length(args)
		if isa(args[i], Number)
			newargs[1] *= args[i]
		else
			push!(newargs, args[i])
		end
	end
	Term([map(string, newargs)...])
end

type Model
	name::ASCIIString
	connection::ASCIIString
	solver::ASCIIString
	workingdir::ASCIIString
	workspace::ASCIIString
	params::Array{Param, 1}
	vars::Array{Var, 1}
	terms::Array{Term, 1}
	bitsolutions::Array{Array{Int32, 1}, 1}
	energies::Array{Float64, 1}
	occurrences::Array{Int32, 1}
	valid::Array{Bool, 1}
	embedding::Dict
	numqubits::Int
end

function Model(name, connection, solver, workingdir, workspace)
	return Model(name, connection, solver, workingdir, workspace, Array(Symbol, 0), Array(Symbol, 0), Array(Expr, 0), Array{Int32, 1}[], Float64[], Int32[], Bool[], Dict(), -1)
end
