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
	Term([map(string, args)...])
end

type Model
	name::ASCIIString
	connection::ASCIIString
	solver::ASCIIString
	workingdir::ASCIIString
	params::Array{Param, 1}
	vars::Array{Var, 1}
	terms::Array{Term, 1}
end

function Model(name, connection, solver, workingdir)
	return Model(name, connection, solver, workingdir, Array(Symbol, 0), Array(Symbol, 0), Array(Expr, 0))
end
