mutable struct Param
	name::Symbol
	value::Float64
end

function Param(s::AbstractString)
	Param(Symbol(s), NaN)
end

mutable struct Var{T}
	name::Symbol
	value::T
end

function Var(s::AbstractString, value=NaN)
	Var(Symbol(s), value)
end

mutable struct VarRef
	v::Var
	args::Tuple
end

abstract type Term end
abstract type AbstractLinearTerm <: Term end
abstract type AbstractQuadraticTerm <: Term end

mutable struct ConstantTerm <: Term
	realcoeff::Float64
end

mutable struct LinearTerm <: AbstractLinearTerm
	realcoeff::Float64
	var::Union{Var, VarRef}
end

mutable struct ParamLinearTerm <: AbstractLinearTerm
	realcoeff::Float64
	param::Param
	var::Union{Var, VarRef}
end

mutable struct QuadraticTerm <: AbstractQuadraticTerm
	realcoeff::Float64
	var1::Union{Var, VarRef}
	var2::Union{Var, VarRef}
	function QuadraticTerm(r, v1, v2)
		if v1 == v2
			return LinearTerm(r, v1)
		else
			return new(r, v1, v2)
		end
	end
end

mutable struct ParamQuadraticTerm <: AbstractQuadraticTerm
	realcoeff::Float64
	param::Param
	var1::Union{Var, VarRef}
	var2::Union{Var, VarRef}
	function ParamQuadraticTerm(r, p, v1, v2)
		if v1 == v2
			return ParamLinearTerm(r, p, v1)
		else
			return new(r, p, v1, v2)
		end
	end
end

function Term(args...)
	realcoeff = 1.
	numvars = 0
	vars = Union{Var, VarRef}[]
	params = Param[]
	for i = 1:length(args)
		if isa(args[i], Number)
			realcoeff *= args[i]
		elseif isa(args[i], Var)
			push!(vars, args[i])
		elseif isa(args[i], VarRef)
			push!(vars, args[i])
		elseif isa(args[i], Param)
			push!(params, args[i])
		else
			error("productant $(args[i]) of unknown type $(typeof(args[i]))")
		end
	end
	if length(params) == 0 && length(vars) == 0
		return ConstantTerm(realcoeff)
	elseif length(params) == 0 && length(vars) == 1
		return LinearTerm(realcoeff, vars[1])
	elseif length(params) == 1 && length(vars) == 1
		return ParamLinearTerm(realcoeff, params[1], vars[1])
	elseif length(params) == 0 && length(vars) == 2
		return QuadraticTerm(realcoeff, vars[1], vars[2])
	elseif length(params) == 1 && length(vars) == 2
		return ParamQuadraticTerm(realcoeff, params[1], vars[1], vars[2])
	else
		error("invalid term $(args)")
	end
end

mutable struct Model
	name::String
	connection::String
	solver::String
	workingdir::String
	workspace::String
	params::Array{Param, 1}
	vars::Array{Var, 1}
	terms::Array{Term, 1}
	bitsolutions::Array{Array{Int32, 1}, 1}
	energies::Array{Float64, 1}
	occurrences::Array{Int32, 1}
	valid::Array{Bool, 1}
	embedding::Dict
end

function Model(name, connection, solver, workingdir, workspace)
	return Model(name, connection, solver, workingdir, workspace, Array{Symbol}(undef, 0), Array{Symbol}(undef, 0), Array{Expr}(undef, 0), Array{Int32, 1}[], Float64[], Int32[], Bool[], Dict())
end
