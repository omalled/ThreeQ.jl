macro defparam(model, x)
	@assert typeof(x) == Symbol
	stringx = string(x)
	code = :($(esc(x)) = Param($stringx))
	code = :($code; addparam!($(esc(model)), $(esc(x))))
	return code
end

macro defvar(model, x)
	if typeof(x) == Symbol
		stringx = string(x)
		code = :($(esc(x)) = Var($stringx))
		code = :($code; addvar!($(esc(model)), $(esc(x))))
		return code
	elseif typeof(x) == Expr && x.head == :ref
		xsym = x.args[1]
		stringx = string(xsym)
		dims = []
		for i = 2:length(x.args)
			@assert typeof(x.args[i]) == Expr
			@assert x.args[i].head == :call
			@assert x.args[i].args[1] == :(:)
			@assert x.args[i].args[2] == :(1)
			@assert length(x.args[i].args) == 3
			push!(dims, :($(esc(x.args[i].args[3]))))
		end
		code = :($(esc(xsym)) = Var($stringx, fill(NaN, $(dims...))))
		code = :($code; addvar!($(esc(model)), $(esc(xsym))))
		return code
	else
		error("variable should be something like x, x[1:7], x[1:3, 1:5], etc")
	end
end

macro addterm(model, termexpr)
	termexpr = macroexpand(Main, termexpr)
	if typeof(termexpr) == Expr
		@assert termexpr.head == :call
		@assert termexpr.args[1] == :*
		code = :(addterm!($(esc(model)), Term($(map(esc, termexpr.args[2:end])...))))
		return code
	elseif typeof(termexpr) == Symbol
		code = :(addterm!($(esc(model)), Term($(esc(termexpr)))))
		return code
	else
		error("incorrect term: $termexpr")
	end
end

macro addquadratic(model, quadexpr)
	quadexpr = macroexpand(Main, quadexpr)
	if typeof(quadexpr) != Expr || quadexpr.head != :call || quadexpr.args[1] != :^ || quadexpr.args[end] != 2
		error("$quadexpr is not a quadratic expression")
	else
		innerexpr = quadexpr.args[2]
		if typeof(innerexpr) == Symbol || isa(innerexpr, Number)
			code = :(addterm!($(esc(model)), Term($(esc(innerexpr)), $(esc(innerexpr)))))
			return code
		elseif typeof(innerexpr) == Expr && innerexpr.head == :call && innerexpr.args[1] == :+
			code = :()
			productants = Any[]
			for i = 2:length(innerexpr.args)
				if (typeof(innerexpr.args[i]) == Expr && innerexpr.args[i].head == :ref) || isa(innerexpr.args[i], Number) || isa(innerexpr.args[i], Symbol)
					push!(productants, Any[innerexpr.args[i]])
				elseif innerexpr.args[i].head == :call && innerexpr.args[i].args[1] == :*
					push!(productants, innerexpr.args[i].args[2:end])
				else
					error("unsupported term in quadratic: $(innerexpr.args[i])")
				end
			end
			for i = 1:length(productants)
				code = :($code; addterm!($(esc(model)), Term($(map(esc, [productants[i]; productants[i]])...))))
				for j = 1:i - 1
					code = :($code; addterm!($(esc(model)), Term($(map(esc, [2; productants[i]; productants[j]])...))))
				end
			end
			return code
		elseif typeof(innerexpr) == Expr && innerexpr.head == :call && innerexpr.args[1] == :*
			code = :(addterm!($(esc(model)), Term($(map(esc, [innerexpr.args[2:end]; innerexpr.args[2:end]])...))))
			return code
		else
			error("cannot handle inner part of quadratic expression: $innerexpr")
		end
	end
end

macro loadsolution(modelexpr, energy, occurrences, valid, solutionnumexpr)
	code = quote
		model = $(esc(modelexpr))
		solutionnum = $(esc(solutionnumexpr))
		$(esc(energy)) = model.energies[solutionnum]
		$(esc(occurrences)) = model.occurrences[solutionnum]
		$(esc(valid)) = model.valid[solutionnum]
		for fullvarname in keys(model.embedding)
			splitvarname = split(fullvarname, "___")
			for j = 1:length(model.vars)
				if splitvarname[1] == string(model.vars[j].name)
					value = model.bitsolutions[solutionnum][rand(model.embedding[fullvarname])]
					if length(splitvarname) == 1
						model.vars[j].value = value
					else
						indices = map(i->parse(Int, i), splitvarname[2:end])
						model.vars[j].value[indices...] = value
					end
				end
			end
		end
	end
	return code
end
