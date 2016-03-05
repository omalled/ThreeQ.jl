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
			@assert x.args[i].head == :(:)
			@assert length(x.args[i].args) == 2
			@assert x.args[i].args[1] == 1
			push!(dims, :($(esc(x.args[i].args[2]))))
		end
		code = :($(esc(xsym)) = Var($stringx, fill(NaN, $(dims...))))
		code = :($code; addvar!($(esc(model)), $(esc(xsym))))
		return code
	else
		error("variable should be something like x, x[1:7], x[1:3, 1:5], etc")
	end
end

macro addterm(model, termexpr)
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

macro loadsolution(model, energy, occurrences, solutionnum)
	code = quote
		filename = string($(esc(model)).name, ".sol_", $(esc(solutionnum)))
		f = open(filename)
		lines = readlines(f)
		$(esc(energy)) = parse(Float64, split(lines[4], " = ")[2])
		$(esc(occurrences)) = parse(Int, split(lines[5], " = ")[2])
		for i = 6:length(lines) - 1
			splitline = split(lines[i], " <== ")
			fullvarname = splitline[1]
			splitvarname = split(fullvarname, "___")
			modelvars = $(esc(model)).vars
			for j = 1:length(modelvars)
				if splitvarname[1] == string(modelvars[j].name)
					if length(splitvarname) == 1
						modelvars[j].value = parse(Float64, splitline[2])
					else
						indices = map(i->parse(Int, i), splitvarname[2:end])
						modelvars[j].value[indices...] = parse(Float64, splitline[2])
					end
				end
			end
		end
		close(f)
	end
	return code
end
