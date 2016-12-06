using ThreeQ

function ntobits(N, numbits=ceil(Int, log2(N)))
	bits = Array(Int, numbits)
	for i = 1:numbits
		bits[i] = div(N, 2 ^ (numbits - i))
		N -= bits[i] * 2 ^ (numbits - i)
	end
	return reverse(bits)
end

function factor3bit(N; connection="laptop", solver="c4-sw_sample", workspace="c4", ancillaryval=1)
	Nbits = ntobits(N, 7)
	m = ThreeQ.Model("factor_model", connection, solver, "workingdir", workspace)
	numbits = 3
	@defparam m ancillary
	@defvar m n1[1:numbits]
	@defvar m n2[1:numbits]
	@defvar m z[1:2 * numbits, 1:2]
	@defvar m s[1:numbits, 1:numbits]
	for i = 1:numbits
		for j = 1:numbits
			#make s[i, j] = n1[i] * n2[j]
			@addterm m ancillary * n1[i] * n2[j]
			@addterm m ancillary * -2 * n1[i] * s[i, j]
			@addterm m ancillary * -2 * n2[j] * s[i, j]
			@addterm m ancillary * 3 * s[i, j]
		end
	end
	@addquadratic m (s[1, 1] + -1 * Nbits[1]) ^ 2
	@addquadratic m (s[2, 1] + s[1, 2] + -(2 ^ 1) * z[2, 1] + -1 * Nbits[2]) ^ 2
	@addquadratic m (s[3, 1] + s[2, 2] + s[1, 3] + (2 ^ 0) * z[2, 1] + -(2 ^ 1) * z[3, 1] + -(2 ^ 2) * z[3, 2] + -1 * Nbits[3]) ^ 2
	@addquadratic m (s[3, 2] + s[2, 3] + (2 ^ 0) * z[3, 1] + (2 ^ 1) * z[3, 2] + -(2 ^ 1) * z[4, 1] + -(2 ^ 2) * z[4, 2] + -1 * Nbits[4]) ^ 2
	@addquadratic m (s[3, 3] + (2 ^ 0) * z[4, 1] + (2 ^ 1) * z[4, 2] + -(2 ^ 1) * z[5, 1] + -(2 ^ 2) * z[5, 2] + -1 * Nbits[5]) ^ 2
	@addquadratic m ((2 ^ 0) * z[5, 1] + (2 ^ 1) * z[5, 2] + -(2 ^ 1) * z[6, 1] + -1 * Nbits[6]) ^ 2
	@addquadratic m ((2 ^ 0) * z[6, 1] + -1 * Nbits[7]) ^ 2
	ThreeQ.qbsolv!(m; minval=-sum(Nbits), ancillary=ancillaryval)
	n1val = 0
	n2val = 0
	for i = 1:numbits
		n1val += 2 ^ (i - 1) * Int(n1.value[i])
		n2val += 2 ^ (i - 1) * Int(n2.value[i])
	end
	return n1val, n2val
end

@generated function factorgeneral{numbits}(N, ::Type{Val{numbits}}; connection="laptop", solver="c4-sw_sample", workingdir="workingdir", workspace="c4", ancillaryval=1, S=0)
	startcode = quote
		Nbits = ntobits(N, 2 * numbits)
		m = ThreeQ.Model("factor_model", connection, solver, workingdir, workspace)
		@defparam m ancillary
		@defvar m n1[1:numbits]
		@defvar m n2[1:numbits]
		@defvar m z[1:3 * numbits, 1:numbits]
		@defvar m s[1:numbits, 1:numbits]
		for i = 1:numbits
			for j = 1:numbits
				#make s[i, j] = n1[i] * n2[j]
				@addterm m ancillary * n1[i] * n2[j]
				@addterm m ancillary * -2 * n1[i] * s[i, j]
				@addterm m ancillary * -2 * n2[j] * s[i, j]
				@addterm m ancillary * 3 * s[i, j]
			end
		end
	end
	finishcode = quote
		ThreeQ.qbsolv!(m; minval=-sum(Nbits), ancillary=ancillaryval, S=S)
		n1val = 0
		n2val = 0
		for i = 1:numbits
			n1val += 2 ^ (i - 1) * Int(n1.value[i])
			n2val += 2 ^ (i - 1) * Int(n2.value[i])
		end
		return n1val, n2val
	end
	middlecode = quote begin end end
	numcarrybits = 0
	for i = 1:2 * numbits
		innerexpr = :(0 + 0)
		numsterms = 0
		for j = 1:numbits
			for k = 1:numbits
				if j + k == i + 1
					push!(innerexpr.args, :(s[$k, $j]))
					numsterms += 1
				end
			end
		end
		for j = 1:numcarrybits
			push!(innerexpr.args, :((2 ^ $(j - 1)) * z[$(i - 1), $j]))
		end
		numcarrybits = floor(Int, log2(numsterms + 2 ^ numcarrybits - 1))
		for j = 1:numcarrybits
			push!(innerexpr.args, :(-(2 ^ $j) * z[$i, $j]))
		end
		push!(innerexpr.args, :(-1 * Nbits[$i]))
		innerexpr.args = [innerexpr.args[1]; innerexpr.args[4:end]]#get rid of the 0 + 0 at the beginning
		push!(middlecode.args, :(@addquadratic m ($innerexpr) ^ 2))
	end
	i = 2 * numbits + 1
	while numcarrybits > 0
		numsterms = 0
		innerexpr = :(0 + 0)
		for j = 1:numcarrybits
			push!(innerexpr.args, :((2 ^ $(j - 1)) * z[$(i - 1), $j]))
		end
		numcarrybits = floor(Int, log2(numsterms + 2 ^ numcarrybits - 1))
		for j = 1:numcarrybits
			push!(innerexpr.args, :(-(2 ^ $j) * z[$i, $j]))
		end
		innerexpr.args = [innerexpr.args[1]; innerexpr.args[4:end]]#get rid of the 0 + 0 at the beginning
		push!(middlecode.args, :(@addquadratic m ($innerexpr) ^ 2))
		i += 1
	end
	code = :($startcode; $middlecode; $finishcode)
	#@show code
	return code
end

function test(numbits, S=0)
	ns = Int[]
	for truen1 = 1:2 ^ numbits - 1
		for truen2 = 1:truen1
			push!(ns, truen1 * truen2)
		end
	end
	for n in sort(unique(ns))
		n1, n2 = factorgeneral(n, Val{numbits})
		if n1 * n2 == n
			println("success: $n1 * $n2 == $n ($numbits bit factorization)")
		else
			println("failure: $n1 * $n2 != $n ($numbits bit factorization)")
		end
	end
end
function test3()
	for truen1 = 1:2 ^ 3 - 1
		for truen2 = 1:truen1
			n1, n2 = factor3bit(truen1 * truen2)
			if n1 * n2 == truen1 * truen2
				println("success: $n1 * $n2 == $(truen1 * truen2)")
			else
				println("failure: $n1 * $n2 != $(truen1 * truen2)")
			end
		end
	end
end

for i = 2:5
	println("Factoring all products of $i bit numbers")
	test(i, false)
end
