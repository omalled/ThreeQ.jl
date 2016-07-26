using ToQ

import Base.factor

function factor(N, numbits; factor_weightval=1, ancillaryval=1, S=20, showoutput=false, connection="laptop", solver="c4-sw_sample", workspace="c4")
	#factor N = n1 * n2 where
	#n1 = ∑2^(i-1) * n1[i] and
	#n2 = ∑2^(i-1) * n2[i]
	#with the sums going from 1 to numbits

	#set up the model
	model = ToQ.Model("factor_model", connection, solver, "workingdir", workspace)

	#set up the params and vars
	@defparam model factor_weight
	@defparam model ancillary
	@defvar model n1[1:numbits]
	@defvar model n2[1:numbits]
	@defvar model s[1:numbits, 1:numbits]

	#assemble the QUBO
	for i = 1:numbits
		for j = 1:numbits
			#make s[i, j] = n1[j] * n2[i]
			@addterm model ancillary * n1[j] * n2[i]
			@addterm model ancillary * -2 * n1[j] * s[i, j]
			@addterm model ancillary * -2 * n2[i] * s[i, j]
			@addterm model ancillary * 3 * s[i, j]
			#make n1 * n2 = N
			for i2 = 1:numbits
				for j2 = 1:numbits
					@addterm model factor_weight * 2 ^ (i + j + i2 + j2 - 4) * s[i, j] * s[i2, j2]
				end
			end
			@addterm model factor_weight * -2 * 2 ^ (i + j - 2) * N * s[i, j]
		end
	end

	#solve the QUBO with qbsolv
	ToQ.qbsolv!(model; minval=-factor_weightval * N ^ 2, S=S, factor_weight=factor_weightval, ancillary=ancillaryval, showoutput=showoutput)

	#translate the binary result into a couple julia Ints
	n1val = 0
	n2val = 0
	for i = 1:numbits
		n1val += 2 ^ (i - 1) * Int(n1.value[i])
		n2val += 2 ^ (i - 1) * Int(n2.value[i])
	end
	return n1val, n2val
end

Ns = [35, 143]#, 899]#, 3599]
for i = 1:length(Ns)
	N = Ns[i]
	numbits = i + 2
	@show N
	@time n1, n2 = factor(N, numbits; S=0, showoutput=false)
	if n1 * n2 != N
		error("factoring not performed correctly: $N != $n1 * $n2")
	else
		println("Success: $N == $n1 * $n2")
	end
end
