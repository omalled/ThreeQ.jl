module HydroToQ

using ToQ

function uq1d(h, klow, khigh, rmax, numbits::Int; modelargs=("uq1d_model", "laptop", "c4-sw_sample", "workingdir", "c4"), tol=1e-6, S=S)
	if numbits < 1
		error("numbits must be at least one, but numbits=$numbits")
	end
	return uq1d(h, klow, khigh, rmax, Val{numbits}; modelargs=modelargs, tol=tol, S=S)
end

@generated function uq1d{numbits}(h, klow, khigh, rmax, ::Type{Val{numbits}}; modelargs=("uq1d_model", "laptop", "c4-sw_sample", "workingdir", "c4"), tol=1e-6, S=0)
	@assert typeof(numbits) == Int
	@assert numbits >= 1
	innerexpr = :((klow * h[i + 1] + (khigh - klow) * q[i] * h[i + 1] +
					-1 * klow * h[i] + -1 * (khigh - klow) * q[i] * h[i] +
					-1 * klow * h[i] + -1 * (khigh - klow) * q[i - 1] * h[i] +
					klow * h[i - 1] + (khigh - klow) * q[i - 1] * h[i - 1]))
	for i = 1:numbits
		push!(innerexpr.args, :(rfactor * (2. ^ $(-i)) * r[$i, i - 1]))
	end
	q = quote
		rfactor = rmax / (1 - 2. ^ -numbits)
		model = ToQ.Model(modelargs...)
		@defvar model q[1:length(h)-1]#1 if perm is high, 0 if perm is low
		@defvar model r[1:numbits, 1:length(h) - 2]#the bits of the groundwater source
		minval = 0
		for i = 2:length(h) - 1
			@addquadratic model ($innerexpr) ^ 2
			minval -= (klow * h[i + 1] - klow * h[i] - klow * h[i] + klow * h[i - 1]) ^ 2
		end
		@time ToQ.qbsolv!(model; showoutput=true, minval=minval + tol, S=S)
		ks = Array(Float64, length(h) - 1)
		rs = Array(Float64, length(h) - 2)
		for i = 1:length(h) - 1
			ks[i] = klow + q.value[i] * (khigh - klow)
		end
		for i = 1:length(h) - 2
			rs[i] = 0.
			for j = 1:numbits
				rs[i] -= rfactor * (2. ^ -j) * r.value[j, i]
			end
		end
		return ks, rs, model
	end
	return q
end

function linearizediteration(uobs, u0, p0, J, p1, p2; modelargs=("lin_model", "laptop", "c4-sw_sample", "workingdir", "c4"))
	return linearizediteration(uobs, u0, p0, J, p1, p2, Val{size(J)}; modelargs=modelargs)
end

@generated function linearizediteration{sizeJ}(uobs, u0, p0, J, p1, p2, ::Type{Val{sizeJ}}; modelargs=("lin_model", "laptop", "c4-sw_sample", "workingdir", "c4"))
	numu = sizeJ[1]
	nump = sizeJ[2]
	innerexpr = :(u0[i] + -1 * uobs[i])
	for i = 1:nump
		push!(innerexpr.args, :(J[i, $i] * (pother[$i] - p0[$i]) * q[$i]))
	end
	q = quote
		model = ToQ.Model(modelargs...)
		@defvar model q[1:$nump]
		pother = map(p_i->p_i == p1 ? p2 : p1, p0)
		for i = 1:$numu
			@addquadratic model ($innerexpr) ^ 2
		end
		ToQ.qbsolv!(model; showoutput=false)
		p = Array(Float64, $nump)
		for i = 1:$nump
			if q.value[i] == 0
				p[i] = p0[i]
			else
				p[i] = pother[i]
			end
		end
		return p, model, q.value
	end
	return q
end

function linearizedinverse(f, g, uobs, p0, p1, p2; modelargs=("lin_model", "laptop", "c4-sw_sample", "workingdir", "c4"), tol=0, showof=false, callback=p->nothing)
	lastp = copy(p0)
	p = similar(p0)
	done = false
	local m
	local q
	bestof = Inf
	bestp = copy(p0)
	callback(p0)
	while !done
		u0 = f(lastp)
		J = g(lastp)
		p, m, q = linearizediteration(uobs, u0, lastp, J, p1, p2; modelargs=modelargs)
		newof = norm(f(p) - uobs) ^ 2
		if bestof - newof <= tol
			done = true
		end
		if newof < bestof
			bestp = copy(p)
			bestof = newof
		end
		if showof
			@show bestof
		end
		callback(bestp)
		lastp = copy(p)
	end
	return bestp, m, q
end

end
