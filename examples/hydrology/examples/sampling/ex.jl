using ThreeQ
import DataFrames
import Gadfly
import PyPlot

@generated function sample{numbits}(h, klow, khigh, rmax, beta, ::Type{Val{numbits}}; modelargs=("uq1d_model", "laptop", "c4-sw_sample", "workingdir", "c4"), num_reads=100, kwargs...)
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
		model = ThreeQ.Model(modelargs...)
		@defvar model q[1:length(h)-1]#1 if perm is high, 0 if perm is low
		@defvar model r[1:numbits, 1:length(h) - 2]#the bits of the groundwater source
		minval = 0
		for i = 2:length(h) - 1
			@addquadratic model ($innerexpr) ^ 2
			minval -= (klow * h[i + 1] - klow * h[i] - klow * h[i] + klow * h[i - 1]) ^ 2
		end
		@time ThreeQ.solvesapi!(model; num_reads=num_reads, kwargs...)
		ks = Array(Float64, length(h) - 1, num_reads)
		rs = Array(Float64, length(h) - 2, num_reads)
		llhoods = Array(Float64, length(model.bitsolutions))
		for j = 1:length(model.bitsolutions)
			@loadsolution(model, energy, occurrences, valid, j)
			llhoods[j] = -(log(occurrences) + beta * energy)
			for i = 1:length(h) - 1
				ks[i, j] = klow + q.value[i] * (khigh - klow)
			end
			for i = 1:length(h) - 2
				rs[i, j] = 0.
				for k = 1:numbits
					rs[i, j] -= rfactor * (2. ^ -k) * r.value[k, i]
				end
			end
		end
		return ks, rs, llhoods, model
	end
	return q
end


#solve 1d groundwater equation assuming the left and right boundary have 0 head, with a spacing of 1 between nodes
function solver(ks::Array{Float64, 1}, rs::Array{Float64, 1})
	I = Int[]
	J = Int[]
	V = Float64[]
	b = Array(Float64, length(ks) + 1)
	#set up left and right boundary condition
	push!(I, 1); push!(J, 1); push!(V, 1.); b[1] = 0.
	push!(I, length(ks) + 1); push!(J, length(ks) + 1); push!(V, 1.); b[end] = 0.
	#set up interior
	for i = 2:length(ks)
		b[i] = rs[i - 1]
		push!(I, i)
		push!(J, i - 1)
		push!(V, ks[i - 1])
		push!(I, i)
		push!(J, i + 1)
		push!(V, ks[i])
		push!(I, i)
		push!(J, i)
		push!(V, -(ks[i - 1] + ks[i]))
	end
	A = sparse(I, J, V)
	return A \ b
end

function b2f(kb, rb, klow, khigh, rmax)
	local ks = Array(Float64, length(kb))
	local rs = Array(Float64, size(rb, 2))
	local numbits = size(rb, 1)
	for i = 1:length(kb)
		if kb[i]
			ks[i] = khigh
		else
			ks[i] = klow
		end
	end
	rfactor = rmax / (1 - 2. ^ -numbits)
	for i = 1:length(rs)
		rs[i] = 0.
		for j = 1:numbits
			rs[i] -= rfactor * (2. ^ -j) * (rb[j, i] ? 1 : 0)
		end
	end
	return ks, rs
end

function solver(kb, rb, klow, khigh, rmax)
	ks, rs = b2f(kb, rb, klow, khigh, rmax)
	return solver(ks, rs)
end

function targetllhood(h0, ks, rs)
	h = solver(ks, rs)
	return -norm(h - h0) ^ 2
end
dwsolver = DWQMI.getdw2xsys4(mytoken)
beta = 6.8258
#=
dwsolver = DWQMI.defaultsolver
beta = 3.
=#
rmax = 1.
#S = 0
S = false
srand(0)
numhs = 10
numbits = 1
kb = bitrand(numhs - 1)
rb = bitrand(numbits, numhs - 2)
klow = 1
khigh = 2
h = solver(kb, rb, klow, khigh, rmax)
ks, rs, llhoods, model = sample(h + 0 * randn(length(h)), klow, khigh, rmax, beta, Val{numbits}; num_reads=10000, dwsolver=dwsolver)
reweightedllhoods = map(i->targetllhood(h, ks[:, i], rs[:, i]) - llhoods[i], 1:length(llhoods))
numshow = length(llhoods)
PyPlot.plot(1:length(llhoods[1:numshow]), llhoods[1:numshow], "k.")
PyPlot.xlabel("Sample Number")
PyPlot.ylabel("Boltzmann Log-Likelihood")
PyPlot.display(PyPlot.gcf()); println()
PyPlot.savefig("boltzmann.pdf")
PyPlot.close(PyPlot.gcf())
PyPlot.plot(1:length(llhoods[1:numshow]), reweightedllhoods[1:numshow], "k.")
PyPlot.xlabel("Sample Number")
PyPlot.ylabel("Target Log-Likelihood")
PyPlot.savefig("target.pdf")
PyPlot.display(PyPlot.gcf()); println()
PyPlot.close(PyPlot.gcf())
#=
df = DataFrames.DataFrame(Sample_Number=1:length(llhoods), Boltzmann_Likelihood=llhoods, Target_Likelihood=reweightedllhoods)
display(Gadfly.plot(df, x=:Sample_Number, y=:Boltzmann_Likelihood, Gadfly.Geom.point)); println()
display(Gadfly.plot(df, x=:Sample_Number, y=:Target_Likelihood, Gadfly.Geom.point)); println()
display(Gadfly.plot(x=1:length(llhoods), y=llhoods, Gadfly.Geom.line)); println()
display(Gadfly.plot(x=1:length(reweightedllhoods), y=reweightedllhoods, Gadfly.Geom.line)); println()
=#
rlhoods = exp(reweightedllhoods) / sum(exp(reweightedllhoods))
#@show sort(rlhoods)
#=
ks2, rs2 = b2f(kb, rb, klow, khigh, rmax)
h2 = solver(ks, rs)
display(Gadfly.plot(Gadfly.layer(x=0:numhs-1, y=h, Gadfly.Geom.point), Gadfly.layer(x=0:numhs-1, y=h2, Gadfly.Geom.line))); println()
display(Gadfly.plot(Gadfly.layer(x=0.5:1:numhs-1.5, y=ks, Gadfly.Geom.point), Gadfly.layer(x=0.5:1:numhs-1.5, y=ks2, Gadfly.Geom.line))); println()
display(Gadfly.plot(Gadfly.layer(x=1:1:numhs-2, y=rs, Gadfly.Geom.point), Gadfly.layer(x=1:1:numhs-2, y=rs2, Gadfly.Geom.line))); println()
=#
#=
@show norm(h-h2)
@show norm(ks-ks2)
@show norm(rs-rs2)
@assert norm(h-h2) == 0.
@assert norm(ks-ks2) == 0.
@assert norm(rs-rs2) == 0.
=#
nothing
