import HydroToQ
import Gadfly

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

modelargs=("uq1d_model", "online", "DW2X_SYS4", "workingdir", "sys4")
rmax = 1.
#S = 0
S = false
srand(0)
numhs = 32
numbits = 1
kb = bitrand(numhs - 1)
rb = bitrand(numbits, numhs - 2)
klow = 1
khigh = 2
h = solver(kb, rb, klow, khigh, rmax)
ks, rs, model = HydroToQ.uq1d(h + 0 * randn(length(h)), klow, khigh, rmax, numbits; modelargs=modelargs, S=S, tol=1e-6)
ks2, rs2 = b2f(kb, rb, klow, khigh, rmax)
h2 = solver(ks, rs)
display(Gadfly.plot(Gadfly.layer(x=0:numhs-1, y=h, Gadfly.Geom.point), Gadfly.layer(x=0:numhs-1, y=h2, Gadfly.Geom.line))); println()
display(Gadfly.plot(Gadfly.layer(x=0.5:1:numhs-1.5, y=ks, Gadfly.Geom.point), Gadfly.layer(x=0.5:1:numhs-1.5, y=ks2, Gadfly.Geom.line))); println()
display(Gadfly.plot(Gadfly.layer(x=1:1:numhs-2, y=rs, Gadfly.Geom.point), Gadfly.layer(x=1:1:numhs-2, y=rs2, Gadfly.Geom.line))); println()
#=
@show norm(h-h2)
@show norm(ks-ks2)
@show norm(rs-rs2)
@assert norm(h-h2) == 0.
@assert norm(ks-ks2) == 0.
@assert norm(rs-rs2) == 0.
=#
