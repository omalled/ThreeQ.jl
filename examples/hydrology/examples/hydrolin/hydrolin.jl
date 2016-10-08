import HydroToQ
import FiniteDifference2D
import PyPlot

srand(0)
numobs = 25
nfd = 1000
obsI = rand(1:nfd, numobs)
obsJ = rand(1:nfd, numobs)
weights = ones(numobs)
u_n(x, y) = 0
u_d(x, y) = y
f(x, y) = 0
a, b, c, d = 0, 1, 0, 1
@generated function fd{T}(x, style::Type{Val{T}}=Val{0})
	if T == :J
		retcode = quote
			return jacobian
		end
	else
		retcode = :(return map((i, j)->u[i, j], obsI, obsJ))
	end
	q = quote
		kx, ky = FiniteDifference2D.x2ks(x, nfd, nfd)
		u, jacobian = FiniteDifference2D.solvedsresiduals(obsI, obsJ, weights, zeros(numobs), kx, ky, u_d, u_n, f, nfd, nfd, a, b, c, d)
		$retcode
	end
	return q
end
#=
@show fd(ones(2 * (nfd + 1) * nfd))
@show maximum(fd(ones(2 * (nfd + 1) * nfd), Val{:J}))
=#

klow = 1.
khigh = 1.01
sqrtnump = 7
truep = map(x->x > 0 ? klow : khigh, randn(sqrtnump * sqrtnump))
@generated function ks{T}(p, style::Type{Val{T}}=Val{1})
	if T == :J
		kdef = :(ks = zeros(Float64, 2 * nfd * (nfd + 1), length(p)))
		setkscode = quote
			ks[i + (j - 1) * (nfd + 1), pix +  (pjx - 1) * sqrtnump] = 1.
			ks[(nfd + 1) * nfd + j + (i - 1) * nfd, piy + (pjy - 1) * sqrtnump] = 1.
		end
	else
		kdef = :(ks = Array(Float64, 2 * nfd * (nfd + 1)))
		setkscode = quote
			ks[i + (j - 1) * (nfd + 1)] = newp[pix, pjx]
			ks[(nfd + 1) * nfd + j + (i - 1) * nfd] = newp[piy, pjy]
		end
	end
	q = quote
		#ks = Array(Float64, 2 * nfd * (nfd + 1))
		$kdef
		newp = reshape(p, sqrtnump, sqrtnump)
		for i = 1:nfd + 1
			for j = 1:nfd
				pix = ceil(Int, i * sqrtnump / (nfd + 1))
				pjx = ceil(Int, j * sqrtnump / nfd)
				pjy = ceil(Int, i * sqrtnump / (nfd + 1))
				piy = ceil(Int, j * sqrtnump / nfd)
				$setkscode
			end
		end
		return ks
	end
end
#=
PyPlot.clf(); PyPlot.imshow(reshape(truep, sqrtnump, sqrtnump), interpolation="nearest"); display(PyPlot.gcf()); println()
PyPlot.clf(); PyPlot.imshow(reshape(ks(truep)[1:(nfd + 1)*nfd], nfd + 1, nfd), interpolation="nearest"); display(PyPlot.gcf()); println()
PyPlot.clf(); PyPlot.imshow(reshape(ks(truep)[(nfd + 1)*nfd+1:end], nfd, nfd + 1), interpolation="nearest"); display(PyPlot.gcf()); println()
J = ks(truep, Val{:J})
for i = 1:9
	PyPlot.clf(); PyPlot.imshow(reshape(J[1:(nfd + 1)*nfd, i], nfd + 1, nfd), extent=[0, 1, 0, 1], interpolation="nearest"); display(PyPlot.gcf()); println()
	PyPlot.clf(); PyPlot.imshow(reshape(J[(nfd + 1)*nfd + 1:end, i], nfd, nfd + 1), extent=[0, 1, 0, 1], interpolation="nearest"); display(PyPlot.gcf()); println()
end
=#

function f(p)
	return fd(ks(p))
end
function g(p)
	x = ks(p)
	Jk = fd(x, Val{:J})
	dkdp = ks(p, Val{:J})
	return Jk * dkdp
end
#=
@show f(truep)
@show g(truep)
=#
uobs = f(truep)
p0 = fill(klow, sqrtnump * sqrtnump)

function callback(p)
	fig, axes = PyPlot.subplots(nrows=1, ncols=2)
	axes[1][:imshow](reshape(truep, sqrtnump, sqrtnump), extent=[a, b, c, d], interpolation="nearest")
	axes[1][:set_title]("true k")
	img = axes[2][:imshow](reshape(p, sqrtnump, sqrtnump), extent=[a, b, c, d], interpolation="nearest")
	axes[2][:set_title]("inverse k")
	fig[:subplots_adjust](top=0.9)
	cbar_ax = fig[:add_axes]([0.05, 0.85, 0.9, 0.05])
	fig[:colorbar](img, cax=cbar_ax, orientation="horizontal")
	for i = 1:numobs
		axes[2][:plot]((nfd - obsJ[i]) / nfd, obsI[i] / nfd, "m.")
	end
	display(PyPlot.gcf()); println()
	PyPlot.close(fig)
end
#p, m, q = HydroToQ.linearizedinverse(f, g, uobs, p0, klow, khigh; showof=true, callback=callback)
nothing
