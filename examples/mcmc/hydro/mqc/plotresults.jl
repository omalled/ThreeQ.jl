import PyPlot

fig, ax = PyPlot.subplots()
barwidth = 0.35
indices = [1, 2, 3]
ax[:bar](indices, results[1, :], barwidth, color="b", label="D-Wave+MQC")
ax[:bar](barwidth .+ indices, results[2, :], barwidth, color="g", label="D-Wave")
ax[:set_ylabel]("Probability of finding optimum")
ax[:set_xlabel]("Permeability contrast")
ax[:set_xticklabels](["", "2", "", "10", "", "100"])
ax[:legend](loc=3)
display(fig)
println()
fig[:savefig]("result.pdf")
PyPlot.close(fig)
