@Requires.require SCIP function bqp(Q, c)
	m = JuMP.Model(solver=SCIP.SCIPSolver())
	@JuMP.variable(m, x[1:size(Q, 1)], Bin)
	@JuMP.objective(m, Min, dot(x, Q * x) + dot(c, x))
	return JuMP.solve(m)
end
