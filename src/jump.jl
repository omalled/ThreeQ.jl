function bqp(Q, c, solver)
	m = JuMP.Model(solver=solver)
	@JuMP.variable(m, x[1:size(Q, 1)], Bin)
	@JuMP.objective(m, Min, dot(x, Q * x) + dot(c, x))
	JuMP.solve(m)
	return m
end
