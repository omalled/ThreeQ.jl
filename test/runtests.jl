using ThreeQ
using Base.Test

model = ThreeQ.Model("test_model", "laptop", "c4-sw_sample", "testdir", "c4")

#test @defvar
@defvar model q
@defvar model qs[1:10]
@defvar model qmat[1:3, 1:3]
@test length(model.vars) == 3
@test model.vars[1].name == :q
@test typeof(model.vars[1].value) == Float64
@test typeof(model.vars[2].value) == Array{Float64, 1}
@test typeof(model.vars[3].value) == Array{Float64, 2}

#test @defparam
@defparam model p
@defparam model p2
@test length(model.params) == 2
@test model.params[1].name == :p
@test model.params[2].name == :p2

#test @addterm
@addterm model (1 + exp(3)) * p * q * qs[2] * 2 * pi
@test isa(model.terms[1], ThreeQ.ParamQuadraticTerm)
@test model.terms[1].realcoeff == 2 * (1 + exp(3)) * pi
@test model.terms[1].param == p
@test model.terms[1].var1 == q
@test model.terms[1].var2 == qs[2]

#test @addquadratic
@addquadratic model q ^ 2
@test isa(model.terms[2], ThreeQ.LinearTerm)
@test model.terms[2].realcoeff == 1.
@test model.terms[2].var == q
@addquadratic model (qs[1] + -1 * qs[2]) ^ 2
@test model.terms[3].realcoeff == 1.
@test model.terms[3].var == qs[1]
@test model.terms[4].realcoeff == 1.
@test model.terms[4].var == qs[2]
@test model.terms[5].realcoeff == -2.
@test model.terms[5].var1 == qs[2]
@test model.terms[5].var2 == qs[1]
@addquadratic model (-1 * qs[2]) ^ 2
@test model.terms[6].realcoeff == 1.
@test model.terms[6].var == qs[2]
macro blah(x)
	return :($(esc(:q)))
end

#test that you can use macros inside the expressions for terms and quadratics
@addterm model @blah 1
@test isa(model.terms[7], ThreeQ.LinearTerm)
@test model.terms[7].realcoeff == 1.
@test model.terms[7].var == q
@addquadratic model (@blah 1) ^ 2
@test model.terms[8].realcoeff == 1.
@test model.terms[8].var == q

x = [1, 2]
@test_throws ErrorException eval(macroexpand(:(@addquadratic model (q + -x[1]) ^ 2)))#the -x[1] isn't supported, and it used to silently do the wrong thing
:passed
