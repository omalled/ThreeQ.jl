using ToQ
using Base.Test

model = ToQ.Model("test_model", "laptop", "c4-sw_sample", "testdir", "c4")

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
@test model.terms[1].prodstrings[1] == string(2 * (1 + exp(3)) * pi)
@test model.terms[1].prodstrings[2] == "p"
@test model.terms[1].prodstrings[3] == "q"
@test model.terms[1].prodstrings[4] == "qs___2"

#test @addquadratic
@addquadratic model q ^ 2
@test model.terms[2].prodstrings[1] == string(1.)
@test model.terms[2].prodstrings[2] == "q"
@test model.terms[2].prodstrings[3] == "q"
@addquadratic model (qs[1] + -1 * qs[2]) ^ 2
@test model.terms[3].prodstrings[1] == string(1.)
@test model.terms[3].prodstrings[2] == "qs___1"
@test model.terms[3].prodstrings[3] == "qs___1"
@test model.terms[4].prodstrings[1] == string(1.)
@test model.terms[4].prodstrings[2] == "qs___2"
@test model.terms[4].prodstrings[3] == "qs___2"
@test model.terms[5].prodstrings[1] == string(-2.)
@test model.terms[5].prodstrings[2] == "qs___2"
@test model.terms[5].prodstrings[3] == "qs___1"
@addquadratic model (-1 * qs[2]) ^ 2
@test model.terms[6].prodstrings[1] == string(1.)
@test model.terms[6].prodstrings[2] == "qs___2"
@test model.terms[6].prodstrings[3] == "qs___2"
