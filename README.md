ToQ.jl: Julia interface to D-Wave's ToQ
===============================

Description
-----------

ToQ.jl is a [Julia](http://julialang.org/) module that includes several macros which effectively extend the syntax of Julia to include the main features of the ToQ programming language developed by [D-Wave](http://dwavesys.com). ToQ.jl uses D-Wave's ToQ as an intermediate representation. It is also capable of using D-Wave's qbsolv as a backend.

Installation
------------

First, Julia and D-Wave's qOp must be installed. Julia binaries can be obtained [here](http://julialang.org/downloads/). QOp must be obtained from D-Wave.

Second, this repo must be cloned into a place where Julia can find it. You can tell Julia where to look for it by adding a line to your .juliarc.jl file like
```
push!(LOAD_PATH, "/dir/where/this/repo/is")
```
For instance, my .juliarc.jl file contains
```
push!(LOAD_PATH, "$(homedir())/codes")
```
and within ~/codes, there is a directory containing this repo called ToQ.jl.

License
-------

ToQ.jl is provided under a BSD-ish license with a "modifications must be indicated" clause.  See LICENSE.md file for the full text.

This package is part of the Hybrid Quantum-Classical Computing suite, known internally as LA-CC-16-032.

Author
------

Daniel O'Malley, <omalled@lanl.gov>
