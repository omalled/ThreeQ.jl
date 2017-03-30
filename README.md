ThreeQ.jl: Julia interface to D-Wave's quantum annealing hardware
===============================

Description
-----------

ThreeQ.jl is a [Julia](http://julialang.org/) module that includes several macros which effectively extend the syntax of Julia to enable rapid prototyping of codes for [D-Wave](https://www.dwavesys.com/) hardware. ThreeQ.jl uses various tools from D-Wave's software suite (including the Python SAPI, [qbsolv](https://github.com/dwavesystems/qbsolv), and dw) as backends to access D-Wave hardware.

Installation
------------

First, Julia  must be installed. Julia binaries can be obtained [here](http://julialang.org/downloads/). In order to actually use the D-Wave, D-Wave hardware must be available and at least part of D-Wave's software stack must also be installed. D-Wave's [qbsolv](https://github.com/dwavesystems/qbsolv) tool can be used without D-Wave hardware. Other components (the Python SAPI and dw) must be obtained from D-Wave, if needed.

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
