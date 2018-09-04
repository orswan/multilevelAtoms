# multilevelAtoms.jl
# Computes properties of multilevel atoms (Sr, in particular) under the influence of near resonant light

using Plots#, #DifferentialEquations

module multilevelAtoms
# constants
h = 6.626e-34
c = 299792458
# Indexing levels by 1=1S0, 2 = 3P0, 3=3P1, 4=3P2, 5=3S1
Gamma = [0 0 0 0 0;0 0 0 0 0;.0465 0 0 0 0;0 0 0 0 0;0 8.9 27 42 0] * 1e6
tau = 1 ./ Gamma
lambda = [0 0 0 0 0;0 0 0 0 0;689 0 0 0 0;0 0 0 0 0;0 679 688 707 0] * 1e-9
freq = c./lambda
ISAT = pi/3 * h*freq ./ (lambda.^2 .* tau)

mutable struct Level
	spin::Int
	L::Int
	J::Int
	name::String
end

mutable struct Transition
	lower::Level
	upper::Level
	lambda::Float64	# In nm
	linewidth::Float64	# In Hz
	ISat::Float64
end
Transition(l::Level,u::Level,lam::Float64,lin::Float64) = Transition(l,u,lam,lin, pi/3 * h*c/lam^3 * lin)

mutable struct Atom
	transitions::Array{Transition,1}
	levels::Array{Level,1}
end
Atom(x::Array{Transition,1}) = Atom(x,collect(union(Set([j.upper for j in x]),Set([j.lower for j in x]))))

mutable struct Laser
	line::Transition
	detuning::Float64
	polarization::Array{Int64,1}
	intensity::Float64
	saturation::Float64
end

#---------------- Sr -----------------
_1s0 = Level(0,0,0,"1S0")
_3p0 = Level(1,1,0,"3P0")
_3p1 = Level(1,1,1,"3P1")
_3p2 = Level(1,1,2,"3P2")
_3s1 = Level(1,0,1,"3S1")
_1p1 = Level(1,1,1,"1P1")
_689 = Transition(_1s0,_3p1,689.449,4.69e4)
_707 = Transition(_3p2,_3s1,707.202,4.2e7)
_688 = Transition(_3p1,_3s1,688.021,2.7e7)
_679 = Transition(_3p0,_3s1,679.290,8.9e6)
_461 = Transition(_1s0,_1p1,460.862,2.01e8)
Sr = Atom([_689,_707,_688,_679,_461 ])

#-------------------------------------

function rateEquations(atom,lasers,init)

end

function BlochEquations(atom,lasers,init)

end


end