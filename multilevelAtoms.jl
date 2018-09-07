# multilevelAtoms.jl
# Computes properties of multilevel atoms (Sr, in particular) under the influence of near resonant light
module multilevelAtoms
import Base.^
import Base.*
using Plots
using LinearAlgebra
#using DifferentialEquations

# constants
const h = 6.626e-34
const c = 299792458
# Indexing levels by 1=1S0, 2 = 3P0, 3=3P1, 4=3P2, 5=3S1
Gamma = [0 0 0 0 0;0 0 0 0 0;.0465 0 0 0 0;0 0 0 0 0;0 8.9 27 42 0] * 1e6
tau = 1 ./ Gamma
lambda = [0 0 0 0 0;0 0 0 0 0;689 0 0 0 0;0 0 0 0 0;0 679 688 707 0] * 1e-9
freq = c./lambda
ISAT = pi/3 * h*freq ./ (lambda.^2 .* tau)

struct Level
	spin::Int
	L::Int
	J::Int
	name::String
end
Base.show(io::IO,l::Level) = print(io,l.name)

struct Transition
	lower::Level
	upper::Level
	lambda::Float64	# In nm
	linewidth::Float64	# In Hz
	ISat::Float64
end
Transition(l::Level,u::Level,lam::Float64,lin::Float64) = Transition(l,u,lam,lin, pi/3 * h*c/(1.0e-9*lam)^3 * lin)
Base.show(io::IO, t::Transition) = print(io,t.lower.name*" -> "*t.upper.name)

struct Atom
	transitions::Array{Transition,1}
	levels::Array{Level,1}
	idx
	Atom(t::Array{Transition,1},l::Array{Level,1}) = new(t,l,Dict(l[j].name=>j for j in 1:length(l)))
end
Atom(x::Array{Transition,1}) = Atom(x,collect(union(Set([j.upper for j in x]),Set([j.lower for j in x]))))
^(a::Atom,t::Transition) = a.idx[t.lower.name],a.idx[t.upper.name]

mutable struct Laser
	line::Transition
	saturation::Float64
	detuning::Float64
	polarization::Array{Float64,1}
	intensity::Float64
end
Laser(l::Transition,s::Number) = Laser(l,convert(Float64,s),0,[0,1,0],s*l.ISat)
*(l::Laser,n::Number) = Laser(l.line,l.saturation*n,l.detuning,l.polarization,l.intensity*n)

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
L689 = Laser(_689,1)
L707 = Laser(_707,1)
L688 = Laser(_688,1)
L679 = Laser(_679,1)
L461 = Laser(_461,1)

Sr = Atom([_689,_707,_688,_679,_461 ])
SrTriplets = Atom([_689,_707,_688,_679])

#-------------------------------------

function rate(l::Laser)
	s,d,y = l.saturation,l.detuning,l.line.linewidth
	#return .5*s*y / (1+s+(2*d/y)^2)
	#return (1/6)*s*y	# Not fully correct...
	return sqrt(s/2)*y/(2*pi)
end

function rateMatrix(atom::Atom,lasers::Array{Laser,1})
	N = length(atom.levels)
	gamma = zeros(N,N)
	pump = zeros(N,N)
	decay = zeros(N,N)
	
	for t in atom.transitions
		gamma[atom^t...] = t.linewidth
	end
	for l in lasers
		pump[atom^l.line...] = rate(l)
	end
	pump += pump'		# Symmetrize
	
	for i=1:N
		decay[i,i] = sum(pump[:,i]) + sum(gamma[:,i])
	end
	
	return gamma + pump - decay
end

function plotRepump(vals,lasers=[L689,L688,L679,L707])
	N = length(vals)
	ev = zeros(N,5)
	for i=1:N
		M = rateMatrix(SrTriplets,lasers*vals[i])
		E = eigen(M)
		mx = findmax(E.values)[2]
		v = E.vectors[:,mx]
		v /= sum(v)
		ev[i,:] = v
	end
	#plot(ev)
	return ev
end


function BlochEquations(atom,lasers,init)
	
end


end