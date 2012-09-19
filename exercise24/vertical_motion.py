
# INF5620 - Exercise 24 - Gunnar Raaen
# A program for vertical motion in a fluid
# This module is also used for exercises 25 and 26
from matplotlib.pyplot import *
from numpy import *
import nose.tools as nt
import argparse
g = 9.81
class Problem:
	def __init__(self, v0, rho, rho_b, V, mu, d, C_D, A):
		self.v0 = v0
		self.rho, self.rho_b, self.V, self.mu, self.d, self.C_D, self.A = rho, rho_b, V, mu, d, C_D, A
		self.b = g*(rho/rho_b -1)
		# Stokes' drag model coeff
		self.a_S = 3*pi*d*mu/(rho_b*V)
		# Quadratic drag model coeff
		self.a_q = 0.5*C_D*rho*A/(rho_b*V)
	def define_command_line_options(self, parser=None):
		if parser is None:
			parser = argparse.ArgumentParser()
		parser.add_argument('--v0', type=float, default=self.v0, help='initial condition, v(0)', metavar='v0');

		return parser
	def init_from_command_line(self, args):
		self.v0 = args.v0
class Solver:
	def __init__(self, problem, dt, T):
		self.dt, self.T, self.problem = dt, T, problem
		self.Fd, self.Fb = 0,0
	def define_command_line_options(self, parser):
		if parser is None:
			parser = argparse.ArgumentParser();
		parser.add_argument('--dt', type=float, default=self.dt, help='Delta t',metavar='dt')
		parser.add_argument('--T', type=float, default=self.T, help='Total time', metavar='T');
		return parser
	def init_from_command_line(self, args):
		self.dt, self.T = args.dt, args.T
	def solve(self):
		T, a_S, a_q, b, dt = self.T, self.problem.a_S, self.problem.a_q, self.problem.b, self.dt
		rho, V, mu, C_D, A = self.problem.rho, self.problem.V, self.problem.mu, self.problem.C_D, self.problem.A
		N = int(T/dt)
		T = N*dt
		self.Fd = zeros(N+1) # relevant drag force
		self.Fb = rho*g*V*(zeros(N+1)+1) # buoyancy force (constant in this problem)
		v = zeros(N+1)	# velocity
		t = linspace(0,T, N+1)
		v[0] = self.problem.v0
		temp = rho*self.problem.d/mu
		# CN scheme
		for n in range(0,N):
			if temp*abs(v[n]) < 1: # Stokes' drag
				self.Fd[n] = -3*pi*mu*v[n]
				v[n+1] = (v[n]*(2-a_S*dt) + 2*b*dt)/(2+a_S*dt)
			else:	# Quad drag
				self.Fd[n] = -0.5*C_D*rho*A*abs(v[n])*v[n]
				v[n+1] = (v[n] + dt*b)/(1 + dt*a_q*abs(v[n]))
		# calc. force at end
		if temp*abs(v[-1]) < 1:	
			self.Fd[-1] = -3*pi*mu*v[-1]
		else:
			self.Fd[-1] = -0.5*C_D*rho*A*abs(v[-1])*v[-1]

		return v, t
	def getDragForce(self):
		return self.Fd
	def getBouyancyForce(self):
		return self.Fb
# Test
if __name__ == '__main__':
	problem = Problem(v0=0.0,rho=0.79 ,rho_b=1003, mu=1.8*10**-5, C_D=1.2, A=0.9, d=0.5, V=0.08)
	solver = Solver(problem, dt=0.01, T=50)
	
	parser = problem.define_command_line_options()
	parser = solver.define_command_line_options(parser)
	args = parser.parse_args()
	problem.init_from_command_line(args)
	solver.init_from_command_line(args);
	v,t = solver.solve()
		
	plot(t,v)
	xlabel('time[s]')
	ylabel('z')
	show()


