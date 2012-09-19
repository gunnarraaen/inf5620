# INF5620 - Exercise 25 - Gunnar Raaen
# Plot forces acting in vertical motion in a fluid using module vertical_motion from ex. 24
from vertical_motion import *
if __name__ == '__main__':
	problem = Problem(v0=1.0,rho=0.01 ,rho_b=503, mu=1.8*10**-5, C_D=1.5, A=0.9, d=0.005, V=0.2)
	solver = Solver(problem, dt=0.01, T=170)
	
	parser = problem.define_command_line_options()
	parser = solver.define_command_line_options(parser)
	args = parser.parse_args()
	problem.init_from_command_line(args)
	solver.init_from_command_line(args);
	v,t = solver.solve()

	subplot(2,1,1)
	plot(t,solver.getBouyancyForce())
	xlabel('time[s]')
	ylabel('Force[N]')
	title('Bouyancy force')
	subplot(2,1,2)
	# getDragForce() returns an array where each element is the relevant drag force model
	plot(t,solver.getDragForce())
	xlabel('time[s]')
	ylabel('Force[N]')
	title('Drag force')
	show()
