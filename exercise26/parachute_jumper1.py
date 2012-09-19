# INF5620 - Exercise 26 - Gunnar Raaen
# Simulates a free fall of a parachute jumper using the vertical_motion module from ex. 24

problem = Problem(v0=0.0,rho=0.79 ,rho_b=1003, mu=1.8*10**-5, C_D=1.2, A=0.9, d=0.5, V=0.08)
solver = Solver(problem, dt=0.01, T=50)
	
parser = problem.define_command_line_options()
parser = solver.define_command_line_options(parser)
args = parser.parse_args()
problem.init_from_command_line(args)
solver.init_from_command_line(args)
v,t = solver.solve()

plot(t,v)
xlabel('time[s]')
ylabel('z')
title('Parachute jumper')
show()


