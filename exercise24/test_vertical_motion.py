# INF5620 - Exercise 24 - Gunnar Raaen
# Nosetests for module vertical_motion
import nose.tools as nt
from vertical_motion import *
def exact_solution(v0,t):
	return v0 - g*t
def test_exact_vs_num_soltion():
	v0 = 0.1;
	# rho=d=0  gives a_S=a_d=0 and b=g
	# and so v'(t) = -g => v(t) = v0 - gt
	problem = Problem(v0=v0,rho=0,rho_b=1.0, mu=1.8*10**(-5),C_D=0.45, A=0.1, d=0, V=0.1)
	solver = Solver(problem, dt=0.001, T=2)
	v,t = solver.solve()
	v_e = np.array([exact_solution(v0,t[n]) for n in range(0,len(t))])
	# choosing the maxval of the diff between the exact and numerical sol.
	diff = np.abs(v_e - v).max()
	nt.assert_almost_equal(diff, 0, delta=1E-10);
	
def test_terminal_velocity():
	# parameters arbitrarly choosen as in exercise 26
	problem = Problem(v0=0,rho=0.79,rho_b=1003, mu=1.8*10**(-5),C_D=1.2, A=0.9, d=0.5, V=0.08)
	# Total time (T=100) has been found by trial and error to ensure max terminal velocity
	solver = Solver(problem, dt=0.01, T=100)
	v,t = solver.solve()
	# numerical terminal velocity
	v_T = abs(v[-1]);
	# exact terminal velocity
	v_Te = sqrt(abs(problem.b/problem.a_q));
	nt.assert_almost_equal(v_T, v_Te, delta=1E-10)
