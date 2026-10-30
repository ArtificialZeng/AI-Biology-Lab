"""Independent analytic/actual-backend observation for GillesPy2 ODE wrapper.

No candidate test code is imported. The spy subclasses the real SciPy backend;
all integrate calls still execute SciPy. Source directory is a command argument.
"""
import contextlib
import importlib
import json
import math
import pathlib
import sys
import warnings

source = pathlib.Path(sys.argv[1]).resolve()
sys.path.insert(0, str(source))
import gillespy2 as gp
from scipy.integrate import ode as RealODE
wrapper = importlib.import_module('gillespy2.solvers.numpy.ode_solver')
observations = []

class TrackingODE(RealODE):
    def integrate(self, target, *args, **kwargs):
        value = super().integrate(target, *args, **kwargs)
        observations.append({'requested': float(target), 'actual': float(self.t),
                             'successful': bool(self.successful()),
                             'backend_return_code': int(self.get_return_code())})
        return value

wrapper.ode = TrackingODE

def make_model():
    model = gp.Model(name='independent_transfer')
    a = gp.Species(name='A', initial_value=1, mode='continuous')
    b = gp.Species(name='B', initial_value=0, mode='continuous')
    k = gp.Parameter(name='k', expression=1)
    model.add_species([a, b]); model.add_parameter(k)
    model.add_reaction(gp.Reaction(name='transfer', reactants={a:1}, products={b:1}, rate=k))
    model.timespan(gp.TimeSpan.linspace(t=2, num_points=3))
    return model

def trial(label, integrator, opts, solver=None):
    observations.clear()
    solver = gp.ODESolver(model=make_model()) if solver is None else solver
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        try:
            result = solver.run(integrator=integrator, integrator_options=opts)
            answer = {'returned_results': True, 'time': result['time'].tolist(),
                      'A': result['A'].tolist(), 'B': result['B'].tolist(),
                      'max_analytic_error': max(abs(float(a)-math.exp(-float(t))) for t,a in zip(result['time'], result['A'])),
                      'max_conservation_error': max(abs(float(a+b)-1) for a,b in zip(result['A'], result['B']))}
        except Exception as exc:
            answer = {'returned_results': False, 'error_type': type(exc).__name__, 'error': str(exc),
                      'cause_type': type(exc.__cause__).__name__ if exc.__cause__ else None,
                      'cause': str(exc.__cause__) if exc.__cause__ else None}
    answer.update({'label':label, 'backend_calls':list(observations),
                   'solver_return_code':solver.rc, 'warnings':[str(x.message) for x in caught]})
    return answer, solver

normal={'nsteps':1000,'rtol':1e-9,'atol':1e-11}
results=[]
for name, integrator, opts in [
    ('vode_step_exhaustion','vode',{'nsteps':1,'rtol':1e-9,'atol':1e-11}),
    ('lsoda_step_exhaustion','lsoda',{'nsteps':1,'rtol':1e-9,'atol':1e-11}),
    ('vode_impossible_accuracy','vode',{'nsteps':1000,'rtol':1e-30,'atol':1e-30}),
    ('vode_success','vode',normal),
    ('dopri5_success','dopri5',normal),
]:
    item, _ = trial(name, integrator, opts); results.append(item)
# A failed attempt followed by a successful retry on the same solver must not
# retain an earlier has_raised_exception marker once the new run succeeds.
failed, reused = trial('reuse_first_failure','vode',{'nsteps':1,'rtol':1e-9,'atol':1e-11})
retried, _ = trial('reuse_successful_retry','vode',normal,solver=reused)
results += [failed,retried]
print(json.dumps({'source':str(source),'oracle':'A=e^(-t), B=1-e^(-t), A+B=1',
                  'results':results},indent=2))
