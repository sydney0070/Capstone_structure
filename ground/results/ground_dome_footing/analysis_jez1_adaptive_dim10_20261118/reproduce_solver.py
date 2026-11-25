"""Same model and force tolerance; early failed-step subdivision for JEZ-1.

Run with exactly the CLI of ground_dome_footing.py. Only solver controls change:
KrylovNewton alone, 48 iterations per attempted increment, 1000 attempts.
No wall-clock timeout. Failed increments halve in the existing continuation.
"""
import sys
from pathlib import Path
sys.path.insert(0, '/Users/sydney/Desktop/capstone/Capstone_structure/ground/ground_modeling')
import ground_dome_footing as model
original_ramp = model.static_ramp
original_test = model.ops.test
original_algorithm = model.ops.algorithm

def bounded_iteration_test(test_name, tolerance, max_iterations, *args):
    return original_test(test_name, tolerance, min(max_iterations, 48), *args)

def adaptive_ramp(*args, **kwargs):
    kwargs['algorithms'] = ('KrylovNewton',)
    kwargs['max_attempts'] = 1000
    return original_ramp(*args, **kwargs)

def expanded_krylov(name, *args):
    if name == 'KrylovNewton':
        return original_algorithm(name, '-maxDim', 10)
    return original_algorithm(name, *args)

model.ops.algorithm = expanded_krylov
model.ops.test = bounded_iteration_test
model.static_ramp = adaptive_ramp
model.main()
