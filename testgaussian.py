from bayespy.nodes import Categorical, Mixture, Gaussian
from bayespy.inference import VB
import numpy as np
import bayespy.plot as bpplt

hidden2 = Categorical((0.7, 0.3))
hidden1 = Mixture(hidden2, Categorical, ((0.6, 0.4), (0.1, 0.9)))

observed1 = Mixture(hidden1, Gaussian, ([-1.0], [0.9]), ([[1.3]], [[0.8]]))
observed2 = Mixture(hidden2, Gaussian, ([-0.9], [1.1]), ([[1.2]], [[0.7]]))

observed1.observe([-1.2])
observed2.observe([1.2])
Q = VB(hidden1, hidden2, observed1, observed2, tol=1e-10)
Q.update(repeat=100)
print(hidden1.get_moments())
print(hidden2.get_moments())


observed1.observe([-0.2])
observed2.observe([1.2])
Q = VB(hidden1, hidden2, observed1, observed2)
Q.update(repeat=100)
print(hidden1.get_moments())
print(hidden2.get_moments())


observed1.observe([1.0])
observed2.observe([1.2])
Q = VB(hidden1, hidden2, observed1, observed2)
Q.update(repeat=10)
print(hidden1.get_moments())
print(hidden2.get_moments())

