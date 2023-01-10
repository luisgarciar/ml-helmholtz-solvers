# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

from IPython.display import display
from sympy import *
from sympy.interactive import printing
printing.init_printing()
import sympy as sym

# <markdowncell>

# In this notebook we compute an exact solution to the Helmholtz boundary value problem
# \begin{align}
# -u''-k^2u&=f \text{ in } (0,1),\\
# u(0)&=0,\\
# u'(1)-iku(1)&=0.
# \end{align}
# 
# The Green function $G(x,s)$ for this problem (see Ihlenburg & Babuska, 1997) is given by
# $$G(x,s)= \frac{1}{k} \begin{cases} \sin(kx)e^{iks}, &0 \leq x \leq s,\\ 
# \sin(ks)e^{ikx} & s \leq x \leq 1.
# \end{cases} $$
# 
# We compute an exact solution to this problem for the right hand side $\tilde{f}=1$. The solution is given by
# \begin{align*}
# \tilde{u}(x)= \int_{0}^{1}G(x,s)\, ds&= \int_{0}^{x}G(x,s)\, ds+ \int_{x}^{1}G(x,s)\, ds \\
# &= \frac{e^{iks}}{k}\int_{0}^{x}\sin(ks) \, ds + \frac{\sin(kx)}{k}\int_{x}^{1}e^{iks}\, ds
# \end{align*}

# <codecell>

u = Symbol("\tilde{u}",complex=True)
k = Symbol("k",real=True,positive=True)
s = Symbol("s",real=True)
x = Symbol("x",real=True)

# <codecell>

u1 = (exp(I*k*x)/k)*integrate(sin(k*s),(s,0,x))
u2 = (sin(k*x)/k)*integrate(exp(I*k*s),(s,x,1))
u = u1+u2

# <markdowncell>

# A closed expression for the solution $\tilde{u}$ is

# <codecell>

collect(simplify(expand(u)),exp(I*k*x))

# <markdowncell>

# We check now that $\tilde{u}$ satisfies the equation $\tilde{u}''-k^2\tilde{u}=1$ and the boundary conditions: 

# <codecell>

simplify(-diff(u,x,x)-(k**2)*u)

# <codecell>

u.subs(x,0)

# <codecell>

diff(u,x).subs(x,1)-I*k*u.subs(x,1)

# <codecell>

N(u.subs([(k, 10), (x, 0.5)]))

# <codecell>


