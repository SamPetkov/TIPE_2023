from math import *
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Initial conditions
R_c = 3.40e-3
R_b = 3.00e-3
l = 14.5e-3
m_c = 0.100e-3
m_b = 0.893e-3
mu = 0.35
alpha = 0.244
g = 9.81


tf = 0.2
N = 1000000

# Moments
Mom_b = m_b*g*R_b*cos(alpha)
Mom_c = m_c*g*(R_c*sin(alpha) - (l/2)*cos(alpha))

def momC(phi):
    return m_c*g*(R_c*sin(alpha) - (l/2)*cos(alpha+phi))

I_b = (7/5)*m_b*R_b**2


def dist2(phi):
    return (l**2)/4+R_c**2 + l*R_c*sin(phi)

def inertiaC(phi):
    return m_c*((1/(2+l/R_c))*((l*(R_c**2+2*l**2))/(2*R_c) + (5*R_c**2+3*l**2)/(12)) + dist2(phi))

def phiddot(phi):
    return (Mom_b-momC(phi))/(I_b+inertiaC(phi)*(R_b/R_c))

def Rot(P):
    [phi,omega] = P
    epsl = phiddot(phi)
    return [omega,epsl]

def Roll(X):
    [x,v] = X
    gamma = (5/7)*g*sin(alpha)
    return [v,gamma]


def Collision(omega):
    return omega*((I_b)/(inertiaC(0) + I_b*(R_c/R_b)))


t=[0]
x=[0]
v=[0]

phi0 = 0

phi=[phi0]
phid=[0]

k=0
# Capsule/ball
phase=0

error = 1e-10
# Euler method
for i in range(N):
    xl=x[len(x)-1]
    vl = v[len(v)-1]
    phil=phi[len(phi)-1]
    phidl=phid[len(phid)-1]

    Xl = [xl,vl]
    Pl = [phil,phidl]

    x0 = k*(pi*R_c+l)

    dt = tf/N

    t.append(i*dt)

    if(xl - x0 >= pi*R_c +l - error):

        k+=1
        phase=0

        x0 = k*(pi*R_c+l) + phase*pi*R_c

        phi.append(0)
        omega = Collision(vl/R_b)
        phid.append(omega)

        x.append(x0+error)
        v.append(omega*R_c)

    elif( xl - x0 > pi*R_c ):

        [vl,gamma] = Roll(Xl)

        xnew = xl+vl*dt
        vnew = vl + gamma*dt

        x.append(xnew)
        v.append(vnew)

    elif( phil >= pi - error):
        phase = 1

        x0 = k*(pi*R_c+l) + phase*pi*R_c
        vnew = 0

        phinew = 0
        phidnew = 0

        phi.append(phinew)
        phid.append(phidnew)

        x.append(x0+error)
        v.append(vnew)
    else:
        [omega,epsl] = Rot(Pl)

        phinew = phil + omega*dt
        phidnew = phidl + epsl*dt

        xnew = x0 + phinew*R_c
        vnew = vl + R_c*epsl*dt

        phi.append(phinew)
        phid.append(phidnew)

        x.append(xnew)
        v.append(vnew)





plt.figure("Temps",facecolor = "white")
plt.title("Evolution temporelle")
plt.plot(t,x,'b-',lw=2)
plt.xlabel("Temps en unite de $T_0$")
plt.ylabel("Position angulaire en radian")
plt.show()
