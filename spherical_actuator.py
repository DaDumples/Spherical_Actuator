from numpy import *
from scipy.integrate import ode
import sys
sys.path.append('../Aero_Funcs')
from Aero_Funcs import *
from Controls_Funcs import *
import matplotlib.pyplot as pltimport
from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import tplquad
from pyquaternion import Quaternion


def force_torque_per_input(theta, r, phi, density, mu, mag_m, x_C):
    R = r*array([sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)])
    J = density*array([-sin(phi), cos(phi), 0])
    B = mu/(4*pi)*((3*R*dot(x_C, R)*mag_m)/norm(R)**5 - mag_m*x_C/norm(R)**3)

    
    Fk = cross(J, B)
    Tk = cross(R, cross(J, B))
    return hstack([Fk, Tk])


def propagate(t, state, mu, Inertia, I_rotor, num_coils, N, Ra, Rb, theta_a,
                    theta_b, rotor_r, rotor_mass,
                    actuator_center, Pgain, Dgain, Kp, Kd, mag_m, As, A_dipole,
                    integrator_t, torques, command_torques,inputs):

    omega = state[0:3]
    eps = state[3:6]
    eta = state[6]
    w_rotor[7:10]
    E_rotor = state[10:13]
    n_rotor = state[13]
    actuator_pos = state[14:17]
    actuator_vel = state[17:]

    eci2body = quat2dcm(E, n)
    eci2rotor = quat2dcm(E_rotor, n_rotor)

    Tc = -Pgain*eps - Dgain@omega
    Fc = -Kp*(actuator_pos - eci2body.T@actuator_center) - Kd*(actuator_vel - eci2body.T@cross(omega, actuator_center))

    

    dipole_axis = eci2body@eci2rotor.T@A_dipole


    Kf = []
    Kt = []

    dtheta = (theta_b - theta_a)/2
    dr = (Rb - Ra)/2
    dphi = 2*pi/2

    density = 2*N/((Rb**2 - Ra**2)*(theta_b - theta_a))
    x_C = array([1,0,0])

    #For each actuator calculate the volumetric integral for the force and torque vectors

    for A in As:
        # y_b2c = cross(A/norm(A), dipole_axis)
        # C_coil2body = vstack([dipole_axis, y_b2c, A/norm(A)]).T
        z_b2c = cross(dipole_axis, A/norm(A))
        C_coil2body = vstack([dipole_axis, A/norm(A), z_b2c]).T
        Fk = zeros((3,))
        Tk = zeros((3,))

        for theta in linspace(theta_a, theta_b, 3)[:-1]:
            for r in linspace(Ra, Rb, 3)[:-1]:
                for phi in linspace(-pi, pi, 3)[:-1]:

                    R = r*array([sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)])
                    J = density*array([-sin(phi), cos(phi), 0])

                    B = mu/(4*pi)*((3*R*dot(x_C, R)*mag_m)/norm(R)**5 - mag_m*x_C/norm(R)**3)
                    #print(R, J, B)
                    scalar = (r**2)*sin(theta)*dtheta*dr*dphi
                    Fk += crux(J)@B*scalar
                    Tk += crux(R)@crux(J)@B*scalar

        


        Kf.append(C_coil2body@Fk)
        Kt.append(C_coil2body@Tk)

    Kf = vstack(Kf).T
    Kt = vstack(Kt).T

    big_K = vstack([Kt, Kf])
    big_C = hstack([Tc, Fc])

    # for row in big_K:
    #     print(row)

    #Calculate the pseudo inverse of the "Big_K" matrix to find the smallest possible input vector that yeilds desired response
    ic = linalg.pinv(big_K)@big_C

    Tc_actual = Kt@ic
    Fc_actual = Kf@ic

    # ic = linalg.pinv(Kt)@Tc
    # Tc_actual = Kt@ic
    # Fc_actual = array([0,0,0])

    integrator_t.append(t)
    torques.append(Tc_actual)
    command_torques.append(Tc)
    inputs.append(ic)

    dw = inv(Inertia)@(Tc_actual - crux(omega)@Inertia@omega)
    dE = .5*(eta*identity(3) + crux(eps))@omega
    dn = -.5*dot(eps, omega)

    dw_rotor = inv(I_rotor)@(-Tc_actual - crux(w_rotor)@I_rotor@w_rotor)
    dE_rotor = .5*(n_rotor*identity(3) + crux(E_rotor))@w_rotor
    dn_rotor = -.5*dot(E_rotor, w_rotor)

    dpos = actuator_vel
    dvel = Fc_actual/rotor_mass
    #hstack is how you concatenate matrices horizontally in python
    #print(hstack([dw, dE, dn, dw_rotor, dE_rotor, dn_rotor, dpos, dvel]))
    return hstack([dw, dE, dn, dw_rotor, dE_rotor, dn_rotor, dpos, dvel])

#End of the propagation function


Inertia = diag([.02, .02, .08])
num_coils = 6
N = 334
Ra = .02
Rb = .024
theta_a = radians(2.5)
theta_b = radians(20)
rotor_r = .018
rotor_mass = .182
actuator_center = array([.1, .1, .1])
actuator_pos = actuator_center + array([1e-2]*3)
actuator_vel = array([0,0,0])

Pgain = .008
Dgain = diag([.001, .001, .003])

# damping_ratio = .65
# settling_time = 60
# w_n = 4.4/(damping_ratio*settling_time)
# Pgain = 2*(w_n**2)*Inertia
# Dgain = 2*damping_ratio*w_n*Inertia

print(Pgain)
print(Dgain)

Kp = 550
Kd = 110
mag_m = 25.66
I_rotor = diag([2/5*rotor_mass*rotor_r**2]*3)
mu = 4*pi*1e-7

psy = (sqrt(5)-1)/2
a = 1/sqrt(3)
b = a/psy
c = a*psy

As = []
for i in [-1, 1]:
    for j in [-1, 1]:
        As.append(array([0, i*c, j*b]))
        As.append(array([i*c, j*b, 0]))
        As.append(array([i*b, 0, j*c]))
        for k in [-1, 1]:
            As.append([i*a, j*a, k*a])

As = vstack(As)
print(As)
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# ax.plot(As[:,0], As[:,1], As[:,2],'.')
# plt.show()

A_dipole = array([0,0,1])

w = array([radians(3)]*3)
q = Quaternion.random()
E = array([q[1], q[2], q[3]])
n = q[0]
w_rotor = array([0,0,0])
E_rotor = array([0,0,0])
n_rotor = 1

integrator_t = []
torques = []
command_torques = []
inputs = []
state = hstack([w, E, n, w_rotor, E_rotor, n_rotor, actuator_pos, actuator_vel])
solver = ode(propagate)
solver.set_integrator('lsoda', rtol = 1e-8, atol = 1e-8)
solver.set_initial_value(state, 0)
solver.set_f_params(mu, Inertia, I_rotor, num_coils, N, Ra, Rb, theta_a,
                    theta_b, rotor_r, rotor_mass,
                    actuator_center, Pgain, Dgain, Kp, Kd, mag_m, As, A_dipole,
                    integrator_t, torques, command_torques, inputs)

dt = .1
t = []
newstate = []
tspan = 60*5
percent = 1
while solver.successful() and solver.t < tspan:
    t.append(solver.t)
    newstate.append(solver.y)
    solver.integrate(solver.t + dt)

    if solver.t/tspan*100 > percent:
        print(percent,'%')
        percent += 1

t = hstack(t)
newstate = vstack(newstate)

#hstack([w, E, n, w_rotor, E_rotor, n_rotor, actuator_pos, actuator_vel])

savetxt('newstate.txt', newstate)
savetxt('t.txt', t)
savetxt('tc_actual.txt', torques)
savetxt('tc_control.txt', command_torques)
savetxt('integrator_t.txt', integrator_t)
savetxt('inputs.txt',inputs)

plt.figure()
plt.plot(t, newstate[:,3:7])
plt.title('quaternion')

plt.figure()
plt.plot(t, newstate[:, 0:3])
plt.title('angular velocity')

plt.figure()
plt.plot(t, newstate[:,14:17])
plt.title('position')

plt.figure()
plt.plot(t, newstate[:,17:])
plt.title('velocity')

plt.figure()
plt.plot(t, newstate[:,7:10])
plt.title('rotor angular velocity')

plt.figure()
plt.plot(integrator_t, torques, 'b')
plt.plot(integrator_t, command_torques, 'r')
plt.title('Torques')
plt.show()