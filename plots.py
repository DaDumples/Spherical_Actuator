from numpy import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from numpy.linalg import *
import sys
sys.path.append('../../Aero_Funcs')
#from Aero_Funcs import *
#from Controls_Funcs import *

def quat2dcm(n, E):
    """
    generates the active rotation matrix from a quaternion.
    :param n: scalar part
    :param E: vector part

    """
    
    frame_rotation = (2*n**2 - 1)*identity(3) + 2*outer(E,E) - 2*n*crux(E) 

    return frame_rotation

def crux(A):
	return array([[0, -A[2], A[1]],
	              [A[2], 0, -A[0]],
	              [-A[1], A[0], 0]])

newstate = loadtxt('newstate.txt')
t = loadtxt('t.txt')
integrator_t = loadtxt('integrator_t.txt')
tc_actual = loadtxt('tc_actual.txt')
tc_control = loadtxt('tc_control.txt')
inputs = loadtxt('inputs.txt')


E1 = newstate[:,10]
E2 = newstate[:,11]
E3 = newstate[:,12]
n = newstate[:,13]
E1 = interp(integrator_t, t, E1)
E2 = interp(integrator_t, t, E2)
E3 = interp(integrator_t, t, E3)
n = interp(integrator_t, t, n)

E = vstack([E1, E2, E3]).T

torque_angle = []
for E, n, torque in zip(E, n, tc_control):
    axis = quat2dcm(n,E)@array([0,0,1])
    torque_angle.append(arccos(dot(axis, torque)/norm(torque)))




fig, (ax1, ax2) = plt.subplots(2, 1, sharex = True)
ax1.plot(integrator_t, tc_control)
ax2.plot(integrator_t, tc_actual)
fig.suptitle('Controller Command vs Actuator Output')
ax1.set_xlabel('Time [s]')
ax1.set_ylabel('Torque [Nm]')
ax2.set_ylabel('Torque [Nm]')
ax1.set_title('Controller Command')
ax2.set_title('Actuator Output')
plt.subplots_adjust(top=0.881,
                    bottom=0.081,
                    left=0.149,
                    right=0.971,
                    hspace=0.316,
                    wspace=0.2)
ax1.legend(['X', 'Y', 'Z'])
ax2.set_ylim(ax1.get_ylim())
plt.savefig('Controller Command vs Actuator Output.png', bbox_inches = 'tight')


plt.figure()
plt.plot(t, newstate[:,0:3])
plt.xlabel('Time [s]')
plt.ylabel('Angular Velocity')
plt.title('SC Angular Velocity vs Time')
plt.savefig('Angular Velocity.png', bbox_inches = 'tight')

plt.figure()
plt.plot(t, newstate[:,3:7])
plt.xlabel('Time [s]')
plt.ylabel('Quaternions')
plt.title('SC Quaternion vs Time')
plt.savefig('Quaternions.png', bbox_inches = 'tight')

plt.figure()
plt.plot(t, newstate[:,7:10])
plt.xlabel('Time [s]')
plt.ylabel('Angular Velocity [rads/s]')
plt.title('Rotor Angular Velocity vs Time')
plt.savefig('rotor angular velocity.png', bbox_inches = 'tight')

plt.figure()
plt.plot(t, norm(newstate[:,7:10], axis = 1)/2*pi*60)
plt.xlabel('Time [s]')
plt.ylabel('Angular Velocity [RPM]')
plt.title('Rotor Angular Velocity vs Time')
plt.savefig('rotor angular velocity magnitude.png', bbox_inches = 'tight')

plt.figure()
plt.plot(integrator_t, torque_angle)
plt.title('Angle from Torque Command to Dipole Axis')
plt.xlabel('Time [s]')
plt.ylabel('Angle [Radians]')
plt.savefig('Angle from Torque Command to Dipole Axis.png',bbox_inches = 'tight')


fig, (ax1, ax2) = plt.subplots(2, 1, sharex = True)
ax1.plot(t, newstate[:,0:3])
ax1.set_ylabel('Angular Rate [rads/s]')
ax2.plot(t, newstate[:,3:7])
ax2.set_ylabel('Quaternion')
ax2.set_xlabel('Time [s]')
# ax3.plot(integrator_t, torque_angle)
# ax3.set_ylabel('Angle between Command Torque and Dipole')
# ax3.set_xlabel('Time [s]')
fig.suptitle('Spacecraft State vs Time')
plt.savefig('States.png',bbox_inches = 'tight')


plt.figure()
plt.plot(integrator_t, inputs)

plt.figure()
plt.plot(integrator_t, norm(inputs, axis = 1))
plt.xlabel('Time [s]')
plt.ylabel('Current [A]')
plt.title('Controller Current vs Time')
plt.savefig('Controller Current vs Time.png', bbox_inches = 'tight')


fig, (ax1, ax2) = plt.subplots(2, 1, sharex = True)
ax1.plot(t, newstate[:,14:17])
ax2.set_xlabel('Time [s]')
ax1.set_ylabel('Position [m]')
ax1.set_title('Position')
ax2.plot(t, newstate[:,17:])
ax2.set_ylabel('Velocity [m/s]')
ax2.set_title('Velocity')
fig.suptitle('Rotor Position and Velocity')
plt.savefig('Actuator Position and Velocity.png',bbox_inches = 'tight')

fig = plt.figure()
ax = fig.add_subplot(111, projection = '3d')
ax.plot(newstate[:,14],newstate[:,15],newstate[:,16])

plt.show()