"""Pure-rigid reference for the 7-DOF grasped-beam comparison.

The arm, joints, actuator torques, beam dimensions, mass, center of mass,
and grasp constraint match exudyn_7dof_grasp_cb_beam.py. The flexible beam
is replaced by a uniform rigid cuboid with the same mass and inertia.
"""

from pathlib import Path

import exudyn as exu
from exudyn.utilities import *
import exudyn.graphics as graphics
import numpy as np


SC = exu.SystemContainer()
mbs = SC.AddSystem()
gravity = [0.0, 0.0, 0.0]

n_arm = 7
arm_lengths = np.array([0.28, 0.30, 0.26, 0.22, 0.18, 0.14, 0.10])
joint_axes = np.array([
    [0.0, 0.0, 1.0],
    [0.0, 1.0, 0.0],
    [0.0, 1.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
], dtype=float)

o_ground = mbs.CreateGround(referencePosition=[0.0, 0.0, 0.0])
arm_bodies = []
arm_position = 0.0

for i, length in enumerate(arm_lengths):
    mass = 0.8 + 0.1 * (i + 1)
    inertia_diag = np.array([
        0.02 + 0.003 * (i + 1),
        0.03 + 0.004 * (i + 1),
        0.04 + 0.005 * (i + 1),
    ])
    inertia = RigidBodyInertia(
        mass=mass,
        inertiaTensor=np.diag(inertia_diag),
        inertiaTensorAtCOM=True,
    )
    graphics_link = graphics.RigidLink(
        p0=[-0.5 * length, 0.0, 0.0],
        p1=[0.5 * length, 0.0, 0.0],
        axis0=[0.0, 0.0, 1.0],
        axis1=[0.0, 0.0, 0.0],
        radius=[0.025, 0.02],
        thickness=0.04,
        width=[0.05, 0.05],
    )
    body = mbs.CreateRigidBody(
        inertia=inertia,
        referencePosition=[arm_position + 0.5 * length, 0.0, 0.0],
        gravity=gravity,
        graphicsDataList=[graphics_link],
    )
    arm_bodies.append(body)
    arm_position += length

joint_position = 0.0
for i, length in enumerate(arm_lengths):
    parent_body = o_ground if i == 0 else arm_bodies[i - 1]
    mbs.CreateRevoluteJoint(
        bodyNumbers=[parent_body, arm_bodies[i]],
        position=[joint_position, 0.0, 0.0],
        axis=joint_axes[i].tolist(),
        axisRadius=0.018,
        axisLength=0.10,
    )
    joint_position += length


# Uniform rigid replacement of the flexible rectangular beam.
beam_length = 1.0
beam_width = 0.04
beam_thickness = 0.025
beam_density = 2700.0
beam_mass = beam_density * beam_length * beam_width * beam_thickness
beam_inertia = np.diag([
    beam_mass * (beam_width**2 + beam_thickness**2) / 12.0,
    beam_mass * (beam_length**2 + beam_thickness**2) / 12.0,
    beam_mass * (beam_length**2 + beam_width**2) / 12.0,
])

rigid_beam = mbs.CreateRigidBody(
    inertia=RigidBodyInertia(
        mass=beam_mass,
        inertiaTensor=beam_inertia,
        inertiaTensorAtCOM=True,
    ),
    referencePosition=[arm_position + 0.5 * beam_length, 0.0, 0.0],
    gravity=gravity,
    graphicsDataList=[graphics.RigidLink(
        p0=[-0.5 * beam_length, 0.0, 0.0],
        p1=[0.5 * beam_length, 0.0, 0.0],
        axis0=[0.0, 0.0, 1.0],
        axis1=[0.0, 0.0, 0.0],
        radius=[0.025, 0.018],
        thickness=beam_thickness,
        width=[beam_width, beam_width],
    )],
)

marker_arm_tip = mbs.AddMarker(MarkerBodyRigid(
    bodyNumber=arm_bodies[-1],
    localPosition=[0.5 * arm_lengths[-1], 0.0, 0.0],
))
marker_beam_base = mbs.AddMarker(MarkerBodyRigid(
    bodyNumber=rigid_beam,
    localPosition=[-0.5 * beam_length, 0.0, 0.0],
))
marker_beam_tip = mbs.AddMarker(MarkerBodyRigid(
    bodyNumber=rigid_beam,
    localPosition=[0.5 * beam_length, 0.0, 0.0],
))
mbs.AddObject(GenericJoint(
    markerNumbers=[marker_arm_tip, marker_beam_base],
    constrainedAxes=[1, 1, 1, 1, 1, 1],
    visualization=VGenericJoint(show=False),
))


def torque_amplitude(i, t):
    amplitudes = np.array([0.35, 0.25, 0.20, 0.16, 0.12, 0.10, 0.08])
    frequencies = np.array([0.8, 0.6, 0.9, 0.7, 0.5, 0.8, 0.4])
    use_sine = np.array([True, False, True, False, True, False, True])
    trig = np.sin if use_sine[i] else np.cos
    return float(amplitudes[i] * trig(2.0 * np.pi * frequencies[i] * t))


def make_torque_user_function(i, sign):
    def user_torque(mbs_local, t, static_load):
        return (sign * joint_axes[i] * torque_amplitude(i, t)).tolist()

    return user_torque


for i, body in enumerate(arm_bodies):
    mbs.CreateTorque(
        bodyNumber=body,
        localPosition=[0.0, 0.0, 0.0],
        bodyFixed=True,
        loadVector=[0.0, 0.0, 0.0],
        loadVectorUserFunction=make_torque_user_function(i, 1.0),
        show=False,
    )
    if i > 0:
        mbs.CreateTorque(
            bodyNumber=arm_bodies[i - 1],
            localPosition=[0.0, 0.0, 0.0],
            bodyFixed=True,
            loadVector=[0.0, 0.0, 0.0],
            loadVectorUserFunction=make_torque_user_function(i, -1.0),
            show=False,
        )


log_data = []


def record_positions(mbs_local, t):
    tip_position = np.asarray(mbs_local.GetMarkerOutput(
        marker_beam_tip, exu.OutputVariableType.Position
    ))
    grasp_position = np.asarray(mbs_local.GetMarkerOutput(
        marker_arm_tip, exu.OutputVariableType.Position
    ))
    row = [float(t), float(tip_position[1]), float(grasp_position[1])]
    if log_data and abs(log_data[-1][0] - row[0]) <= 1.0e-12:
        log_data[-1] = row
    else:
        log_data.append(row)


def pre_step_user(mbs_local, t):
    record_positions(mbs_local, t)
    return True


mbs.SetPreStepUserFunction(pre_step_user)
mbs.Assemble()
record_positions(mbs, 0.0)

simulation_time = 1.0
time_step = 1.0e-3
settings = exu.SimulationSettings()
settings.timeIntegration.numberOfSteps = int(round(simulation_time / time_step))
settings.timeIntegration.endTime = simulation_time
settings.timeIntegration.generalizedAlpha.spectralRadius = 1.0
settings.timeIntegration.newton.useModifiedNewton = True
settings.timeIntegration.newton.relativeTolerance = 1.0e-6
settings.timeIntegration.verboseMode = 1
settings.linearSolverType = exu.LinearSolverType.EigenDense

mbs.SolveDynamic(
    settings,
    solverType=exu.DynamicSolverType.GeneralizedAlpha,
)

if not log_data or abs(log_data[-1][0] - simulation_time) > 1.0e-12:
    record_positions(mbs, simulation_time)

log_array = np.asarray(log_data)
project_dir = Path(__file__).resolve().parents[1]
data_dir = project_dir / "data" / "arm_grasp_cb_beam"
data_dir.mkdir(parents=True, exist_ok=True)
output_path = data_dir / "exudyn_rigid_beam_tip_y.csv"
np.savetxt(
    output_path,
    log_array,
    delimiter=",",
    header="time,exudyn_rigid_tip_y,exudyn_rigid_grasp_y",
    comments="",
)

print("Exudyn 7-DOF rigid-beam grasping simulation completed.")
print(f"Rigid beam mass: {beam_mass:.9g} kg")
print(f"Rigid beam inertia diagonal: {np.diag(beam_inertia)} kg m^2")
print(f"Data saved to: {output_path}")
