"""Exudyn reference model for test_FD_flex_7dof_grasp_cb_beam.m.

The seven-link arm and the 3D geometrically exact beam use the same joint locations,
center-of-mass data, section properties, rigid grasp, and actuator torques as
the MATLAB model.  The exported file contains ``time,tip_y``.
"""

from pathlib import Path

import exudyn as exu
from exudyn.utilities import *
import exudyn.graphics as graphics
import numpy as np


SC = exu.SystemContainer()
mbs = SC.AddSystem()

g = [0.0, 0.0, 0.0]

# Same seven-link arm parameters as the MATLAB test.
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
    mass = 0.8
    inertia_diag = np.array([
        0.02,
        0.03,
        0.04,
    ])
    inertia = RigidBodyInertia(
        mass=mass,
        inertiaTensor=np.diag(inertia_diag),
        inertiaTensorAtCOM=True,
    )
    p_com = [arm_position + 0.5 * length, 0.0, 0.0]
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
        referencePosition=p_com,
        gravity=g,
        graphicsDataList=[graphics_link],
    )
    arm_bodies.append(body)
    arm_position += length

# Revolute joints are defined in the reference configuration by global
# position and axis; Exudyn computes the corresponding local marker frames.
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


# Beam section shared with the MATLAB Euler-Bernoulli model.
beam_length = 1.0
beam_width = 0.04
beam_thickness = 0.025
beam_density = 2700.0
beam_youngs_modulus = 7.0e7
beam_poisson_ratio = 0.33
beam_area = beam_width * beam_thickness
beam_Iyy = beam_thickness * beam_width**3 / 12.0
beam_Izz = beam_width * beam_thickness**3 / 12.0
beam_Ip = beam_Iyy + beam_Izz
beam_G = beam_youngs_modulus / (2.0 * (1.0 + beam_poisson_ratio))

a_rect = max(beam_width, beam_thickness)
b_rect = min(beam_width, beam_thickness)
beam_J = a_rect * b_rect**3 * (
    1.0 / 3.0
    - 0.21 * (b_rect / a_rect) * (1.0 - b_rect**4 / (12.0 * a_rect**4))
)

section_data = exu.BeamSection()
section_data.stiffnessMatrix = np.diag([
    beam_youngs_modulus * beam_area,
    1.0e4 * beam_G * beam_area,
    1.0e4 * beam_G * beam_area,
    beam_G * beam_J,
    beam_youngs_modulus * beam_Iyy,
    beam_youngs_modulus * beam_Izz,
])
section_data.massPerLength = beam_density * beam_area
section_data.inertia = beam_density * np.diag([
    beam_Ip, beam_Iyy, beam_Izz
])

arm_tip = float(np.sum(arm_lengths))
n_beam_elements = 16
element_length = beam_length / n_beam_elements
node_beam_base = None
node_beam_tip = None
previous_node = None

for element_index in range(n_beam_elements + 1):
    x = arm_tip + element_index * element_length
    node = mbs.AddNode(NodeRigidBodyEP(
        referenceCoordinates=[
            x, 0.0, 0.0,
            1.0, 0.0, 0.0, 0.0,
        ]
    ))
    if element_index == 0:
        node_beam_base = node
    if previous_node is not None:
        mbs.AddObject(ObjectBeamGeometricallyExact(
            nodeNumbers=[previous_node, node],
            physicsLength=element_length,
            sectionData=section_data,
            visualization=VBeam3D(),
        ))
    previous_node = node
    node_beam_tip = node

# A six-axis joint gives the same rigid grasp as MATLAB's fixed link 8.
marker_arm_tip = mbs.AddMarker(MarkerBodyRigid(
    bodyNumber=arm_bodies[-1],
    localPosition=[0.5 * arm_lengths[-1], 0.0, 0.0],
))
marker_beam_base = mbs.AddMarker(MarkerNodeRigid(nodeNumber=node_beam_base))
mbs.AddObject(GenericJoint(
    markerNumbers=[marker_arm_tip, marker_beam_base],
    constrainedAxes=[1, 1, 1, 1, 1, 1],
    visualization=VGenericJoint(show=False),
))


# MATLAB torque law.  The torque is applied to the child body in body-fixed
# coordinates; at q=0 the body frames coincide with the specified joint axes.
def torque_amplitude(i, t):
    amplitudes = np.array([0.35, 0.25, 0.20, 0.16, 0.12, 0.10, 0.08])
    frequencies = np.array([0.8, 0.6, 0.9, 0.7, 0.5, 0.8, 0.4])
    use_sine = np.array([True, False, True, False, True, False, True])
    value = amplitudes[i] * (np.sin if use_sine[i] else np.cos)(2.0 * np.pi * frequencies[i] * t)
    return float(value)


def make_torque_user_function(i, sign):
    def user_torque(mbs_local, t, static_load):
        return (sign * joint_axes[i] * torque_amplitude(i, t)).tolist()

    return user_torque


for i, body in enumerate(arm_bodies):
    # Internal actuator torque: +tau on the child and -tau on the parent.
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


# Record only time and the ANCF rod endpoint y-coordinate.
log_data = []


def record_tip(mbs_local, t):
    tip_position = mbs_local.GetNodeOutput(
        node_beam_tip, exu.OutputVariableType.Position
    )
    grasp_position = np.asarray(mbs_local.GetMarkerOutput(
        marker_arm_tip, exu.OutputVariableType.Position
    ))
    log_data.append([
        float(t), float(tip_position[1]), float(grasp_position[1])
    ])


def post_step_user(mbs_local, t):
    record_tip(mbs_local, t)
    return True


# Record the converged state at t; a pre-step callback sees the previous state.
mbs.SetPostStepUserFunction(post_step_user)
mbs.Assemble()
record_tip(mbs, 0.0)

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

if len(log_data) == 0 or abs(log_data[-1][0] - simulation_time) > 1.0e-12:
    record_tip(mbs, simulation_time)

log_array = np.asarray(log_data)
project_dir = Path(__file__).resolve().parents[1]
data_dir = project_dir / "data" / "arm_grasp_cb_beam"
data_dir.mkdir(parents=True, exist_ok=True)
output_path = data_dir / "exudyn_beam_tip_y.csv"
np.savetxt(
    output_path,
    log_array,
    delimiter=",",
    header="time,exudyn_tip_y,exudyn_grasp_y",
    comments="",
)

print("Exudyn 7-DOF 3D-beam grasping simulation completed.")
print(f"Data saved to: {output_path}")
