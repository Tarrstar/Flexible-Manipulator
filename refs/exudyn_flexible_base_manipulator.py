"""Planar Exudyn reference for flexible_base_manipulator.m.

A 5 kg floating rigid base is fixed to the beam root. The beam tip is fixed
to a four-body chain and the three revolute joints receive sinusoidal internal
torques. Running without arguments solves and displays the default animation.
"""

from pathlib import Path
import time

import exudyn as exu
from exudyn.utilities import *
import exudyn.graphics as graphics
import numpy as np


SC = exu.SystemContainer()
mbs = SC.AddSystem()

base_mass = 5.0
base_length = 0.20
base_height = 0.12
base_depth = 0.08
base_com_x = -0.5 * base_length
base_Izz_com = base_mass * (base_length**2 + base_height**2) / 12.0
base_Izz_reference = base_Izz_com + base_mass * base_com_x**2

base_node = mbs.AddNode(NodeRigidBody2D(referenceCoordinates=[0.0, 0.0, 0.0]))
base_graphics = graphics.Brick(
    centerPoint=[base_com_x, 0.0, 0.0],
    size=[base_length, base_height, base_depth],
    color=[0.18, 0.22, 0.28, 1.0],
    addEdges=True,
)
base_body = mbs.AddObject(ObjectRigidBody2D(
    nodeNumber=base_node,
    physicsMass=base_mass,
    physicsInertia=base_Izz_reference,
    physicsCenterOfMass=[base_com_x, 0.0],
    visualization=VObjectRigidBody2D(graphicsData=[base_graphics]),
))
base_marker = mbs.AddMarker(MarkerBodyRigid(
    bodyNumber=base_body,
    localPosition=[0.0, 0.0, 0.0],
))

beam_length = 1.0
beam_width = 0.032
beam_thickness = 0.0015
beam_density = 7.542e3
beam_youngs_modulus = 1.983e11
beam_elements = 5
beam_zeta = 0.0
beam_omega_n = 7.808

beam_area = beam_width * beam_thickness
beam_Izz = beam_width * beam_thickness**3 / 12.0
beam_EA = beam_youngs_modulus * beam_area
beam_EI = beam_youngs_modulus * beam_Izz
# Exudyn Cable2D uses constitutive damping.  This stiffness-proportional
# coefficient has the same damping ratio as MATLAB at beam_omega_n.
stiffness_damping = 2.0 * beam_zeta / beam_omega_n

beam_nodes = []
element_length = beam_length / beam_elements
for i in range(beam_elements + 1):
    node = mbs.AddNode(NodePoint2DSlope1(referenceCoordinates=[
        i * element_length, 0.0, 1.0, 0.0,
    ]))
    beam_nodes.append(node)
    if i:
        mbs.AddObject(ObjectANCFCable2D(
            nodeNumbers=[beam_nodes[i - 1], node],
            physicsLength=element_length,
            physicsMassPerLength=beam_density * beam_area,
            physicsBendingStiffness=beam_EI,
            physicsAxialStiffness=beam_EA,
            physicsBendingDamping=stiffness_damping * beam_EI,
            physicsAxialDamping=stiffness_damping * beam_EA,
            visualization=VObjectANCFCable2D(drawHeight=beam_thickness),
        ))

# Body reference frames are located at their left joints.  This is the
# same convention as model.Xtree and RigidBody.com in MATLAB.
link_lengths = np.array([0.2, 0.1418, 0.1448, 2.0 * 0.0504])
link_masses = np.array([1.5093, 0.1029, 0.1139, 0.1009])
link_com_x = np.array([0.1, 0.0709, 0.0723, 0.0504])
# MATLAB RigidBody receives inertia about the center of mass and shifts it
# to the left-joint reference frame. ObjectRigidBody2D expects inertia about
# that reference point when physicsCenterOfMass is nonzero.
link_Izz_com = np.array([0.015155, 0.002069, 0.002387, 0.001107])
link_Izz_reference = link_Izz_com + link_masses * link_com_x**2

link_nodes = []
link_bodies = []
link_frame_markers = []
x_left = beam_length
for i in range(4):
    node = mbs.AddNode(NodeRigidBody2D(referenceCoordinates=[
        x_left, 0.0, 0.0,
    ]))
    body_graphics = graphics.RigidLink(
        p0=[0.0, 0.0, 0.0],
        p1=[link_lengths[i], 0.0, 0.0],
        axis0=[0.0, 0.0, 1.0],
        axis1=[0.0, 0.0, 0.0],
        radius=[0.009, 0.007],
        thickness=0.012,
        width=[0.018, 0.018],
        color=graphics.color.steelblue,
    )
    body = mbs.AddObject(ObjectRigidBody2D(
        nodeNumber=node,
        physicsMass=float(link_masses[i]),
        physicsInertia=float(link_Izz_reference[i]),
        physicsCenterOfMass=[float(link_com_x[i]), 0.0],
        visualization=VObjectRigidBody2D(graphicsData=[body_graphics]),
    ))
    frame_marker = mbs.AddMarker(MarkerBodyRigid(
        bodyNumber=body,
        localPosition=[0.0, 0.0, 0.0],
    ))
    link_nodes.append(node)
    link_bodies.append(body)
    link_frame_markers.append(frame_marker)
    x_left += link_lengths[i]

# Rigidly attach the beam root to the floating base and body 1 to the tip.
marker_beam_root = mbs.AddMarker(MarkerNodeRigid(nodeNumber=beam_nodes[0]))
marker_beam_tip = mbs.AddMarker(MarkerNodeRigid(nodeNumber=beam_nodes[-1]))
mbs.AddObject(GenericJoint(
    markerNumbers=[base_marker, marker_beam_root],
    constrainedAxes=[1, 1, 0, 0, 0, 1],
    visualization=VGenericJoint(show=False),
))
mbs.AddObject(GenericJoint(
    markerNumbers=[marker_beam_tip, link_frame_markers[0]],
    constrainedAxes=[1, 1, 0, 0, 0, 1],
    visualization=VGenericJoint(show=False),
))

# Three planar revolute joints.  Dedicated 2D joints remain regular during
# repeated rotations caused by the constant torques.
for i in range(1, 4):
    parent_marker = mbs.AddMarker(MarkerBodyRigid(
        bodyNumber=link_bodies[i - 1],
        localPosition=[float(link_lengths[i - 1]), 0.0, 0.0],
    ))
    child_marker = link_frame_markers[i]
    mbs.AddObject(ObjectJointRevolute2D(
        markerNumbers=[parent_marker, child_marker],
        visualization=VObjectJointRevolute2D(drawSize=0.014),
    ))

# Equal and opposite generalized coordinate loads are internal joint torques.
joint_torque_amplitude = 0.02
joint_torque_frequency = 1.0


def make_joint_torque_user_function(sign):
    def joint_torque_user_function(mbs_local, t, static_load):
        del mbs_local, static_load
        return sign * joint_torque_amplitude * np.sin(
            2.0 * np.pi * joint_torque_frequency * t
        )

    return joint_torque_user_function


for i in range(1, 4):
    child_coordinate = mbs.AddMarker(MarkerNodeCoordinate(
        nodeNumber=link_nodes[i], coordinate=2,
    ))
    parent_coordinate = mbs.AddMarker(MarkerNodeCoordinate(
        nodeNumber=link_nodes[i - 1], coordinate=2,
    ))
    mbs.AddLoad(LoadCoordinate(
        markerNumber=child_coordinate,
        load=0.0,
        loadUserFunction=make_joint_torque_user_function(1.0),
    ))
    mbs.AddLoad(LoadCoordinate(
        markerNumber=parent_coordinate,
        load=0.0,
        loadUserFunction=make_joint_torque_user_function(-1.0),
    ))


def marker_rotation(marker_number):
    values = mbs.GetMarkerOutput(marker_number, exu.OutputVariableType.RotationMatrix)
    return np.asarray(values).reshape(3, 3)


log_data = []


def record_outputs(mbs_local, t):
    time_value = float(t)
    if log_data and abs(log_data[-1][0] - time_value) <= 1.0e-12:
        return

    root_position = np.asarray(mbs_local.GetNodeOutput(
        base_node, exu.OutputVariableType.Position
    ))
    tip_position = np.asarray(mbs_local.GetNodeOutput(
        beam_nodes[-1], exu.OutputVariableType.Position
    ))
    root_rotation_values = mbs_local.GetMarkerOutput(
        base_marker, exu.OutputVariableType.RotationMatrix
    )
    root_rotation = np.asarray(root_rotation_values).reshape(3, 3)

    joint_angles = []
    for i in range(1, 4):
        parent_rotation = marker_rotation(link_frame_markers[i - 1])
        child_rotation = marker_rotation(link_frame_markers[i])
        relative_rotation = parent_rotation.T @ child_rotation
        joint_angles.append(float(np.arctan2(
            relative_rotation[1, 0], relative_rotation[0, 0]
        )))

    local_tip = root_rotation.T @ (tip_position - root_position)
    base_theta = float(np.arctan2(root_rotation[1, 0], root_rotation[0, 0]))
    log_data.append([
        time_value,
        *joint_angles,
        float(local_tip[1]),
        float(root_position[0]),
        float(root_position[1]),
        base_theta,
        float(tip_position[0]),
        float(tip_position[1]),
    ])


def post_step_user_function(mbs_local, t):
    record_outputs(mbs_local, t)
    return True


mbs.SetPostStepUserFunction(post_step_user_function)
mbs.Assemble()
record_outputs(mbs, 0.0)

project_dir = Path(__file__).resolve().parents[1]
data_dir = project_dir / "data" / "flexible_base_manipulator"
data_dir.mkdir(parents=True, exist_ok=True)

simulation_time = 5.0
time_step = 1.0e-3
settings = exu.SimulationSettings()
settings.timeIntegration.numberOfSteps = int(round(simulation_time / time_step))
settings.timeIntegration.endTime = simulation_time
settings.timeIntegration.generalizedAlpha.spectralRadius = 0.9
settings.timeIntegration.newton.useModifiedNewton = False
settings.timeIntegration.newton.relativeTolerance = 1.0e-6
settings.timeIntegration.newton.absoluteTolerance = 1.0e-8
settings.timeIntegration.verboseMode = 1
settings.solutionSettings.writeSolutionToFile = True
settings.solutionSettings.solutionWritePeriod = 0.01
settings.solutionSettings.coordinatesSolutionFileName = str(
    data_dir / "exudyn_animation_solution.txt"
)
settings.linearSolverType = exu.LinearSolverType.EigenDense

start_time = time.perf_counter()
mbs.SolveDynamic(settings, solverType=exu.DynamicSolverType.GeneralizedAlpha)
elapsed = time.perf_counter() - start_time
record_outputs(mbs, simulation_time)

output_path = data_dir / "exudyn_floating_base.csv"
log_array = np.asarray(log_data)
log_array[:, 1:4] = np.unwrap(log_array[:, 1:4], axis=0)
np.savetxt(
    output_path,
    log_array,
    delimiter=",",
    header=(
        "time,joint_1_angle,joint_2_angle,joint_3_angle,"
        "beam_tip_elastic_y,base_x,base_y,base_theta,"
        "beam_tip_global_x,beam_tip_global_y"
    ),
    comments="",
)

print(f"Exudyn simulation completed in {elapsed:.3f} s.")
print(f"Data saved to: {output_path}")

SC.visualizationSettings.general.showHelpOnStartup = 1
SC.visualizationSettings.dialogs.openTreeView = False
SC.visualizationSettings.interactive.pauseWithSpacebar = False
SC.visualizationSettings.general.graphicsUpdateInterval = 0.02
SC.visualizationSettings.openGL.multiSampling = 4
SC.visualizationSettings.openGL.initialCenterPoint = [0.65, 0.0, 0.0]
SC.visualizationSettings.openGL.initialMaxSceneSize = 1.6
SC.visualizationSettings.window.renderWindowSize = [1280, 720]
SC.visualizationSettings.nodes.show = False
SC.visualizationSettings.markers.show = False
mbs.SolutionViewer()
