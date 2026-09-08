"""Exudyn reference for the fixed-base flexible-beam manipulator."""

from pathlib import Path
import sys
import time

import exudyn as exu
from exudyn.utilities import *
import exudyn.graphics as graphics
import numpy as np


SC = exu.SystemContainer()
mbs = SC.AddSystem()

ground = mbs.AddObject(ObjectGround())
ground_marker = mbs.AddMarker(MarkerBodyRigid(
    bodyNumber=ground, localPosition=[0.0, 0.0, 0.0]
))

beam_length = 1.0
beam_width = 0.032
beam_thickness = 0.0015
beam_density = 7.542e3
beam_youngs_modulus = 1.983e11
beam_elements = 3
beam_zeta = 0.0
beam_omega_n = 7.808

beam_area = beam_width * beam_thickness
beam_Izz = beam_width * beam_thickness**3 / 12.0
beam_EA = beam_youngs_modulus * beam_area
beam_EI = beam_youngs_modulus * beam_Izz
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

marker_beam_root = mbs.AddMarker(MarkerNodeRigid(nodeNumber=beam_nodes[0]))
marker_beam_tip = mbs.AddMarker(MarkerNodeRigid(nodeNumber=beam_nodes[-1]))
mbs.AddObject(GenericJoint(
    markerNumbers=[ground_marker, marker_beam_root],
    constrainedAxes=[1, 1, 0, 0, 0, 1],
    visualization=VGenericJoint(show=False),
))

link_lengths = np.array([0.2, 0.1418, 0.1448, 2.0 * 0.0504])
link_masses = np.array([1.5093, 0.1029, 0.1139, 0.1009])
link_com_x = np.array([0.1, 0.0709, 0.0723, 0.0504])
link_Izz_com = np.array([0.015155, 0.002069, 0.002387, 0.001107])
link_Izz_reference = link_Izz_com + link_masses * link_com_x**2

link_nodes = []
link_bodies = []
link_frame_markers = []
x_left = beam_length
for i in range(4):
    node = mbs.AddNode(NodeRigidBody2D(referenceCoordinates=[x_left, 0.0, 0.0]))
    body_graphics = graphics.RigidLink(
        p0=[0.0, 0.0, 0.0], p1=[link_lengths[i], 0.0, 0.0],
        axis0=[0.0, 0.0, 1.0], axis1=[0.0, 0.0, 0.0],
        radius=[0.009, 0.007], thickness=0.012, width=[0.018, 0.018],
        color=graphics.color.steelblue,
    )
    body = mbs.AddObject(ObjectRigidBody2D(
        nodeNumber=node,
        physicsMass=float(link_masses[i]),
        physicsInertia=float(link_Izz_reference[i]),
        physicsCenterOfMass=[float(link_com_x[i]), 0.0],
        visualization=VObjectRigidBody2D(graphicsData=[body_graphics]),
    ))
    marker = mbs.AddMarker(MarkerBodyRigid(
        bodyNumber=body, localPosition=[0.0, 0.0, 0.0]
    ))
    link_nodes.append(node)
    link_bodies.append(body)
    link_frame_markers.append(marker)
    x_left += link_lengths[i]

mbs.AddObject(GenericJoint(
    markerNumbers=[marker_beam_tip, link_frame_markers[0]],
    constrainedAxes=[1, 1, 0, 0, 0, 1],
    visualization=VGenericJoint(show=False),
))

for i in range(1, 4):
    parent_marker = mbs.AddMarker(MarkerBodyRigid(
        bodyNumber=link_bodies[i - 1],
        localPosition=[float(link_lengths[i - 1]), 0.0, 0.0],
    ))
    mbs.AddObject(ObjectJointRevolute2D(
        markerNumbers=[parent_marker, link_frame_markers[i]],
        visualization=VObjectJointRevolute2D(drawSize=0.014),
    ))

joint_torque = 0.02
for i in range(1, 4):
    child_coordinate = mbs.AddMarker(MarkerNodeCoordinate(
        nodeNumber=link_nodes[i], coordinate=2
    ))
    parent_coordinate = mbs.AddMarker(MarkerNodeCoordinate(
        nodeNumber=link_nodes[i - 1], coordinate=2
    ))
    mbs.AddLoad(LoadCoordinate(markerNumber=child_coordinate, load=joint_torque))
    mbs.AddLoad(LoadCoordinate(markerNumber=parent_coordinate, load=-joint_torque))


def marker_rotation(marker_number):
    values = mbs.GetMarkerOutput(marker_number, exu.OutputVariableType.RotationMatrix)
    return np.asarray(values).reshape(3, 3)


rows = []


def record_outputs(mbs_local, t):
    time_value = float(t)
    if rows and abs(rows[-1][0] - time_value) <= 1.0e-12:
        return
    root = np.asarray(mbs_local.GetNodeOutput(
        beam_nodes[0], exu.OutputVariableType.Position
    ))
    tip = np.asarray(mbs_local.GetNodeOutput(
        beam_nodes[-1], exu.OutputVariableType.Position
    ))
    joint_angles = []
    for i in range(1, 4):
        relative = marker_rotation(link_frame_markers[i - 1]).T @ marker_rotation(
            link_frame_markers[i]
        )
        joint_angles.append(float(np.arctan2(relative[1, 0], relative[0, 0])))
    elastic_tip = tip - root - np.array([beam_length, 0.0, 0.0])
    rows.append([
        time_value, *joint_angles, float(elastic_tip[0]), float(elastic_tip[1]),
        float(tip[0]), float(tip[1]),
    ])


def post_step_user_function(mbs_local, t):
    record_outputs(mbs_local, t)
    return True


mbs.SetPostStepUserFunction(post_step_user_function)
mbs.Assemble()
record_outputs(mbs, 0.0)

simulation_time = 5.0
time_step = 1.0e-3
settings = exu.SimulationSettings()
settings.timeIntegration.numberOfSteps = int(round(simulation_time / time_step))
settings.timeIntegration.endTime = simulation_time
settings.timeIntegration.generalizedAlpha.spectralRadius = 0.9
settings.timeIntegration.newton.useModifiedNewton = False
settings.timeIntegration.newton.relativeTolerance = 1.0e-6
settings.timeIntegration.verboseMode = 1
settings.solutionSettings.writeSolutionToFile = False
settings.linearSolverType = exu.LinearSolverType.EigenDense

start_time = time.perf_counter()
mbs.SolveDynamic(settings, solverType=exu.DynamicSolverType.GeneralizedAlpha)
elapsed = time.perf_counter() - start_time
record_outputs(mbs, simulation_time)

data = np.asarray(rows)
data[:, 1:4] = np.unwrap(data[:, 1:4], axis=0)
project_dir = Path(__file__).resolve().parents[1]
data_dir = project_dir / "data" / "fixed_FBM"
data_dir.mkdir(parents=True, exist_ok=True)
output_path = data_dir / "exudyn_fixed_FBM.csv"
np.savetxt(
    output_path, data, delimiter=",",
    header=(
        "time,joint_1_angle,joint_2_angle,joint_3_angle,"
        "beam_tip_elastic_x,beam_tip_elastic_y,"
        "beam_tip_global_x,beam_tip_global_y"
    ),
    comments="",
)

print(f"Exudyn fixed-FBM simulation completed in {elapsed:.3f} s.")
print(f"Data saved to: {output_path}")

if "--show" in sys.argv:
    SC.visualizationSettings.openGL.initialCenterPoint = [0.8, 0.0, 0.0]
    SC.visualizationSettings.openGL.initialMaxSceneSize = 1.8
    SC.visualizationSettings.window.renderWindowSize = [1280, 720]
    SC.visualizationSettings.nodes.show = False
    SC.visualizationSettings.markers.show = False
    mbs.SolutionViewer()
