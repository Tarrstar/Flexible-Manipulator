"""Exudyn reference for a floating rigid base hinged to a flexible beam."""

from pathlib import Path
import sys
import time

import exudyn as exu
from exudyn.utilities import *
import exudyn.graphics as graphics
import numpy as np


SC = exu.SystemContainer()
mbs = SC.AddSystem()

# The body reference frame is at the beam hinge; its center of mass is left of it.
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
    color=[0.18, 0.22, 0.28, 1.0], addEdges=True,
)
base_body = mbs.AddObject(ObjectRigidBody2D(
    nodeNumber=base_node,
    physicsMass=base_mass,
    physicsInertia=base_Izz_reference,
    physicsCenterOfMass=[base_com_x, 0.0],
    visualization=VObjectRigidBody2D(graphicsData=[base_graphics]),
))
base_marker = mbs.AddMarker(MarkerBodyRigid(
    bodyNumber=base_body, localPosition=[0.0, 0.0, 0.0]
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
mbs.AddObject(ObjectJointRevolute2D(
    markerNumbers=[base_marker, marker_beam_root],
    visualization=VObjectJointRevolute2D(drawSize=0.018),
))

# Equal and opposite torques form one internal hinge actuation.
hinge_torque = 0.02
mbs.AddLoad(LoadTorqueVector(
    markerNumber=marker_beam_root,
    loadVector=[0.0, 0.0, hinge_torque], bodyFixed=False,
))
mbs.AddLoad(LoadTorqueVector(
    markerNumber=base_marker,
    loadVector=[0.0, 0.0, -hinge_torque], bodyFixed=False,
))


def marker_rotation(marker_number):
    values = mbs.GetMarkerOutput(marker_number, exu.OutputVariableType.RotationMatrix)
    return np.asarray(values).reshape(3, 3)


rows = []


def record_outputs(mbs_local, t):
    time_value = float(t)
    if rows and abs(rows[-1][0] - time_value) <= 1.0e-12:
        return
    base_position = np.asarray(mbs_local.GetNodeOutput(
        base_node, exu.OutputVariableType.Position
    ))
    root_position = np.asarray(mbs_local.GetNodeOutput(
        beam_nodes[0], exu.OutputVariableType.Position
    ))
    tip_position = np.asarray(mbs_local.GetNodeOutput(
        beam_nodes[-1], exu.OutputVariableType.Position
    ))
    base_rotation = marker_rotation(base_marker)
    beam_rotation = marker_rotation(marker_beam_root)
    relative_rotation = base_rotation.T @ beam_rotation
    hinge_angle = float(np.arctan2(relative_rotation[1, 0], relative_rotation[0, 0]))
    local_tip = beam_rotation.T @ (tip_position - root_position)
    base_theta = float(np.arctan2(base_rotation[1, 0], base_rotation[0, 0]))
    rows.append([
        time_value, hinge_angle,
        float(local_tip[0] - beam_length), float(local_tip[1]),
        float(base_position[0]), float(base_position[1]), base_theta,
        float(tip_position[0]), float(tip_position[1]),
    ])


def post_step_user_function(mbs_local, t):
    record_outputs(mbs_local, t)
    return True


mbs.SetPostStepUserFunction(post_step_user_function)
mbs.Assemble()
record_outputs(mbs, 0.0)

simulation_time = 2.0
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
data[:, 1] = np.unwrap(data[:, 1])
data[:, 6] = np.unwrap(data[:, 6])
project_dir = Path(__file__).resolve().parents[1]
data_dir = project_dir / "data" / "hinged_floating_beam"
data_dir.mkdir(parents=True, exist_ok=True)
output_path = data_dir / "exudyn_hinged_floating_beam.csv"
np.savetxt(
    output_path, data, delimiter=",",
    header=(
        "time,hinge_angle,beam_tip_elastic_x,beam_tip_elastic_y,"
        "base_x,base_y,base_theta,beam_tip_global_x,beam_tip_global_y"
    ),
    comments="",
)

print(f"Exudyn hinged floating-beam simulation completed in {elapsed:.3f} s.")
print(f"Data saved to: {output_path}")

if "--show" in sys.argv:
    SC.visualizationSettings.openGL.initialCenterPoint = [0.4, 0.0, 0.0]
    SC.visualizationSettings.openGL.initialMaxSceneSize = 1.4
    SC.visualizationSettings.window.renderWindowSize = [1280, 720]
    SC.visualizationSettings.nodes.show = False
    SC.visualizationSettings.markers.show = False
    mbs.SolutionViewer()
