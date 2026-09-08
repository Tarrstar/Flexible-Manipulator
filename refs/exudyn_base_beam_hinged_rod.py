"""Exudyn model of a floating base, flexible beam, and hinged rigid rod."""

from pathlib import Path
import sys
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

beam_root_marker = mbs.AddMarker(MarkerNodeRigid(nodeNumber=beam_nodes[0]))
beam_tip_marker = mbs.AddMarker(MarkerNodeRigid(nodeNumber=beam_nodes[-1]))
mbs.AddObject(GenericJoint(
    markerNumbers=[base_marker, beam_root_marker],
    constrainedAxes=[1, 1, 0, 0, 0, 1],
    visualization=VGenericJoint(show=False),
))

rod_length = 0.20
rod_mass = 1.5093
rod_com_x = 0.10
rod_Izz_com = 0.015155
rod_Izz_reference = rod_Izz_com + rod_mass * rod_com_x**2
rod_node = mbs.AddNode(NodeRigidBody2D(
    referenceCoordinates=[beam_length, 0.0, 0.0]
))
rod_graphics = graphics.RigidLink(
    p0=[0.0, 0.0, 0.0],
    p1=[rod_length, 0.0, 0.0],
    axis0=[0.0, 0.0, 1.0],
    axis1=[0.0, 0.0, 0.0],
    radius=[0.009, 0.007],
    thickness=0.012,
    width=[0.018, 0.018],
    color=graphics.color.steelblue,
)
rod_body = mbs.AddObject(ObjectRigidBody2D(
    nodeNumber=rod_node,
    physicsMass=rod_mass,
    physicsInertia=rod_Izz_reference,
    physicsCenterOfMass=[rod_com_x, 0.0],
    visualization=VObjectRigidBody2D(graphicsData=[rod_graphics]),
))
rod_joint_marker = mbs.AddMarker(MarkerBodyRigid(
    bodyNumber=rod_body, localPosition=[0.0, 0.0, 0.0]
))
rod_tip_marker = mbs.AddMarker(MarkerBodyPosition(
    bodyNumber=rod_body, localPosition=[rod_length, 0.0, 0.0]
))
mbs.AddObject(ObjectJointRevolute2D(
    markerNumbers=[beam_tip_marker, rod_joint_marker],
    visualization=VObjectJointRevolute2D(drawSize=0.014),
))

# Equal and opposite torques make the hinge actuation internal.
hinge_torque = 0.02
mbs.AddLoad(LoadTorqueVector(
    markerNumber=rod_joint_marker,
    loadVector=[0.0, 0.0, hinge_torque],
    bodyFixed=False,
))
mbs.AddLoad(LoadTorqueVector(
    markerNumber=beam_tip_marker,
    loadVector=[0.0, 0.0, -hinge_torque],
    bodyFixed=False,
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
    beam_tip_position = np.asarray(mbs_local.GetNodeOutput(
        beam_nodes[-1], exu.OutputVariableType.Position
    ))
    rod_tip_position = np.asarray(mbs_local.GetMarkerOutput(
        rod_tip_marker, exu.OutputVariableType.Position
    ))
    base_rotation = marker_rotation(base_marker)
    beam_tip_rotation = marker_rotation(beam_tip_marker)
    rod_rotation = marker_rotation(rod_joint_marker)
    relative_rotation = beam_tip_rotation.T @ rod_rotation
    hinge_angle = float(np.arctan2(
        relative_rotation[1, 0], relative_rotation[0, 0]
    ))
    base_theta = float(np.arctan2(base_rotation[1, 0], base_rotation[0, 0]))
    local_beam_tip = base_rotation.T @ (beam_tip_position - base_position)
    rows.append([
        time_value,
        hinge_angle,
        base_theta,
        float(beam_tip_position[1]),
        float(rod_tip_position[1]),
        float(local_beam_tip[1]),
        float(base_position[0]),
        float(base_position[1]),
        float(beam_tip_position[0]),
        float(rod_tip_position[0]),
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
data[:, 1:3] = np.unwrap(data[:, 1:3], axis=0)
project_dir = Path(__file__).resolve().parents[1]
data_dir = project_dir / "data" / "base_beam_hinged_rod"
data_dir.mkdir(parents=True, exist_ok=True)
output_path = data_dir / "exudyn_base_beam_hinged_rod.csv"
np.savetxt(
    output_path,
    data,
    delimiter=",",
    header=(
        "time,hinge_angle,base_theta,beam_tip_global_y,rod_tip_global_y,"
        "beam_tip_elastic_y,base_x,base_y,beam_tip_global_x,rod_tip_global_x"
    ),
    comments="",
)

print(f"Exudyn base-beam-hinged-rod simulation completed in {elapsed:.3f} s.")
print(f"Data saved to: {output_path}")

if "--show" in sys.argv:
    SC.visualizationSettings.openGL.initialCenterPoint = [0.5, 0.0, 0.0]
    SC.visualizationSettings.openGL.initialMaxSceneSize = 1.5
    SC.visualizationSettings.window.renderWindowSize = [1280, 720]
    SC.visualizationSettings.nodes.show = False
    SC.visualizationSettings.markers.show = False
    mbs.SolutionViewer()
