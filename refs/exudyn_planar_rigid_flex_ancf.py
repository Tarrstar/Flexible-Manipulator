"""Planar rigid-flexible double pendulum for comparison with FD_flex.

The first link is rigid and the second link is a planar ANCF Cable2D.  The
flexible-link root is position-constrained to the rigid-link tip while its
out-of-plane rotation is left free, representing the second revolute joint.
The output is written as ``time,exudyn_tip_y``.
"""

from pathlib import Path

import exudyn as exu
from exudyn.utilities import *
import exudyn.graphics as graphics
import numpy as np


SC = exu.SystemContainer()
mbs = SC.AddSystem()

g = [0.0, 0.0, 0.0]
L1 = 0.50
L2 = 0.50

# Same torque functions as the MATLAB script.
def tau1(t):
    return 0.1


def tau2(t):
    return 0.0


o_ground = mbs.CreateGround(referencePosition=[0.0, 0.0, 0.0])

# Rigid first link: mass, COM and principal inertia match the MATLAB model.
rigid_mass = 1.0
rigid_inertia = RigidBodyInertia(
    mass=rigid_mass,
    inertiaTensor=np.diag([0.01, 0.01, rigid_mass * L1**2 / 12.0]),
    # The Exudyn rigid-body reference position is at the COM below.
    # Keeping com=0 avoids counting the L1/2 offset twice.
    com=[0.0, 0.0, 0.0],
    inertiaTensorAtCOM=True,
)
rigid_graphics = graphics.RigidLink(
    p0=[-L1 / 2.0, 0.0, 0.0],
    p1=[L1 / 2.0, 0.0, 0.0],
    axis0=[0.0, 0.0, 1.0],
    axis1=[0.0, 0.0, 0.0],
    radius=[0.035, 0.025],
    thickness=0.05,
    width=[0.06, 0.06],
)
body_rigid = mbs.CreateRigidBody(
    inertia=rigid_inertia,
    referencePosition=[L1 / 2.0, 0.0, 0.0],
    gravity=g,
    graphicsDataList=[rigid_graphics],
)

# First revolute joint, around the global z-axis.
mbs.CreateRevoluteJoint(
    bodyNumbers=[o_ground, body_rigid],
    position=[0.0, 0.0, 0.0],
    axis=[0.0, 0.0, 1.0],
    axisRadius=0.025,
    axisLength=0.12,
)

# Planar ANCF flexible link.
beam_width = 0.04
beam_thickness = 0.025
beam_density = 2700.0
beam_youngs_modulus = 7.0e7
beam_area = beam_width * beam_thickness
beam_I = beam_width * beam_thickness**3 / 12.0
beam_template = Cable2D(
    physicsMassPerLength=beam_density * beam_area,
    physicsAxialStiffness=beam_youngs_modulus * beam_area,
    physicsAxialDamping=0.0,
    physicsBendingStiffness=beam_youngs_modulus * beam_I,
    physicsBendingDamping=0.0,
    visualization=VCable2D(drawHeight=beam_thickness),
)

beam_data = GenerateStraightBeam(
    mbs,
    positionOfNode0=[L1, 0.0, 0.0],
    positionOfNode1=[L1 + L2, 0.0, 0.0],
    numberOfElements=4,
    beamTemplate=beam_template,
    gravity=g,
    fixedConstraintsNode0=None,
    fixedConstraintsNode1=None,
)
node_beam_root = beam_data[0][0]
node_beam_tip = beam_data[0][-1]

# The second joint locks x/y translation and leaves the planar z rotation free.
marker_rigid_tip = mbs.AddMarker(
    MarkerBodyRigid(bodyNumber=body_rigid, localPosition=[L1 / 2.0, 0.0, 0.0])
)
marker_beam_root = mbs.AddMarker(MarkerNodeRigid(nodeNumber=node_beam_root))
mbs.AddObject(GenericJoint(
    markerNumbers=[marker_rigid_tip, marker_beam_root],
    constrainedAxes=[1, 1, 0, 0, 0, 0],
    visualization=VGenericJoint(axesRadius=0.018, axesLength=0.08),
))


def uf_tau1(mbs_local, t, static_load):
    # FD_flex uses spatial coordinate transforms where positive Rz(q)
    # corresponds to a negative geometric rotation in Xpt; reverse the
    # global Exudyn torque sign to use the same convention.
    return [0.0, 0.0, float(tau1(t))]


def uf_tau2_positive(mbs_local, t, static_load):
    return [0.0, 0.0, float(tau2(t))]


def uf_tau2_negative(mbs_local, t, static_load):
    return [0.0, 0.0, float(tau2(t))]


# Joint 1 torque: the ground reaction is implicit because the base is fixed.
mbs.CreateTorque(
    bodyNumber=body_rigid,
    localPosition=[0.0, 0.0, 0.0],
    bodyFixed=True,
    loadVector=[0.0, 0.0, 0.0],
    loadVectorUserFunction=uf_tau1,
    show=False,
)

# Joint 2 torque: equal and opposite torques on the flexible root and rigid
# link implement an internal generalized joint torque.
mbs.AddLoad(Torque(
    markerNumber=marker_beam_root,
    loadVector=[0.0, 0.0, 0.0],
    bodyFixed=False,
    loadVectorUserFunction=uf_tau2_positive,
))
mbs.AddLoad(Torque(
    markerNumber=marker_rigid_tip,
    loadVector=[0.0, 0.0, 0.0],
    bodyFixed=False,
    loadVectorUserFunction=uf_tau2_negative,
))


log_data = []


def record_tip(mbs_local, t):
    p_tip = mbs_local.GetNodeOutput(node_beam_tip, exu.OutputVariableType.Position)
    log_data.append([float(t), float(p_tip[1])])


def pre_step_user(mbs_local, t):
    record_tip(mbs_local, t)
    return True


mbs.SetPreStepUserFunction(pre_step_user)
mbs.Assemble()
record_tip(mbs, 0.0)

simulation_time = 1.0
time_step = 1.0e-3
settings = exu.SimulationSettings()
settings.timeIntegration.numberOfSteps = int(round(simulation_time / time_step))
settings.timeIntegration.endTime = simulation_time
settings.timeIntegration.newton.useModifiedNewton = True
settings.timeIntegration.newton.relativeTolerance = 1.0e-7
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
data_dir = project_dir / "data" / "planar_rigid_flex"
data_dir.mkdir(parents=True, exist_ok=True)
output_path = data_dir / "exudyn_endpoint_y.csv"
np.savetxt(
    output_path,
    log_array,
    delimiter=",",
    header="time,exudyn_tip_y",
    comments="",
)

print("Exudyn planar rigid-flexible double pendulum completed.")
print(f"Data saved to: {output_path}")
