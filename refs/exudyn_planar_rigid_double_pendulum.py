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

beam_width = 0.04
beam_thickness = 0.025
beam_density = 2700.0

# --- 1. 修正 Link 1 惯量（满足物理三角不等式） ---
W1, T1 = 0.06, 0.05
m1 = 1.0
I1_x = (1.0 / 12.0) * m1 * (W1**2 + T1**2)
I1_y = (1.0 / 12.0) * m1 * (L1**2 + T1**2)
I1_z = (1.0 / 12.0) * m1 * (L1**2 + W1**2)
rigid_inertia = RigidBodyInertia(
    mass=m1,
    inertiaTensor=np.diag([I1_x, I1_y, I1_z]),
    inertiaTensorAtCOM=True,
)

# --- 2. Link 2 惯量 ---
rod_mass = beam_density * (beam_width * beam_thickness) * L2
I2_x = (1.0 / 12.0) * rod_mass * (beam_width**2 + beam_thickness**2)
I2_y = (1.0 / 12.0) * rod_mass * (L2**2 + beam_thickness**2)
I2_z = (1.0 / 12.0) * rod_mass * (L2**2 + beam_width**2)
rod_inertia = RigidBodyInertia(
    mass=rod_mass,
    inertiaTensor=np.diag([I2_x, I2_y, I2_z]),
    inertiaTensorAtCOM=True,
)

o_ground = mbs.CreateGround(referencePosition=[0.0, 0.0, 0.0])

link1_graphics = graphics.RigidLink(
    p0=[-L1 / 2.0, 0.0, 0.0],
    p1=[L1 / 2.0, 0.0, 0.0],
    axis0=[0.0, 0.0, 1.0],
    axis1=[0.0, 0.0, 0.0],
    radius=[0.035, 0.025],
    thickness=T1,
    width=[W1, W1],
)
link2_graphics = graphics.RigidLink(
    p0=[-L2 / 2.0, 0.0, 0.0],
    p1=[L2 / 2.0, 0.0, 0.0],
    axis0=[0.0, 0.0, 1.0],
    axis1=[0.0, 0.0, 0.0],
    radius=[0.025, 0.018],
    thickness=beam_thickness,
    width=[beam_width, beam_width],
)

body1 = mbs.CreateRigidBody(
    inertia=rigid_inertia,
    referencePosition=[L1 / 2.0, 0.0, 0.0],
    gravity=g,
    graphicsDataList=[link1_graphics],
)
body2 = mbs.CreateRigidBody(
    inertia=rod_inertia,
    referencePosition=[L1 + L2 / 2.0, 0.0, 0.0],
    gravity=g,
    graphicsDataList=[link2_graphics],
)

mbs.CreateRevoluteJoint(
    bodyNumbers=[o_ground, body1],
    position=[0.0, 0.0, 0.0],
    axis=[0.0, 0.0, 1.0],
    axisRadius=0.025,
    axisLength=0.12,
)
mbs.CreateRevoluteJoint(
    bodyNumbers=[body1, body2],
    position=[L1, 0.0, 0.0],
    axis=[0.0, 0.0, 1.0],
    axisRadius=0.02,
    axisLength=0.10,
)

def tau1(t):
    return 0.1

def uf_tau1(mbs_local, t, static_load):
    return [0.0, 0.0, float(tau1(t))]

mbs.CreateTorque(
    bodyNumber=body1,
    localPosition=[0.0, 0.0, 0.0],
    bodyFixed=True,
    loadVector=[0.0, 0.0, 0.0],
    loadVectorUserFunction=uf_tau1,
    show=False,
)

log_data = []

def record_tip(mbs_local, t):
    # --- 核心修正：使用 GetObjectOutput 直接提取 ObjectIndex 的输出 ---
    position = mbs_local.GetObjectOutput(body2, exu.OutputVariableType.Position)
    rotation = mbs_local.GetObjectOutput(body2, exu.OutputVariableType.RotationMatrix)
    tip = np.asarray(position) + np.asarray(rotation).reshape(3, 3).dot(
        np.array([L2 / 2.0, 0.0, 0.0]))
    
    if len(log_data) == 0 or abs(log_data[-1][0] - float(t)) > 1.0e-12:
        log_data.append([float(t), float(tip[1])])

def pre_step_user(mbs_local, t):
    record_tip(mbs_local, t)
    return True

mbs.SetPreStepUserFunction(pre_step_user)
mbs.Assemble()

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

record_tip(mbs, simulation_time)

log_array = np.asarray(log_data)
project_dir = Path(__file__).resolve().parents[1]
data_dir = project_dir / "data" / "planar_rigid_flex"
data_dir.mkdir(parents=True, exist_ok=True)
output_path = data_dir / "exudyn_endpoint_y_rigid.csv"

np.savetxt(
    output_path,
    log_array,
    delimiter=",",
    header="time,exudyn_rigid_tip_y",
    comments="",
)

print("Exudyn pure-rigid planar double pendulum completed successfully.")
print(f"Data saved to: {output_path}")