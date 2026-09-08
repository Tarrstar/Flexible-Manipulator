import exudyn as exu
from exudyn.utilities import *
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


SC = exu.SystemContainer()
mbs = SC.AddSystem()

# ================== 1. Beam parameters ==================

Lbeam = 1.0
W_beam = 0.032
T_beam = 0.0015
rho = 7543
E = 1.983e11

A = W_beam * T_beam
Izz = W_beam * T_beam**3 / 12.0
rhoA = rho * A
EA = E * A
EI = E * Izz

# 3D ANCF cable:
# Node coordinates are [x, y, z, slope_x, slope_y, slope_z].
# The cable uses one scalar bending stiffness for both transverse
# directions, so this is a 3D reference for force directions rather than
# a fully anisotropic rectangular-section beam.
beam = Cable(
    physicsMassPerLength=rhoA,
    physicsAxialStiffness=EA,
    physicsAxialDamping=0,
    physicsBendingStiffness=EI,
    physicsBendingDamping=0,
    visualization=VCable(radius=0.5 * T_beam)
)

n_elements = 5
beamData = GenerateStraightBeam(
    mbs,
    positionOfNode0=[0, 0, 0],
    positionOfNode1=[Lbeam, 0, 0],
    numberOfElements=n_elements,
    beamTemplate=beam,
    gravity=[0, 0, 0]
)

node_base = beamData[0][0]
node_tip = beamData[0][-1]

mBase = mbs.AddMarker(MarkerNodePosition(nodeNumber=node_base))

# Same load convention as rhs_flexbody:
# [Fx, Fy, Fz] = [0, 1, 1], no externally applied torque.
force_global = [0.0, 1.0, 2.0]
mbs.AddLoad(LoadForceVector(markerNumber=mBase, loadVector=force_global))

# ================== 2. Data logging ==================

log_data = []


def RecordState(mbs, t):
    pos0 = mbs.GetNodeOutput(node_base, exu.OutputVariableType.Position)
    coords0 = mbs.GetNodeOutput(node_base, exu.OutputVariableType.Coordinates)

    slope_x = float(coords0[3])
    slope_y = float(coords0[4])
    base_theta = np.arctan2(slope_y, slope_x)

    pos_tip = mbs.GetNodeOutput(node_tip, exu.OutputVariableType.Position)

    log_data.append([
        t,
        float(pos0[0]),
        float(pos0[1]),
        float(pos0[2]),
        base_theta,
        float(pos_tip[0]),
        float(pos_tip[1]),
        float(pos_tip[2]),
    ])


def PreStepUser(mbs, t):
    RecordState(mbs, t)
    return True


mbs.SetPreStepUserFunction(PreStepUser)
mbs.Assemble()

# ================== 3. Simulation settings ==================

simulationTime = 2.0
timeStep = 2e-3

sim = exu.SimulationSettings()
sim.timeIntegration.numberOfSteps = int(round(simulationTime / timeStep))
sim.timeIntegration.endTime = simulationTime
sim.linearSolverType = exu.LinearSolverType.EigenDense
sim.timeIntegration.newton.useModifiedNewton = True
sim.timeIntegration.newton.relativeTolerance = 1e-6
sim.timeIntegration.verboseMode = 1

mbs.SolveDynamic(sim, solverType=exu.DynamicSolverType.GeneralizedAlpha)

if len(log_data) == 0 or abs(log_data[-1][0] - simulationTime) > 1e-12:
    RecordState(mbs, simulationTime)

# ================== 4. Save and plot ==================

log_array = np.array(log_data)
project_dir = Path(__file__).resolve().parents[1]
data_dir = project_dir / "data" / "floating_beam"
data_dir.mkdir(parents=True, exist_ok=True)

output_path = data_dir / "floating_beam_3d_exudyn.csv"
figure_path = data_dir / "floating_beam_3d_exudyn.png"
np.savetxt(
    output_path,
    log_array,
    delimiter=",",
    header="time,base_x,base_y,base_z,base_theta,tip_x,tip_y,tip_z",
    comments=""
)

plt.figure(figsize=(10, 8))

plt.subplot(2, 1, 1)
plt.plot(log_array[:, 0], log_array[:, 1], label="Base X")
plt.plot(log_array[:, 0], log_array[:, 2], label="Base Y")
plt.plot(log_array[:, 0], log_array[:, 3], label="Base Z")
plt.ylabel("Base position")
plt.legend()
plt.grid(True)

plt.subplot(2, 1, 2)
plt.plot(log_array[:, 0], log_array[:, 5], label="Tip X")
plt.plot(log_array[:, 0], log_array[:, 6], label="Tip Y")
plt.plot(log_array[:, 0], log_array[:, 7], label="Tip Z")
plt.xlabel("Time / s")
plt.ylabel("Tip position")
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.savefig(figure_path, dpi=160)
plt.close()

print("3D Exudyn beam simulation completed.")
print(f"Data saved to:\n{output_path}")
print(f"Figure saved to:\n{figure_path}")
