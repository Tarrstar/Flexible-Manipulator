"""PyChrono reference for the fixed-base flexible-beam manipulator."""

from pathlib import Path
import math
import time

import numpy as np
import pychrono as chrono
import pychrono.fea as fea


def add_box_visual(body, lengths, offset, color):
    shape = chrono.ChVisualShapeBox(*lengths)
    shape.SetColor(chrono.ChColor(*color))
    body.AddVisualShape(shape, chrono.ChFramed(chrono.ChVector3d(*offset)))


def yaw(body):
    q = body.GetRot()
    return 2.0 * math.atan2(q.e3, q.e0)


system = chrono.ChSystemSMC()
system.SetGravitationalAcceleration(chrono.ChVector3d(0.0, 0.0, 0.0))

beam_length = 1.0
beam_width = 0.032
beam_thickness = 0.0015
beam_density = 7.542e3
beam_youngs_modulus = 1.983e11
beam_elements = 3
beam_zeta = 0.0
beam_omega_n = 7.808

section = fea.ChBeamSectionCable()
section.SetArea(beam_width * beam_thickness)
section.SetInertia(beam_width * beam_thickness**3 / 12.0)
section.SetYoungModulus(beam_youngs_modulus)
section.SetDensity(beam_density)
section.SetRayleighDamping(2.0 * beam_zeta / beam_omega_n)

mesh = fea.ChMesh()
beam_builder = fea.ChBuilderCableANCF()
beam_builder.BuildBeam(
    mesh, section, beam_elements,
    chrono.ChVector3d(0.0, 0.0, 0.0),
    chrono.ChVector3d(beam_length, 0.0, 0.0),
)
beam_nodes = list(beam_builder.GetLastBeamNodes())
beam_nodes[0].SetFixed(True)
system.Add(mesh)

link_lengths = np.array([0.2, 0.1418, 0.1448, 2.0 * 0.0504])
link_masses = np.array([1.5093, 0.1029, 0.1139, 0.1009])
link_com_x = np.array([0.1, 0.0709, 0.0723, 0.0504])
link_Izz_com = np.array([0.015155, 0.002069, 0.002387, 0.001107])

links = []
x_left = beam_length
for length, mass, com_x, inertia in zip(
        link_lengths, link_masses, link_com_x, link_Izz_com):
    body = chrono.ChBody()
    body.SetMass(float(mass))
    body.SetInertiaXX(chrono.ChVector3d(float(inertia), float(inertia), float(inertia)))
    body.SetPos(chrono.ChVector3d(float(x_left + com_x), 0.0, 0.0))
    add_box_visual(
        body, (float(length), 0.018, 0.012),
        (float(0.5 * length - com_x), 0.0, 0.0), (0.25, 0.48, 0.72),
    )
    system.Add(body)
    links.append(body)
    x_left += length

tip_position_link = fea.ChLinkNodeFrame()
tip_position_link.Initialize(beam_nodes[-1], links[0])
system.Add(tip_position_link)
tip_slope_link = fea.ChLinkNodeSlopeFrame()
tip_slope_link.Initialize(beam_nodes[-1], links[0])
system.Add(tip_slope_link)

x_joint = beam_length + link_lengths[0]
for i in range(3):
    motor = chrono.ChLinkMotorRotationTorque()
    motor.Initialize(
        links[i + 1], links[i],
        chrono.ChFramed(chrono.ChVector3d(float(x_joint), 0.0, 0.0)),
    )
    motor.SetTorqueFunction(chrono.ChFunctionConst(0.02))
    system.Add(motor)
    x_joint += link_lengths[i + 1]

solver = chrono.ChSolverSparseQR()
system.SetSolver(solver)
timestepper = chrono.ChTimestepperHHT(system)
timestepper.SetAlpha((0.9 - 1.0) / (0.9 + 1.0))
timestepper.SetMaxIters(30)
timestepper.SetAbsTolerances(1.0e-9)
timestepper.SetModifiedNewton(False)
timestepper.SetStepControl(False)
system.SetTimestepper(timestepper)

simulation_time = 5.0
time_step = 1.0e-3
n_steps = int(round(simulation_time / time_step))
rows = []


def record(time_value):
    link_angles = np.array([yaw(body) for body in links])
    joint_angles = np.diff(link_angles)
    root = beam_nodes[0].GetPos()
    tip = beam_nodes[-1].GetPos()
    rows.append([
        time_value, *joint_angles,
        tip.x - root.x - beam_length, tip.y - root.y,
        tip.x, tip.y,
    ])


record(0.0)
start_time = time.perf_counter()
for step in range(1, n_steps + 1):
    system.DoStepDynamics(time_step)
    record(step * time_step)
elapsed = time.perf_counter() - start_time

data = np.asarray(rows)
data[:, 1:4] = np.unwrap(data[:, 1:4], axis=0)
project_dir = Path(__file__).resolve().parents[1]
data_dir = project_dir / "data" / "fixed_FBM"
data_dir.mkdir(parents=True, exist_ok=True)
output_path = data_dir / "chrono_fixed_FBM.csv"
np.savetxt(
    output_path, data, delimiter=",",
    header=(
        "time,joint_1_angle,joint_2_angle,joint_3_angle,"
        "beam_tip_elastic_x,beam_tip_elastic_y,"
        "beam_tip_global_x,beam_tip_global_y"
    ),
    comments="",
)

print(f"Chrono fixed-FBM simulation completed in {elapsed:.3f} s.")
print(f"Data saved to: {output_path}")
