"""PyChrono reference for the floating flexible-base manipulator."""

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

# Five-kilogram floating base. Its right-face center is the beam attachment.
base_mass = 5.0
base_length = 0.20
base_height = 0.12
base_depth = 0.08
base_com_x = -0.5 * base_length
base_inertia = base_mass / 12.0 * np.array([
    base_height**2 + base_depth**2,
    base_length**2 + base_depth**2,
    base_length**2 + base_height**2,
])

base = chrono.ChBody()
base.SetMass(base_mass)
base.SetInertiaXX(chrono.ChVector3d(*base_inertia))
base.SetPos(chrono.ChVector3d(base_com_x, 0.0, 0.0))
add_box_visual(
    base,
    (base_length, base_height, base_depth),
    (0.0, 0.0, 0.0),
    (0.18, 0.22, 0.28),
)
system.Add(base)

# ANCF cable with the same section, mesh, stiffness, mass, and damping data.
beam_length = 1.0
beam_width = 0.032
beam_thickness = 0.0015
beam_density = 7.542e3
beam_youngs_modulus = 1.983e11
beam_elements = 5
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
    mesh,
    section,
    beam_elements,
    chrono.ChVector3d(0.0, 0.0, 0.0),
    chrono.ChVector3d(beam_length, 0.0, 0.0),
)
beam_nodes = list(beam_builder.GetLastBeamNodes())
system.Add(mesh)

# Four rigid arm links. Chrono body frames are located at each center of mass.
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
        body,
        (float(length), 0.018, 0.012),
        (float(0.5 * length - com_x), 0.0, 0.0),
        (0.25, 0.48, 0.72),
    )
    system.Add(body)
    links.append(body)
    x_left += length

# Rigid attachment of beam root to base: position plus slope direction.
root_position_link = fea.ChLinkNodeFrame()
root_position_link.Initialize(beam_nodes[0], base)
system.Add(root_position_link)
root_slope_link = fea.ChLinkNodeSlopeFrame()
root_slope_link.Initialize(beam_nodes[0], base)
system.Add(root_slope_link)

# Rigid attachment of beam tip to the first arm body.
tip_position_link = fea.ChLinkNodeFrame()
tip_position_link.Initialize(beam_nodes[-1], links[0])
system.Add(tip_position_link)
tip_slope_link = fea.ChLinkNodeSlopeFrame()
tip_slope_link.Initialize(beam_nodes[-1], links[0])
system.Add(tip_slope_link)

# Three revolute torque motors; positive torque acts on each child link.
joint_torque_amplitude = 0.02
joint_torque_frequency = 1.0
motors = []
x_joint = beam_length + link_lengths[0]
for i in range(3):
    motor = chrono.ChLinkMotorRotationTorque()
    motor.Initialize(
        links[i + 1],
        links[i],
        chrono.ChFramed(chrono.ChVector3d(float(x_joint), 0.0, 0.0)),
    )
    motor.SetTorqueFunction(chrono.ChFunctionSine(
        joint_torque_amplitude, joint_torque_frequency
    ))
    system.Add(motor)
    motors.append(motor)
    x_joint += link_lengths[i + 1]

# Sparse direct solve and generalized-alpha-equivalent HHT integration.
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
    base_angle = yaw(base)
    link_angles = np.array([yaw(body) for body in links])
    joint_angles = np.diff(link_angles)

    root = beam_nodes[0].GetPos()
    tip = beam_nodes[-1].GetPos()
    relative_tip = base.GetRot().GetConjugate().Rotate(tip - root)

    rows.append([
        time_value,
        *joint_angles,
        relative_tip.y,
        root.x,
        root.y,
        base_angle,
        tip.x,
        tip.y,
    ])


record(0.0)
start_time = time.perf_counter()
for step in range(1, n_steps + 1):
    system.DoStepDynamics(time_step)
    record(step * time_step)
elapsed = time.perf_counter() - start_time

data = np.asarray(rows)
data[:, 1:4] = np.unwrap(data[:, 1:4], axis=0)
data[:, 7] = np.unwrap(data[:, 7])

project_dir = Path(__file__).resolve().parents[1]
data_dir = project_dir / "data" / "flexible_base_manipulator"
data_dir.mkdir(parents=True, exist_ok=True)
output_path = data_dir / "chrono_floating_base.csv"
np.savetxt(
    output_path,
    data,
    delimiter=",",
    header=(
        "time,joint_1_angle,joint_2_angle,joint_3_angle,"
        "beam_tip_elastic_y,base_x,base_y,base_theta,"
        "beam_tip_global_x,beam_tip_global_y"
    ),
    comments="",
)

print(f"Chrono simulation completed in {elapsed:.3f} s.")
print(f"Data saved to: {output_path}")
