"""PyChrono reference for a floating base rigidly connected to a flexible beam."""

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
    quaternion = body.GetRot()
    return 2.0 * math.atan2(quaternion.e3, quaternion.e0)


system = chrono.ChSystemSMC()
system.SetGravitationalAcceleration(chrono.ChVector3d(0.0, 0.0, 0.0))

# Chrono locates the base body frame at its center of mass.
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
    mesh,
    section,
    beam_elements,
    chrono.ChVector3d(0.0, 0.0, 0.0),
    chrono.ChVector3d(beam_length, 0.0, 0.0),
)
beam_nodes = list(beam_builder.GetLastBeamNodes())
system.Add(mesh)

# Fix both the root position and slope to the base frame.
root_position_link = fea.ChLinkNodeFrame()
root_position_link.Initialize(beam_nodes[0], base)
system.Add(root_position_link)
root_slope_link = fea.ChLinkNodeSlopeFrame()
root_slope_link.Initialize(beam_nodes[0], base)
system.Add(root_slope_link)

solver = chrono.ChSolverSparseQR()
system.SetSolver(solver)
timestepper = chrono.ChTimestepperHHT(system)
timestepper.SetAlpha((0.9 - 1.0) / (0.9 + 1.0))
timestepper.SetMaxIters(30)
timestepper.SetAbsTolerances(1.0e-9)
timestepper.SetModifiedNewton(False)
timestepper.SetStepControl(False)
system.SetTimestepper(timestepper)

force_amplitude = 0.1
force_frequency = 1.0
simulation_time = 10.0
time_step = 1.0e-3
n_steps = int(round(simulation_time / time_step))
rows = []


def force_y(time_value):
    return force_amplitude * math.sin(2.0 * math.pi * force_frequency * time_value)


def record(time_value):
    root = beam_nodes[0].GetPos()
    tip = beam_nodes[-1].GetPos()
    relative_tip = base.GetRot().GetConjugate().Rotate(tip - root)
    rows.append([
        time_value,
        yaw(base),
        tip.y,
        relative_tip.y,
        root.x,
        root.y,
        tip.x,
        force_y(time_value),
    ])


record(0.0)
start_time = time.perf_counter()
for step in range(1, n_steps + 1):
    # SetForce is persistent for the following integration step.
    load_time = (step - 1) * time_step
    beam_nodes[-1].SetForce(chrono.ChVector3d(0.0, force_y(load_time), 0.0))
    system.DoStepDynamics(time_step)
    record(step * time_step)
elapsed = time.perf_counter() - start_time

data = np.asarray(rows)
data[:, 1] = np.unwrap(data[:, 1])
project_dir = Path(__file__).resolve().parents[1]
data_dir = project_dir / "data" / "base_beam"
data_dir.mkdir(parents=True, exist_ok=True)
output_path = data_dir / "chrono_base_beam.csv"
np.savetxt(
    output_path,
    data,
    delimiter=",",
    header=(
        "time,base_theta,beam_tip_global_y,beam_tip_elastic_y,"
        "base_x,base_y,beam_tip_global_x,force_y"
    ),
    comments="",
)

print(f"Chrono base-beam simulation completed in {elapsed:.3f} s.")
print(f"Data saved to: {output_path}")
