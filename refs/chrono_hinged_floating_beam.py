"""PyChrono reference for a floating rigid base hinged to a flexible beam."""

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

# Five-kilogram base; the hinge is at x=0 and its center of mass at x=-0.1.
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
    base, (base_length, base_height, base_depth),
    (0.0, 0.0, 0.0), (0.18, 0.22, 0.28),
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
    mesh, section, beam_elements,
    chrono.ChVector3d(0.0, 0.0, 0.0),
    chrono.ChVector3d(beam_length, 0.0, 0.0),
)
beam_nodes = list(beam_builder.GetLastBeamNodes())
system.Add(mesh)

# xyzD cable nodes have a slope coordinate rather than a rigid rotation.
# A negligible-inertia adapter transmits the revolute motor torque to both
# the root position and slope without changing the physical model mass.
adapter = chrono.ChBody()
adapter.SetMass(1.0e-12)
adapter.SetInertiaXX(chrono.ChVector3d(1.0e-12, 1.0e-12, 1.0e-12))
adapter.SetPos(chrono.ChVector3d(0.0, 0.0, 0.0))
system.Add(adapter)

root_position_link = fea.ChLinkNodeFrame()
root_position_link.Initialize(beam_nodes[0], adapter)
system.Add(root_position_link)
root_slope_link = fea.ChLinkNodeSlopeFrame()
root_slope_link.Initialize(beam_nodes[0], adapter)
system.Add(root_slope_link)

hinge_motor = chrono.ChLinkMotorRotationTorque()
hinge_motor.Initialize(
    adapter, base,
    chrono.ChFramed(chrono.ChVector3d(0.0, 0.0, 0.0)),
)
hinge_motor.SetTorqueFunction(chrono.ChFunctionConst(0.02))
system.Add(hinge_motor)

solver = chrono.ChSolverSparseQR()
system.SetSolver(solver)
timestepper = chrono.ChTimestepperHHT(system)
timestepper.SetAlpha((0.9 - 1.0) / (0.9 + 1.0))
timestepper.SetMaxIters(30)
timestepper.SetAbsTolerances(1.0e-9)
timestepper.SetModifiedNewton(False)
timestepper.SetStepControl(False)
system.SetTimestepper(timestepper)

simulation_time = 2.0
time_step = 1.0e-3
n_steps = int(round(simulation_time / time_step))
rows = []


def record(time_value):
    base_angle = yaw(base)
    beam_angle = yaw(adapter)
    root = beam_nodes[0].GetPos()
    tip = beam_nodes[-1].GetPos()
    relative_tip = adapter.GetRot().GetConjugate().Rotate(tip - root)
    rows.append([
        time_value, beam_angle - base_angle,
        relative_tip.x - beam_length, relative_tip.y,
        root.x, root.y, base_angle, tip.x, tip.y,
    ])


record(0.0)
start_time = time.perf_counter()
for step in range(1, n_steps + 1):
    system.DoStepDynamics(time_step)
    record(step * time_step)
elapsed = time.perf_counter() - start_time

data = np.asarray(rows)
data[:, 1] = np.unwrap(data[:, 1])
data[:, 6] = np.unwrap(data[:, 6])
project_dir = Path(__file__).resolve().parents[1]
data_dir = project_dir / "data" / "hinged_floating_beam"
data_dir.mkdir(parents=True, exist_ok=True)
output_path = data_dir / "chrono_hinged_floating_beam.csv"
np.savetxt(
    output_path, data, delimiter=",",
    header=(
        "time,hinge_angle,beam_tip_elastic_x,beam_tip_elastic_y,"
        "base_x,base_y,base_theta,beam_tip_global_x,beam_tip_global_y"
    ),
    comments="",
)

print(f"Chrono hinged floating-beam simulation completed in {elapsed:.3f} s.")
print(f"Data saved to: {output_path}")
