"""PyChrono model of a floating base, flexible beam, and hinged rigid rod."""

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

root_position_link = fea.ChLinkNodeFrame()
root_position_link.Initialize(beam_nodes[0], base)
system.Add(root_position_link)
root_slope_link = fea.ChLinkNodeSlopeFrame()
root_slope_link.Initialize(beam_nodes[0], base)
system.Add(root_slope_link)

# A negligible-inertia adapter converts the cable tip position/slope into a
# frame that can be connected to a revolute torque motor.
adapter = chrono.ChBody()
adapter.SetMass(1.0e-12)
adapter.SetInertiaXX(chrono.ChVector3d(1.0e-12, 1.0e-12, 1.0e-12))
adapter.SetPos(chrono.ChVector3d(beam_length, 0.0, 0.0))
system.Add(adapter)
tip_position_link = fea.ChLinkNodeFrame()
tip_position_link.Initialize(beam_nodes[-1], adapter)
system.Add(tip_position_link)
tip_slope_link = fea.ChLinkNodeSlopeFrame()
tip_slope_link.Initialize(beam_nodes[-1], adapter)
system.Add(tip_slope_link)

rod_length = 0.20
rod_mass = 1.5093
rod_com_x = 0.10
rod_Izz_com = 0.015155
rod = chrono.ChBody()
rod.SetMass(rod_mass)
rod.SetInertiaXX(chrono.ChVector3d(
    rod_Izz_com, rod_Izz_com, rod_Izz_com
))
rod.SetPos(chrono.ChVector3d(beam_length + rod_com_x, 0.0, 0.0))
add_box_visual(
    rod, (rod_length, 0.018, 0.012),
    (0.0, 0.0, 0.0), (0.25, 0.48, 0.72),
)
system.Add(rod)

hinge_torque = 0.02
hinge_motor = chrono.ChLinkMotorRotationTorque()
hinge_motor.Initialize(
    rod, adapter,
    chrono.ChFramed(chrono.ChVector3d(beam_length, 0.0, 0.0)),
)
hinge_motor.SetTorqueFunction(chrono.ChFunctionConst(hinge_torque))
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

simulation_time = 5.0
time_step = 1.0e-3
n_steps = int(round(simulation_time / time_step))
rows = []


def record(time_value):
    base_angle = yaw(base)
    beam_tip_angle = yaw(adapter)
    hinge_angle = yaw(rod) - beam_tip_angle
    root = beam_nodes[0].GetPos()
    beam_tip = beam_nodes[-1].GetPos()
    rod_tip = rod.TransformPointLocalToParent(
        chrono.ChVector3d(rod_length - rod_com_x, 0.0, 0.0)
    )
    relative_beam_tip = base.GetRot().GetConjugate().Rotate(beam_tip - root)
    rows.append([
        time_value,
        hinge_angle,
        base_angle,
        beam_tip.y,
        rod_tip.y,
        relative_beam_tip.y,
        root.x,
        root.y,
        beam_tip.x,
        rod_tip.x,
    ])


record(0.0)
start_time = time.perf_counter()
for step in range(1, n_steps + 1):
    system.DoStepDynamics(time_step)
    record(step * time_step)
elapsed = time.perf_counter() - start_time

data = np.asarray(rows)
data[:, 1:3] = np.unwrap(data[:, 1:3], axis=0)
project_dir = Path(__file__).resolve().parents[1]
data_dir = project_dir / "data" / "base_beam_hinged_rod"
data_dir.mkdir(parents=True, exist_ok=True)
output_path = data_dir / "chrono_base_beam_hinged_rod.csv"
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

print(f"Chrono base-beam-hinged-rod simulation completed in {elapsed:.3f} s.")
print(f"Data saved to: {output_path}")
