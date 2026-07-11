import exudyn as exu
from exudyn.utilities import *
import exudyn.graphics as graphics
import numpy as np
import matplotlib.pyplot as plt
import os
from pathlib import Path

SC = exu.SystemContainer()
mbs = SC.AddSystem()

# ============================================================
# 1. 原梁参数：用于构造刚性等效杆
# ============================================================
Lbeam = 2.0
W_beam = 0.1
T_beam = 0.05
rho = 2700.0
E = 7.0e7  # 刚性杆仿真中不使用，保留仅用于参数对应

A = W_beam * T_beam
rhoA = rho * A

# 原 Cable2D 梁的总质量
m_beam = rhoA * Lbeam

# 采用 Cable2D 线质量分布对应的质心转动惯量
I_beam_COM = m_beam * Lbeam**2 / 12.0

# 若你想采用实体矩形杆的 z 轴惯量，可改用：
# I_beam_COM = m_beam * (Lbeam**2 + T_beam**2) / 12.0
# 但这不再与 Cable2D 的纯线质量模型严格一致。

print('Rigid equivalent beam parameters:')
print(f'mass = {m_beam:.8f} kg')
print(f'I_COM = {I_beam_COM:.8f} kg*m^2')

# ============================================================
# 2. 刚性杆建模：完全悬空，无约束，无重力
# ============================================================
c_beam = 0.5 * Lbeam

# 刚体节点位于质心；初始时杆从 x=0 延伸到 x=Lbeam
nBeam = mbs.AddNode(
    NodeRigidBody2D(
        referenceCoordinates=[c_beam, 0, 0]
    )
)

bBeam = mbs.AddObject(
    ObjectRigidBody2D(
        nodeNumber=nBeam,
        physicsMass=m_beam,
        physicsInertia=I_beam_COM,
        visualization=VObjectRigidBody2D(
            graphicsData=[
                graphics.Brick(
                    centerPoint=[0, 0, 0],
                    size=[Lbeam, T_beam, W_beam],
                    color=graphics.color.blue
                )
            ]
        )
    )
)

# 左端和右端 marker
mBase = mbs.AddMarker(
    MarkerBodyRigid(
        bodyNumber=bBeam,
        localPosition=[-c_beam, 0, 0]
    )
)

mTip = mbs.AddMarker(
    MarkerBodyRigid(
        bodyNumber=bBeam,
        localPosition=[c_beam, 0, 0]
    )
)

# ============================================================
# 3. 外载荷：与原 ANCF 代码一致
# ============================================================
# 你的原代码实际使用 Fy = 1.0 N
Fy = 0.0
Mz = 1.0

# 力作用在左端，惯性系方向
mbs.AddLoad(
    LoadForceVector(
        markerNumber=mBase,
        loadVector=[0.0, 1.0, 1.0],
        bodyFixed=False
    )
)

# 力矩作用在左端；当前为 0
mbs.AddLoad(
    LoadTorqueVector(
        markerNumber=mBase,
        loadVector=[0.0, 0.0, Mz],
        bodyFixed=False
    )
)

# ============================================================
# 4. 数据记录
# ============================================================
log_data = []

def RecordState(mbs, t):
    # 刚体节点质心的绝对位置
    pos_com = mbs.GetNodeOutput(
        nBeam,
        exu.OutputVariableType.Position
    )

    # NodeRigidBody2D 的第三坐标为刚体绝对转角
    coords = mbs.GetNodeOutput(
        nBeam,
        exu.OutputVariableType.Coordinates
    )

    x_com = float(pos_com[0])
    y_com = float(pos_com[1])
    theta = float(coords[2])

    # 左端点位置
    base_x = x_com - c_beam * np.cos(theta)
    base_y = y_com - c_beam * np.sin(theta)

    # 右端点位置
    tip_x = x_com + c_beam * np.cos(theta)
    tip_y = y_com + c_beam * np.sin(theta)

    log_data.append([
        t,
        base_x,
        base_y,
        theta,
        tip_x,
        tip_y
    ])

def PreStepUser(mbs, t):
    RecordState(mbs, t)
    return True

mbs.SetPreStepUserFunction(PreStepUser)
mbs.Assemble()

# ============================================================
# 5. 仿真设置：保持与原 ANCF 代码一致
# ============================================================
simulationTime = 2.0
timeStep = 2e-3

sim = exu.SimulationSettings()
sim.timeIntegration.numberOfSteps = int(round(simulationTime / timeStep))
sim.timeIntegration.endTime = simulationTime

sim.linearSolverType = exu.LinearSolverType.EigenDense
sim.timeIntegration.newton.useModifiedNewton = True
sim.timeIntegration.newton.relativeTolerance = 1e-6
sim.timeIntegration.verboseMode = 1

# 刚性杆模型不应忽略奇异约束；这里没有约束，因此应为 False
sim.linearSolverSettings.ignoreSingularJacobian = False

exu.StartRenderer()
mbs.SolveDynamic(
    sim,
    solverType=exu.DynamicSolverType.GeneralizedAlpha
)
exu.StopRenderer()

if len(log_data) == 0 or abs(log_data[-1][0] - simulationTime) > 1e-12:
    RecordState(mbs, simulationTime)

# ============================================================
# 6. 导出与绘图
# ============================================================
log_array = np.array(log_data)
project_dir = Path(__file__).resolve().parents[1]
data_dir = project_dir / "data" / "floating_beam"
data_dir.mkdir(parents=True, exist_ok=True)
save_dir = data_dir
os.makedirs(save_dir, exist_ok=True)

save_path = os.path.join(save_dir, "floating_rigid_beam_exudyn.csv")

np.savetxt(
    save_path,
    log_array,
    delimiter=",",
    header="time,base_x,base_y,base_theta,tip_x,tip_y",
    comments=""
)

# 单独导出末端坐标，便于 MATLAB 对比
tip_path = os.path.join(save_dir, "floating_rigid_tip_exudyn.csv")

np.savetxt(
    tip_path,
    log_array[:, [0, 4, 5]],
    delimiter=",",
    header="time,tip_x,tip_y",
    comments=""
)

print(f"Rigid-beam simulation completed.")
print(f"Full data saved to:\n{save_path}")
print(f"Tip data saved to:\n{tip_path}")

plt.figure(figsize=(10, 8))

plt.subplot(3, 1, 1)
plt.plot(log_array[:, 0], log_array[:, 1], label="Base X")
plt.plot(log_array[:, 0], log_array[:, 2], label="Base Y")
plt.ylabel("Base position (m)")
plt.legend()
plt.grid(True)

plt.subplot(3, 1, 2)
plt.plot(log_array[:, 0], log_array[:, 3], "g", label="Rigid rod angle")
plt.ylabel("Angle (rad)")
plt.legend()
plt.grid(True)

plt.subplot(3, 1, 3)
plt.plot(log_array[:, 0], log_array[:, 4], label="Tip X")
plt.plot(log_array[:, 0], log_array[:, 5], label="Tip Y")
plt.xlabel("Time (s)")
plt.ylabel("Tip position (m)")
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.show()