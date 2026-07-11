import exudyn as exu
from exudyn.utilities import *
import exudyn.graphics as graphics
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

SC = exu.SystemContainer()
mbs = SC.AddSystem()

# ================== 1. 梁参数与物理属性 ==================
Lbeam = 2.0
W_beam = 0.1
T_beam = 0.05
rho = 2700.0
E = 7.0e7

A = W_beam * T_beam
I = W_beam * (T_beam**3) / 12.0
rhoA = rho * A
EA = E * A
EI = E * I

# 创建 ANCF 梁 (无阻尼, 无重力)
beam = Cable2D(
    physicsMassPerLength=rhoA,
    physicsAxialStiffness=EA,
    physicsAxialDamping=0,
    physicsBendingStiffness=EI,
    physicsBendingDamping=0,
    visualization=VCable2D(drawHeight=T_beam)
)

# 生成直线梁，去掉了 fixedConstraintsNode0，使其完全悬空
n_elements = 5
beamData = GenerateStraightLineANCFCable2D(
    mbs, 
    positionOfNode0=[0, 0, 0], 
    positionOfNode1=[Lbeam, 0, 0], 
    numberOfElements=n_elements, 
    cableTemplate=beam,
    massProportionalLoad=[0, 0, 0] # 无重力
)

# ================== 2. 施加外力 (0.1N, 0.1Nm) ==================
# 在梁的最左端节点 (Node 0) 施加载荷
node_base = beamData[0][0]
node_tip = beamData[0][-1]

mBase = mbs.AddMarker(MarkerNodeRigid(nodeNumber=node_base))
mTip = mbs.AddMarker(MarkerNodeRigid(nodeNumber=node_tip))
# Fx = 0.1 N, Fy = 0.1 N, Torque = 0.1 Nm
mbs.AddLoad(LoadForceVector(markerNumber=mBase, loadVector=[0.0, 1.0, 1.0]))
mbs.AddLoad(LoadTorqueVector(markerNumber=mBase, loadVector=[0.0, 0.0, 0.0]))

# ================== 3. 提取数据日志 ==================
log_data = []

def PreStepUser(mbs, t):
    # --- 修改这里：使用 Position 而不是 Coordinates ---
    # Position 会返回节点在惯性系下的绝对 [X, Y, Z] 坐标
    pos0 = mbs.GetNodeOutput(node_base, exu.OutputVariableType.Position)
    base_x, base_y = pos0[0], pos0[1]
    
    # 对于角度，Coordinates 里的 [2, 3] 位是 dx/ds, dy/ds，这是正确的
    coords0 = mbs.GetNodeOutput(node_base, exu.OutputVariableType.Coordinates)
    dx, dy = coords0[2], coords0[3] 
    base_theta = np.arctan2(dy, dx)
    
    # --- 修改这里：获取末端的绝对位置 ---
    posTip = mbs.GetNodeOutput(node_tip, exu.OutputVariableType.Position)
    tip_x, tip_y = posTip[0], posTip[1]
    
    log_data.append([t, base_x, base_y, base_theta, tip_x, tip_y])
    return True

mbs.SetPreStepUserFunction(PreStepUser)
mbs.Assemble()

# ================== 4. 仿真设置 ==================
simulationTime = 2.0
timeStep = 2e-3

sim = exu.SimulationSettings()
sim.timeIntegration.numberOfSteps = int(simulationTime / timeStep)
sim.timeIntegration.endTime = simulationTime
sim.linearSolverType = exu.LinearSolverType.EigenDense
sim.timeIntegration.newton.useModifiedNewton = True
sim.timeIntegration.newton.relativeTolerance = 1e-6
sim.timeIntegration.verboseMode = 1

# 运行仿真
exu.StartRenderer()
mbs.SolveDynamic(sim, solverType=exu.DynamicSolverType.GeneralizedAlpha)
exu.StopRenderer()

# ================== 5. 保存并绘图 ==================
log_array = np.array(log_data)
# Export CSV for MATLAB comparison.
project_dir = Path(__file__).resolve().parents[1]
data_dir = project_dir / "data" / "floating_beam"
data_dir.mkdir(parents=True, exist_ok=True)
output_path = data_dir / "floating_beam_exudyn.csv"
np.savetxt(output_path, log_array, delimiter=",",
           header="time,base_x,base_y,base_theta,tip_x,tip_y", comments="")

plt.figure(figsize=(10, 8))
plt.subplot(3, 1, 1)
plt.plot(log_array[:, 0], log_array[:, 1], label='Base X')
plt.plot(log_array[:, 0], log_array[:, 2], label='Base Y')
plt.legend(); plt.grid(True)

plt.subplot(3, 1, 2)
plt.plot(log_array[:, 0], log_array[:, 3], 'g', label='Base Theta (rad)')
plt.legend(); plt.grid(True)

plt.subplot(3, 1, 3)
plt.plot(log_array[:, 0], log_array[:, 4], label='Tip X')
plt.plot(log_array[:, 0], log_array[:, 5], label='Tip Y')
plt.legend(); plt.grid(True)
plt.tight_layout()
plt.show()
print("仿真完成，数据已保存至 floating_beam_exudyn.csv")