# 空间七关节、七杆全柔性机械臂：MATLAB / Exudyn 末端位置对比

数据需通过下面的步骤重新生成；历史结果与收敛汇总脚本已删除。

## 模型和驱动

- 固定基座，7 根连杆均为柔性梁，没有额外的第八根末端梁。
- 长度为 `[0.28,0.30,0.26,0.22,0.18,0.14,0.10] m`；关节轴为
  `Rz,Ry,Ry,Rx,Ry,Rx,Ry`。沿用原七关节演示的长度和关节轴。
- 示例梁截面为 `0.025 m × 0.025 m`，密度 `2700 kg/m³`，杨氏模量
  `7e8 Pa`，泊松比 `0.33`。这些是本算例的示例参数，不是实物辨识数据。
- 零重力、零物理阻尼；初始关节角、弹性坐标及其速度全部为零。
- MATLAB 使用当前 `ID_flex`、`FlexBody`、`FEM_Beam_CB_Reduction` 和
  `generalized_alpha_id_integrator`。每杆 4 个 Euler-Bernoulli 单元，保留
  6 个末端 Craig-Bampton 约束坐标，不保留内部固定界面振型，合计
  7 个关节坐标＋42 个弹性坐标。6 个约束坐标不等于 6 个固有振型。
- Exudyn 1.10.0 独立建立 `ObjectBeamGeometricallyExact` 梁和转动副；最终
  每杆 128 个单元。剪切刚度取 `100 GA`，逼近 MATLAB 的 Euler-Bernoulli
  假设。细网格用 Newton 坐标增量判据，绝对/相对容差均为 `1e-10`，
  避免大刚度下力残差的浮点噪声；32 单元下已复核更换判据的轨迹影响。
- 全刚性基准保持每杆质量、质心和惯量一致；MATLAB 在未变形状态逐杆
  核对刚柔空间惯量相等。
- 输入是同一组预先给定的关节力矩，不是位移约束或反馈控制：根据等惯量
  刚性臂的五次多项式参考运动计算前馈力矩，1 s 后力矩为零，仿真至 2 s。
  两套软件独立积分，没有向 Exudyn 输入 MATLAB 的运动轨迹或加速度。

由于是开环力矩驱动，柔性也会改变实际关节运动。因此，柔刚末端位置差
包含关节运动变化与各杆弹性变形，不能当成某根杆的局部挠度。

## 复现

MATLAB（从项目目录启动）：

```matlab
addpath('D:/project/Flexible_manipulator/src/tests/dynamics');
out = 'D:/project/Flexible_manipulator/data/space_7dof_all_flexible/reproduction';
test_7dof_all_flexible_exudyn('prepare',2.5e-4,2,out);
test_7dof_all_flexible_exudyn('simulate',2.5e-4,2,out);
```

PowerShell（Exudyn 环境和画图）：

```powershell
& 'D:/tools/anaconda/envs/Exudyn/python.exe' refs/exudyn_7dof_all_flexible.py --directory data/space_7dof_all_flexible/reproduction --case flexible --elements 128 --step 5e-4
& 'D:/tools/anaconda/envs/Exudyn/python.exe' refs/exudyn_7dof_all_flexible.py --directory data/space_7dof_all_flexible/reproduction --case rigid --step 1.25e-4
& 'D:/tools/anaconda/envs/Exudyn/python.exe' src/tests/dynamics/plot_7dof_all_flexible_exudyn.py --directory data/space_7dof_all_flexible/reproduction
```

MATLAB 的刚性基准自动使用不大于 1.25e-4 s 的步长。

MATLAB 小变形 FFRF 与 Exudyn 几何精确梁并非同一离散模型，二者差值包含模型近似及离散误差。
