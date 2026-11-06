"""Independent Exudyn spatial seven-flexible-link model and rigid baseline.

Loads/physical data come from MATLAB's exported model.json and torques.csv;
neither MATLAB accelerations nor trajectories are used by this solver.
"""
import argparse
import json
from pathlib import Path
import time

import exudyn as exu
from exudyn.itemInterface import (ObjectGround, MarkerBodyRigid, MarkerNodeRigid,
    NodeRigidBodyEP, ObjectBeamGeometricallyExact, GenericJoint, LoadTorqueVector)
from exudyn.utilities import RigidBodyInertia, RotationMatrixX, RotationMatrixY
import numpy as np


def simulate(cfg,torque,rigid,elements,step,end_time,out):
    sc=exu.SystemContainer(); mbs=sc.AddSystem()
    ground=mbs.AddObject(ObjectGround())
    parent_tip=mbs.AddMarker(MarkerBodyRigid(bodyNumber=ground,localPosition=[0,0,0]))
    axes=np.asarray(cfg['jointAxes'],dtype=float)
    lengths=np.asarray(cfg['lengths'],dtype=float).ravel()
    width=cfg['width']; thickness=cfg['thickness']; density=cfg['density']
    E=cfg['youngsModulus']; nu=cfg['poissonRatio']; G=E/(2*(1+nu))
    area=width*thickness
    Iyy=width*thickness**3/12; Izz=thickness*width**3/12
    aa=max(width,thickness); bb=min(width,thickness)
    J=aa*bb**3*(1/3-.21*(bb/aa)*(1-bb**4/(12*aa**4)))
    section=exu.BeamSection()
    # Large shear stiffness approaches the Euler-Bernoulli assumption used
    # in MATLAB. Square sections avoid an axis-label ambiguity in old cases.
    section.stiffnessMatrix=np.diag([E*area,100*G*area,100*G*area,G*J,E*Iyy,E*Izz])
    section.massPerLength=density*area
    section.inertia=density*np.diag([Iyy+Izz,Iyy,Izz])
    x_left=0.; last_tip_node=None; last_body=None
    keep_callbacks=[]
    def load_function(i,sign):
        def law(system,t,static):
            return (sign*axes[i]*np.interp(t,torque[:,0],torque[:,i+1])).tolist()
        keep_callbacks.append(law)
        return law
    for i,L in enumerate(lengths):
        if rigid:
            mass=density*area*L
            inertia=mass/12*np.diag([width**2+thickness**2,L**2+thickness**2,L**2+width**2])
            body=mbs.CreateRigidBody(inertia=RigidBodyInertia(mass,inertia,[L/2,0,0],True),
                                      referencePosition=[x_left,0,0],gravity=[0,0,0])
            root=mbs.AddMarker(MarkerBodyRigid(bodyNumber=body,localPosition=[0,0,0]))
            tip=mbs.AddMarker(MarkerBodyRigid(bodyNumber=body,localPosition=[L,0,0]))
            last_body=body
        else:
            nodes=[]
            for j in range(elements+1):
                node=mbs.AddNode(NodeRigidBodyEP(referenceCoordinates=
                    [x_left+L*j/elements,0,0,1,0,0,0]))
                nodes.append(node)
                if j:
                    mbs.AddObject(ObjectBeamGeometricallyExact(nodeNumbers=[nodes[j-1],node],
                        physicsLength=L/elements,sectionData=section))
            root=mbs.AddMarker(MarkerNodeRigid(nodeNumber=nodes[0]))
            tip=mbs.AddMarker(MarkerNodeRigid(nodeNumber=nodes[-1]))
            last_tip_node=nodes[-1]
        axis=axes[i]
        if axis[2]: rotation=np.eye(3)
        elif axis[1]: rotation=RotationMatrixX(-np.pi/2)
        else: rotation=RotationMatrixY(np.pi/2)
        mbs.AddObject(GenericJoint(markerNumbers=[parent_tip,root],
            constrainedAxes=[1,1,1,1,1,0],rotationMarker0=rotation,rotationMarker1=rotation))
        mbs.AddLoad(LoadTorqueVector(markerNumber=root,bodyFixed=True,
                    loadVectorUserFunction=load_function(i,1.)))
        if i:
            mbs.AddLoad(LoadTorqueVector(markerNumber=parent_tip,bodyFixed=True,
                        loadVectorUserFunction=load_function(i,-1.)))
        parent_tip=tip; x_left+=L
    rows=[]
    def record(system,t):
        if rows and abs(rows[-1][0]-t)<1e-12: return
        if rigid:
            p=system.GetObjectOutputBody(last_body,exu.OutputVariableType.Position,
                                          localPosition=[lengths[-1],0,0])
        else:
            p=system.GetNodeOutput(last_tip_node,exu.OutputVariableType.Position)
        rows.append([float(t),*np.asarray(p).ravel().tolist()])
    def post(system,t):
        record(system,t); return True
    mbs.SetPostStepUserFunction(post)
    mbs.Assemble(); record(mbs,0.)
    settings=exu.SimulationSettings()
    settings.timeIntegration.endTime=end_time
    settings.timeIntegration.numberOfSteps=round(end_time/step)
    settings.timeIntegration.generalizedAlpha.spectralRadius=cfg['spectralRadius']
    # The stiff, fine mesh can hit a force-residual roundoff floor. Check
    # Newton coordinate increments instead, at a fixed tight tolerance.
    settings.timeIntegration.newton.newtonResidualMode=1
    settings.timeIntegration.newton.relativeTolerance=1e-10
    settings.timeIntegration.newton.absoluteTolerance=1e-10
    settings.timeIntegration.newton.useModifiedNewton=False
    settings.timeIntegration.verboseMode=1
    settings.linearSolverType=exu.LinearSolverType.EigenSparse
    settings.solutionSettings.writeSolutionToFile=False
    start=time.perf_counter()
    mbs.SolveDynamic(settings,solverType=exu.DynamicSolverType.GeneralizedAlpha)
    elapsed=time.perf_counter()-start
    record(mbs,end_time)
    values=np.asarray(rows)
    assert np.isfinite(values).all() and abs(values[-1,0]-end_time)<1e-10
    np.savetxt(out,values,delimiter=',',header='time,x,y,z',comments='')
    return dict(Exudyn=exu.__version__,rigid=rigid,elements_per_link=0 if rigid else elements,
                step=step,endTime=end_time,elapsedSeconds=elapsed,output=str(out),samples=len(values),
                shear_stiffness_factor=100,newton_residual_mode=1,
                relative_newton_tolerance=1e-10,absolute_newton_tolerance=1e-10)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory',type=Path,default=Path(__file__).resolve().parents[1]/'data/space_7dof_all_flexible')
    parser.add_argument('--elements',type=int,default=128)
    parser.add_argument('--step',type=float)
    parser.add_argument('--end',type=float)
    parser.add_argument('--case',choices=['both','flexible','rigid'],default='both')
    parser.add_argument('--suffix',default='')
    args=parser.parse_args()
    cfg=json.loads((args.directory/'model.json').read_text())
    torque=np.loadtxt(args.directory/'torques.csv',delimiter=',',skiprows=1)
    reports=[]
    for rigid in (False,True):
        name='rigid' if rigid else 'flexible'
        if args.case not in ('both',name): continue
        print(f'Exudyn all-{name}, elements/link={args.elements}',flush=True)
        case_step=args.step or (min(cfg['step'],1.25e-4) if rigid else cfg['step'])
        reports.append(simulate(cfg,torque,rigid,args.elements,case_step,
            args.end or cfg['endTime'],args.directory/f'exudyn_{name}_tip{args.suffix}.csv'))
    (args.directory/f'exudyn_report{args.suffix}.json').write_text(json.dumps(reports,indent=2))


if __name__=='__main__':
    main()
