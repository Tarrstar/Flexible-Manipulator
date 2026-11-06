function report = test_7dof_all_flexible_exudyn(action, step, endTime, outputDir)
% Seven spatial revolute joints, seven flexible links, and a rigid baseline.
% prepare: export one common physical model and prescribed torque history.
% simulate: integrate the current ID_flex model and the equal-inertia rigid arm.
if nargin<1, action='simulate'; end
if nargin<2, step=2.5e-4; end
if nargin<3, endTime=2; end
projectDir=fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
if nargin<4, outputDir=fullfile(projectDir,'data','space_7dof_all_flexible'); end
if ~exist(outputDir,'dir'), mkdir(outputDir); end
srcDir=fullfile(projectDir,'src');
addpath(fullfile(srcDir,'classes'),fullfile(srcDir,'dynamics'), ...
    fullfile(srcDir,'kinematics'),fullfile(srcDir,'utils'), ...
    fullfile(srcDir,'integrators'),fullfile(srcDir,'tests','dynamics','support'));
addpath('D:/libraries/casadi/casadi-windows-matlabR2016a-v3.5.5');

cfg.lengths=[.28;.30;.26;.22;.18;.14;.10];
cfg.jointAxes=[0 0 1;0 1 0;0 1 0;1 0 0;0 1 0;1 0 0;0 1 0];
cfg.jointTypes={'Rz';'Ry';'Ry';'Rx';'Ry';'Rx';'Ry'};
cfg.width=.025; cfg.thickness=.025;
cfg.density=2700; cfg.youngsModulus=7e8; cfg.poissonRatio=.33;
cfg.gravity=[0;0;0]; cfg.step=step; cfg.endTime=endTime;
cfg.rigidStep=min(step,1.25e-4);
cfg.driveTime=1; cfg.goalAngles=[.25;-.35;.30;.20;-.25;.18;.22];
cfg.matlabElements=4; cfg.interiorModes=0; cfg.spectralRadius=.9;
cfg.initialConditions='All joint/modal coordinates and velocities are zero';
cfg.flexibleModel='Linear Euler-Bernoulli Craig-Bampton FFRF, 6 interface coordinates per link';
cfg.referenceModel='Exudyn geometrically exact 3D beams';
[flexModel,rigidModel,tipShapes]=build_models(cfg);
assert(all(strcmp(flexModel.link_type,'flexible')) && flexModel.NB==7);
assert(flexModel.modal_num==42);

% Feed-forward torques from an equal-mass, equal-inertia rigid arm. The same
% table is supplied to every simulation; there is no trajectory constraint.
if strcmpi(action,'prepare') || ~exist(fullfile(outputDir,'torques.csv'),'file')
    torqueTime=linspace(0,endTime,round(endTime/(step/4))+1).';
    torques=zeros(numel(torqueTime),7);
    for k=1:numel(torqueTime)
        [q,v,a]=reference_motion(torqueTime(k),cfg);
        torques(k,:)=ID_flex(rigidModel,q,v,a,[],[],[]).';
    end
    write_csv(fullfile(outputDir,'torques.csv'),[torqueTime,torques], ...
        {'time','tau1','tau2','tau3','tau4','tau5','tau6','tau7'});
    fid=fopen(fullfile(outputDir,'model.json'),'w');
    cleanup=onCleanup(@() fclose(fid));
    fprintf(fid,'%s',jsonencode(cfg,PrettyPrint=true)); clear cleanup;
else
    torqueData=readmatrix(fullfile(outputDir,'torques.csv'));
    torqueTime=torqueData(:,1); torques=torqueData(:,2:8);
end
fprintf('7-DOF model: seven flexible links, %d elastic coordinates; gravity=0.\n',flexModel.modal_num);
report=cfg;
if strcmpi(action,'prepare'), return; end

names={'flexible','rigid'};
models={flexModel,rigidModel};
cases=1:2;
if strcmpi(action,'rigid'), cases=2; end
for c=cases
    model=models{c}; n=model.ndof;
    caseStep=step;
    if c==2, caseStep=cfg.rigidStep; end
    time=(0:round(endTime/caseStep)).'*caseStep;
    fprintf('Building ID_flex AD integration function: %s (%d DOFs)...\n',names{c},n);
    setup=tic; dynamics=ID_AD_func(model); setupSeconds=toc(setup);
    force=@(t) [interp1(torqueTime,torques,t,'linear').';zeros(model.modal_num,1)];
    options.MaxStep=caseStep; options.AbsTol=1e-10; options.RelTol=1e-9;
    options.MaxIterations=20; options.SpectralRadius=cfg.spectralRadius;
    fprintf('Integrating MATLAB %s, h=%g s, T=%g s...\n',names{c},caseStep,endTime);
    [t,x,info]=generalized_alpha_id_integrator(dynamics,force,time,zeros(2*n,1),options);
    tip=zeros(numel(t),3); maxElastic=0; maxRotation=0;
    for k=1:numel(t)
        tip(k,:)=endpoint_position(model,cfg,x(k,1:n).',tipShapes).';
        if model.modal_num>0
            for i=1:7
                e=x(k,7+model.modal_idx{i}).';
                maxElastic=max(maxElastic,norm(tipShapes{i}*e));
                maxRotation=max(maxRotation,norm(model.Phi{i}(1:3,:)*e));
            end
        end
    end
    assert(all(isfinite(x),'all') && all(isfinite(tip),'all'));
    write_csv(fullfile(outputDir,['matlab_',names{c},'_tip.csv']),[t,tip],{'time','x','y','z'});
    report.(names{c})=struct('step',caseStep,'setupSeconds',setupSeconds, ...
        'integrationSeconds',info.integrationSeconds,'maxTipElasticDisplacement',maxElastic, ...
        'maxInterfaceElasticRotation',maxRotation,'maxIterations',info.maximumIterations);
    save(fullfile(outputDir,['matlab_',names{c},'_state.mat']),'t','x','info','cfg');
    fprintf('MATLAB %s completed: %.3f s integration, max local elastic displacement %.3e m.\n', ...
        names{c},info.integrationSeconds,maxElastic);
end
fid=fopen(fullfile(outputDir,'matlab_report.json'),'w');
cleanup=onCleanup(@() fclose(fid)); fprintf(fid,'%s',jsonencode(report,PrettyPrint=true));
end

function [model,rigid,tipShapes]=build_models(cfg)
N=7; model.NB=N; model.parent=(0:N-1).'; model.jtype=cfg.jointTypes;
model.Xtree=cell(N,1); model.Phi=cell(N,1); model.modal_idx=cell(N,1);
model.flexbody=cell(N,1); model.rigidbody=cell(N,1); tipShapes=cell(N,1);
model.Xtree{1}=eye(6); rigidBodies=cell(N,1); offset=0;
for i=1:N
    L=cfg.lengths(i);
    if i>1, model.Xtree{i}=xlt([cfg.lengths(i-1),0,0]); end
    beam=struct('length',L,'width',cfg.width,'thickness',cfg.thickness, ...
        'density',cfg.density,'youngs_modulus',cfg.youngsModulus, ...
        'poisson_ratio',cfg.poissonRatio,'elements_num',cfg.matlabElements, ...
        'fixed_root',true,'nInteriorModes',cfg.interiorModes, ...
        'boundaryNodeIds',[1,cfg.matlabElements+1]);
    fem=FEM_Beam_CB_Reduction(beam);
    [u0,S,weights,shape,K]=fem.build_ffrf_data(cfg.interiorModes,4,3,3);
    m=size(shape,2);
    p=struct('density',cfg.density,'u0',u0,'S',S,'weights',weights, ...
        'Kff',K,'Dff',zeros(m));
    body=FlexBody(p); model.flexbody{i}=body;
    nodeDofs=fem.get_node_dofs(cfg.matlabElements+1);
    tipNode=shape(nodeDofs,:);
    model.Phi{i}=[tipNode(4:6,:);tipNode(1:3,:)];
    tipShapes{i}=fem.evaluate_displacement(shape,L,0,0);
    model.modal_idx{i}=offset+(1:m); offset=offset+m;
    mass=cfg.density*cfg.width*cfg.thickness*L;
    inertia=mass/12*diag([cfg.width^2+cfg.thickness^2, ...
        L^2+cfg.thickness^2,L^2+cfg.width^2]);
    rigidBodies{i}=RigidBody(struct('mass',mass,'inertia',inertia,'com',[L/2;0;0]));
    [M,~,~,~,~]=body.update_matrices(zeros(6,1),zeros(m,1),zeros(m,1));
    assert(norm(M-rigidBodies{i}.get_spatial_inertia(),'fro')<1e-12);
end
model.link_type=repmat({'flexible'},N,1); model.gravity=cfg.gravity;
model.active_joints=(1:N).'; model.joint_num=N; model.modal_num=offset; model.ndof=N+offset;
rigid=model; rigid.flexbody=cell(N,1); rigid.rigidbody=rigidBodies;
rigid.link_type=repmat({'rigid'},N,1); rigid.Phi=repmat({zeros(6,0)},N,1);
rigid.modal_idx=repmat({[]},N,1); rigid.modal_num=0; rigid.ndof=N;
end

function [q,v,a]=reference_motion(t,cfg)
s=min(max(t/cfg.driveTime,0),1);
q=cfg.goalAngles*(10*s^3-15*s^4+6*s^5);
v=cfg.goalAngles*(30*s^2-60*s^3+30*s^4)/cfg.driveTime;
a=cfg.goalAngles*(60*s-180*s^2+120*s^3)/cfg.driveTime^2;
end

function tip=endpoint_position(model,cfg,coordinates,tipShapes)
R=eye(3); p=zeros(3,1);
for i=1:7
    if i>1
        delta=zeros(6,1); Rf=eye(3);
        if model.modal_num>0
            e=coordinates(7+model.modal_idx{i-1});
            [Rf,~,delta]=get_deformation_transform(model.Phi{i-1},e);
        end
        p=p+R*([cfg.lengths(i-1);0;0]+delta(4:6)); R=R*Rf;
    end
    XJ=jcalc(model.jtype{i},coordinates(i)); R=R*XJ(1:3,1:3).';
end
local=[cfg.lengths(7);0;0];
if model.modal_num>0, local=local+tipShapes{7}*coordinates(7+model.modal_idx{7}); end
tip=p+R*local;
end

function write_csv(path,data,headers)
writetable(array2table(data,'VariableNames',headers),path);
end
