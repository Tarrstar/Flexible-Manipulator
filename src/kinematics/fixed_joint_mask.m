function fixed = fixed_joint_mask(model)

fixed=false(model.NB,1);
if ~isfield(model,'fixed_joints') || isempty(model.fixed_joints)
    return;
end
spec=model.fixed_joints;
if islogical(spec)
    assert(numel(spec)==model.NB, 'Fixed-joint mask must have NB entries.');
    fixed=spec(:);
else
    assert(all(isfinite(spec(:))) && all(spec(:)==fix(spec(:))) && ...
        all(spec(:)>=1 & spec(:)<=model.NB), 'Invalid fixed-joint indices.');
    fixed(spec(:))=true;
end
end
