function a_grav = get_gravity(model)
% Test double for the spatial-vector toolbox helper.
if isfield(model, 'gravity')
    a_grav = model.gravity(:);
else
    a_grav = zeros(6, 1);
end
end
