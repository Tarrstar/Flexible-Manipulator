function Phi = endpoint_shape(model, child, parent, defaultPhi, nEta)
% Use an edge-specific endpoint map for branched flexible structures.
% Existing serial models continue to use model.Phi{parent}.
Phi = defaultPhi;
if isfield(model, 'Phi_to_child') && numel(model.Phi_to_child) >= child && ...
        ~isempty(model.Phi_to_child{child})
    Phi = model.Phi_to_child{child};
end

if ~is_flexible_link(model, parent)
    Phi = zeros(6, 0);
    return;
end
assert(size(Phi,1)==6 && size(Phi,2)==nEta, ...
    'Endpoint map must be 6-by-number-of-parent-modes.');
end
