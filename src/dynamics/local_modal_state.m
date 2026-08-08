function [Phi, eta_i, eta_d_i] = local_modal_state(model, i, eta, eta_d)
Phi = zeros(6, 0);
if isfield(model, 'Phi') && numel(model.Phi) >= i && ...
        ~isempty(model.Phi{i}) && is_flexible_link(model, i)
    Phi = model.Phi{i};
else
    body = [];
    if is_flexible_link(model, i)
        body = local_body(model, i, true);
    end
    if ~isempty(body) && body.nEta > 0
        Phi = zeros(6, body.nEta);
    end
end

n = size(Phi, 2);
if n == 0
    eta_i = zeros(0, 1); eta_d_i = zeros(0, 1); 
    return;
end
eta_i = local_link_vector(model, i, eta, n);
eta_d_i = local_link_vector(model, i, eta_d, n);
end