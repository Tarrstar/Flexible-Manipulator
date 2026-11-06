function Phi = local_modal_state(model, i)
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

end