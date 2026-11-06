function x = local_link_vector(model, i, values, n)
if iscell(values)
    x = values{i}(:);
elseif isfield(model, 'modal_idx') && numel(model.modal_idx) >= i
    x = values(model.modal_idx{i}); 
    x = x(:);
else
    x = values(:);
end
if numel(x) ~= n
    error('Modal coordinate length does not match size(model.Phi{%d},2).', i);
end
end