function body = local_body(model, i, isFlexible)
body = [];
if isFlexible
    names = {'flexbody', 'flex_body', 'FlexBody'};
    className = 'FlexBody';
else
    names = {'rigidbody', 'rigid_body', 'RigidBody'};
    className = 'RigidBody';
end
for name = names
    name = name{1};
    if ~isfield(model, name) || numel(model.(name)) < i
        continue
    end
    if iscell(model.(name))
        candidate = model.(name){i};
    else
        candidate = model.(name)(i);
    end
    if ~isempty(candidate)
        if isa(candidate, className)
            body = candidate;
            return
        end
    end
end
if isempty(body)
    error('Missing %s object for link %d.', className, i);
end
end
