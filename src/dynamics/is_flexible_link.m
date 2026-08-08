function tf = is_flexible_link(model, i)
types = model.link_type;
if iscell(types)
    type = types{i};
elseif isstring(types)
    if isscalar(types)
        type = types;
    else
        type = types(i);
    end
elseif ischar(types)
    type = types;
else
    error('model.link_type must be a character or string array/cell array.');
end
tf = strcmpi(strtrim(char(type)), 'flexible');
if ~tf && ~strcmpi(strtrim(char(type)), 'rigid')
    error('model.link_type{%d} must be ''rigid'' or ''flexible''.', i);
end
end
