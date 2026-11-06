function tf = is_casadi_type(value)
% Numeric validation is deferred while building a CasADi computation graph.
    tf = isa(value, 'casadi.MX') || isa(value, 'casadi.SX') || isa(value, 'casadi.DM');
end
