function gamma = flex_connection_bias(v, vRelative, vFlex, XJ, vFlexPreJoint)
%FLEX_CONNECTION_BIAS Body-component derivative of the joint/elastic twist.
% Constant S and Phi. Including the component derivative of the endpoint
% mapping avoids treating Map*Phi*eta_dot as a constant motion subspace.
if nargin<4
    W = get_skew_sym_matrix(vRelative(1:3));
    mapDerivative = -[W*vFlex(1:3); W*vFlex(4:6)];
else
    % A prismatic joint shifts the child origin from the attachment point.
    % Map = XJ*Rpre; dot(XJ)=-crm(vJoint)*XJ. The rotation-only algorithm
    % is recovered when XJ contains no translational offset.
    W = get_skew_sym_matrix(vFlexPreJoint(1:3));
    mapDerivative = -crm(vRelative-vFlex)*vFlex ...
        - XJ*[W*vFlexPreJoint(1:3); W*vFlexPreJoint(4:6)];
end
gamma = crm(v)*vRelative + mapDerivative;
end
