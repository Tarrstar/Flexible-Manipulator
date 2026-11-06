function J=inertia_from_moment(T)
J=trace(T)*eye(3)-T; J=(J+J.')/2;
end

