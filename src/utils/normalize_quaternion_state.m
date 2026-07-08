function [q, qdot] = normalize_quaternion_state(q, qdot)

    q = q(:);
    qdot = qdot(:);

    theta = q(1:4);

    nTheta = norm(theta);

    if nTheta < eps
        error('Quaternion norm is too small.');
    end

    theta = theta / nTheta;

    thetaDot = qdot(1:4);

    % Project thetaDot onto the tangent space of unit quaternion.
    thetaDot = thetaDot - theta * (theta.' * thetaDot);

    q(1:4) = theta;
    qdot(1:4) = thetaDot;
end