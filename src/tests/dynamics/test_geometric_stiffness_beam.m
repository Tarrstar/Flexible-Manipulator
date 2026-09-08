%% Verify the optional first-order geometric-stiffness correction
clear; clc;

thisFile = mfilename('fullpath');
testDir = fileparts(thisFile);
srcDir = fileparts(fileparts(testDir));
addpath(fullfile(srcDir, 'classes'), fullfile(srcDir, 'utils'));

beam.length = 1.0;
beam.width = 0.032;
beam.thickness = 0.0015;
beam.density = 7.542e3;
beam.youngs_modulus = 1.983e11;
beam.poisson_ratio = 0.30;
beam.elements_num = 3;
beam.fixed_root = true;
beam.geometric_stiffness_plane = 'xy';

femBeam = FEM_Beam_GeometricStiffness(beam);
[~, Kfull] = femBeam.get_matrices_3d();
[shapeMatrix, stiffnessDofs] = femBeam.build_direct_fe_shape_matrix();
Kff = full(Kfull(stiffnessDofs, stiffnessDofs));
[u0, S, weights] = femBeam.build_ffrf_data(shapeMatrix, 4, 3, 3);
[C, D] = femBeam.get_geometric_stiffness_coefficients(shapeMatrix);

assert(norm(C - C.', 'fro') < 1e-12);
assert(norm(D - D.', 'fro') < 1e-12);
assert(min(eig(C)) > -1e-10 * max(1, norm(C, 2)));
assert(min(eig(D)) > -1e-10 * max(1, norm(D, 2)));

nEta = size(shapeMatrix, 2);
common.density = beam.density;
common.u0 = u0;
common.S = S;
common.weights = weights;
common.Kff = Kff;
common.Dff = zeros(nEta);
common.geometric_C = C;
common.geometric_D = D;
common.geometric_rotation_axis = [0; 0; 1];

disabledParams = common;
disabledParams.geometric_stiffening = false;
disabledBody = FlexBody(disabledParams);

legacyParams = rmfield(common, {'geometric_C', 'geometric_D', ...
    'geometric_rotation_axis'});
legacyBody = FlexBody(legacyParams);

enabledParams = common;
enabledParams.geometric_stiffening = true;
enabledBody = FlexBody(enabledParams);

eta = 1e-4 * (1:nEta).' / nEta;
q = [1; 0; 0; 0; zeros(3, 1); eta];
qdot = zeros(7 + nEta, 1);
[~, Gbar] = get_quaternion_G_matrices(q(1:4));
omega = 2.0;
qdot(1:4) = 0.25 * Gbar.' * [0; 0; omega];

forceDisabled = disabledBody.get_elastic_force(q, qdot);
forceLegacy = legacyBody.get_elastic_force(q, qdot);
assert_relative_error(forceLegacy, forceDisabled, 1e-12, ...
    'legacy FlexBody behavior');
forceEnabled = enabledBody.get_elastic_force(q, qdot);
expectedIncrement = omega^2 * D * eta;
assert_relative_error(forceEnabled(8:end) - forceDisabled(8:end), ...
    expectedIncrement, 1e-12, 'rotational geometric stiffness');

axialAcceleration = 0.3;
enabledBody.set_axial_base_acceleration(axialAcceleration);
forceWithAcceleration = enabledBody.get_elastic_force(q, qdot);
expectedIncrement = (omega^2 * D - axialAcceleration * C) * eta;
assert_relative_error(forceWithAcceleration(8:end) - ...
    forceDisabled(8:end), expectedIncrement, 1e-12, ...
    'axial-acceleration geometric stiffness');

enabledBody.set_geometric_stiffening(false);
forceSwitchedOff = enabledBody.get_elastic_force(q, qdot);
assert_relative_error(forceSwitchedOff, forceDisabled, 1e-12, ...
    'disabled geometric stiffness');

fprintf('Geometric-stiffness beam tests passed.\n');

function assert_relative_error(actual, expected, tolerance, label)
scale = max(1, norm(expected));
relativeError = norm(actual - expected) / scale;
assert(relativeError < tolerance, ...
    '%s error %.3e exceeds tolerance %.3e.', ...
    label, relativeError, tolerance);
end
