function [M_beam, K_beam] = get_beam_matrices_3d(beam)
    femBeam = FEM_Beam(beam);
    [M_beam, K_beam] = femBeam.get_matrices_3d();
end
