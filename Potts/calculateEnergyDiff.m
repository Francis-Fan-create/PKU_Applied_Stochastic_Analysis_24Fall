% Function to calculate energy difference for a spin flip
function delta_E = calculateEnergyDiff(lattice, i, j, new_spin, J, h)
    N = size(lattice, 1);
    old_spin = lattice(i,j);
    
    % Get neighboring spins (periodic boundary conditions)
    up = lattice(mod(i-2,N)+1, j);
    down = lattice(mod(i,N)+1, j);
    left = lattice(i, mod(j-2,N)+1);
    right = lattice(i, mod(j,N)+1);
    
    % Calculate energy difference
    old_energy = -J * sum([up,down,left,right] == old_spin) - h*old_spin;
    new_energy = -J * sum([up,down,left,right] == new_spin) - h*new_spin;
    delta_E = new_energy - old_energy;
end