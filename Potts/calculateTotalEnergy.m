function total_energy = calculateTotalEnergy(lattice, J, h)
    N = size(lattice, 1);
    total_energy = 0;
    
    % Loop over all sites
    for i = 1:N
        for j = 1:N
            % Get current spin
            spin = lattice(i,j);
            
            % Get right neighbor (periodic boundary)
            right_neighbor = lattice(i, mod(j,N)+1);
            % Get down neighbor (periodic boundary)
            down_neighbor = lattice(mod(i,N)+1, j);
            
            % Add interaction energies
            % Only count right and down to avoid double counting
            total_energy = total_energy - J * (double(spin == right_neighbor) + ...
                                             double(spin == down_neighbor));
            
            % Add external field contribution
            total_energy = total_energy - h * spin;
        end
    end
end