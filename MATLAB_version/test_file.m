%% test 1
% FEM solver for EM modes in dielectric waveguides
% Original formulation is based on: https://doi.org/10.1080/02726340290084012
%
% Usage:
%   Run the script directly for a demo 
%   at 1550 nm, or call compute_modes() from your own script.
%
% Requirements: MATLAB R2020b+ (for complex sparse eigs)

%%
clear; close all;

wavelength = 1.55; % µm
% ---- Geometry (all lengths in micrometres) ----
w_core   = 1.6;   % core width
h_core   = 0.7;   % core height
h_clad   = 2.7;   % cladding height above core
h_box    = 2.00;   % buried-oxide thickness
w_sim    = 6.00;   % total simulation width
% --- materials
n_core   = get_refractive_index('Si3N4', wavelength); % refractive index of the core
n_clad   = get_refractive_index('SiO2', wavelength);  % refractive index of the cladding
n_box    = n_clad;  % refractive index of the substrate

num_modes  = 6;         % number of modes to search for
mesh_res   = round(w_sim/wavelength*100);   % approximate triangles along the longest side

compute_overlaps = 0; % set it to 1, if you want to calculate the mode overlap matrix
plot_mesh = 0;
% end of inputs


fprintf('Building mesh ...\n');
[nodes, elems, epsilon_r, regions] = build_soi_mesh( ...
    w_core, h_core, h_clad, h_box, w_sim, ...
    n_core, n_clad, n_box, mesh_res);

fprintf('Nodes: %d,  Elements: %d\n', size(nodes,1), size(elems,1));


if plot_mesh ==1
    figure('Name', 'Mesh Visualization (Permittivity)', 'Color', 'w');
    ph = patch('Faces', elems, 'Vertices', nodes, ...
               'FaceVertexCData', epsilon_r, ...
               'FaceColor', 'flat', ...
               'EdgeColor', 'k', ...
               'LineWidth', 0.1);
    axis equal tight;
    xlabel('x [\mum]'); ylabel('y [\mum]');
    title('Zoomed in Mesh Material Distribution (\epsilon_r)');
    cb = colorbar;
    ylabel(cb, 'Relative Permittivity \epsilon_r');
    xlim([-w_core*1.2, w_core*1.2]); 
    ylim([-h_core*0.2, h_core*1.2]);
    % print -dpng figure_mesh
end

tic
fprintf('Assembling & solving eigenvalue problem ...\n');
modes = compute_modes(nodes, elems, epsilon_r, wavelength, ...
    'num_modes', num_modes, 'mu_r', 1.0);
toc
%% ---- Print results ----
fprintf('\n--- Guided modes ---\n');
for m = 1:length(modes)
    fprintf('Mode %d:  n_eff = %.6f + %.2ei,  TE-frac = %.3f\n', ...
        m, real(modes(m).n_eff), imag(modes(m).n_eff), modes(m).te_fraction);
end

%% ---- Plot dominant mode ----
plot_mode_fields(modes(1), nodes, elems, 'Mode 1 (fundamental TE)');
plot_mode_fields(modes(2), nodes, elems, 'Mode 2 ');
plot_mode_fields(modes(3), nodes, elems, 'Mode 3');

if compute_overlaps ==1
%% ---- Overlap matrix ----
    fprintf('\nComputing overlap matrix ...\n');
    N = length(modes);
    OL = zeros(N);
    for i = 1:N
        for j = 1:N
            OL(i,j) = calculate_overlap(modes(i), modes(j));
        end
    end
    figure('Name','Overlap matrix');
    imagesc(real(OL)); colorbar; axis square;
    title('Re\{Overlap integrals\}');
    xlabel('Mode j'); ylabel('Mode i');
    set(gca,'XTick',1:N,'YTick',1:N);
end


