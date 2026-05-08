function a()
    mc = 1.4399;      
    alpha_s = 0.4827; 
    b = 0.1488;       
    sigma = 1.2819;   

    fm_to_inv_GeV = 1 / 0.197326; 
    a_fm = 10.0; 
    a = a_fm * fm_to_inv_GeV; 
    
    N_max = 300;     
    N_points = 8000; 
   
    l = 0; 
    S = 1; 
    spin_factor = S*(S + 1) - 1.5;

    r = linspace(1e-6, a, N_points);

    V_cent = (l * (l + 1)) ./ (mc * r.^2); 
    V_coul = - (4 * alpha_s) ./ (3 * r);
    V_lin = b * r;
    
    delta_sigma = (sigma^3 / pi^1.5) * exp(-(sigma * r).^2);
    V_spin = (16 * pi * alpha_s / (9 * mc^2)) * spin_factor * delta_sigma;
    V_eff = V_cent + V_coul + V_lin + V_spin;

   
    H = zeros(N_max, N_max);
    mu = mc / 2; % equal-mass quarks (mc * mc / (mc + mc))

    fprintf('Assembling the %dx%d Hamiltonian matrix...\n', N_max, N_max);
    
    for n = 1:N_max
        psi_n = sqrt(2/a) * sin(n * pi * r / a);
        H(n,n) = (n^2 * pi^2) / (2 * mu * a^2);
        
        for m = n:N_max
            psi_m = sqrt(2/a) * sin(m * pi * r / a);
            integrand = psi_n .* V_eff .* psi_m;
            H_pot = trapz(r, integrand);
            H(n,m) = H(n,m) + H_pot;
            if n ~= m
                H(m,n) = H(n,m);
            end
        end
    end

    
    eigenvalues = eig(H);
    eigenvalues = sort(eigenvalues); % Ensure lowest energy states are first

    masses = eigenvalues + 2 * mc;

    fprintf('\n--- Charmonium Mass Spectrum (l=%d, S=%d) ---\n', l, S);
    fprintf('State 1 (e.g., J/psi or eta_c) : %.4f GeV\n', masses(1));
    fprintf('State 2 (e.g., psi(2S))        : %.4f GeV\n', masses(2));
    fprintf('State 3 (e.g., psi(3S))        : %.4f GeV\n', masses(3));
    fprintf('State 4 (e.g., psi(4S))        : %.4f GeV\n', masses(4));
end