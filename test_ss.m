% Find K_p that gives L = 0.33 in steady state
alpha = 0.35; sigma_L = 2.0; omega_E = 0.045; sigma_E = 0.6;
delta_p = 0.025; delta_b = 0.03; delta_g = 0.0125;
u = 0.97; E = 0.15;
I_bat = delta_b * 0.19;
I_grid = delta_g * 1.14;
L_target = 0.33;

fprintf('Finding K_p for L = %.2f with sigma_L = %.1f\n\n', L_target, sigma_L);

% Solve for K_p
% W_demand = (1-alpha) * (u*K_p)^alpha * L^(-alpha) = (1-alpha) * u^alpha * K_p^alpha * L^(-alpha)
% W_supply = C * L^sigma_L
% C = Y - delta_p*K_p - I_bat - I_grid
% Y = CES(V, E), V = (u*K_p)^alpha * L^(1-alpha)

res_fun = @(Kp) kp_residual(Kp, L_target, u, alpha, sigma_L, omega_E, sigma_E, E, delta_p, I_bat, I_grid);

% Scan
for Kp = [0.5, 1, 2, 3, 5, 7, 10, 15]
    r = res_fun(Kp);
    fprintf('  K_p=%.1f: residual=%.6f\n', Kp, r);
end

% Solve
[Kp_sol, fval, exitflag] = fzero(res_fun, [0.1, 20]);
fprintf('\nSolution: K_p=%.6f, residual=%.12e\n', Kp_sol, fval);

% Compute full SS
K_p = Kp_sol;
uk = u*K_p;
L = L_target;
V = uk^alpha * L^(1-alpha);
Y = ((1-omega_E)*V^((sigma_E-1)/sigma_E) + omega_E*E^((sigma_E-1)/sigma_E))^(sigma_E/(sigma_E-1));
I_p = delta_p * K_p;
C = Y - I_p - I_bat - I_grid;
W = (1-alpha)*V/L;
R = alpha*V/(u*K_p);

fprintf('\nFull Steady State:\n');
fprintf('  K_p   = %.6f\n', K_p);
fprintf('  L     = %.6f\n', L);
fprintf('  V     = %.6f\n', V);
fprintf('  Y     = %.6f\n', Y);
fprintf('  C     = %.6f\n', C);
fprintf('  W     = %.6f\n', W);
fprintf('  R     = %.6f\n', R);
fprintf('  I_p   = %.6f\n', I_p);
fprintf('  C/Y   = %.4f\n', C/Y);
fprintf('  I_p/Y = %.4f\n', I_p/Y);
fprintf('  K_p/Y = %.4f\n', K_p/Y);

% Also try with sigma_L = 1
sigma_L2 = 1.0;
res_fun2 = @(Kp) kp_residual(Kp, L_target, u, alpha, sigma_L2, omega_E, sigma_E, E, delta_p, I_bat, I_grid);
[Kp_sol2, ~] = fzero(res_fun2, [0.1, 20]);

K_p2 = Kp_sol2;
V2 = (u*K_p2)^alpha * L_target^(1-alpha);
Y2 = ((1-omega_E)*V2^((sigma_E-1)/sigma_E) + omega_E*E^((sigma_E-1)/sigma_E))^(sigma_E/(sigma_E-1));
C2 = Y2 - delta_p*K_p2 - I_bat - I_grid;
fprintf('\nWith sigma_L=1.0: K_p=%.4f, Y=%.4f, C=%.4f, C/Y=%.4f\n', K_p2, Y2, C2, C2/Y2);

function res = kp_residual(K_p, L, u, alpha, sigma_L, omega_E, sigma_E, E, delta_p, I_bat, I_grid)
    uk = u*K_p;
    V = uk^alpha * L^(1-alpha);
    Y = ((1-omega_E)*V^((sigma_E-1)/sigma_E) + omega_E*E^((sigma_E-1)/sigma_E))^(sigma_E/(sigma_E-1));
    I_p = delta_p * K_p;
    C = Y - I_p - I_bat - I_grid;
    if C <= 0
        res = 1e10; return;
    end
    W_demand = (1-alpha) * V / L;
    W_supply = C * L^sigma_L;
    res = W_demand - W_supply;
end
