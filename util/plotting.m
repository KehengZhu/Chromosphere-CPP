% Define constants
CNI = 1;        % Ion density
CNN = 2;        % Neutral density
CNV = 3;        % Ion velocity
CNU = 4;        % Neutral velocity
CEI = 5;        % Ion energy
CEN = 6;        % Neutral energy
% Primitive Variable Indices
PNI = 1;        % Ion density
PNN = 2;        % Neutral density
PV  = 3;        % Ion velocity
PU  = 4;        % Neutral velocity
PPI = 5;        % Ion pressure
PPN = 6;        % Neutral pressure

% Physical constants
gammamono = 5.0 / 3.0;  % Ratio of specific heats (adiabatic index for a monatomic ideal gas)
alpha_p = 0.18;  % Divergence error propagation parameter
pi_const = 3.14159265358979323846;  % Value of π (pi)
m_i = 1.6726219e-27;  % Mass of proton in kilograms
m_n = m_i;  % Mass of neutron equals mass of proton
m_e = 9.10938356e-31;  % Mass of electron in kilograms
g = 0.27395e3;  % Gravitational acceleration in normalized units (km/s^2)
k_b = 1.380649e-23;  % Boltzmann constant in SI units
e_ = 1.602176634e-19;  % Elementary charge in Coulombs

% Read data
[xx, xn_III] = read_XIII('../chromo_build/output.txt');

% Process data
nn = xn_III(:, CNN) / m_n;
ni = xn_III(:, CNI) / m_i;
uu = xn_III(:, CNU) ./ (m_n * nn);
vv = xn_III(:, CNV) ./ (m_n * ni);
pi = 2.0/3.0 * xn_III(:, CEI) - 1.0/3.0 * m_i * ni .* vv .* vv;
pn = 2.0/3.0 * xn_III(:, CEN) - 1.0/3.0 * m_n * nn .* uu .* uu;
Ti = pi ./ (2.0 * ni * k_b);
Tn = pn ./ (nn * k_b);

% Plotting
figure(1);
sgtitle('Chromo');

subplot(2,4,1);
plot(xx, nn);
title('nn (m^{-3})');

subplot(2,4,5);
plot(xx, ni);
title('ni (m^{-3})');
xlabel('height (m)');

subplot(2,4,2);
plot(xx, uu);
title('u (m/s)');

subplot(2,4,6);
plot(xx, vv);
title('v (m/s)');
xlabel('height (m)');

subplot(2,4,3);
plot(xx, pn);
title('Pn (Pa)');

subplot(2,4,7);
plot(xx, pi);
title('Pi (Pa)');
xlabel('height (m)');

subplot(2,4,4);
plot(xx, Tn);
title('Tn (K)');

subplot(2,4,8);
plot(xx, Ti);
title('Ti (K)');
xlabel('height (m)');

% Function to read data
function [xx, xn_III] = read_XIII(filename)
    % Read the entire file into cell array of lines
    fid = fopen(filename, 'r');
    lines = {};
    tline = fgetl(fid);
    while ischar(tline)
        lines{end+1} = tline;
        tline = fgetl(fid);
    end
    fclose(fid);

    % Read the first line and store it in xx as an array
    xx = sscanf(lines{1}, '%f')';

    % Determine ns and n_eq from the shape of the data
    ns = length(xx);  % Number of elements in the first line
    n_eq = numel(sscanf(lines{2}, '%f'));  % Number of columns in the rest of the data

    % Initialize an array for ns rows and n_eq columns
    xn_III = zeros(ns, n_eq);

    % Populate xn_III with the remaining data
    for i = 1:ns
        xn_III(i, :) = sscanf(lines{i+1}, '%f')';
    end
end
