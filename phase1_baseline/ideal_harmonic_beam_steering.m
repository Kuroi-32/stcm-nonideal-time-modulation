%% ========================================================================
%  Array-Level Far-Field Reproduction of Fig. 2 (Zhang et al. 2018)
%  Harmonic beam steering via a time-gradient sequence across 8 elements.
%  Validation target: theta_m = arcsin(m/4) -> 0, +-14.48, +-30, +-48.59 deg
% =========================================================================
clear; clc; close all;

%% ---------------------- 1. INPUTS ----------------------
Gamma0 = 0.911*exp(1j*deg2rad(-102));   % CST-extracted OFF state
Gamma1 = 0.952*exp(1j*deg2rad(94));     % CST-extracted ON state

seq = [1 0 0 0 0 0 0 0];   % L=8 time-coding sequence (one full period)
L   = length(seq);

f0  = 1e6;                 % modulation frequency [Hz] (placeholder)
T0  = 1/f0;

Nq   = 8;                  % number of elements along y (matches paper)
dy_over_lambda = 0.5;      % d_y = lambda_c/2, as in the paper

Mmax = 3;                  % harmonics to evaluate: m = -3..+3
m_idx = -Mmax:Mmax;

theta_deg = -90:0.05:90;
theta = deg2rad(theta_deg);

%% ---------------------- 2. SINGLE-ELEMENT FOURIER COEFFICIENTS a^m ----------------------
% Same analytical engine as before (Zhang et al. Eq. 6), now only need m = -Mmax..Mmax
GammaSeq = zeros(1,L);
for k = 1:L
    GammaSeq(k) = (seq(k)==0)*Gamma0 + (seq(k)==1)*Gamma1;
end

a_m = zeros(size(m_idx));
for idx = 1:length(m_idx)
    m = m_idx(idx);
    s = 0;
    for k = 1:L
        s = s + (GammaSeq(k)/L) * sinc(m/L) * exp(-1j*pi*m*(2*k-1)/L);
    end
    a_m(idx) = s;
end

%% ---------------------- 3. TIME-GRADIENT SEQUENCE ACROSS ELEMENTS ----------------------
% Element q (q = 1..Nq) uses the SAME Gamma(t), just delayed by t_q = (q-1)*T0/Nq.
% Supplementary Eq. (4): a^m(shifted) = a^m * exp(-j*2*pi*m*f0*t_q)
t_q = (0:Nq-1) * T0/Nq;

a_qm = zeros(Nq, length(m_idx));   % rows = elements, cols = harmonics
for q = 1:Nq
    phase_shift = exp(-1j*2*pi*m_idx*f0*t_q(q));
    a_qm(q,:) = a_m .* phase_shift;
end

%% ---------------------- 4. ARRAY FACTOR / FAR-FIELD (Eq. 7, isotropic Epq=1) ----------------------
% F_m(theta) = sum_q exp[j*2*pi*(q-1)*dy_over_lambda*sin(theta)] * a_qm(q,m)
F = zeros(length(m_idx), length(theta));
for mi = 1:length(m_idx)
    for ti = 1:length(theta)
        array_phase = exp(1j*2*pi*(0:Nq-1)*dy_over_lambda*sin(theta(ti)));
        F(mi,ti) = sum(array_phase .* a_qm(:,mi).');
    end
end

F_dB = 20*log10(abs(F) / max(abs(F(:))));

%% ---------------------- 5. EXTRACT PEAK STEERING ANGLES ----------------------
peak_angle = zeros(size(m_idx));
for mi = 1:length(m_idx)
    [~, pk] = max(abs(F(mi,:)));
    peak_angle(mi) = theta_deg(pk);
end

target_angle = asind(m_idx/4);   % theoretical: theta_m = arcsin(m/4)

fprintf('%6s %12s %12s %10s\n','m','Simulated','Theoretical','Diff(deg)');
for mi = 1:length(m_idx)
    fprintf('%6d %12.2f %12.2f %10.2f\n', m_idx(mi), peak_angle(mi), target_angle(mi), ...
             peak_angle(mi)-target_angle(mi));
end

%% ---------------------- 6. PLOTS ----------------------
figure('Name','Harmonic beam steering (Fig. 2 reproduction)');
hold on;
colors = lines(length(m_idx));
legend_entries = cell(1,length(m_idx));
for mi = 1:length(m_idx)
    plot(theta_deg, F_dB(mi,:), 'Color', colors(mi,:), 'LineWidth',1.3);
    legend_entries{mi} = sprintf('m = %+d', m_idx(mi));
end
xlabel('\theta (degrees)'); ylabel('Normalized scattering pattern (dB)');
title('Harmonic beam steering via time-gradient sequence');
legend(legend_entries, 'Location','eastoutside');
ylim([-40 0]); xlim([-90 90]); grid on;

figure('Name','Simulated vs theoretical steering angles');
plot(m_idx, target_angle, 'bo-', 'DisplayName','Theoretical: arcsin(m/4)'); hold on;
plot(m_idx, peak_angle, 'rx', 'MarkerSize',10, 'DisplayName','Simulated peak');
xlabel('Harmonic index m'); ylabel('Steering angle \theta_m (deg)');
title('Validation: simulated vs. theoretical beam-steering angles');
legend('Location','best'); grid on;