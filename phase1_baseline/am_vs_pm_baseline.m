%% ========================================================================
%  Baseline Reproduction — Part 2
%  (1) Clean ideal-case harmonic efficiency reference (tr = 0, no non-ideality)
%  (2) AM vs PM comparison, reproducing Fig. 2c-f of Zhang et al. 2018
%  Reuses the SAME validated array engine as fig2_array_beam_steering.m
% =========================================================================
clear; clc; close all;

%% ---------------------- SHARED SETTINGS ----------------------
% NOTE: this must be a LOW-DUTY-CYCLE pulse (single "1" in the sequence),
% matching Fig. 2a/b's diagonal identity-matrix coding matrix. A symmetric
% 50%-duty sequence (like [0 0 0 0 1 1 1 1]) cancels the DC/carrier term
% almost entirely, which is correct physics but the WRONG sequence for
% reproducing the AM/PM carrier-dominance comparison in Fig. 2e,f.
seq = [1 0 0 0 0 0 0 0];   % single-pulse sequence -> matches Fig. 2's coding matrix
L   = length(seq);
f0  = 1e6;
T0  = 1/f0;
Nq  = 8;
dy_over_lambda = 0.5;
Mmax = 3;
m_idx = -Mmax:Mmax;
theta_deg = -90:0.05:90;
theta = deg2rad(theta_deg);
t_q = (0:Nq-1) * T0/Nq;

%% ---------------------- PART 1: IDEAL BASELINE HARMONIC EFFICIENCY ----------------------
Gamma0_PM = 0.911*exp(1j*deg2rad(-102));   % CST-extracted OFF state
Gamma1_PM = 0.952*exp(1j*deg2rad(94));     % CST-extracted ON state

a_m_PM = compute_am(seq, L, Gamma0_PM, Gamma1_PM, m_idx);

P_n   = abs(a_m_PM).^2;
eta_n = P_n / sum(P_n);

figure('Name','Baseline harmonic efficiency (ideal, tr=0)');
bar(m_idx, eta_n);
xlabel('Harmonic index m'); ylabel('\eta_m (efficiency)');
title('Baseline (ideal) harmonic efficiency — reference for all Phase 2 comparisons');
grid on;

fprintf('=== Baseline ideal harmonic efficiency table ===\n');
fprintf('%6s %12s\n','m','eta_m');
for i = 1:length(m_idx)
    fprintf('%6d %12.4f\n', m_idx(i), eta_n(i));
end
fprintf('Carrier efficiency eta_0 = %.4f\n\n', eta_n(m_idx==0));

%% ---------------------- PART 2: AM vs PM COMPARISON (Fig. 2c-f) ----------------------
% PM (phase modulation): both states have |Gamma|~1, 196 deg apart (your CST values)
% AM (amplitude modulation): one state fully reflective, other near-zero reflection
Gamma0_AM = 0;                  % "off" -> no reflection
Gamma1_AM = 1*exp(1j*0);        % "on"  -> full reflection, 0 deg reference phase

a_m_AM = compute_am(seq, L, Gamma0_AM, Gamma1_AM, m_idx);

% --- equivalent amplitude/phase per element (Fig 2c,d style) ---
% Apply per-element time shift to get a_qm for both schemes
a_qm_PM = apply_time_shift(a_m_PM, m_idx, f0, t_q);
a_qm_AM = apply_time_shift(a_m_AM, m_idx, f0, t_q);

figure('Name','Equivalent amplitude/phase: AM vs PM (Fig. 2c,d)');
subplot(2,2,1); imagesc(m_idx, 1:Nq, abs(a_qm_AM)); colorbar;
xlabel('Harmonic m'); ylabel('Element q'); title('AM: equivalent amplitude'); axis xy;
subplot(2,2,2); imagesc(m_idx, 1:Nq, rad2deg(angle(a_qm_AM))); colorbar;
xlabel('Harmonic m'); ylabel('Element q'); title('AM: equivalent phase (deg)'); axis xy;
subplot(2,2,3); imagesc(m_idx, 1:Nq, abs(a_qm_PM)); colorbar;
xlabel('Harmonic m'); ylabel('Element q'); title('PM: equivalent amplitude'); axis xy;
subplot(2,2,4); imagesc(m_idx, 1:Nq, rad2deg(angle(a_qm_PM))); colorbar;
xlabel('Harmonic m'); ylabel('Element q'); title('PM: equivalent phase (deg)'); axis xy;

% --- far-field patterns for both schemes (Fig 2e,f style) ---
F_AM = array_farfield(a_qm_AM, m_idx, Nq, dy_over_lambda, theta);
F_PM = array_farfield(a_qm_PM, m_idx, Nq, dy_over_lambda, theta);

% Normalize PM relative to AM peak, as the paper does (Fig 2f is normalized to Fig 2e)
global_ref = max(abs(F_AM(:)));
F_AM_dB = 20*log10(abs(F_AM)/global_ref);
F_PM_dB = 20*log10(abs(F_PM)/global_ref);

figure('Name','AM vs PM scattering patterns (Fig. 2e,f)');
subplot(1,2,1);
plot(theta_deg, F_AM_dB', 'LineWidth',1.2);
xlabel('\theta (deg)'); ylabel('Normalized scattering pattern (dB)');
title('Amplitude Modulation (AM)'); ylim([-50 15]); xlim([-90 90]); grid on;
legend(arrayfun(@(m) sprintf('m=%+d',m), m_idx, 'UniformOutput',false), 'Location','eastoutside');

subplot(1,2,2);
plot(theta_deg, F_PM_dB', 'LineWidth',1.2);
xlabel('\theta (deg)'); ylabel('Normalized scattering pattern (dB)');
title('Phase Modulation (PM) — your baseline'); ylim([-50 15]); xlim([-90 90]); grid on;
legend(arrayfun(@(m) sprintf('m=%+d',m), m_idx, 'UniformOutput',false), 'Location','eastoutside');

% --- quantify the PM gain over AM at each harmonic peak, as the paper reports ---
peak_AM_dB = max(F_AM_dB, [], 2);
peak_PM_dB = max(F_PM_dB, [], 2);
gain_dB = peak_PM_dB - peak_AM_dB;

fprintf('=== PM vs AM peak gain per harmonic (paper reports +15.56 dB at fc) ===\n');
fprintf('%6s %14s\n','m','PM-AM gain(dB)');
for i = 1:length(m_idx)
    fprintf('%6d %14.2f\n', m_idx(i), gain_dB(i));
end

%% ---------------------- LOCAL FUNCTIONS ----------------------
function a_m = compute_am(seq, L, Gamma0, Gamma1, m_idx)
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
end

function a_qm = apply_time_shift(a_m, m_idx, f0, t_q)
    Nq = length(t_q);
    a_qm = zeros(Nq, length(m_idx));
    for q = 1:Nq
        phase_shift = exp(-1j*2*pi*m_idx*f0*t_q(q));
        a_qm(q,:) = a_m .* phase_shift;
    end
end

function F = array_farfield(a_qm, m_idx, Nq, dy_over_lambda, theta)
    F = zeros(length(m_idx), length(theta));
    for mi = 1:length(m_idx)
        for ti = 1:length(theta)
            array_phase = exp(1j*2*pi*(0:Nq-1)*dy_over_lambda*sin(theta(ti)));
            F(mi,ti) = sum(array_phase .* a_qm(:,mi).');
        end
    end
end