%% ========================================================================
%  PHASE 2 — Composite Simulation: All 4 Non-Idealities Simultaneously
%  v2: Architecture faithful to individual experiment scripts.
%
%  Key fix: keeps temporal (a_m from waveform) and spatial (apply_time_shift
%  via t_q) strictly separated, exactly as in jitter.m / risefalltime.m /
%  temporalmismatch.m / waveform_distortion.m.
%
%  Rise/fall time  → applied as sinc filter on a_m  (from risefalltime.m)
%  Droop           → applied via numerical Fourier on build_droop waveform
%                    (from waveform_distortion.m)
%  Jitter          → perturbs t_q per element       (from jitter.m)
%  Phase mismatch  → perturbs Gamma0/Gamma1 phases  (from temporalmismatch.m)
%  Amplitude mism. → perturbs Gamma magnitudes      (from temporalmismatch.m)
%
%  Harmonic efficiency computed from single-element a_m BEFORE spatial
%  shift, so it matches waveform_distortion.m output exactly.
% =========================================================================
clear; clc; close all; rng(42);

%% ─── 1. SHARED BASELINE (identical to all 4 experiment scripts) ─────────
Gamma0_nom = 0.911*exp(1j*deg2rad(-102));
Gamma1_nom = 0.952*exp(1j*deg2rad(94));
seq        = [1 0 0 0 0 0 0 0];
L          = length(seq);
f0         = 1e6;  T0 = 1/f0;
Nq         = 8;    dy_over_lambda = 0.5;
Mmax       = 3;    m_idx = -Mmax:Mmax;
t_q_ideal  = (0:Nq-1)*T0/Nq;          % ideal spatial time offsets
theta_deg  = -90:0.02:90;
theta      = deg2rad(theta_deg);
Nt         = 8192;                     % time samples (from waveform_distortion.m)
t_fine     = linspace(0,T0,Nt+1); t_fine = t_fine(1:end-1);
N_mc       = 500;

%% ─── 2. IDEAL BASELINE ──────────────────────────────────────────────────
% Waveform-level Fourier coefficients (no spatial shift yet)
a_m_ideal = compute_am(seq, L, Gamma0_nom, Gamma1_nom, m_idx);

% Harmonic efficiency from single-element waveform
P_ideal    = abs(a_m_ideal).^2;
eta_ideal  = P_ideal / sum(P_ideal);

% Spatial shift applied: a_qm(q,:) = a_m .* exp(-j*2*pi*m*f0*t_q(q))
a_qm_ideal = apply_time_shift(a_m_ideal, m_idx, f0, t_q_ideal);
F_ideal    = array_farfield(a_qm_ideal, m_idx, Nq, dy_over_lambda, theta);
[pa_ideal, ~] = analyze_pattern(F_ideal, theta_deg, m_idx);

fprintf('=== Ideal baseline ===\n');
fprintf('%-8s %-14s %-10s\n','m','Steer (deg)','eta_m');
for i = 1:length(m_idx)
    fprintf('%-8d %-14.3f %-10.4f\n', m_idx(i), pa_ideal(i), eta_ideal(i));
end
fprintf('\nExpected: steering = arcsin(m/4), eta_0=0.6490, eta_+1=0.0702\n\n');

%% ─── 3. NON-IDEALITY PARAMETER SETS ─────────────────────────────────────
% Physical operating point — Ch 6.6, Zhang et al. hardware
phys.tr_frac   = 0.010;   % MACOM MADP-000907: 2-3 ns / 250 ns slot
phys.sigma_j   = 0.0002;  % Altera Cyclone IV PLL: ~50 ps / 250 ns slot
phys.sigma_phi = 4;        % deg — F4B PCB standard fabrication
phys.sigma_A   = 0.02;    % 2% amplitude variation
phys.alpha     = 3.4e-5;  % droop: negligible at 0.5 MHz

% Worst-case — outer edges of individual experiment sweep ranges
worst.tr_frac   = 0.30;   % risefalltime.m: tr_fracs(end)
worst.sigma_j   = 0.20;   % jitter.m: sigma_fracs(end)
worst.sigma_phi = 20;      % temporalmismatch.m: sigma_phi_vals(end) [deg]
worst.sigma_A   = 0.20;   % temporalmismatch.m: sigma_A_vals(end)
worst.alpha     = 0.50;   % waveform_distortion.m: alpha_droop(end)


%% ─── 4. MONTE CARLO FUNCTION ─────────────────────────────────────────────
function [pe_rms, pe_mean, overlap_mean, overlap_p5, eta_mean, eta_p5] = ...
        run_mc(p, a_m_ideal, m_idx, f0, T0, L, Nq, t_q_ideal, ...
               Gamma0_nom, Gamma1_nom, seq, t_fine, dy_over_lambda, ...
               theta, theta_deg, F_ideal, pa_ideal, N_mc)

    pe_mc      = zeros(N_mc, length(m_idx));
    overlap_mc = zeros(N_mc, length(m_idx));
    eta_mc     = zeros(N_mc, length(m_idx));   % per-trial, per-harmonic

    for mc = 1:N_mc
        % ── Draw stochastic realisations ──────────────────────────────────
        jitter_q = p.sigma_j   * randn(1, Nq);   % normalised to T0/L
        dphi_q   = p.sigma_phi * randn(1, Nq);   % common phase offset (deg)
        dA_q     = p.sigma_A   * randn(1, Nq);   % fractional amplitude offset

        a_qm_pert = zeros(Nq, length(m_idx));
        a_m_trial = zeros(Nq, length(m_idx));     % store per-element a_m

        for q = 1:Nq
            % ── Step 1: perturbed Gamma values (temporalmismatch.m style) ─
            G0_q = abs(Gamma0_nom)*(1+dA_q(q)) * ...
                   exp(1j*(angle(Gamma0_nom) + deg2rad(dphi_q(q))));
            G1_q = abs(Gamma1_nom)*(1+dA_q(q)) * ...
                   exp(1j*(angle(Gamma1_nom) + deg2rad(dphi_q(q))));

            % ── Step 2: waveform-level a_m ────────────────────────────────
            if p.alpha > 1e-6
                % Droop: numerical Fourier (waveform_distortion.m style)
                Gamma_t_q = build_droop(t_fine, T0, seq, G0_q, G1_q, p.alpha);
                a_m_q     = fourier_coeffs(Gamma_t_q, t_fine, T0, m_idx);
            else
                % No droop: analytical (compute_am style)
                a_m_q = compute_am(seq, L, G0_q, G1_q, m_idx);
            end

            % ── Step 3: rise/fall time as sinc filter (risefalltime.m) ───
            if p.tr_frac > 0
                a_m_q = a_m_q .* sinc(m_idx * f0 * p.tr_frac * (T0/L));
            end

            % ── Step 4: jitter perturbs spatial time offset (jitter.m) ───
            t_q_pert = t_q_ideal(q) + jitter_q(q) * (T0/L);

            % ── Step 5: apply spatial phase gradient ──────────────────────
            a_qm_pert(q,:) = a_m_q .* exp(-1j*2*pi*m_idx*f0*t_q_pert);
            a_m_trial(q,:) = a_m_q;   % save pre-shift for efficiency
        end

        % ── Far-field and metrics ─────────────────────────────────────────
        F_pert = array_farfield(a_qm_pert, m_idx, Nq, dy_over_lambda, theta);
        [pa_pert, ~] = analyze_pattern(F_pert, theta_deg, m_idx);

        for mi = 1:length(m_idx)
            Fi = F_ideal(mi,:);  Fp = F_pert(mi,:);
            pe_mc(mc,mi)      = pa_pert(mi) - pa_ideal(mi);
            overlap_mc(mc,mi) = abs(dot(Fi,Fp))^2 / ...
                                 (sum(abs(Fi).^2) * sum(abs(Fp).^2));
        end

        % Efficiency: mean power per harmonic across elements, then normalise.
        % Average |a_m_q|^2 across elements (power average, not amplitude
        % average) — matches waveform_distortion.m: P_n = abs(a_m).^2
        P_trial = mean(abs(a_m_trial).^2, 1);   % 1 x length(m_idx)
        eta_mc(mc,:) = P_trial / sum(P_trial);
    end

    pe_rms      = std(pe_mc, 0, 1);
    pe_mean     = mean(pe_mc, 1);
    overlap_mean = mean(overlap_mc, 1);
    overlap_p5  = prctile(overlap_mc, 5, 1);
    eta_mean    = mean(eta_mc, 1);
    eta_p5      = prctile(eta_mc, 5, 1);
end

%% ─── 5. SINGLE-POINT EVALUATION ─────────────────────────────────────────
param_sets = {phys, worst};
labels     = {'Physical operating point (Ch 6.6)', 'Worst-case combined'};
res        = cell(2,1);

for ps = 1:2
    fprintf('Running MC: %s ...\n', labels{ps});
    [pe_rms, pe_mean, ov_mean, ov_p5, eta_mean, eta_p5] = run_mc( ...
        param_sets{ps}, a_m_ideal, m_idx, f0, T0, L, Nq, t_q_ideal, ...
        Gamma0_nom, Gamma1_nom, seq, t_fine, dy_over_lambda, ...
        theta, theta_deg, F_ideal, pa_ideal, N_mc);

    res{ps}.pe_rms   = pe_rms;
    res{ps}.pe_mean  = pe_mean;
    res{ps}.ov_mean  = ov_mean;
    res{ps}.ov_p5    = ov_p5;
    res{ps}.eta_mean = eta_mean;
    res{ps}.eta_p5   = eta_p5;

    fprintf('\n--- %s ---\n', labels{ps});
    fprintf('%-8s %-14s %-14s %-14s %-10s\n', ...
            'm','PE RMS (deg)','PE mean (deg)','Overlap mean','Overlap 5pct');
    for i = 1:length(m_idx)
        fprintf('%-8d %-14.4f %-14.4f %-14.4f %-10.4f\n', ...
            m_idx(i), pe_rms(i), pe_mean(i), ov_mean(i), ov_p5(i));
    end
    fprintf('\nHarmonic efficiency (mean across MC trials):\n');
    fprintf('%-8s %-12s %-12s %-10s\n','m','eta_mean','eta_5pct','eta_ideal');
    for i = 1:length(m_idx)
        fprintf('%-8d %-12.4f %-12.4f %-10.4f\n', ...
            m_idx(i), eta_mean(i), eta_p5(i), eta_ideal(i));
    end
    fprintf('\n');
end

%% ─── 6. SEVERITY SWEEP ───────────────────────────────────────────────────
fprintf('Running composite severity sweep (this takes a few minutes)...\n');
scale_vals = [0, 0.05, 0.10, 0.20, 0.30, 0.50, 0.70, 1.0];
n_sc  = length(scale_vals);
m1i   = find(m_idx==1);
m3i   = find(m_idx==3);

ov_sweep   = zeros(n_sc, length(m_idx));
pe_sweep   = zeros(n_sc, length(m_idx));
eta_sweep  = zeros(n_sc, length(m_idx));

for si = 1:n_sc
    s = scale_vals(si);
    p_sc.tr_frac   = s * worst.tr_frac;
    p_sc.sigma_j   = s * worst.sigma_j;
    p_sc.sigma_phi = s * worst.sigma_phi;
    p_sc.sigma_A   = s * worst.sigma_A;
    p_sc.alpha     = s * worst.alpha;

    [pe_r, ~, ov_m, ~, eta_m, ~] = run_mc( ...
        p_sc, a_m_ideal, m_idx, f0, T0, L, Nq, t_q_ideal, ...
        Gamma0_nom, Gamma1_nom, seq, t_fine, dy_over_lambda, ...
        theta, theta_deg, F_ideal, pa_ideal, N_mc);

    ov_sweep(si,:)  = ov_m;
    pe_sweep(si,:)  = pe_r;
    eta_sweep(si,:) = eta_m;
    fprintf('  scale=%.2f done\n', s);
end

%% ─── 7. FIGURES ─────────────────────────────────────────────────────────

% Fig 1: Pattern overlap vs severity, all harmonics
figure('Name','Composite: Pattern overlap vs severity');
plot(scale_vals, ov_sweep, '-o', 'LineWidth', 1.4);
xlabel('Composite severity scale  (0 = ideal,  1 = worst-case limits)');
ylabel('Mean pattern overlap  O_m');
title('Composite non-ideality: Pattern overlap across harmonics');
legend(arrayfun(@(m) sprintf('m=%+d',m), m_idx, 'UniformOutput',false), ...
    'Location','southwest');
h1 = yline(0.95,'r--','LineWidth',1.2); h1.HandleVisibility='off';
h2 = yline(0.99,'b:','LineWidth',1.2); h2.HandleVisibility='off';
text(0.01,0.951,'0.95 threshold','Color','r','FontSize',9);
text(0.01,0.991,'0.99 threshold','Color','b','FontSize',9);
ylim([0.85 1.005]); xlim([0 1]); grid on;

% Fig 2: PE RMS vs severity, m=+1 and m=+3
figure('Name','Composite: Pointing error vs severity');
plot(scale_vals, pe_sweep(:,m1i), '-o','LineWidth',1.4,'DisplayName','m=+1'); hold on;
plot(scale_vals, pe_sweep(:,m3i), '-s','LineWidth',1.4,'DisplayName','m=+3');
xlabel('Composite severity scale  (0 = ideal,  1 = worst-case limits)');
ylabel('RMS pointing error (deg)');
title('Composite non-ideality: RMS beam-pointing error');
legend('Location','northwest'); grid on; xlim([0 1]);

% Fig 3: Harmonic efficiency vs severity (carrier + sideband)
m0i = find(m_idx==0);
figure('Name','Composite: Harmonic efficiency vs severity');
plot(scale_vals, eta_sweep(:,m0i), '-o','LineWidth',1.4,'DisplayName','m=0 (carrier)'); hold on;
plot(scale_vals, eta_sweep(:,m1i), '-s','LineWidth',1.4,'DisplayName','m=+1 (sideband)');
h3 = yline(eta_ideal(m0i),'b:','LineWidth',1.2); h3.HandleVisibility='off';
h4 = yline(eta_ideal(m1i),'r:','LineWidth',1.2); h4.HandleVisibility='off';
text(0.82, eta_ideal(m0i)+0.01, 'Ideal \eta_0',  'Color','b','FontSize',9);
text(0.82, eta_ideal(m1i)+0.01, 'Ideal \eta_{+1}','Color','r','FontSize',9);
xlabel('Composite severity scale');
ylabel('\eta_m');
title('Composite non-ideality: Harmonic efficiency');
legend('Location','east'); grid on; xlim([0 1]);

% Fig 4: Worst-case far-field spread at m=+1, 20 realisations
figure('Name','Composite: Far-field spread m=+1, worst-case');
hold on;
for k = 1:20
    jq = worst.sigma_j   * randn(1,Nq);
    dq = worst.sigma_phi * randn(1,Nq);
    aq = worst.sigma_A   * randn(1,Nq);
    a_qm_k = zeros(Nq, length(m_idx));
    for q = 1:Nq
        G0q = abs(Gamma0_nom)*(1+aq(q))*exp(1j*(angle(Gamma0_nom)+deg2rad(dq(q))));
        G1q = abs(Gamma1_nom)*(1+aq(q))*exp(1j*(angle(Gamma1_nom)+deg2rad(dq(q))));
        if worst.alpha > 1e-6
            Gt  = build_droop(t_fine,T0,seq,G0q,G1q,worst.alpha);
            amq = fourier_coeffs(Gt,t_fine,T0,m_idx);
        else
            amq = compute_am(seq,L,G0q,G1q,m_idx);
        end
        amq = amq .* sinc(m_idx*f0*worst.tr_frac*(T0/L));
        tqp = t_q_ideal(q) + jq(q)*(T0/L);
        a_qm_k(q,:) = amq .* exp(-1j*2*pi*m_idx*f0*tqp);
    end
    Fk = array_farfield(a_qm_k, m_idx, Nq, dy_over_lambda, theta);
    plot(theta_deg, 20*log10(abs(Fk(m1i,:))/max(abs(F_ideal(m1i,:)))), ...
        'Color',[0.6 0.6 0.6 0.4],'LineWidth',0.8,'HandleVisibility','off');
end
plot(theta_deg, 20*log10(abs(F_ideal(m1i,:))/max(abs(F_ideal(m1i,:)))), ...
    'b','LineWidth',2,'DisplayName','Ideal (m=+1)');
xlabel('\theta (deg)'); ylabel('Normalised power (dB)');
title('Composite worst-case: m=+1 far-field spread (20 realisations)');
legend; ylim([-40 5]); xlim([-90 90]); grid on;

%% ─── 8. SUMMARY TABLE ───────────────────────────────────────────────────
fprintf('\n╔══════════════════════════════════════════════════════════════╗\n');
fprintf('║         COMPOSITE NON-IDEALITY: KEY RESULTS SUMMARY         ║\n');
fprintf('╠══════════════════════════════════════════════════════════════╣\n');
fprintf('║ Physical operating point (Zhang et al. hardware)            ║\n');
fprintf('║   Overlap (mean):  %.4f  (5th pct: %.4f)               ║\n', ...
    mean(res{1}.ov_mean), mean(res{1}.ov_p5));
fprintf('║   PE RMS m=+1:     %.4f°   m=+3: %.4f°                ║\n', ...
    res{1}.pe_rms(m1i), res{1}.pe_rms(m3i));
fprintf('║   eta_0:           %.4f    eta_+1: %.4f (ideal: %.4f)║\n', ...
    res{1}.eta_mean(m0i), res{1}.eta_mean(m1i), eta_ideal(m1i));
fprintf('╠══════════════════════════════════════════════════════════════╣\n');
fprintf('║ Worst-case combined (all 4 non-idealities at sweep limits)  ║\n');
fprintf('║   Overlap (mean):  %.4f  (5th pct: %.4f)               ║\n', ...
    mean(res{2}.ov_mean), mean(res{2}.ov_p5));
fprintf('║   PE RMS m=+1:     %.4f°   m=+3: %.4f°                ║\n', ...
    res{2}.pe_rms(m1i), res{2}.pe_rms(m3i));
fprintf('║   eta_0:           %.4f    eta_+1: %.4f (ideal: %.4f)║\n', ...
    res{2}.eta_mean(m0i), res{2}.eta_mean(m1i), eta_ideal(m1i));
fprintf('╠══════════════════════════════════════════════════════════════╣\n');
fprintf('║ Overlap drops below 0.95 at composite severity scale ≈      ║\n');
% find where overlap at m=+3 crosses 0.95
m3_ov = ov_sweep(:, m3i);
cross_idx = find(m3_ov < 0.95, 1, 'first');
if ~isempty(cross_idx) && cross_idx > 1
    s_cross = interp1(m3_ov(cross_idx-1:cross_idx), ...
                      scale_vals(cross_idx-1:cross_idx)', 0.95);
    fprintf('║   s ≈ %.2f  (m=±3, most sensitive harmonic)              ║\n', s_cross);
else
    fprintf('║   Not crossed within sweep range                            ║\n');
end
fprintf('╚══════════════════════════════════════════════════════════════╝\n');

%% ─── LOCAL FUNCTIONS (faithful copies from individual scripts) ───────────

function a_m = compute_am(seq, L, Gamma0, Gamma1, m_idx)
    % Analytical Fourier coefficients — from jitter.m / risefalltime.m
    GammaSeq = zeros(1,L);
    for k = 1:L
        GammaSeq(k) = (seq(k)==0)*Gamma0 + (seq(k)==1)*Gamma1;
    end
    a_m = zeros(size(m_idx));
    for idx = 1:length(m_idx)
        m = m_idx(idx); s = 0;
        for k = 1:L
            s = s + (GammaSeq(k)/L)*sinc(m/L)*exp(-1j*pi*m*(2*k-1)/L);
        end
        a_m(idx) = s;
    end
end

function Gamma_t = build_droop(t, T0, seq, G0, G1, alpha)
    % Linear droop during active slots — from waveform_distortion.m
    L = length(seq); tau = T0/L; Gamma_t = zeros(size(t));
    for i = 1:length(t)
        tt   = mod(t(i),T0);
        k    = min(floor(tt/tau)+1, L);
        t_in = tt - (k-1)*tau;
        G    = (seq(k)==0)*G0 + (seq(k)==1)*G1;
        if seq(k)==1
            G = G + alpha*(t_in/tau)*(G0-G1);
        end
        Gamma_t(i) = G;
    end
end

function a_m = fourier_coeffs(Gamma_t, t, T0, m_idx)
    % Numerical DFT — from waveform_distortion.m
    dt = t(2)-t(1); a_m = zeros(size(m_idx));
    for idx = 1:length(m_idx)
        a_m(idx) = sum(Gamma_t .* exp(-1j*2*pi*m_idx(idx)*t/T0))*dt/T0;
    end
end

function a_qm = apply_time_shift(a_m, m_idx, f0, t_q)
    % Spatial phase gradient — from all 4 experiment scripts
    Nq = length(t_q); a_qm = zeros(Nq, length(m_idx));
    for q = 1:Nq
        a_qm(q,:) = a_m .* exp(-1j*2*pi*m_idx*f0*t_q(q));
    end
end

function F = array_farfield(a_qm, m_idx, Nq, dy_over_lambda, theta)
    % Array factor — from all 4 experiment scripts
    F = zeros(length(m_idx), length(theta));
    for mi = 1:length(m_idx)
        for ti = 1:length(theta)
            ap = exp(1j*2*pi*(0:Nq-1)*dy_over_lambda*sin(theta(ti)));
            F(mi,ti) = sum(ap .* a_qm(:,mi).');
        end
    end
end

function [peak_angle, peak_sidelobe_dB] = analyze_pattern(F, theta_deg, m_idx)
    % Pattern analysis — from all 4 experiment scripts
    peak_angle = zeros(size(m_idx)); peak_sidelobe_dB = zeros(size(m_idx));
    for mi = 1:length(m_idx)
        mag = abs(F(mi,:)); [main_val, main_idx] = max(mag);
        peak_angle(mi) = theta_deg(main_idx);
        win  = round(3/mean(diff(theta_deg))); mask = true(size(mag));
        mask(max(1,main_idx-win):min(length(mag),main_idx+win)) = false;
        if any(mask)
            peak_sidelobe_dB(mi) = 20*log10(max(mag(mask))/main_val);
        else
            peak_sidelobe_dB(mi) = -Inf;
        end
    end
end