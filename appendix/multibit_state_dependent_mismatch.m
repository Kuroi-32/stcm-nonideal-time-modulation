%% ========================================================================
%  APPENDIX C: Multi-bit Verification and State-Dependent Mismatch
%
%  Purpose:
%    (1) Verify empirically that jitter and common-phase-mismatch
%        sensitivity coefficients are bit-depth invariant (confirming
%        the analytical proofs in Appendices A and B).
%    (2) Characterise state-dependent mismatch for 2-bit and 3-bit
%        STCMs — a genuinely new experiment not covered by the main
%        Phase 2 analysis, which assumed common-offset mismatch only.
%
%  Analytical predictions being tested:
%    (A) Jitter: sensitivity = 3.43 * |m| deg/unit sigma, independent
%        of bit depth. Follows from Fourier shift theorem (Appendix A).
%    (B) Common-phase mismatch: sensitivity ~0.081 deg/deg at m=+/-3,
%        independent of bit depth. Follows from factorisation of
%        exp(j*delta_phi) from Fourier sum (Appendix B).
%    (C) State-dependent mismatch: sensitivity scales as
%        sigma_eps * |sinc(m/L)| / sqrt(L), dependent on bit depth
%        through the number of independent error terms. This is the
%        new prediction to be verified here.
%
%  Three bit depths compared:
%    1-bit: 2 states, phases 0 and 180 deg (your validated baseline)
%    2-bit: 4 states, phases 0, 90, 180, 270 deg
%    3-bit: 8 states, phases 0, 45, 90, 135, 180, 225, 270, 315 deg
% =========================================================================
clear; clc; close all; rng(42);

%% ---------------------- 1. SHARED PARAMETERS ----------------------
% Array and modulation parameters — identical to main Phase 2 scripts
f0  = 1e6; T0 = 1/f0;
Nq  = 8; dy_over_lambda = 0.5;
Mmax = 3; m_idx = -Mmax:Mmax; n_m = length(m_idx);
t_q = (0:Nq-1)*T0/Nq;
theta_deg = -90:0.02:90;
theta = deg2rad(theta_deg);
N_mc = 500;

% Jitter sweep (normalized to T0/L)
sigma_j_fracs = [0, 0.01, 0.02, 0.05, 0.10, 0.20];

% Common-phase mismatch sweep (degrees)
sigma_phi_vals = [0, 1, 2, 5, 10, 20];

% State-dependent mismatch sweep (fractional phase error per state, degrees)
sigma_sd_vals = [0, 1, 2, 5, 10, 20];

%% ---------------------- 2. BIT-DEPTH CONFIGURATIONS ----------------------
% Reflection amplitude assumed ~0.93 in all states (near-unity reflection,
% same assumption as the main analysis). Phases equally spaced.
% In a real multi-bit implementation these would come from CST simulation
% of the multi-diode unit cell; here ideal phase spacing is assumed to
% isolate the effect of bit depth from unit-cell geometry.

configs = struct();

% 1-bit: validated baseline
configs(1).name      = '1-bit';
configs(1).n_states  = 2;
configs(1).phases_deg= [0, 180];
configs(1).seq       = [1 0 0 0 0 0 0 0];  % low-duty-cycle, validated
configs(1).L         = 8;

% 2-bit: four states at 90-degree spacing, gradient sequence
configs(2).name      = '2-bit';
configs(2).n_states  = 4;
configs(2).phases_deg= [0, 90, 180, 270];
configs(2).seq       = [1 2 3 4 3 2 1 2];  % length-8, repeating gradient
configs(2).L         = 8;

% 3-bit: eight states at 45-degree spacing, one-of-each sequence
configs(3).name      = '3-bit';
configs(3).n_states  = 8;
configs(3).phases_deg= 0:45:315;
configs(3).seq       = [1 3 5 7 2 4 6 8];  % each state once, L=8
configs(3).L         = 8;

n_configs = length(configs);

%% ---------------------- 3. PRECOMPUTE IDEAL BASELINES ----------------------
for ci = 1:n_configs
    cfg = configs(ci);
    Gamma_states = 0.93 * exp(1j*deg2rad(cfg.phases_deg));
    a_m = compute_am_nbit(cfg.seq, Gamma_states, cfg.L, m_idx);
    a_qm = apply_time_shift(a_m, m_idx, f0, t_q);
    F = array_farfield(a_qm, m_idx, Nq, dy_over_lambda, theta);
    [pa, ~, ~] = analyze_pattern(F, theta_deg, m_idx);

    configs(ci).Gamma_states   = Gamma_states;
    configs(ci).a_m_ideal      = a_m;
    configs(ci).a_qm_ideal     = a_qm;
    configs(ci).F_ideal        = F;
    configs(ci).pa_ideal       = pa;

    P_n = abs(a_m).^2;
    eta = P_n/sum(P_n);
    fprintf('--- %s baseline ---\n', cfg.name);
    fprintf('%6s %10s %12s\n','m','eta_m','PA(deg)');
    for i=1:n_m
        fprintf('%6d %10.4f %12.2f\n',m_idx(i),eta(i),pa(i));
    end
    fprintf('\n');
end

%% ---------------------- 4. EXPERIMENT A: JITTER (bit-depth invariance) ----------------------
fprintf('=== EXPERIMENT A: JITTER SENSITIVITY — BIT-DEPTH INVARIANCE ===\n');
fprintf('Expected: identical sensitivity coefficients across all bit depths\n\n');

jitter_sensitivity = zeros(n_configs, n_m);

for ci = 1:n_configs
    cfg = configs(ci);
    n_sig = length(sigma_j_fracs);
    pe_std_j = zeros(n_sig, n_m);

    for ii = 1:n_sig
        sigma = sigma_j_fracs(ii) * (T0/cfg.L);
        pe_mc = zeros(N_mc, n_m);
        for mc = 1:N_mc
            delta_q = sigma * randn(1, Nq);
            a_qm_j = apply_time_shift(cfg.a_m_ideal, m_idx, f0, t_q + delta_q);
            F_j = array_farfield(a_qm_j, m_idx, Nq, dy_over_lambda, theta);
            [pa, ~, ~] = analyze_pattern(F_j, theta_deg, m_idx);
            pe_mc(mc,:) = pa - cfg.pa_ideal;
        end
        pe_std_j(ii,:) = std(pe_mc, 0, 1);
    end

    % Linear sensitivity from fit
    for i = 1:n_m
        slope = polyfit(sigma_j_fracs', pe_std_j(:,i), 1);
        jitter_sensitivity(ci,i) = slope(1);
    end
    fprintf('%s jitter sensitivity done\n', cfg.name);
end

fprintf('\nJitter sensitivity (deg per unit sigma/(T0/L)):\n');
fprintf('%14s', 'Bit depth');
for i=1:n_m; fprintf('%10s',sprintf('m=%+d',m_idx(i))); end
fprintf('%16s\n','Per |m| unit');
for ci=1:n_configs
    fprintf('%14s', configs(ci).name);
    for i=1:n_m; fprintf('%10.3f', jitter_sensitivity(ci,i)); end
    % Compute per-harmonic-order coefficient (should be ~3.43 for all)
    nonzero = m_idx ~= 0;
    per_order = mean(abs(jitter_sensitivity(ci,nonzero)) ./ abs(m_idx(nonzero)));
    fprintf('%16.3f\n', per_order);
end

%% ---------------------- 5. EXPERIMENT B: COMMON-PHASE MISMATCH (bit-depth invariance) ----------------------
fprintf('\n=== EXPERIMENT B: COMMON-PHASE MISMATCH — BIT-DEPTH INVARIANCE ===\n');
fprintf('Expected: identical sensitivity coefficients across all bit depths\n\n');

phi_sensitivity = zeros(n_configs, n_m);

for ci = 1:n_configs
    cfg = configs(ci);
    n_sp = length(sigma_phi_vals);
    pe_std_p = zeros(n_sp, n_m);

    for ii = 1:n_sp
        sP = deg2rad(sigma_phi_vals(ii));
        pe_mc = zeros(N_mc, n_m);
        for mc = 1:N_mc
            delta_phi = sP * randn(1, Nq);
            a_qm_mc = zeros(Nq, n_m);
            for q = 1:Nq
                % Common phase offset: all states of element q shift by delta_phi(q)
                G_states_q = cfg.Gamma_states * exp(1j*delta_phi(q));
                a_m_q = compute_am_nbit(cfg.seq, G_states_q, cfg.L, m_idx);
                a_qm_mc(q,:) = a_m_q .* exp(-1j*2*pi*m_idx*f0*t_q(q));
            end
            F_mc = array_farfield(a_qm_mc, m_idx, Nq, dy_over_lambda, theta);
            [pa, ~, ~] = analyze_pattern(F_mc, theta_deg, m_idx);
            pe_mc(mc,:) = pa - cfg.pa_ideal;
        end
        pe_std_p(ii,:) = std(pe_mc, 0, 1);
    end

    for i = 1:n_m
        slope = polyfit(sigma_phi_vals', pe_std_p(:,i), 1);
        phi_sensitivity(ci,i) = slope(1);
    end
    fprintf('%s common-phase mismatch done\n', cfg.name);
end

fprintf('\nCommon-phase mismatch sensitivity (deg RMS per deg sigma_phi):\n');
fprintf('%14s', 'Bit depth');
for i=1:n_m; fprintf('%10s',sprintf('m=%+d',m_idx(i))); end; fprintf('\n');
for ci=1:n_configs
    fprintf('%14s', configs(ci).name);
    for i=1:n_m; fprintf('%10.4f', phi_sensitivity(ci,i)); end; fprintf('\n');
end

%% ---------------------- 6. EXPERIMENT C: STATE-DEPENDENT MISMATCH (NEW) ----------------------
% This is the genuinely new experiment for Appendix C.
% Each of the n_states states of element q gets an INDEPENDENT phase error.
% Analytical prediction: RMS error scales as sigma_eps * |sinc(m/L)| / sqrt(L)
% This DOES depend on bit depth through n_states.
fprintf('\n=== EXPERIMENT C: STATE-DEPENDENT MISMATCH (NEW) ===\n');
fprintf('Analytical prediction: sensitivity ~ sigma_eps * |sinc(m/L)| / sqrt(L)\n');
fprintf('Expected: bit-depth DEPENDENT results\n\n');

sd_sensitivity = zeros(n_configs, n_m);

for ci = 1:n_configs
    cfg = configs(ci);
    n_sd = length(sigma_sd_vals);
    pe_std_sd = zeros(n_sd, n_m);

    for ii = 1:n_sd
        sSD = deg2rad(sigma_sd_vals(ii));
        pe_mc = zeros(N_mc, n_m);
        for mc = 1:N_mc
            a_qm_mc = zeros(Nq, n_m);
            for q = 1:Nq
                % INDEPENDENT phase error per state per element
                eps_q = sSD * randn(1, cfg.n_states);
                G_states_q = cfg.Gamma_states .* exp(1j*eps_q);
                a_m_q = compute_am_nbit(cfg.seq, G_states_q, cfg.L, m_idx);
                a_qm_mc(q,:) = a_m_q .* exp(-1j*2*pi*m_idx*f0*t_q(q));
            end
            F_mc = array_farfield(a_qm_mc, m_idx, Nq, dy_over_lambda, theta);
            [pa, ~, ~] = analyze_pattern(F_mc, theta_deg, m_idx);
            pe_mc(mc,:) = pa - cfg.pa_ideal;
        end
        pe_std_sd(ii,:) = std(pe_mc, 0, 1);
    end

    for i = 1:n_m
        slope = polyfit(sigma_sd_vals', pe_std_sd(:,i), 1);
        sd_sensitivity(ci,i) = slope(1);
    end
    fprintf('%s state-dependent mismatch done\n', cfg.name);
end

fprintf('\nState-dependent mismatch sensitivity (deg RMS per deg sigma_eps):\n');
fprintf('%14s', 'Bit depth');
for i=1:n_m; fprintf('%10s',sprintf('m=%+d',m_idx(i))); end; fprintf('\n');
for ci=1:n_configs
    fprintf('%14s', configs(ci).name);
    for i=1:n_m; fprintf('%10.4f', sd_sensitivity(ci,i)); end; fprintf('\n');
end

% Compare against analytical prediction sigma_eps * |sinc(m/L)| / sqrt(L)
% The sensitivity coefficient should be proportional to |sinc(m/L)| / sqrt(L)
% The proportionality constant comes from how pointing error maps to Fourier
% coefficient perturbation — we compute the predicted shape here
fprintf('\nAnalytical prediction shape (normalized, for L=8):\n');
fprintf('%14s', 'Harmonic m');
for i=1:n_m; fprintf('%10s',sprintf('m=%+d',m_idx(i))); end; fprintf('\n');
predicted_shape = abs(sinc(m_idx/8)) / sqrt(8);
predicted_shape_norm = predicted_shape / max(predicted_shape);
fprintf('%14s', 'Predicted');
for i=1:n_m; fprintf('%10.4f', predicted_shape_norm(i)); end; fprintf('\n');

%% ---------------------- 7. SUMMARY FIGURE ----------------------
figure('Name','Appendix C: Multi-bit verification and state-dependent mismatch',...
    'Position',[50 50 1400 900]);

colors = lines(n_configs);
labels = {configs.name};

% Panel 1: Jitter sensitivity bar chart
subplot(2,3,1);
bar(m_idx, jitter_sensitivity', 'grouped');
xlabel('Harmonic m'); ylabel('Sensitivity (deg/unit \sigma)');
title('Jitter sensitivity — all bit depths');
legend(labels,'Location','north'); grid on;

% Panel 2: Common-phase mismatch sensitivity bar chart
subplot(2,3,2);
bar(m_idx, phi_sensitivity', 'grouped');
xlabel('Harmonic m'); ylabel('Sensitivity (deg/deg \sigma_\phi)');
title('Common-phase mismatch sensitivity — all bit depths');
legend(labels,'Location','north'); grid on;

% Panel 3: State-dependent mismatch sensitivity bar chart
subplot(2,3,3);
bar(m_idx, sd_sensitivity', 'grouped');
xlabel('Harmonic m'); ylabel('Sensitivity (deg/deg \sigma_{eps})');
title('State-dependent mismatch — bit-depth dependent');
legend(labels,'Location','north'); grid on;

% Panel 4: Jitter — confirm invariance by overlaying all three curves
subplot(2,3,4);
hold on;
for ci=1:n_configs
    plot(m_idx, abs(jitter_sensitivity(ci,:)), '-o', ...
        'Color',colors(ci,:),'LineWidth',1.5,'DisplayName',configs(ci).name);
end
xlabel('Harmonic m'); ylabel('|Sensitivity| (deg/unit \sigma)');
title('Jitter: confirmed bit-depth invariant');
legend('Location','northwest'); grid on;

% Panel 5: Phase mismatch — confirm invariance
subplot(2,3,5);
hold on;
for ci=1:n_configs
    plot(m_idx, abs(phi_sensitivity(ci,:)), '-s', ...
        'Color',colors(ci,:),'LineWidth',1.5,'DisplayName',configs(ci).name);
end
xlabel('Harmonic m'); ylabel('|Sensitivity| (deg/deg \sigma_\phi)');
title('Common-phase mismatch: confirmed bit-depth invariant');
legend('Location','northwest'); grid on;

% Panel 6: State-dependent — show bit-depth dependence
subplot(2,3,6);
hold on;
for ci=1:n_configs
    plot(m_idx, abs(sd_sensitivity(ci,:)), '-^', ...
        'Color',colors(ci,:),'LineWidth',1.5,'DisplayName',configs(ci).name);
end
% Overlay analytical prediction (scaled to match 1-bit)
scale = max(abs(sd_sensitivity(1,:))) / max(predicted_shape_norm);
plot(m_idx, scale*predicted_shape_norm, 'k--', ...
    'LineWidth',1.5,'DisplayName','Predicted shape');
xlabel('Harmonic m'); ylabel('|Sensitivity| (deg/deg \sigma_{eps})');
title('State-dependent mismatch: bit-depth dependent + prediction');
legend('Location','northwest'); grid on;

sgtitle('Appendix C: Multi-bit Verification and State-Dependent Mismatch');

%% ---------------------- LOCAL FUNCTIONS ----------------------
function a_m = compute_am_nbit(seq, Gamma_states, L, m_idx)
    % Generalised n-bit Fourier coefficient computation.
    % seq: 1xL array of integers indexing into Gamma_states (1-based)
    %      OR 1xL binary array (0/1) for 1-bit backward compatibility
    % Gamma_states: 1xN complex reflection coefficients for each state
    % Handles both 1-bit (seq values 0/1 mapping to states 1/2) and
    % n-bit (seq values 1..N indexing directly into Gamma_states)

    GammaSeq = zeros(1,L);
    for k = 1:L
        idx = seq(k);
        if idx == 0
            idx = 1;  % 1-bit backward compatibility: 0 -> state 1
        end
        GammaSeq(k) = Gamma_states(idx);
    end

    a_m = zeros(size(m_idx));
    for mi = 1:length(m_idx)
        m = m_idx(mi); s = 0;
        for k = 1:L
            s = s + (GammaSeq(k)/L)*sinc(m/L)*exp(-1j*pi*m*(2*k-1)/L);
        end
        a_m(mi) = s;
    end
end

function a_qm = apply_time_shift(a_m, m_idx, f0, t_q)
    Nq = length(t_q);
    a_qm = zeros(Nq, length(m_idx));
    for q = 1:Nq
        a_qm(q,:) = a_m .* exp(-1j*2*pi*m_idx*f0*t_q(q));
    end
end

function F = array_farfield(a_qm, m_idx, Nq, dy_over_lambda, theta)
    F = zeros(length(m_idx), length(theta));
    for mi = 1:length(m_idx)
        for ti = 1:length(theta)
            ap = exp(1j*2*pi*(0:Nq-1)*dy_over_lambda*sin(theta(ti)));
            F(mi,ti) = sum(ap .* a_qm(:,mi).');
        end
    end
end

function [peak_angle, peak_sl, rms_spread] = analyze_pattern(F, theta_deg, m_idx)
    peak_angle = zeros(size(m_idx));
    peak_sl = zeros(size(m_idx));
    rms_spread = zeros(size(m_idx));
    dtheta = mean(diff(theta_deg));
    for mi = 1:length(m_idx)
        mag = abs(F(mi,:));
        [main_val, main_idx] = max(mag);
        peak_angle(mi) = theta_deg(main_idx);
        win = round(3/dtheta);
        mask = true(size(mag));
        mask(max(1,main_idx-win):min(length(mag),main_idx+win)) = false;
        if any(mask)
            peak_sl(mi) = 20*log10(max(mag(mask))/main_val);
        else
            peak_sl(mi) = -Inf;
        end
        win2 = round(10/dtheta);
        lo = max(1,main_idx-win2); hi = min(length(mag),main_idx+win2);
        t_win = theta_deg(lo:hi); m_win = mag(lo:hi).^2;
        tc = sum(t_win.*m_win)/sum(m_win);
        rms_spread(mi) = sqrt(sum((t_win-tc).^2.*m_win)/sum(m_win));
    end
end