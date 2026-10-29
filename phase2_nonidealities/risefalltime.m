%% ========================================================================
%  PHASE 2 — Experiment 1: Finite Rise/Fall Time
%  Key finding: rise/fall time introduces zero pointing error and negligible
%  sidelobe degradation. Only trivial carrier efficiency change observed.
%  ONE figure produced: pointing error confirmation + pattern overlay.
% =========================================================================
clear; clc; close all;

%% ---------------------- 1. FIXED BASELINE ----------------------
Gamma0 = 0.911*exp(1j*deg2rad(-102));
Gamma1 = 0.952*exp(1j*deg2rad(94));
seq = [1 0 0 0 0 0 0 0];
L   = length(seq);
f0  = 1e6; T0 = 1/f0;
Nq  = 8; dy_over_lambda = 0.5;
Mmax = 3; m_idx = -Mmax:Mmax;
t_q = (0:Nq-1)*T0/Nq;
theta_deg = -90:0.02:90;
theta = deg2rad(theta_deg);
tr_fracs = [0, 0.02, 0.05, 0.10, 0.20, 0.30];
tr_vals  = tr_fracs*(T0/L);

%% ---------------------- 2. IDEAL BASELINE ----------------------
a_m_ideal = compute_am(seq, L, Gamma0, Gamma1, m_idx);
a_qm_ideal = apply_time_shift(a_m_ideal, m_idx, f0, t_q);
F_ideal = array_farfield(a_qm_ideal, m_idx, Nq, dy_over_lambda, theta);
[peak_angle_ideal, ~] = analyze_pattern(F_ideal, theta_deg, m_idx);
target_angle = asind(m_idx/4);

fprintf('=== Baseline (tr=0) steering angle check ===\n');
fprintf('%6s %12s %12s\n','m','Simulated','Theoretical');
for i=1:length(m_idx)
    fprintf('%6d %12.2f %12.2f\n',m_idx(i),peak_angle_ideal(i),target_angle(i));
end

%% ---------------------- 3. SWEEP ----------------------
n_tr = length(tr_vals);
pointing_error = zeros(n_tr, length(m_idx));
carrier_eff    = zeros(n_tr,1);

for ii = 1:n_tr
    tr = tr_vals(ii);
    filt = sinc(m_idx*f0*tr);
    a_m_tr = a_m_ideal.*filt;
    a_qm_tr = apply_time_shift(a_m_tr, m_idx, f0, t_q);
    F_tr = array_farfield(a_qm_tr, m_idx, Nq, dy_over_lambda, theta);
    [peak_angle_tr, ~] = analyze_pattern(F_tr, theta_deg, m_idx);
    pointing_error(ii,:) = peak_angle_tr - peak_angle_ideal;
    P_n = abs(a_m_tr).^2;
    carrier_eff(ii) = P_n(m_idx==0)/sum(P_n);
end

fprintf('\n=== Carrier efficiency vs rise/fall time ===\n');
fprintf('%10s %12s\n','tr/T0_bit','eta_0');
for ii=1:n_tr
    fprintf('%10.2f %12.4f\n',tr_fracs(ii),carrier_eff(ii));
end

%% ---------------------- 4. FIGURE: pointing error + pattern overlay ----------------------
figure('Name','Exp1: Rise/fall time summary');
subplot(1,2,1);
plot(tr_fracs, pointing_error, '-o', 'LineWidth',1.2);
xlabel('t_r / (T_0/L)'); ylabel('Pointing error (deg)');
title('Pointing error vs rise/fall time');
legend(arrayfun(@(m) sprintf('m=%+d',m),m_idx,'UniformOutput',false),'Location','eastoutside');
grid on;

subplot(1,2,2);
filt_worst = sinc(m_idx*f0*tr_vals(end));
a_m_worst  = a_m_ideal.*filt_worst;
a_qm_worst = apply_time_shift(a_m_worst, m_idx, f0, t_q);
F_worst    = array_farfield(a_qm_worst, m_idx, Nq, dy_over_lambda, theta);
m0 = find(m_idx==0);
plot(theta_deg, 20*log10(abs(F_ideal(m0,:))/max(abs(F_ideal(m0,:)))),'LineWidth',1.5); hold on;
plot(theta_deg, 20*log10(abs(F_worst(m0,:))/max(abs(F_ideal(m0,:)))),'--','LineWidth',1.2);
xlabel('\theta (deg)'); ylabel('Normalized pattern (dB)');
title(sprintf('Carrier pattern: ideal vs t_r=%.2f T_0/L',tr_fracs(end)));
legend('Ideal (t_r=0)','Worst-case t_r'); ylim([-40 5]); grid on;
sgtitle('Experiment 1: Finite Rise/Fall Time');

%% ---------------------- LOCAL FUNCTIONS ----------------------
function a_m = compute_am(seq, L, Gamma0, Gamma1, m_idx)
    GammaSeq = zeros(1,L);
    for k=1:L; GammaSeq(k)=(seq(k)==0)*Gamma0+(seq(k)==1)*Gamma1; end
    a_m = zeros(size(m_idx));
    for idx=1:length(m_idx)
        m=m_idx(idx); s=0;
        for k=1:L; s=s+(GammaSeq(k)/L)*sinc(m/L)*exp(-1j*pi*m*(2*k-1)/L); end
        a_m(idx)=s;
    end
end
function a_qm = apply_time_shift(a_m, m_idx, f0, t_q)
    Nq=length(t_q); a_qm=zeros(Nq,length(m_idx));
    for q=1:Nq; a_qm(q,:)=a_m.*exp(-1j*2*pi*m_idx*f0*t_q(q)); end
end
function F = array_farfield(a_qm, m_idx, Nq, dy_over_lambda, theta)
    F=zeros(length(m_idx),length(theta));
    for mi=1:length(m_idx)
        for ti=1:length(theta)
            ap=exp(1j*2*pi*(0:Nq-1)*dy_over_lambda*sin(theta(ti)));
            F(mi,ti)=sum(ap.*a_qm(:,mi).');
        end
    end
end
function [peak_angle, peak_sidelobe_dB] = analyze_pattern(F, theta_deg, m_idx)
    peak_angle=zeros(size(m_idx)); peak_sidelobe_dB=zeros(size(m_idx));
    for mi=1:length(m_idx)
        mag=abs(F(mi,:)); [main_val,main_idx]=max(mag);
        peak_angle(mi)=theta_deg(main_idx);
        win=round(3/mean(diff(theta_deg))); mask=true(size(mag));
        mask(max(1,main_idx-win):min(length(mag),main_idx+win))=false;
        if any(mask); peak_sidelobe_dB(mi)=20*log10(max(mag(mask))/main_val);
        else; peak_sidelobe_dB(mi)=-Inf; end
    end
end
