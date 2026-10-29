%% ========================================================================
%  PHASE 2 — Experiment 2: Timing Jitter
%  Key finding: zero mean pointing error, non-zero RMS that grows with sigma
%  and scales with harmonic order |m|. Sidelobe unaffected (not plotted).
%  TWO figures: pointing error RMS vs sigma + pattern spread visualization.
% =========================================================================
clear; clc; close all; rng(42);

%% ---------------------- 1. FIXED BASELINE ----------------------
Gamma0 = 0.911*exp(1j*deg2rad(-102));
Gamma1 = 0.952*exp(1j*deg2rad(94));
seq = [1 0 0 0 0 0 0 0];
L   = length(seq);
f0  = 1e6; T0 = 1/f0;
Nq  = 8; dy_over_lambda = 0.5;
Mmax = 3; m_idx = -Mmax:Mmax;
t_q_ideal = (0:Nq-1)*T0/Nq;
theta_deg = -90:0.02:90;
theta = deg2rad(theta_deg);
sigma_fracs = [0, 0.01, 0.02, 0.05, 0.10, 0.20];
sigma_vals  = sigma_fracs*(T0/L);
N_mc = 500;

%% ---------------------- 2. IDEAL BASELINE ----------------------
a_m_ideal = compute_am(seq, L, Gamma0, Gamma1, m_idx);
a_qm_ideal = apply_time_shift(a_m_ideal, m_idx, f0, t_q_ideal);
F_ideal = array_farfield(a_qm_ideal, m_idx, Nq, dy_over_lambda, theta);
[peak_angle_ideal, ~] = analyze_pattern(F_ideal, theta_deg, m_idx);
target_angle = asind(m_idx/4);

fprintf('=== Baseline steering angle check ===\n');
fprintf('%6s %12s %12s\n','m','Simulated','Theoretical');
for i=1:length(m_idx)
    fprintf('%6d %12.2f %12.2f\n',m_idx(i),peak_angle_ideal(i),target_angle(i));
end

%% ---------------------- 3. MONTE CARLO SWEEP ----------------------
n_sig=length(sigma_vals); n_m=length(m_idx);
pe_mean=zeros(n_sig,n_m); pe_std=zeros(n_sig,n_m);

for ii=1:n_sig
    sigma=sigma_vals(ii);
    pe_mc=zeros(N_mc,n_m);
    for mc=1:N_mc
        delta_q=sigma*randn(1,Nq);
        a_qm_j=apply_time_shift(a_m_ideal,m_idx,f0,t_q_ideal+delta_q);
        F_j=array_farfield(a_qm_j,m_idx,Nq,dy_over_lambda,theta);
        [pa,~]=analyze_pattern(F_j,theta_deg,m_idx);
        pe_mc(mc,:)=pa-peak_angle_ideal;
    end
    pe_mean(ii,:)=mean(pe_mc,1);
    pe_std(ii,:) =std(pe_mc,0,1);
    fprintf('sigma=%.2f T0/L done\n',sigma_fracs(ii));
end

% Save worst-case std for master summary (m=+1, worst sigma)
m1 = find(m_idx==1);
pe_std_worst_jitter = pe_std(end,:);   %#ok<NASGU>
fprintf('\n=== RMS pointing error (deg) vs jitter sigma ===\n');
fprintf('%12s','sigma/T0_L');
for i=1:n_m; fprintf('%10s',sprintf('m=%+d',m_idx(i))); end; fprintf('\n');
for ii=1:n_sig
    fprintf('%12.2f',sigma_fracs(ii));
    for i=1:n_m; fprintf('%10.3f',pe_std(ii,i)); end; fprintf('\n');
end

%% ---------------------- 4. FIGURES ----------------------
% Fig 1: RMS pointing error vs sigma
figure('Name','Exp2: Pointing error RMS vs jitter');
errorbar(repmat(sigma_fracs',1,n_m), pe_mean, pe_std, '-o','LineWidth',1.2);
xlabel('\sigma / (T_0/L)'); ylabel('Pointing error (deg)   mean \pm std');
title('Experiment 2: Beam-pointing error vs timing jitter (N=500)');
legend(arrayfun(@(m) sprintf('m=%+d',m),m_idx,'UniformOutput',false),'Location','eastoutside');
grid on;

% Fig 2: Pattern spread at worst-case sigma (m=+1)
figure('Name','Exp2: Pattern spread at worst-case jitter');
sigma_worst=sigma_vals(end);
hold on;
for mc=1:20
    delta_q=sigma_worst*randn(1,Nq);
    a_qm_j=apply_time_shift(a_m_ideal,m_idx,f0,t_q_ideal+delta_q);
    F_j=array_farfield(a_qm_j,m_idx,Nq,dy_over_lambda,theta);
    plot(theta_deg,20*log10(abs(F_j(m1,:))/max(abs(F_ideal(m1,:)))),...
        'Color',[0.7 0.7 0.7 0.5],'LineWidth',0.8);
end
plot(theta_deg,20*log10(abs(F_ideal(m1,:))/max(abs(F_ideal(m1,:)))),...
    'b','LineWidth',2);
xlabel('\theta (deg)'); ylabel('Normalized pattern (dB)');
title(sprintf('Experiment 2: Pattern spread, m=+1, \\sigma=%.2f T_0/L, 20 realizations',sigma_fracs(end)));
legend('Jitter realizations','Ideal'); ylim([-40 5]); grid on;

%% ---------------------- LOCAL FUNCTIONS ----------------------
function a_m = compute_am(seq, L, Gamma0, Gamma1, m_idx)
    GammaSeq=zeros(1,L);
    for k=1:L; GammaSeq(k)=(seq(k)==0)*Gamma0+(seq(k)==1)*Gamma1; end
    a_m=zeros(size(m_idx));
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
