%% ========================================================================
%  PHASE 2 — Experiment 3: Element-wise Mismatch
%  Key findings:
%    - Amplitude mismatch: zero pointing error (not plotted, stated in text)
%    - Phase mismatch: non-zero RMS pointing error, linear growth with sigma_phi
%  Sidelobe unaffected in both cases (not plotted, stated in text).
%  ONE figure produced: phase mismatch pointing error vs sigma_phi.
% =========================================================================
clear; clc; close all; rng(42);

%% ---------------------- 1. FIXED BASELINE ----------------------
Gamma0_nom = 0.911*exp(1j*deg2rad(-102));
Gamma1_nom = 0.952*exp(1j*deg2rad(94));
seq = [1 0 0 0 0 0 0 0];
L   = length(seq);
f0  = 1e6; T0 = 1/f0;
Nq  = 8; dy_over_lambda = 0.5;
Mmax = 3; m_idx = -Mmax:Mmax;
t_q = (0:Nq-1)*T0/Nq;
theta_deg = -90:0.02:90;
theta = deg2rad(theta_deg);
N_mc = 500;
sigma_A_vals   = [0, 0.01, 0.02, 0.05, 0.10, 0.20];
sigma_phi_vals = [0, 1, 2, 5, 10, 20];

%% ---------------------- 2. IDEAL BASELINE ----------------------
a_m_ideal = compute_am_nominal(seq, L, Gamma0_nom, Gamma1_nom, m_idx);
a_qm_ideal = apply_time_shift(a_m_ideal, m_idx, f0, t_q);
F_ideal = array_farfield(a_qm_ideal, m_idx, Nq, dy_over_lambda, theta);
[peak_angle_ideal, ~] = analyze_pattern(F_ideal, theta_deg, m_idx);
target_angle = asind(m_idx/4);

fprintf('=== Baseline steering angle check ===\n');
fprintf('%6s %12s %12s\n','m','Simulated','Theoretical');
for i=1:length(m_idx)
    fprintf('%6d %12.2f %12.2f\n',m_idx(i),peak_angle_ideal(i),target_angle(i));
end

%% ---------------------- 3. SWEEP A: AMPLITUDE MISMATCH (results table only) ----------------------
fprintf('\nRunning amplitude mismatch sweep...\n');
n_sA=length(sigma_A_vals); n_m=length(m_idx);
pe_std_A=zeros(n_sA,n_m);

for ii=1:n_sA
    sA=sigma_A_vals(ii); pe_mc=zeros(N_mc,n_m);
    for mc=1:N_mc
        delta_A=sA*randn(1,Nq);
        a_qm_mc=zeros(Nq,n_m);
        for q=1:Nq
            G0_q=abs(Gamma0_nom)*(1+delta_A(q))*exp(1j*angle(Gamma0_nom));
            G1_q=abs(Gamma1_nom)*(1+delta_A(q))*exp(1j*angle(Gamma1_nom));
            a_m_q=compute_am_nominal(seq,L,G0_q,G1_q,m_idx);
            a_qm_mc(q,:)=a_m_q.*exp(-1j*2*pi*m_idx*f0*t_q(q));
        end
        F_mc=array_farfield(a_qm_mc,m_idx,Nq,dy_over_lambda,theta);
        [pa,~]=analyze_pattern(F_mc,theta_deg,m_idx);
        pe_mc(mc,:)=pa-peak_angle_ideal;
    end
    pe_std_A(ii,:)=std(pe_mc,0,1);
    fprintf('  sigma_A=%.2f done\n',sA);
end

fprintf('\n=== RMS pointing error vs amplitude mismatch (all zeros expected) ===\n');
fprintf('%12s','sigma_A');
for i=1:n_m; fprintf('%10s',sprintf('m=%+d',m_idx(i))); end; fprintf('\n');
for ii=1:n_sA
    fprintf('%12.2f',sigma_A_vals(ii));
    for i=1:n_m; fprintf('%10.3f',pe_std_A(ii,i)); end; fprintf('\n');
end

%% ---------------------- 4. SWEEP B: PHASE MISMATCH ----------------------
fprintf('\nRunning phase mismatch sweep...\n');
n_sP=length(sigma_phi_vals);
pe_mean_P=zeros(n_sP,n_m); pe_std_P=zeros(n_sP,n_m);

for ii=1:n_sP
    sP=deg2rad(sigma_phi_vals(ii)); pe_mc=zeros(N_mc,n_m);
    for mc=1:N_mc
        delta_phi=sP*randn(1,Nq);
        a_qm_mc=zeros(Nq,n_m);
        for q=1:Nq
            G0_q=abs(Gamma0_nom)*exp(1j*(angle(Gamma0_nom)+delta_phi(q)));
            G1_q=abs(Gamma1_nom)*exp(1j*(angle(Gamma1_nom)+delta_phi(q)));
            a_m_q=compute_am_nominal(seq,L,G0_q,G1_q,m_idx);
            a_qm_mc(q,:)=a_m_q.*exp(-1j*2*pi*m_idx*f0*t_q(q));
        end
        F_mc=array_farfield(a_qm_mc,m_idx,Nq,dy_over_lambda,theta);
        [pa,~]=analyze_pattern(F_mc,theta_deg,m_idx);
        pe_mc(mc,:)=pa-peak_angle_ideal;
    end
    pe_mean_P(ii,:)=mean(pe_mc,1);
    pe_std_P(ii,:) =std(pe_mc,0,1);
    fprintf('  sigma_phi=%.1f deg done\n',sigma_phi_vals(ii));
end

fprintf('\n=== RMS pointing error vs phase mismatch ===\n');
fprintf('%12s','sigma_phi');
for i=1:n_m; fprintf('%10s',sprintf('m=%+d',m_idx(i))); end; fprintf('\n');
for ii=1:n_sP
    fprintf('%12.1f',sigma_phi_vals(ii));
    for i=1:n_m; fprintf('%10.3f',pe_std_P(ii,i)); end; fprintf('\n');
end

%% ---------------------- 5. FIGURE: phase mismatch pointing error only ----------------------
figure('Name','Exp3: Phase mismatch pointing error');
errorbar(repmat(sigma_phi_vals',1,n_m), pe_mean_P, pe_std_P, '-o','LineWidth',1.2);
xlabel('\sigma_\phi  (phase mismatch std, degrees)');
ylabel('Pointing error (deg)   mean \pm std');
title('Experiment 3: Beam-pointing error vs element phase mismatch (N=500)');
legend(arrayfun(@(m) sprintf('m=%+d',m),m_idx,'UniformOutput',false),'Location','eastoutside');
grid on;

%% ---------------------- LOCAL FUNCTIONS ----------------------
function a_m = compute_am_nominal(seq, L, Gamma0, Gamma1, m_idx)
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
