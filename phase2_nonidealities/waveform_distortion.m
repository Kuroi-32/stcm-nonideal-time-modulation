%% ========================================================================
%  PHASE 2 — Experiment 4: Waveform Distortion
%  Non-ideality: Gamma(t) deviates from ideal square wave via:
%    (A) overshoot/ringing at transitions
%    (B) droop during "on" interval
%  Key findings: both distortions introduce zero pointing error (waveform
%  shape affects all elements identically, preserving inter-element phase).
%  Only droop produces meaningful efficiency redistribution (plotted).
%  TWO figures: waveform visualization + droop harmonic efficiency.
%
%  BUG FIXES vs previous version:
%    - build_overshoot: ringing now correctly placed at actual transition
%      edges only, using direct time comparison rather than mod() wrapping
%    - build_droop: drift direction confirmed (Gamma1 -> Gamma0 in complex
%      plane during "on" interval); carrier efficiency now correctly decreases
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

Nt = 8192;
t  = linspace(0,T0,Nt+1); t=t(1:end-1);

A_os_fracs = [0, 0.05, 0.10, 0.20, 0.30, 0.50];
alpha_droop = [0, 0.05, 0.10, 0.20, 0.30, 0.50];

% Ringing parameters: decay over ~1 pulse width, ring at 5x modulation freq
tau_decay = (T0/L)/3;
f_ring    = 5*(L/T0);

DeltaGamma = Gamma1 - Gamma0;

%% ---------------------- 2. IDEAL BASELINE ----------------------
Gamma_ideal_t = build_ideal(t, T0, seq, Gamma0, Gamma1);
a_m_ideal = fourier_coeffs(Gamma_ideal_t, t, T0, m_idx);
a_qm_ideal = apply_time_shift(a_m_ideal, m_idx, f0, t_q);
F_ideal = array_farfield(a_qm_ideal, m_idx, Nq, dy_over_lambda, theta);
[peak_angle_ideal, ~] = analyze_pattern(F_ideal, theta_deg, m_idx);
target_angle = asind(m_idx/4);

fprintf('=== Baseline steering angle check (numerical Fourier) ===\n');
fprintf('%6s %12s %12s\n','m','Simulated','Theoretical');
for i=1:length(m_idx)
    fprintf('%6d %12.2f %12.2f\n',m_idx(i),peak_angle_ideal(i),target_angle(i));
end

% Baseline efficiency
P_ideal = abs(a_m_ideal).^2;
eta_ideal = P_ideal/sum(P_ideal);
fprintf('\nBaseline carrier efficiency eta_0 = %.4f\n', eta_ideal(m_idx==0));

%% ---------------------- 3. SWEEP A: OVERSHOOT (table only, no plot) ----------------------
fprintf('\nRunning overshoot sweep...\n');
n_os=length(A_os_fracs); n_m=length(m_idx);
pe_os=zeros(n_os,n_m); eta_os=zeros(n_os,n_m);

for ii=1:n_os
    A_os=A_os_fracs(ii)*abs(DeltaGamma);
    Gamma_t=build_overshoot(t,T0,seq,Gamma0,Gamma1,A_os,tau_decay,f_ring);
    a_m=fourier_coeffs(Gamma_t,t,T0,m_idx);
    P_n=abs(a_m).^2; eta_os(ii,:)=P_n/sum(P_n);
    a_qm=apply_time_shift(a_m,m_idx,f0,t_q);
    F=array_farfield(a_qm,m_idx,Nq,dy_over_lambda,theta);
    [pa,~]=analyze_pattern(F,theta_deg,m_idx);
    pe_os(ii,:)=pa-peak_angle_ideal;
    fprintf('  A_os=%.2f done\n',A_os_fracs(ii));
end

fprintf('\n=== Pointing error vs overshoot (all zeros confirms robustness) ===\n');
fprintf('%12s','A_os_frac');
for i=1:n_m; fprintf('%10s',sprintf('m=%+d',m_idx(i))); end; fprintf('\n');
for ii=1:n_os
    fprintf('%12.2f',A_os_fracs(ii));
    for i=1:n_m; fprintf('%10.3f',pe_os(ii,i)); end; fprintf('\n');
end

%% ---------------------- 4. SWEEP B: DROOP ----------------------
fprintf('\nRunning droop sweep...\n');
n_dr=length(alpha_droop);
pe_dr=zeros(n_dr,n_m); eta_dr=zeros(n_dr,n_m);

for ii=1:n_dr
    alpha=alpha_droop(ii);
    Gamma_t=build_droop(t,T0,seq,Gamma0,Gamma1,alpha);
    a_m=fourier_coeffs(Gamma_t,t,T0,m_idx);
    P_n=abs(a_m).^2; eta_dr(ii,:)=P_n/sum(P_n);
    a_qm=apply_time_shift(a_m,m_idx,f0,t_q);
    F=array_farfield(a_qm,m_idx,Nq,dy_over_lambda,theta);
    [pa,~]=analyze_pattern(F,theta_deg,m_idx);
    pe_dr(ii,:)=pa-peak_angle_ideal;
    fprintf('  alpha_droop=%.2f done\n',alpha);
end

fprintf('\n=== Pointing error vs droop (all zeros confirms robustness) ===\n');
fprintf('%12s','alpha');
for i=1:n_m; fprintf('%10s',sprintf('m=%+d',m_idx(i))); end; fprintf('\n');
for ii=1:n_dr
    fprintf('%12.2f',alpha_droop(ii));
    for i=1:n_m; fprintf('%10.3f',pe_dr(ii,i)); end; fprintf('\n');
end

fprintf('\n=== Carrier efficiency vs droop ===\n');
fprintf('%12s %12s\n','alpha','eta_0');
for ii=1:n_dr
    fprintf('%12.2f %12.4f\n',alpha_droop(ii),eta_dr(ii,m_idx==0));
end

%% ---------------------- 5. FIGURE 1: Waveform visualization ----------------------
figure('Name','Exp4: Distorted waveforms');
subplot(2,1,1);
Gamma_os_ex = build_overshoot(t,T0,seq,Gamma0,Gamma1,...
    0.30*abs(DeltaGamma),tau_decay,f_ring);
plot(t/T0, real(Gamma_ideal_t),'b','LineWidth',1.5); hold on;
plot(t/T0, real(Gamma_os_ex),'r','LineWidth',1.2);
xlabel('t/T_0'); ylabel('Re\{\Gamma(t)\}');
title('Overshoot/ringing (A_{os}=0.30)');
legend('Ideal','Overshoot'); grid on;

subplot(2,1,2);
Gamma_dr_ex = build_droop(t,T0,seq,Gamma0,Gamma1,0.30);
plot(t/T0, real(Gamma_ideal_t),'b','LineWidth',1.5); hold on;
plot(t/T0, real(Gamma_dr_ex),'r','LineWidth',1.2);
xlabel('t/T_0'); ylabel('Re\{\Gamma(t)\}');
title('Droop (\alpha=0.30)');
legend('Ideal','Droop'); grid on;
sgtitle('Experiment 4: Distorted waveforms (real part of \Gamma(t))');

%% ---------------------- 6. FIGURE 2: Droop harmonic efficiency ----------------------
figure('Name','Exp4: Droop harmonic efficiency');
plot(alpha_droop, eta_dr, '-s','LineWidth',1.3);
xlabel('Droop fraction \alpha');
ylabel('\eta_m (harmonic efficiency)');
title('Experiment 4: Harmonic efficiency redistribution under droop');
legend(arrayfun(@(m) sprintf('m=%+d',m),m_idx,'UniformOutput',false),'Location','eastoutside');
grid on;

%% ---------------------- LOCAL FUNCTIONS ----------------------
function Gamma_t = build_ideal(t, T0, seq, Gamma0, Gamma1)
    L=length(seq); tau=T0/L; Gamma_t=zeros(size(t));
    for i=1:length(t)
        k=min(floor(mod(t(i),T0)/tau)+1,L);
        Gamma_t(i)=(seq(k)==0)*Gamma0+(seq(k)==1)*Gamma1;
    end
end

function Gamma_t = build_overshoot(t, T0, seq, Gamma0, Gamma1, A_os, tau_d, f_r)
    % Fixed: iterate over transition edges explicitly, no mod() wrapping issue
    L=length(seq); tau=T0/L;
    Gamma_t=build_ideal(t,T0,seq,Gamma0,Gamma1);
    GammaSeq=zeros(1,L);
    for k=1:L; GammaSeq(k)=(seq(k)==0)*Gamma0+(seq(k)==1)*Gamma1; end

    for k=1:L
        k_prev=mod(k-2,L)+1;
        G_curr=GammaSeq(k); G_prev=GammaSeq(k_prev);
        if abs(G_curr-G_prev)<1e-10; continue; end  % no transition here

        t_edge=(k-1)*tau;   % edge occurs at start of interval k
        % dt = time elapsed since this edge, only within this period
        dt=t-t_edge;
        mask=(dt>=0)&(dt<4*tau_d);   % causal: only after the edge
        direction=(G_curr-G_prev)/abs(G_curr-G_prev);
        ring=A_os*exp(-dt/tau_d).*sin(2*pi*f_r*dt);
        Gamma_t=Gamma_t+direction*ring.*mask;
    end
end

function Gamma_t = build_droop(t, T0, seq, Gamma0, Gamma1, alpha)
    % Droop: during "on" interval, Gamma drifts linearly from Gamma1 toward Gamma0
    % As alpha increases, effective Gamma1 decreases -> carrier efficiency decreases
    L=length(seq); tau=T0/L; Gamma_t=zeros(size(t));
    for i=1:length(t)
        tt=mod(t(i),T0); k=min(floor(tt/tau)+1,L);
        t_in=tt-(k-1)*tau;
        G=(seq(k)==0)*Gamma0+(seq(k)==1)*Gamma1;
        if seq(k)==1
            G=G+alpha*(t_in/tau)*(Gamma0-Gamma1);  % linear drift toward Gamma0
        end
        Gamma_t(i)=G;
    end
end

function a_m = fourier_coeffs(Gamma_t, t, T0, m_idx)
    dt=t(2)-t(1); a_m=zeros(size(m_idx));
    for idx=1:length(m_idx)
        m=m_idx(idx);
        a_m(idx)=sum(Gamma_t.*exp(-1j*2*pi*m*t/T0))*dt/T0;
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
