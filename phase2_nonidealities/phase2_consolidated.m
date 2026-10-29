%% ========================================================================
%  PHASE 2 — Consolidated Quantitative Metrics
%  Computes ALL five supervisor-requested metrics for all four experiments:
%    1. Harmonic efficiency (eta_m)
%    2. Beam-pointing error (deg)
%    3. RMS beam spread (deg) — NEW
%    4. Peak sidelobe level (dB)
%    5. Pattern overlap with ideal (0-1) — NEW
%  Outputs: printed sensitivity tables + summary figure
% =========================================================================
clear; clc; close all; rng(42);

%% ---------------------- 1. SHARED BASELINE ----------------------
Gamma0 = 0.911*exp(1j*deg2rad(-102));
Gamma1 = 0.952*exp(1j*deg2rad(94));
seq    = [1 0 0 0 0 0 0 0];
L      = length(seq);
f0     = 1e6; T0 = 1/f0;
Nq     = 8; dy_over_lambda = 0.5;
Mmax   = 3; m_idx = -Mmax:Mmax; n_m = length(m_idx);
t_q    = (0:Nq-1)*T0/Nq;
theta_deg = -90:0.02:90;
theta  = deg2rad(theta_deg);
N_mc   = 500;

% Ideal baseline
a_m_ideal  = compute_am(seq, L, Gamma0, Gamma1, m_idx);
a_qm_ideal = apply_time_shift(a_m_ideal, m_idx, f0, t_q);
F_ideal    = array_farfield(a_qm_ideal, m_idx, Nq, dy_over_lambda, theta);
P_ideal    = abs(a_m_ideal).^2;
eta_ideal  = P_ideal/sum(P_ideal);
[pa_ideal, sl_ideal, spread_ideal] = analyze_pattern(F_ideal, theta_deg, m_idx);

fprintf('=== IDEAL BASELINE METRICS ===\n');
fprintf('%6s %10s %12s %12s %12s\n','m','eta_m','PA(deg)','SLL(dB)','Spread(deg)');
for i=1:n_m
    fprintf('%6d %10.4f %12.2f %12.2f %12.4f\n',...
        m_idx(i),eta_ideal(i),pa_ideal(i),sl_ideal(i),spread_ideal(i));
end

%% ---------------------- 2. EXP 1: RISE/FALL TIME ----------------------
tr_fracs = [0, 0.02, 0.05, 0.10, 0.20, 0.30];
tr_vals  = tr_fracs*(T0/L);
n_tr     = length(tr_vals);

metrics_tr = struct('eta',zeros(n_tr,n_m),'pe',zeros(n_tr,n_m),...
    'sl',zeros(n_tr,n_m),'spread',zeros(n_tr,n_m),'overlap',zeros(n_tr,n_m));

for ii=1:n_tr
    filt = sinc(m_idx*f0*tr_vals(ii));
    a_m  = a_m_ideal.*filt;
    P_n  = abs(a_m).^2; metrics_tr.eta(ii,:) = P_n/sum(P_n);
    a_qm = apply_time_shift(a_m, m_idx, f0, t_q);
    F    = array_farfield(a_qm, m_idx, Nq, dy_over_lambda, theta);
    [pa,sl,sp] = analyze_pattern(F, theta_deg, m_idx);
    metrics_tr.pe(ii,:)     = pa - pa_ideal;
    metrics_tr.sl(ii,:)     = sl;
    metrics_tr.spread(ii,:) = sp;
    metrics_tr.overlap(ii,:)= pattern_overlap(F, F_ideal, m_idx);
end

fprintf('\n=== EXP 1: RISE/FALL TIME — ALL METRICS ===\n');
print_metrics_table(tr_fracs, 'tr/(T0/L)', metrics_tr, m_idx, eta_ideal);

% Quantitative sensitivity: efficiency slope (dB per unit tr fraction)
d_eta0_dtr = (metrics_tr.eta(end,m_idx==0) - eta_ideal(m_idx==0)) / tr_fracs(end);
fprintf('Carrier efficiency sensitivity: d(eta_0)/d(tr_frac) = %.4f per unit\n', d_eta0_dtr);

%% ---------------------- 3. EXP 2: TIMING JITTER ----------------------
sigma_fracs = [0, 0.01, 0.02, 0.05, 0.10, 0.20];
sigma_vals  = sigma_fracs*(T0/L);
n_sig       = length(sigma_vals);

pe_mean_j = zeros(n_sig,n_m); pe_std_j  = zeros(n_sig,n_m);
sl_mean_j = zeros(n_sig,n_m); sl_std_j  = zeros(n_sig,n_m);
sp_mean_j = zeros(n_sig,n_m); sp_std_j  = zeros(n_sig,n_m);
ov_mean_j = zeros(n_sig,n_m); ov_std_j  = zeros(n_sig,n_m);

for ii=1:n_sig
    sigma=sigma_vals(ii);
    pe_mc=zeros(N_mc,n_m); sl_mc=zeros(N_mc,n_m);
    sp_mc=zeros(N_mc,n_m); ov_mc=zeros(N_mc,n_m);
    for mc=1:N_mc
        delta_q=sigma*randn(1,Nq);
        a_qm_j=apply_time_shift(a_m_ideal,m_idx,f0,t_q+delta_q);
        F_j=array_farfield(a_qm_j,m_idx,Nq,dy_over_lambda,theta);
        [pa,sl,sp]=analyze_pattern(F_j,theta_deg,m_idx);
        pe_mc(mc,:)=pa-pa_ideal;
        sl_mc(mc,:)=sl; sp_mc(mc,:)=sp;
        ov_mc(mc,:)=pattern_overlap(F_j,F_ideal,m_idx);
    end
    pe_mean_j(ii,:)=mean(pe_mc,1); pe_std_j(ii,:)=std(pe_mc,0,1);
    sl_mean_j(ii,:)=mean(sl_mc,1); sl_std_j(ii,:)=std(sl_mc,0,1);
    sp_mean_j(ii,:)=mean(sp_mc,1); sp_std_j(ii,:)=std(sp_mc,0,1);
    ov_mean_j(ii,:)=mean(ov_mc,1); ov_std_j(ii,:)=std(ov_mc,0,1);
    fprintf('Jitter sigma=%.2f done\n',sigma_fracs(ii));
end

fprintf('\n=== EXP 2: JITTER — RMS POINTING ERROR (deg) ===\n');
print_stoch_table(sigma_fracs,'sigma/(T0/L)',pe_std_j,m_idx);
fprintf('\n=== EXP 2: JITTER — MEAN PATTERN OVERLAP ===\n');
print_stoch_table(sigma_fracs,'sigma/(T0/L)',ov_mean_j,m_idx);
fprintf('\n=== EXP 2: JITTER — MEAN RMS BEAM SPREAD (deg) ===\n');
print_stoch_table(sigma_fracs,'sigma/(T0/L)',sp_mean_j,m_idx);

% Quantitative sensitivity coefficients
fprintf('\n=== EXP 2: LINEAR SENSITIVITY COEFFICIENTS ===\n');
fprintf('(RMS pointing error per unit sigma/(T0/L), from linear fit)\n');
fprintf('%6s %12s\n','m','deg/unit_sigma');
for i=1:n_m
    slope = polyfit(sigma_fracs', pe_std_j(:,i), 1);
    fprintf('%6d %12.3f\n', m_idx(i), slope(1));
end

%% ---------------------- 4. EXP 3: ELEMENT MISMATCH ----------------------
sigma_phi_vals = [0,1,2,5,10,20];
n_sP = length(sigma_phi_vals);

pe_std_P=zeros(n_sP,n_m); sp_std_P=zeros(n_sP,n_m);
ov_mean_P=zeros(n_sP,n_m); sl_std_P=zeros(n_sP,n_m);

for ii=1:n_sP
    sP=deg2rad(sigma_phi_vals(ii));
    pe_mc=zeros(N_mc,n_m); sp_mc=zeros(N_mc,n_m);
    ov_mc=zeros(N_mc,n_m); sl_mc=zeros(N_mc,n_m);
    for mc=1:N_mc
        delta_phi=sP*randn(1,Nq);
        a_qm_mc=zeros(Nq,n_m);
        for q=1:Nq
            G0_q=abs(Gamma0)*exp(1j*(angle(Gamma0)+delta_phi(q)));
            G1_q=abs(Gamma1)*exp(1j*(angle(Gamma1)+delta_phi(q)));
            a_m_q=compute_am(seq,L,G0_q,G1_q,m_idx);
            a_qm_mc(q,:)=a_m_q.*exp(-1j*2*pi*m_idx*f0*t_q(q));
        end
        F_mc=array_farfield(a_qm_mc,m_idx,Nq,dy_over_lambda,theta);
        [pa,sl,sp]=analyze_pattern(F_mc,theta_deg,m_idx);
        pe_mc(mc,:)=pa-pa_ideal;
        sl_mc(mc,:)=sl; sp_mc(mc,:)=sp;
        ov_mc(mc,:)=pattern_overlap(F_mc,F_ideal,m_idx);
    end
    pe_std_P(ii,:)=std(pe_mc,0,1);
    sp_std_P(ii,:)=std(sp_mc,0,1);
    ov_mean_P(ii,:)=mean(ov_mc,1);
    sl_std_P(ii,:)=std(sl_mc,0,1);
    fprintf('Phase mismatch sigma_phi=%.1f deg done\n',sigma_phi_vals(ii));
end

fprintf('\n=== EXP 3: PHASE MISMATCH — RMS POINTING ERROR (deg) ===\n');
print_stoch_table(sigma_phi_vals,'sigma_phi(deg)',pe_std_P,m_idx);
fprintf('\n=== EXP 3: PHASE MISMATCH — MEAN PATTERN OVERLAP ===\n');
print_stoch_table(sigma_phi_vals,'sigma_phi(deg)',ov_mean_P,m_idx);

fprintf('\n=== EXP 3: LINEAR SENSITIVITY COEFFICIENTS ===\n');
fprintf('(RMS pointing error per degree of phase spread)\n');
fprintf('%6s %14s\n','m','deg/deg_spread');
for i=1:n_m
    slope=polyfit(sigma_phi_vals',pe_std_P(:,i),1);
    fprintf('%6d %14.4f\n',m_idx(i),slope(1));
end

%% ---------------------- 5. EXP 4: WAVEFORM DISTORTION ----------------------
Nt=8192; t_wf=linspace(0,T0,Nt+1); t_wf=t_wf(1:end-1);
alpha_droop=[0,0.05,0.10,0.20,0.30,0.50];
n_dr=length(alpha_droop);
eta_dr=zeros(n_dr,n_m); ov_dr=zeros(n_dr,n_m);
pe_dr=zeros(n_dr,n_m); sp_dr=zeros(n_dr,n_m);

Gamma_ideal_t=build_ideal(t_wf,T0,seq,Gamma0,Gamma1);
a_m_ideal_num=fourier_coeffs(Gamma_ideal_t,t_wf,T0,m_idx);

for ii=1:n_dr
    Gamma_t=build_droop(t_wf,T0,seq,Gamma0,Gamma1,alpha_droop(ii));
    a_m=fourier_coeffs(Gamma_t,t_wf,T0,m_idx);
    P_n=abs(a_m).^2; eta_dr(ii,:)=P_n/sum(P_n);
    a_qm=apply_time_shift(a_m,m_idx,f0,t_q);
    F=array_farfield(a_qm,m_idx,Nq,dy_over_lambda,theta);
    % reference F_ideal_num for overlap
    a_qm_ref=apply_time_shift(a_m_ideal_num,m_idx,f0,t_q);
    F_ref=array_farfield(a_qm_ref,m_idx,Nq,dy_over_lambda,theta);
    [pa,~,sp]=analyze_pattern(F,theta_deg,m_idx);
    pe_dr(ii,:)=pa-pa_ideal;
    sp_dr(ii,:)=sp;
    ov_dr(ii,:)=pattern_overlap(F,F_ref,m_idx);
end

fprintf('\n=== EXP 4: DROOP — CARRIER EFFICIENCY & OVERLAP ===\n');
fprintf('%10s %12s %14s\n','alpha','eta_0','overlap_m0');
m0idx=find(m_idx==0);
for ii=1:n_dr
    fprintf('%10.2f %12.4f %14.4f\n',alpha_droop(ii),eta_dr(ii,m0idx),ov_dr(ii,m0idx));
end

d_eta0_dalpha=(eta_dr(end,m0idx)-eta_dr(1,m0idx))/alpha_droop(end);
fprintf('Carrier efficiency sensitivity: d(eta_0)/d(alpha) = %.4f per unit droop\n',d_eta0_dalpha);

%% ---------------------- 6. SUMMARY FIGURE ----------------------
figure('Name','Phase 2: Quantitative Metrics Summary','Position',[50 50 1400 900]);

% Panel 1: RMS pointing error — jitter vs phase mismatch
subplot(2,3,1);
plot(sigma_fracs, pe_std_j(:,[1,4,7]),'-o','LineWidth',1.3);
xlabel('\sigma/(T_0/L)'); ylabel('RMS pointing error (deg)');
title('Jitter: RMS pointing error');
legend('m=\pm3','m=0','m=\pm3','Location','northwest'); grid on;

subplot(2,3,2);
plot(sigma_phi_vals, pe_std_P(:,[1,4,7]),'-s','LineWidth',1.3);
xlabel('\sigma_\phi (deg)'); ylabel('RMS pointing error (deg)');
title('Phase mismatch: RMS pointing error');
legend('m=\pm3','m=0','Location','northwest'); grid on;

% Panel 2: Pattern overlap
subplot(2,3,3);
plot(sigma_fracs, ov_mean_j(:,[1,4,7]),'-o','LineWidth',1.3); hold on;
plot(sigma_phi_vals/100, ov_mean_P(:,[1,4,7]),'--s','LineWidth',1.3);
xlabel('Normalized perturbation'); ylabel('Pattern overlap');
title('Pattern overlap vs perturbation');
legend('Jitter m=\pm3','Jitter m=0','Jitter m=\pm3',...
    'Mismatch m=\pm3','Mismatch m=0','Location','southwest'); grid on;
ylim([0 1]);

% Panel 3: Beam spread — jitter
subplot(2,3,4);
errorbar(sigma_fracs, sp_mean_j(:,[1,4,7]), sp_std_j(:,[1,4,7]),'-o','LineWidth',1.2);
xlabel('\sigma/(T_0/L)'); ylabel('RMS beam spread (deg)');
title('Jitter: RMS beam spread');
legend('m=\pm3','m=0','Location','northwest'); grid on;

% Panel 4: Droop efficiency
subplot(2,3,5);
plot(alpha_droop, eta_dr,'-s','LineWidth',1.3);
xlabel('Droop \alpha'); ylabel('\eta_m');
title('Droop: harmonic efficiency redistribution');
legend(arrayfun(@(m)sprintf('m=%+d',m),m_idx,'UniformOutput',false),...
    'Location','eastoutside'); grid on;

% Panel 5: Sensitivity bar chart (RMS pointing error at worst case, m=+1)
subplot(2,3,6);
m1=find(m_idx==1);
worst_case = [0, pe_std_j(end,m1), pe_std_P(end,m1), 0];
bar(worst_case);
set(gca,'XTickLabel',{'Rise/fall','Jitter','Phase mismatch','Droop'});
ylabel('RMS pointing error at m=+1 (deg)');
title('Worst-case comparison across experiments');
grid on;

sgtitle('Phase 2 Non-Ideality Study: Quantitative Metrics Summary');

%% ---------------------- LOCAL FUNCTIONS ----------------------
function a_m = compute_am(seq, L, Gamma0, Gamma1, m_idx)
    GammaSeq=zeros(1,L);
    for k=1:L; GammaSeq(k)=(seq(k)==0)*Gamma0+(seq(k)==1)*Gamma1; end
    a_m=zeros(size(m_idx));
    for idx=1:length(m_idx)
        m=m_idx(idx); s=0;
        for k=1:L
            s=s+(GammaSeq(k)/L)*sinc(m/L)*exp(-1j*pi*m*(2*k-1)/L);
        end
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

function [peak_angle, peak_sl, rms_spread] = analyze_pattern(F, theta_deg, m_idx)
    peak_angle=zeros(size(m_idx)); peak_sl=zeros(size(m_idx));
    rms_spread=zeros(size(m_idx));
    dtheta=mean(diff(theta_deg));
    for mi=1:length(m_idx)
        mag=abs(F(mi,:)); [main_val,main_idx]=max(mag);
        peak_angle(mi)=theta_deg(main_idx);
        % sidelobe
        win=round(3/dtheta); mask=true(size(mag));
        mask(max(1,main_idx-win):min(length(mag),main_idx+win))=false;
        if any(mask); peak_sl(mi)=20*log10(max(mag(mask))/main_val);
        else; peak_sl(mi)=-Inf; end
        % RMS beam spread: second moment over main lobe window (±10 deg)
        win2=round(10/dtheta);
        lo=max(1,main_idx-win2); hi=min(length(mag),main_idx+win2);
        theta_win=theta_deg(lo:hi); mag_win=mag(lo:hi).^2;
        theta_c=sum(theta_win.*mag_win)/sum(mag_win);
        rms_spread(mi)=sqrt(sum((theta_win-theta_c).^2.*mag_win)/sum(mag_win));
    end
end

function ov = pattern_overlap(F_pert, F_ref, m_idx)
    % Normalized inner product per harmonic, real part of cosine similarity
    ov=zeros(size(m_idx));
    for mi=1:length(m_idx)
        a=F_pert(mi,:); b=F_ref(mi,:);
        ov(mi)=abs(sum(conj(a).*b))/(norm(a)*norm(b));
    end
end

function print_metrics_table(sweep_vals, sweep_label, metrics, m_idx, eta_ref)
    fprintf('%14s | %s\n', sweep_label, 'eta_0  | PE_m0 | SLL_m0 | Spread_m0 | Overlap_m0');
    m0=find(m_idx==0);
    for ii=1:length(sweep_vals)
        fprintf('%14.3f | %.4f | %+.3f | %6.2f | %9.4f | %10.4f\n',...
            sweep_vals(ii), metrics.eta(ii,m0), metrics.pe(ii,m0),...
            metrics.sl(ii,m0), metrics.spread(ii,m0), metrics.overlap(ii,m0));
    end
end

function print_stoch_table(sweep_vals, sweep_label, data, m_idx)
    fprintf('%14s', sweep_label);
    for i=1:length(m_idx); fprintf('%10s',sprintf('m=%+d',m_idx(i))); end; fprintf('\n');
    for ii=1:length(sweep_vals)
        fprintf('%14.2f', sweep_vals(ii));
        for i=1:length(m_idx); fprintf('%10.4f',data(ii,i)); end; fprintf('\n');
    end
end

function Gamma_t = build_ideal(t, T0, seq, Gamma0, Gamma1)
    L=length(seq); tau=T0/L; Gamma_t=zeros(size(t));
    for i=1:length(t)
        k=min(floor(mod(t(i),T0)/tau)+1,L);
        Gamma_t(i)=(seq(k)==0)*Gamma0+(seq(k)==1)*Gamma1;
    end
end

function Gamma_t = build_droop(t, T0, seq, Gamma0, Gamma1, alpha)
    L=length(seq); tau=T0/L; Gamma_t=zeros(size(t));
    for i=1:length(t)
        tt=mod(t(i),T0); k=min(floor(tt/tau)+1,L);
        t_in=tt-(k-1)*tau;
        G=(seq(k)==0)*Gamma0+(seq(k)==1)*Gamma1;
        if seq(k)==1; G=G+alpha*(t_in/tau)*(Gamma0-Gamma1); end
        Gamma_t(i)=G;
    end
end

function a_m = fourier_coeffs(Gamma_t, t, T0, m_idx)
    dt=t(2)-t(1); a_m=zeros(size(m_idx));
    for idx=1:length(m_idx)
        a_m(idx)=sum(Gamma_t.*exp(-1j*2*pi*m_idx(idx)*t/T0))*dt/T0;
    end
end