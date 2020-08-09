function [Ave_PDC,delta_PDC,theta_PDC,alpha_PDC,beta_PDC,gamma_PDC,gpdc,pdc,dtf,dc,coh,pcoh,pcoh2]= FunctionPDC(A,Q,nfft,fc,ch,p,fr)
%-------------------------------------------------------------------------%
%           Compute Frequency-domain Connectivity Measures                %
%-------------------------------------------------------------------------%
K=1;N=ch;
% GPDC_PDC = zeros(N,N,nfft);
Ave_PDC = zeros(N,N,K);
delta_PDC = zeros(N,N,K); 
theta_PDC = zeros(N,N,K); 
alpha_PDC = zeros(N,N,K); 
beta_PDC = zeros(N,N,K);
gamma_PDC = zeros(N,N,K);

for j=1:K
    
Am = squeeze(A(1:N,1:N*p,j));
Su = squeeze(Q(1:N,1:N,j));

[dc,dtf,pdc,gpdc,coh,pcoh,pcoh2,h,sp,pp,f] = fdMVAR(Am,Su,nfft,fc);
% S=abs(sp); % spectral matrix
% P=abs(pp); % inverse spectral matrix
% DC=abs(dc).^2; % directed coherence
PDC=abs(pdc).^2; % partial directed coherence%gpdc
% COH=abs(coh).^2; %coherence
% PCOH=abs(pcoh).^2; % partial coherence
% PCOH2=abs(pcoh2).^2; % partial coherence, other definition to check equality
    
    
    % All-band
    for fi = 1:nfft
        Ave_PDC(:,:,j) = Ave_PDC(:,:,j) + PDC(:,:,fi); end
    Ave_PDC(:,:,j) = Ave_PDC(:,:,j)/nfft;
%     for i=1:N
%         Ave_PDC(i,i,j) = 0; end

    % Delta-band
    lf = 1; uf = 3;
    for fi = lf:uf*fr
        delta_PDC(:,:,j) = delta_PDC(:,:,j) + PDC(:,:,fi); end
    delta_PDC(:,:,j) = delta_PDC(:,:,j)/(uf*fr - lf*fr);
    for i=1:N
        theta_PDC(i,i,j) = 0; end

    
    % Theta-band
    lf = 4; uf = 7;
    for fi = lf*fr:uf*fr
        theta_PDC(:,:,j) = theta_PDC(:,:,j) + PDC(:,:,fi); end
    theta_PDC(:,:,j) = theta_PDC(:,:,j)/(uf*fr - lf*fr);
    for i=1:N
        theta_PDC(i,i,j) = 0; end

    % Alpha-band
    lf = 8; uf = 12;
    for fi = lf*fr:uf*fr
        alpha_PDC(:,:,j) = alpha_PDC(:,:,j) + PDC(:,:,fi); end
    alpha_PDC(:,:,j) = alpha_PDC(:,:,j)/(uf*fr - lf*fr);
    for i=1:N
        alpha_PDC(i,i,j) = 0; end
    
    % Beta-band
    lf = 13; uf = 30;
    for fi = lf*fr:uf*fr
        beta_PDC(:,:,j) = beta_PDC(:,:,j) + PDC(:,:,fi); end
    beta_PDC(:,:,j) = beta_PDC(:,:,j)/(uf*fr - lf*fr);
    for i=1:N
        beta_PDC(i,i,j) = 0; end


    % Gamma-band
    lf = 31; uf = nfft/2;
    for fi = lf*fr:uf*fr
        gamma_PDC(:,:,j) = gamma_PDC(:,:,j) + PDC(:,:,fi); end
    gamma_PDC(:,:,j) = gamma_PDC(:,:,j)/(uf*fr - lf*fr);
    for i=1:N
        beta_PDC(i,i,j) = 0; end


end