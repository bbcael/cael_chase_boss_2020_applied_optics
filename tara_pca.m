close all; clear all; clc;
% load ac-s data -- contact cael@noc.ac.uk for Tara ac-s data (large file)
load a_interpolated.mat; % spectra interpolated onto a common set of wavelengths

%%

% if using chl residual
%{
load ap_global_fromChl_AB.mat; % avaiable in github repository

al = interp1(wl,Achi2,l);
bl = interp1(wl,Bchi2,l);

for i = 1:size(asi,2);
    ap_a = interp1(l,asi(:,i),[650 676 715],'linear'); 
    acs_lh = (ap_a(2)-(39/65*ap_a(1)+26/65*ap_a(3)));
    chl(i) = 157*acs_lh.^1.22; 
    ap_chl(:,i) = al.*chl(i).^bl;
    i
end

asir = asi - ap_chl;
%}

%%

% normalize
sd_asir = sqrt(var(asir'));
mu_asir = mean(asir');
asirn = (asir-mu_asir')./sd_asir';

%%

[eivecs,~,~,~,pcts,~] = pca(asirn'); % or asir' if not normalizing, or just on normalized spectra if not looking at Chl-residual
raweivecs = (eivecs+mu_asir').*sd_asir';

%%

% normalize eigenvectors for plotting
eivecs = eivecs./max(abs(eivecs));

%plotting
%{
figure
plot(l,eivecs(:,1:3),'linewidth',2);
set(gca,'fontsize',18,'ticklabelinterpreter','latex','ytick',[])
ylabel('normalized absorption','interpreter','latex')
axis([min(l) 700 -1 1.1]);%1.05.*min(min(eivecs(:,1:4))) 1.05.*max(max(eivecs(:,1:4)))])
xlabel('wavelength [nm]','interpreter','latex')
hold on
plot(min(l):max(l),0.*(min(l):max(l)),'k')
box on
title('eigenvectors: residual spectra','interpreter','latex')
%lgd = legend('mode 1 (40.74\%)','mode 2 (16.06\%)','mode 3 (9.97\%)','mode 4 (4.05\%)','mode 5 (2.91\%)');%,'mode 6 (1.69\%)','mode 7')
lgd = legend('mode 1 (91.20\%)','mode 2 (5.37\%)','mode 3 (1.90\%)');%,'mode 4 (4.05\%)','mode 5 (2.91\%)');%,'mode 6 (1.69\%)','mode 7')
set(lgd,'interpreter','latex')
hold on

figure
scatter(1:length(pcts),pcts./100,'k','filled')
set(gca,'fontsize',20,'yscale','log','ticklabelinterpreter','latex')
ylabel('fraction of variance explained','interpreter','latex')
%axis([0 length(l)+1 1e-7 1])
xlabel('mode number','interpreter','latex')
box on
%title('residual, weighted','interpreter','latex')
%text(5,.4,'0.922, 0.00211, 0.00166, 0.00102...','interpreter','latex','fontsize',20)
axis([0 length(pcts)+1 1e-4 1])
%}

% that's it!
