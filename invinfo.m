close all; clear all; clc; % clear workspace
%load data_stramski.mat; % load data

%if including water:
%load abs_h2o_data.mat;
%absh2o = interp1(wvlh2o,ah2o,wvl);
%clear wvlh2o ah2o;

k = spectra; % each column is an absorption spectrum [1/m]
l = wvl'; % column vector of wavelengths [nm]
labels = labels;
%k(:,end+1) = absh2o; % add a_w as needed
%k(:,end+1) = a_NAP, a_NAP + a_CDOM, or a_NAP+CDOM; % add other absorption spectra as needed
n = size(k,2); % number of kernels
clearvars -EXCEPT k l n labels;

% normalize kernels by square norm
% n.b. if error varies with wavelength, easy to incorporate,
% just rescale all kernels by error weights (1/eps^2)
m = (trapz(l,k.^2)).^.5;
m = repmat(m,length(l),1);
k = k./m;

% compute covariance matrix for kernels
Cij = zeros(n);
for i = 1:n;
    for j = 1:n;
        Cij(i,j) = trapz(l,k(:,i).*k(:,j));
    end
end

% find square roots of eigenvectors of covariance matrix
% these give a sense of error tolerance
sli = real(sqrt(flipud(eig(Cij))))';
clear n m j i Cij;
%that's it!

% plotting
%{
plot(l,k,'linewidth',1.25);
axis([-Inf Inf 0 Inf])
set(gca,'ticklabelinterpreter','latex','fontsize',18,'ytick',[])
xlabel('wavelength [nm]','interpreter','latex')
ylabel('normalized absorption','interpreter','latex')
legend(labels,'interpreter','latex')
slis = regexprep(num2str(sli,2),'\s+',', ');
title(['$\sqrt{\lambda_i}$ = ',slis],'interpreter','latex')
axis([min(l) max(l) 0 1.05.*max(k(:))];
clear slis;
%}