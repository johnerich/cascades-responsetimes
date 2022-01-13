clear all; close all

% Choose output to plot - reponse times based on horizontal or vertical MB
% gradients:
load('final_scaling_estimates_dbdx')
%load('final_scaling_estimates_dbdz')

figure(1)

subplot 421; hold on
histogram(hf(idx_big),[0:5:130],'facecolor',0.3*[1 1 1]); alpha(1)
ylabel('count')
xlabel('thickness (m)')
set(gca,'xtick',[0:20:120])
xlim([0 130])
ylim([0 100])

subplot 423;
histogram(bt(idx_big)*rho/1000,[-15:0],'facecolor',0.3*[1 1 1]); alpha(1) % note, converted back to water-equiv.
ylabel('count')
xlabel('b_{term} (m w.e. yr^{-1})')
%ylim([0 150])
xlim([-14 0])
set(gca,'xtick',[-14:2:0],'ytick',[0:50:200])

%%% very inelegant method for creating area-weighted distribution of tau...

% bin EDGES (for making histograms)
bins_tau = [0:5:80];

% bin CENTERS (for plotting in non histogram mode)
centers_tau = [(bins_tau(1)+bins_tau(2))/2:bins_tau(2)-bins_tau(1):(bins_tau(end)+bins_tau(end-1))/2];

[tau_sort,tau_idx] = sort(tau(idx_big));

figure(2); % keep this in, Matlab for some reason needs a separate figure handle to retain histogram data. 
T = histogram(tau_sort,bins_tau);T = T.Values;

close(figure(2))

A = Area(idx_big);          % original order of Volume vector
Atau(:) = A(tau_idx);       % sort Volume by increasing tau
A_at_tau = zeros(length(bins_tau)-1,1);

%%%% find total area in each bin range. T.Values is a vector describing
%%%% the number of glaciers in each bin. Since they are sorted by
%%%% increasing tau, we can then find total area in each
%%%% bin. 
% jj is index for A (ordered list of glacier volumes
% ii is bin index

jj = 1;
for ii = 1:(length(bins_tau) - 1);
    if T(ii) > 0
        A_at_tau(ii) = sum(Atau(jj:(jj+T(ii)-1)));
    else
        A_at_tau(ii) = 0;
    end
    jj = jj + T(ii);
end
jj = 1;

subplot 425; box on; hold on
histogram(tau(idx_big),[0:5:100],'facecolor',0.3*[1 1 1]); alpha(1)
%histogram(2*tau(idx_big),[0:5:130]+0.5,'facecolor','none','edgecolor',0.8*[0.5 0.5 1],'linewidth',1)

ylabel('count')
xlabel('\tau (years)')
set(gca,'xtick',[0:20:100],'ytick',[0:20:100])

subplot 427; box on; hold on
bar(centers_tau,100*A_at_tau/sum(A),1,'facecolor',0.3*[1 1 1])
%bar(centers_tau*2,100*A_at_tau/sum(A),1,'facecolor','none','edgecolor',0.8*[0.5 0.5 1],'linewidth',1)
set(gca,'xtick',bins_tau(1:5:end))
ylabel('% of total Area')
xlabel('\tau (years)')
ylim([0 30])
xlim([0 80]); set(gca,'Xtick',[0:10:80])