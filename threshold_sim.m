% threshold_sim.m

% Program is for n alleles.  
% It calculates exact frequencies of all alleles over a range of times for particular initial
% frequencies and a particular set of fitness effects.

% It then sets to zero all frequencies below a threshold of 0.05. All
% non-zero frequencies are rescaled (divided by their sum), so they sum to
% 1 and these frequencies are also plotted.

clear all
close all
clc

% Parameters
n=10;                    % No. alleles
T=1e3;                  % Final time
s0=0.01;                % Strength selection
Threshold=0.05;

% Useful
F=ones(n,1);

% Initialisation
X=zeros(n,T+1);         % Contains exact frequency trajectories

A= [0.0037   -0.0050    0.0154   -0.0014    0.0004    0.0002   -0.0039    0.0155   -0.0070   -0.0120
   -0.0050    0.0041    0.0088    0.0047    0.0101   -0.0012    0.0129    0.0014    0.0003    0.0046
    0.0154    0.0088   -0.0172    0.0064    0.0027    0.0027   -0.0019   -0.0034   -0.0000   -0.0020
   -0.0014    0.0047    0.0064   -0.0028   -0.0045   -0.0005    0.0032   -0.0036    0.0079   -0.0068
    0.0004    0.0101    0.0027   -0.0045   -0.0029    0.0024   -0.0007   -0.0016    0.0118    0.0021
    0.0002   -0.0012    0.0027   -0.0005    0.0024   -0.0067    0.0054    0.0085   -0.0095   -0.0031
   -0.0039    0.0129   -0.0019    0.0032   -0.0007    0.0054   -0.0009   -0.0030   -0.0024   -0.0012
    0.0155    0.0014   -0.0034   -0.0036   -0.0016    0.0085   -0.0030    0.0054   -0.0046   -0.0067
   -0.0070    0.0003   -0.0000    0.0079    0.0118   -0.0095   -0.0024   -0.0046    0.0011    0.0103
   -0.0120    0.0046   -0.0020   -0.0068    0.0021   -0.0031   -0.0012   -0.0067    0.0103   -0.0232];


X0=[0.0271
    0.1077
    0.1874
    0.0833
    0.0447
    0.0163
    0.1931
    0.1248
    0.0781
    0.1376];

X(:,1)=X0;
x=X0;
for t=1:T
    V=diag(x)-x*x';
    D=V*A*x/(1+x'*A*x);
    xp=x+D;
    xp=xp/sum(xp);                  % Forces normalisation every generation (may not be necessary)
    x=xp;
    X(:,t+1)=x;                     % 
end

I=X>Threshold;
Xth=(X.*I)./(F*sum(X.*I));          % Contains thresholded frequencies

l_width = 3;

plot(0:T,Xth,'color', [0, 0.447, 0.741],'linewidth',l_width)
hold on
plot(0:T,X, ':','color', [0.85, 0.325, 0.098],'linewidth',l_width)
axis([0,1e3,0,0.8])
fs=20;
set(gca,'linewidth',l_width,'fontsize',fs*0.9)
set(gca,'xtick',(0:0.2:1)*T)
set(gca,'ytick',0:0.2:1)
set(gca, 'box', 'off')

xlabel('generations','fontsize',fs);
ylabel('allele frequencies','fontsize',fs);



lgd=legend('{  }5% threshold{  }','','','','','','','','','','{  }exact{  }', ...
           'FontSize', fs, 'Location', 'northwest');

legendLines = findobj(lgd, 'Type', 'Line');
set(legendLines, 'LineWidth', l_width); % Set desired linewidth for legend entries

ounits = get(gcf,'Units');
set(gcf,'Units','normalized','OuterPosition',[0 0 1 1],'Units',ounits)
orient landscape
print -dpdf -bestfit FigThreshold
