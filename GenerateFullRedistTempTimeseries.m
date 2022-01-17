% I need to generate the 362x292x64x240 timeseries of redistributed
% temperature. This is going to be a pain in the ass. To do this, I first
% need to create a drift controlled Cnat change file.

load('/home/ct/MATLAB/BIGMATFILES/DIC_RAD.mat');
load('/home/ct/MATLAB/BIGMATFILES/DIC_CTR.mat');

dCnat = DIC_RAD - DIC_CTR;
clear DIC_CTR
% Now we need to get dCanth too.

load('/home/ct/MATLAB/BIGMATFILES/DIC_COU.mat');
dCanth = DIC_COU - DIC_RAD;
clear DIC_COU DIC_RAD

% Now load kappa_r values in. This .mat file also contains polynomial fits
% which were used as a sanity check.
% % This file is generated in
% /home/ct/MATLAB/2020/November/TestingShortAndLongTermGlobalOcean.m
load('/home/ct/MATLAB/2020/November/KappaRandPolyfitTest.mat','kappa_r','ecc');
kappa_r = kappa_r .* ecc;

% Now, we have our kappa_r table. We need to use dCanth and dCnat
% to determine gamma for every year now. 

%%  The following 5 lines generate the globalGammaVals list. Note everything
%  is turned into a single in order to use a parfor loop: otherwise we run
%  out of memory on my machine and crash matlab.
dCnat = single(dCnat);
dCanth = single(dCanth);
%
gammaValuesTemp = calculateGammaValue(dCnat,dCanth,kappa_r,1);
gammaValuesCarbon = calculateGammaValue(dCnat,dCanth,kappa_r,0);
save globalGammaValsYEARLYCTR.mat gammaValues                     

%% Load it instead of generating it again, if rerunning this script
load globalGammaValsYEARLYCTR.mat

%
load grid_data_areas_masks.mat dV;

dV = repmat(dV,[1 1 1 240]);

gblCanthInv = squeeze(nansum(nansum(nansum(dCanth .* dV,1),2),3));


%
figure(1); clf;
figure(1); hold on
plot(1860:2099,gammaValues,'Color','k');
plot(1860:2099,smooth(gammaValues,10,'rloess'),'Color','k','Linewidth',2);
% Probably also want to plot the global inventory of Canth since the gamma
% value swings about a lot for a while
ylabel('$\gamma$');
title('$\gamma$ Factor: Repartitioning C$_{anth}$ into C$_{sat}$');
ax1= gca;
ax1.YLim = [-0.01 0.25];

legend({'Raw Data', 'Smoothed (Decadal, rloess)'},'interpreter','latex')

yyaxis right; hold on;
%
ax2 = gca;
plot(1860:2099,12*gblCanthInv/1e18,'Color','r','Linewidth',2);
ylabel('Global C$_{anth}$ Inventory / PgC')
ax2.YColor =  'red';

% Save this

smoothedGammaValues = smooth(gammaValues,10,'rloess');
%save globalGammaValsYEARLYCTR.mat gammaValues smoothedGammaValues


% Now build the adjusted Cnat field


AdjCnat = dCnat + permute(repmat(smoothedGammaValues,[1 362 292 64]), [2 3 4 1]) .* dCanth;

%% And save it off

save ADJ_CNAT_YEARLYCTR.mat AdjCnat -v7.3

%% Finally, we can build the temperature fields we are after

clear dCnat dCanth
load ADJ_CNAT_YEARLYCTR.mat
%
dTMPr = AdjCnat .* repmat(kappa_r,[1 1 1 240]);

% We should run some sanity checks: global integral is zero, scale grows with time (probably)

load grid_data_areas_masks.mat dV;


gblTMPrInv = 3850*1025*squeeze(nansum(nansum(nansum(dTMPr .* dV,1),2),3));
gblTMPrInvAb =  3850*1025*squeeze(nansum(nansum(nansum(abs(dTMPr) .* dV,1),2),3));


%
figure(2); clf;
figure(2); hold on; grid on;

plot(1860:2099,gblTMPrInv/1e21,'Color','k','Linewidth',2);
plot(1860:2099,gblTMPrInvAb/1e21,'Color','r','Linewidth',2);

title('Globally Integrated Redistributed Heat / ZJ');
xlabel('Year');
ylabel('OHC (Redistributed) / ZJ');

legend({'$\int \Delta \Theta_r dV$','$\int | \Delta \Theta_r| dV$'},'interpreter','latex');

%% Looks about right if you ask me. Now build the excess temperature field.
clear AdjCnat

load('/home/ct/MATLAB/BIGMATFILES/TMP_COU.mat');
load('/home/ct/MATLAB/BIGMATFILES/TMP_CTR.mat');

dTMP = TMP_COU - TMP_CTR;

clear TMP_COU TMP_CTR

dTMPe = dTMP - dTMPr;

save TEMP_REDIST_YEARLY_CTR.mat dTMPr -v7.3
save TEMP_EXCESS_YEARLY_CTR.mat dTMPe -v7.3
save TEMP_CHANGE.mat dTMP  -v7.3
%clear
%% 
load TEMP_EXCESS_YEARLY_CTR.mat
load TEMP_REDIST_YEARLY_CTR.mat
load grid_data_areas_masks.mat dV;
%%
gblTMPeInv = 3850*1025*squeeze(nansum(nansum(nansum(dTMPe .* dV,1),2),3));


gblTMPrInv = 3850*1025*squeeze(nansum(nansum(nansum(dTMPr .* dV,1),2),3));
gblTMPrInvAb =  3850*1025*squeeze(nansum(nansum(nansum(abs(dTMPr) .* dV,1),2),3));



%%
figure(3); clf;
figure(3); hold on; grid on;

plot(1860:2099,gblTMPeInv/1e21,'Color','b','Linewidth',2);
plot(1860:2099,gblTMPrInv/1e21,'Color','k','Linewidth',2);
plot(1860:2099,gblTMPrInvAb/1e21,'Color','r','Linewidth',2);

title('Globally Integrated Redistributed and Excess Heat / ZJ');
xlabel('Year');
ylabel('OHC (Redistributed) / ZJ');

legend({'$\int \Delta \Theta_e dV$','$\int \Delta \Theta_r dV$','$\int | \Delta \Theta_r| dV$'},'interpreter','latex');




% This all looks like it has worked pretty well. In summary, I'm actually
% pretty happy with how its panned out.




