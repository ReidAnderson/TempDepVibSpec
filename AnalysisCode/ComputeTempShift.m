T = 250:350;
%T = [T 500 750 1000];
allT = zeros(length(T),35000);
% ACN-1milIR-20iter-1run-208
%fileBase = 'ACN/ACN-RamanInt100k-Raman-';
runName = '9-5-2016-IRInt1mil';
molAbbrev = 'ACN';
fileBase = ['../Results/' molAbbrev '/' runName '/TempDepVibSpec/Raman/' runName '-'];
cmap = colormap(jet(length(T)));
allMax = zeros(length(T),1);
allSums = zeros(length(T),1);
figure
for i = 1:length(T)
    i
    hold on
    load([fileBase num2str(T(i)) '.mat']);
    allSums(i) = sum(I_T);
    I_T = I_T./max(I_T);
%     plot(I_T(20000:25000));
plot(I_T,'Color',cmap(i,:));
    
    [v d] = max(I_T(20000:25000));
    allMax(i) = d;
%      pause;
end
hold off;
%set(gca,'XTickLabel',str2num(get(gca,'XTickLabel')).*1000);
figure,plot(allMax);