 dataInP=xlsread('InP_Peak.xlsx','A15:B2061');
 dataInGaP=xlsread('InGaP_Peak.xlsx','A15:B2061');
 
 intime_InP=0.2;
 intime_InGaP=0.3;
 
 wavelengthInP=dataInP(:,1);
 intensityInP=dataInP(:,2)./0.2;
 energyInP=1239513.*wavelengthInP.^-1;
 figure(1)
 plot(wavelengthInP,intensityInP);
 hold on
 figure(2)
 plot(energyInP,intensityInP);
 hold on
 
 [~,maxindxInP]=max(intensityInP);
 peak_nmInP=wavelengthInP(maxindxInP);
 peak_eVInP=1239513/peak_nmInP;
 
 FWHM_nmInP=fwhm(wavelengthInP,intensityInP);
 FWHM_eVInP=abs(fwhm(energyInP,intensityInP));
 
 error_nmInP=FWHM_nmInP/2.3548;
 error_eVInP=FWHM_eVInP/2.3548;
 
 wavelengthInGaP=dataInP(:,1);
 intensityInGaP=dataInGaP(:,2)./0.3;
 energyInGaP=1239513.*wavelengthInGaP.^-1;
 figure(1)
 plot(wavelengthInGaP,intensityInGaP);
 title('Intensity Spectrum in terms of wavelength');
 xlabel('Wavelength in nm');
 ylabel('Intensity in photon counts per second');
 legend('InP Sample','InGaP Sample');
 figure(2)
 plot(energyInGaP,intensityInGaP);
 title('Intensity Spectrum of InGaP in terms of energy');
 xlabel('Energy in meV');
 ylabel('Intensity in photon counts per second');
 legend('InP Sample','InGaP Sample');
 
 [~,maxindxInGaP]=max(intensityInGaP);
 peak_nmInGaP=wavelengthInGaP(maxindxInGaP);
 peak_eVInGaP=1239513/peak_nmInGaP;
 
 FWHM_nmInGaP=fwhm(wavelengthInGaP,intensityInGaP);
 FWHM_eVInGaP=abs(fwhm(energyInGaP,intensityInGaP));
 
 error_nmInGaP=FWHM_nmInGaP/2.3548;
 error_eVInGaP=FWHM_eVInGaP/2.3548;
 
 ev_delta=peak_eVInP-peak_eVInGaP;
 error_ev_delta= sqrt(error_eVInGaP^2+error_eVInP^2);
 
 p=[603 512 ev_delta];
 r=roots(p);
 
 perrorp=[603 512 ev_delta+error_ev_delta];
 perrorr=roots(perrorp);
 
 nerrorp=[603 512 ev_delta-error_ev_delta];
 nerrorr=roots(nerrorp);
 
 Ga_concentration=r(2,1);
 Ga_concentration_range=[perrorr(2,1) nerrorr(2,1)];
 
function width = fwhm(x,y)
% function width = fwhm(x,y)
%
% Full-Width at Half-Maximum (FWHM) of the waveform y(x)
% and its polarity.
% The FWHM result in 'width' will be in units of 'x'
%
%
% Rev 1.2, April 2006 (Patrick Egan)
y = y / max(y);
N = length(y);
lev50 = 0.5;
if y(1) < lev50                  % find index of center (max or min) of pulse
    [garbage,centerindex]=max(y);
    Pol = +1;
    %disp('Pulse Polarity = Positive');
else
    [garbage,centerindex]=min(y);
    Pol = -1;
    %disp('Pulse Polarity = Negative');
end
i = 2;
while sign(y(i)-lev50) == sign(y(i-1)-lev50)
    i = i+1;
end                                   %first crossing is between v(i-1) & v(i)
interp = (lev50-y(i-1)) / (y(i)-y(i-1));
tlead = x(i-1) + interp*(x(i)-x(i-1));
i = centerindex+1;                    %start search for next crossing at center
while ((sign(y(i)-lev50) == sign(y(i-1)-lev50)) & (i <= N-1))
    i = i+1;
end
if i ~= N
    Ptype = 1;  
    %disp('Pulse is Impulse or Rectangular with 2 edges');
    interp = (lev50-y(i-1)) / (y(i)-y(i-1));
    ttrail = x(i-1) + interp*(x(i)-x(i-1));
    width = ttrail - tlead;
else
    Ptype = 2; 
    %disp('Step-Like Pulse, no second edge');
    ttrail = NaN;
    width = NaN;
end
end