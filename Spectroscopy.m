 data=xlsread('GaAs_Peak.xlsx','A15:B2061');
 
 wavelength=data(:,1);
 intensity=data(:,2);
 energy=1239513.*wavelength.^-1;
 figure(1)
 plot(wavelength,intensity);
 title('Intensity Spectrum of GaAs in terms of wavelength');
 xlabel('Wavelength in nm');
 ylabel('Intensity in photon counts');
 hold on
 figure(2)
 plot(energy,intensity);
 title('Intensity Spectrum of GaAs in terms of wavelength');
 xlabel('Wavelength in nm');
 ylabel('Intensity in photon counts');
 hold on
 
 [~,maxindx]=max(intensity);
 peak_nm=wavelength(maxindx);
 peak_eV=1239513/peak_nm;
 
 FWHM_nm=fwhm(wavelength,intensity);
 FWHM_eV=abs(fwhm(energy,intensity));
 
 error_nmInP=FWHM_nm/2.3548;
 error_eVInP=FWHM_eV/2.3548;
 
 
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