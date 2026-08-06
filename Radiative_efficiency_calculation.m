clear all
% Calculate the instantaneous or adjusted radiative efficiency from a cross-section spectrum
% Author: K. Le Bris

%Enter the cross-section spectrum  (first column: wavenumber in cm-1, second column: cross-section in cm^2/molecule) 
   file1 = load('CHL315.dat'); 
% Enter the instantaneous or adjusted Pinnock curve (first column: wavenumber in cm-1, second column: radiative forcing in mW m^-2 cm (per 10e-18 cm^2 per molecule) 
Pinnock= load('NewPinnock_1cm.dat'); 

x=file1(:,1);



index=find(x>0);x=x(index);data =file1(:,2);
xmin=round(min(x));
xmax=round(max(x));
if xmax>3000, xmax=3000; end;
 
for i = xmin:xmax-1
    xlow=i;
    xhigh=i+1;
    xi=find(x>xlow & x<xhigh);
    v(i+1-xmin)=i+0.5;
    
        yj=mean(data(xi));
     
        y(i+1-xmin)=yj;
     
   
end

ind=find(Pinnock(:,1)>=xmin&Pinnock(:,1)<xmax);
    RE=y'.*Pinnock(ind,2)*1E15; % radiative efficiency in mW m^-2 per molecule
    InsRE=trapz(v',RE)

