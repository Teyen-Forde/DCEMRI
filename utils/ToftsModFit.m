function pars =ToftsModFit(SI_aif,SI_toi,HTres,LTres,TR,FA,Ktrans_init,ShowPlot)
%   pars = fitdcemri(Ctoi,Cp,time,x0,lb,ub,'Tofts')
%    /::INPUTS::\
%       Ctoi: a column vector of Contrast Agent(CA) Concentration over time
%           for the Tissue of Interest.
%           NOTE: Must be concentration, not delta R1 (deltaR1 = r1*Conc
%           where r1 is the relaxivity of the CA)
%       Cp: an Arterial Input Function (AIF) for the tissue in question in
%           the form of a column vector equal in size to toi
%       time: a column vector of time in minutes corresponding to toi and
%           Cp
%       x0: a vector of best guesses for parameters from the non-linear
%           least squared fitting (Ex:[Ktrans_guess,kepTOI_guess])
%       lb: vector of lowerbounds for the parameters (same form as above)
%       ub: vector of upperbounds for the parameters (same form as above)
%       'Tofts': This string indicates that you want the non-linear Tofts
%           Model (versus the non-linear RRM)
%    /::OUTPUT::\
%       pars: a vector with the following components
%           pars(1)=  Ktrans
%           pars(2)=  kep for the tissue of interest, kepTOI
%           pars(3)=  rsquare for the fitting
% Authors:
% Joey DeGrandchamp                 Julio Cardenas-Rodriguez
% University of Arizona             University of Arizona
% jdegrandchamp@email.arizona.edu   cardenaj@email.arizona.edu
%
%                       www.cardenaslab.org/resources
% v2.0 09/24/2015

len=length(SI_aif);
if size(SI_aif,1) == 1
    SI_aif = SI_aif';
end
if size(SI_toi,1) == 1
    SI_toi = SI_toi';
end


%% setup some fixed variables in this model
tlist=(1/2:1:len)'*HTres;
% tlist=(0:len-1)'*Tres;
T10_blood=1.4;
T10_kidney=1.2;
r1=4.5;% the same for blood, EV and IV
Hct_large=0.41;


%% Convert signal intensity to concentration
Cb_aif=SI2C(SI_aif,TR,FA,T10_blood,r1);
Cp_aif=Cb_aif/(1-Hct_large);
%tissue concentration
%interpolate into higher resolution
tlow = (1/2:1:length(SI_toi))'*LTres;
SI_toi = interp1(tlow,SI_toi,tlist,'linear','extrap');
Ct=SI2C(SI_toi,TR,FA,T10_kidney,r1);

%% Fitting using curve fit

% Set options for non-linear LSQ fitting
S=optimset; S.Algorithm='trust-region-reflective'; S.Display='off';
S.TolFun=1e-6; S.TolX=1e-6;
S.MaxIter=1000;

% intial and bound
x0 = [Ktrans_init 0.5 4 1];
lb = [1e-6  1e-6 1e-6 1e-6];
ub = [1  1 10  10];
% Perform the fitting (the function "rrm" is in the nested functions below)
[B,res_norm] = lsqcurvefit(@ToftsMod,x0,[tlist,Cp_aif],Ct,lb,ub,S);

% Store each parameter (see reference)
pars(1,1)=B(1);    %Ktrans
pars(1,2)=B(2);    %vb
pars(1,3)=B(3);    %Tg
pars(1,4)=B(4);    %Delta(delay)
pars(1,5)=res_norm;

% pred=tofts(B,[time,Cp]);


if (ShowPlot)
%     figure
    [Ct_fit, Ct_IV,Ct_EV]= ToftsMod(B,[tlist,Cp_aif]);
    %[Cfit, Cart, Ctub]= ThreeCompartment(xFit,xdata1);
    plot(tlist,Ct_fit,'-k');
    hold on
    plot(tlist,Ct_IV,'--r');
    plot(tlist,Ct_EV,'--k');
    plot(tlist,Ct,'*r');
    plot(tlist,Cp_aif,'+b');
    plot(tlist,abs(Ct_fit-Ct),'+m');
    hold off
    legend('Fit','IV','EV','parenchyma','AIF','residual')
    title(['K21 = ' num2str(B(1)) ' ' 'vb = ' num2str(B(2))  ' resnorm = ' num2str(res_norm) ])
    xlabel('Time (s)')
    ylabel('Concentration (mM)')
end


end



function [ C ] = SI2C( SI,TR,FA,T10,r1 )
% This funtion takes a signal SI(t) and returns the concentration C(t)
% The function assumes at t=0 there is no contrast agent.
% SI is the signal intensity
% TR is the TR in sec
% FA is the Flip Angle in radians
% T10 is the T1 of the tissue before contrast agent (i.e. at t=0) in sec
% r1 is the relaxivity of contrast agent in (s^-1)*(mM^-1)
%
% RETURN VALUE => C in mM or C = zeros(size(SI)) if threshold value is
% exceeded

C_limit = 20;   % Threshold value for estimation in mM

R10 = 1/T10;

S0 = SI(1)*(1-exp(-R10*TR)*cos(FA))/((1-exp(-R10*TR))*sin(FA));

R1_limit = R10 + r1 * C_limit;
SI_limit = S0*((1-exp(-R1_limit*TR))*sin(FA))./(1-exp(-R1_limit*TR)*cos(FA));

if (max(SI)>SI_limit) %#ok<BDSCI>
    C = zeros(size(SI));
    warning('exceed maximum SI intensity');
else
    R1 = 1/TR*( log(SI*cos(FA)/(S0*sin(FA))-1) - log(SI/(S0*sin(FA))-1));
    %     figure(100);hold on; plot(1./R1,'-');
    
    C = (R1-R10)/r1;
end

end

function [Ct,Ct_IV,Ct_EV]=ToftsMod(beta,X)
%const
Hct_small=0.0;
Kep=0;

time = X(:,1);
time_res=time(2)-time(1);%in s
% time_new=(time(1):time_res:time(end))';
Cp_aif   = X(:,2);
% Cp_aif= interp1(time,Cp_aif,time_new,'linear','extrap');

Ktrans = beta(1)/60;%tranform into seconds unit
vb     = beta(2);
Tg     = beta(3);
Delta  = beta(4);
% Cd = Ktrans*convolution(Cp,exp(-kepTOI.*time));
%% VIRF, delayed exponential
newtime_res=time_res;
newtime=(time(1):newtime_res:time(end))';
Cp_aif=interp1(time,Cp_aif,newtime-Delta,'linear','extrap');
% figure(2), hold on; plot(Cp_aif); hold off;
VIRF=1/Tg*exp(-(newtime-time(1))/Tg);
% 
%Normalize to ensure having unit area under curve
VIRF=VIRF/(sum(VIRF)*newtime_res);
Cp_toi=convolution(Cp_aif,VIRF)*newtime_res;
Ct_EV=Ktrans*convolution(Cp_toi,exp(-Kep*time))*time_res;
Ct_IV=vb*(1-Hct_small)*Cp_toi;
Ct=Ct_IV+Ct_EV;

% Ct = interp1(newtime,Ct,time,'linear','extrap');
% Ct_IV = interp1(newtime,Ct_IV,time,'linear','extrap');
% Ct_EV = interp1(newtime,Ct_EV,time,'linear','extrap');

end

function c = convolution(a, b, shape)
%CONVOLUTION MODIFIED BY JULIO CARDENAS, MAY OF 2011.
%   SAME THAN CONV BUT RETURN A VECTOR FROM a(1) to a(end), not the central
%   section as described for the usual convolution function.
%

if ~isvector(a) || ~isvector(b)
    error(message('MATLAB:conv:AorBNotVector'));
end

if nargin < 3
    shape = 'full';
end

if ~ischar(shape)
    error(message('MATLAB:conv:unknownShapeParameter'));
end

% compute as if both inputs are column vectors
[rows,~]=size(a);
c = conv(a(:),b(:),shape);
c=c(1:rows);

% restore orientation
if shape(1) == 'f'
    if length(a) > length(b)
        if size(a,1) == 1 %row vector
            c = c.';
        end
    else
        if size(b,1) == 1 %row vector
            c = c.';
        end
    end
else
    if size(a,1) == 1 %row vector
        c = c.';
    end
end

end