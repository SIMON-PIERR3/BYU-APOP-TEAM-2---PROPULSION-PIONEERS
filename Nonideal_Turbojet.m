clc;

% ALTITUDE = 22000 ft, 
mc = 192; %lbm/s
H = 17900; %Btu/lbm
PID = .93;%diffuser pressure ratio
PIC = 17;%compressor pressure ratio
PIB = .93; %burner pressure ratio
PIAB = .97;%a/b pressure ratio
Ma = .88; %mach #
g = 32.17; %lbm/slug
BTU = 778.2; %lbf-ft/BTU
R = 53.35; %imperial...lbf

nC = .89;
nB = .95;
nT = .87;
nS = .997; %shaft
nN = .97;
nAB = .91;

Pa = 6.242; %psia -found by interpolation
Ta = 440.4; %R    - "
ua = sqrt(ga*R*Ta*g)*Ma;

ga = 1.378;
gd = ga;
gc = gd;
gt = gc;
gn = gt;

for i = 1:4
    Tta = Ta*(1+(gd-1)/2*Ma^2);
    cpd = .2269807*exp(.000097247*(Tta+Ta)/2);
    gd = cpd/(cpd-.068559);
end
Pta = Pa*(1+(gd-1)/2*Ma^2)^(gd/(gd-1));

Pt2 = PID*Pta;
Tt2 = PID^((gd-1)/gd)*Tta;

for i = 1:4
    Tt3 = Tt2*((PIC^((gc-1)/gc)-1)/nC+1);
    cpc = .2269807*exp(.000097247*(Tt3+Tt2)/2);
    gc = cpc/(cpc-.068559);
end
Pt3 = Pt2*PIC;

Pt4 = Pt3*PIB;
Tt4 = 2350; %R

cpb = .2269807*exp(.000097247*(Tt3+Tt4)/2);
mf = (mc*cpb*(Tt4-Tt3))/(nB*H-cpb*Tt4);

cpt = cpc;
for i = 1:4
    Tt5 = Tt4-(mc*cpc*(Tt3-Tt2))/(nS*mc*cpt*(1+mf/mc));
    cpt = .2269807*exp(.000097247*(Tt5+Tt4)/2);
    gt = cpt/(cpt-.068559);
end
Pt5 = Pt4*(1-(1-Tt5/Tt4)/nT)^(gt/(gt-1));

Pt6 = PIAB*Pt5;
Tt6 = 2980; %R
cpab = .2269807*exp(.000097247*(Tt5+Tt6)/2);
mfab = (mc+mf)*cpab*(Tt6-Tt5)/(H*nAB-cpab*Tt6);
cpn = .2269807*exp(.000097247*(Tt6));
gn = cpn/(cpn-.068559);
P_star = Pt6*(1+(1-gn)/(nN*(1+gn)))^(gn/(gn-1));

T8 = 0;
if P_star > Pa
    P8 = P_star;
    M8 = 1;
    T8 = 2*Tt6/(1+gn);
else
    P8 = Pa;
end

u8 = sqrt(2*cpn*(Tt6-T8)*g*778.16);
m8 = mc+mf+mfab;
rho = P8*144/(R*T8*g);
A8 = m8*144/(rho*u8*g);
F = m8*u8/g-mc*ua/g+A8*(P8-Pa)