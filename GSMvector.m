clear all
close all
clc

%% ========================================================================
% ECI -> GSM
%
% Input  : satellite position (Altitude, longitide, latitude) [km]
% Output : r_gsm [km]
%
% Example:
% GEO satellite at 42164 km from Earth center
% ========================================================================

%% date and time
utc = datetime(2014,7,19,10,0,0);
JD = juliandate(utc);

%% satellite position (in ECEF and convert to ECEF)

lon = 75;              % deg
lat= 0; 
Re= 6378.137;          % km - Earth radius
altitude= 35786;       % km - Satellite altitude
r = Re+altitude;

% satellite position in ECEF
x_ecef = r*cosd(lat)*cosd(lon);
y_ecef = r*cosd(lat)*sind(lon);
z_ecef = r*sind(lat);

r_ecef = [x_ecef;y_ecef;z_ecef];

gmst = gmstangle(JD);

R_ecef2eci = Recef2eci(gmst);

% satellite position in ECI

r_eci=R_ecef2eci*r_ecef

%% Earth dipole

Dipole_Latitude  = 80.65;       %deg
Dipole_Longitude = -72.68;    %deg

lat = deg2rad(Dipole_Latitude);
lon = deg2rad(Dipole_Longitude);

m_ecef=Mecef(lon, lat);  % Earth dipole in ECEF

m_ecef = m_ecef/norm(m_ecef);

m_eci = R_ecef2eci*m_ecef;  % convert earth dipole in ECI

%% calculate sun position in ECI frame

[n,lambda,eps]=sunangles(JD);
lambda = deg2rad(lambda);
eps = deg2rad(23.439 - 0.0000004*n);    % Obliquity of the ecliptic (Earth's axial tilt) [radians]

% Sun unit vector in the Earth-Centered Inertial (ECI) frame
sun_eci = [ ...
    cos(lambda);
    cos(eps)*sin(lambda);
    sin(eps)*sin(lambda)];

sun_eci = sun_eci/norm(sun_eci);


%% GSM axes:

x_gsm = sun_eci

%y_gsm = cross(m_eci,x_gsm)

y_gsm = cross(m_eci,sun_eci);

y_gsm = y_gsm/norm(y_gsm);

z_gsm = cross(x_gsm,y_gsm);

z_gsm = z_gsm/norm(z_gsm);

%% Direction Cosine Matri


C_gsm_eci = [
    x_gsm';
    y_gsm';
    z_gsm'];
% Satellite position from ECI to GSM

r_gsm = (C_gsm_eci * r_eci)/Re
norm(r_gsm)*Re

%% B field convert from gsm to ECI

Bxgsm= 0.19;  Bygsm= 0.0023; Bzgsm= 20.48; 

B_gsm=[Bxgsm;Bygsm;Bzgsm];
NBgsm=norm(B_gsm)
B_eci = C_gsm_eci' * B_gsm;
NB_eci=norm(B_eci)

%% theoratical B = KB/r^3
KB=8e15;
B_theory=KB/r^3

%% Functions
%% Calculating GMST

function GMST=gmstangle(JD)

T= ( JD - 2451545.0 ) / 36525.0;
GMST= 280.46061837 + 360.98564736629 * (JD - 2451545.0)+ 0.000387933 * T * T- (T * T * T) / 38710000.0;
     GMST = mod(GMST,360);   

GMST = mod(GMST,360);

if GMST < 0
    GMST = GMST + 360;
end

GMST = deg2rad(GMST);

end

%% convert ECEF yo ECI

function R_ecef2eci=Recef2eci(gmst)

R_ecef2eci = [
 cos(gmst) -sin(gmst) 0;
 sin(gmst)  cos(gmst) 0;
 0          0         1];

end

%% Earth dipole axis in ECEF

function m_ecef=Mecef (lon, lat)

m_ecef = [ ...
    cos(lat)*cos(lon);
    cos(lat)*sin(lon);
    sin(lat)];
end

%% calculate sun angles 
function [n,L,g,lambda,eps]=sunangles(JD)

n = JD - 2451545.0;                  % Number of days elapsed since J2000 epoch
L = mod(280.460 + 0.9856474*n,360);  % Mean longitude of the Sun (degrees)

g = mod(357.528 + 0.9856003*n,360);  % Mean anomaly of the Sun (degrees)
g = deg2rad(g);

lambda = L + 1.915*sind(rad2deg(g)) ...  % Apparent ecliptic longitude of the Sun (degrees)
           + 0.020*sind(2*rad2deg(g));

lambda = deg2rad(lambda);
eps = deg2rad(23.439 - 0.0000004*n);    % Obliquity of the ecliptic (Earth's axial tilt) [radians]
end