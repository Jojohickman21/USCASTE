
epoch = toJulianDate('2-Sep-2024',8);
start = epoch;

drawEarthMap;
hold on;
mu = 3.986*10^5;

inc = asind(sqrt(4/5));

semiMajor = (mu *( 43082.1 /(2*pi) )^2)^(1/3);
%Programming the corrected SemiMajor Axis
R0 = 6378;
ecc = 0.745;
J2 = 0.00108263; % Zonal coefficient for Earth
n = sqrt( mu / semiMajor^3 );  % Mean motion (rad/s, uncorrected)

dOmegaByDt = -3*n * J2 * R0^2 * cosd(inc) / ...
    ( 2*semiMajor^2 * (1 - ecc^2)^2 ); % In rad/sec

dOmegaByDtDegDay = dOmegaByDt * 180/pi * 86164; % In degrees/day

EarthRotationRateInertialDegDay = 360; % In degrees per sidereal day
EarthRotationRateApparentDegDay = EarthRotationRateInertialDegDay ...
    - dOmegaByDtDegDay;

% Want spacecraft to orbit twice per apparent Earth rotation
nCorrectedDegDay = 2 * EarthRotationRateApparentDegDay;
nCorrected = nCorrectedDegDay * pi/180 / 86164; % In rad/sec

% n = sqrt( mu / a^3 )
aCorrected = ( mu / nCorrected^2 )^(1/3);



drawGroundTrack1(aCorrected,0.745,-10,inc,-90,0,epoch,start,4);

function drawGroundTrack1( a, ecc, Omega, inc, omega, theta, ...
                                epoch, start, nOrbits )
      %DRAWGROUNDTRACK Draw the ground track of an Earth satellite with
      %specified classical orbital elements as follows:
      % a: Semimajor axis, km
      % ecc: Eccentricity
      % Omega: RAAN (degrees)
      % inc: inclination (degrees)
      % omega: Argument of periapsis (degrees)
      % theta: True anomaly (degrees)
      % epoch: Time at which orbital elements are specified (Julian date)
      % start: Time at which to draw the ground track (Julian date)
      % nOrbits: Number of orbital passes to draw
      
      %Earth Gravitational Parameter defined in km^3 * s^-2
      mu = physicalConstant('muEarth');

      elems = [a,ecc,Omega, inc, omega, theta];
      
      period = 2*pi * sqrt( a^3 / mu );

      orbitSeconds = period*nOrbits;
      
      
      deltaDays = orbitSeconds/86400;
      endDate = start + deltaDays;
     
      
      %Creating an array of times from start to end, with an interval of
      %1 minute in julian date
      time = start:(1/1440):endDate;

      
      rECEFVec = findECEFLocation(elems,epoch,time);
      La = zeros(size(rECEFVec,1),1);
      Lo = zeros(size(rECEFVec,1),1); 

      for i = 1:size(rECEFVec,1)
        La(i) = asind( rECEFVec(i,3) / norm(rECEFVec(i,:)) );
        Lo(i) = atan2d( rECEFVec(i,2), rECEFVec(i,1) );
      end
         
      plot( Lo,La, 'b', 'DisplayName', "Orbital"); 
            

    
end
