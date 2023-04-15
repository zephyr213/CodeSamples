function vizplot(dt,tstep,Temp,elements,varargin)% plot T evolution of individual element
% elements - vector of the elements of interest
% varargin - optional time step start and end (default mode, entire
% simulation)

c2k = 273.15;

if nargin ==4
    tstart = 1;
    tend = tstep;
elseif nargin == 6
    tstart = varargin{1};
    tend = varargin{2};
end
tvec = tstart:tend;
tvec = tvec*dt;
%elements = [1,8,41,48]; % index of elements of interest
plot(tvec',Temp(tstart:tend,elements)-c2k,'-');
xlabel('Time [s]')
ylabel('Temperature [\circC]')
% legends
s1 = strings(1,length(elements));
for i = 1:length(elements)
    s1(i) = "Element: " + num2str(elements(i));
end
legend(s1,'Location','best')
%aaa
grid on

end