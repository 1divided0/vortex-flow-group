%% Hauptprogramm: Transportgleichung_1d

function Transportgleichung_1d(u,c,r,T,options)
arguments
    u double % Transportgeschwindigkeit
    c double % Diffusionskonstante
    r double % Quelle (proportional zur Konzentration)
    T double % Simulationsendzeit
    options.gif_export logical = false % optionaler gif Export
end

% Diskretisierung des eindimensionalen Rechenraumes
[x,lbda] = Raumdiskretisierung();

% Startverteilung
phi0 = Startverteilung(x);

% Definition der in der Zeit zu integrierenden Funktion
% -> Zeitfunktion(phi,u,c,r,lbda)

% Lösung mit Zeitschrittverfahren
[t,phi] = ode45(@(t,phi) Zeitfunktion(phi,u,c,r,lbda),[0,T],phi0);

% Plot der zeitlichen Entwicklung
Plot(u,c,r,T,x,t,phi,gif_export=options.gif_export);
end


%% Subfunktion: Raumdiskretisierung

function [x,lbda] = Raumdiskretisierung()
% Rückgabe
% --------
% x:    äquidistantes 1d-Gitter
% lbda: diskrete Wellenzahlen

N    = 256; % Anzahl der Stützstellen
L    = 30;  % Gesamtlänge

x    = linspace(0,L,N)';
lbda = (2*pi/L)*(-N/2:N/2-1);
lbda = fftshift(lbda)'; % für die FFT benötigte Anordnung
end


%% Subfunktion: Startverteilung

function phi0 = Startverteilung(x)
% Rückgabe
% --------
% phi0: Startverteilung der Strömungsgröße als Hut-Funktion

N = length(x); % Anzahl der Stützstellen
a = 1;         % linke Position
b = 2;         % rechte Position
s = 0.8;       % Skalierungsfaktor
phimax = 0.1;  % Maximalwert

% linke Seite
phi1 = (x-a*ones(N,1)+s*ones(N,1)).^3.*(a*ones(N,1)-x+s*ones(N,1)).^3;
phi1(x<a-s|x>a+s) = 0;

% rechte Seite
phi2 = (x-b*ones(N,1)+s*ones(N,1)).^3.*(b*ones(N,1)-x+s*ones(N,1)).^3;
phi2(x<b-s|x>b+s) = 0;

% Superposition und Normierung
phi0 = (phi1+phi2)/max(phi1+phi2)*phimax;
end


%% Subfunktion: Zeitfunktion

function phi_dot = Zeitfunktion(phi,u,c,r,lbda)
% Rückgabe
% --------
% phi_dot: Zeitableitung der Strömungsgröße

phi_dot = real(ifft((c*(1i*lbda).^2-u*(1i*lbda)).*fft(phi))) + r*phi;
end


%% Subfunktion: Plot

function Plot(u,c,r,T,x,t,phi,options)
arguments
    u double         % Transportgeschwindigkeit
    c double         % Diffusionskonstante
    r double         % Quelle (proportional zur Konzentration)
    T double         % Simulationsendzeit
    x (:,1) double   % äquidistantes 1d-Gitter
    t (1,:) double   % Zeitwerte
    phi (:,:) double % Funktionswerte
    options.gif_export logical = false % optionaler gif Export
end

% Ploteinstellungen
line = plot(nan,nan);
line.XData = x;
ylim([-0.01 .12]);
ylabel('$\Phi$','Interpreter','latex');
xlabel('$x$','Interpreter','latex');
set(gcf,'color','w');
grid on;

% Darstellung von 100 Zeitschritten
for i = 1:round(length(t)/100):length(t)
    line.YData = phi(i,:);
    title([
        '$u=' num2str(u) '$, ' ...
        '$c=' num2str(c) '$, ' ...
        '$r=' num2str(r) '$, ' ...
        '$t=' num2str(t(i),'%.2f') '$' ...
    ],'Interpreter','latex');
    
    % optionaler gif Export
    if options.gif_export
        file = ['u' num2str(u) ...
            '_c' num2str(c) ...
            '_r' num2str(r) ...
            '_T' num2str(T) '.gif'];
        if i==1
            exportgraphics(gcf,file);
        else
            exportgraphics(gcf,file,'Append',true);
        end
    else
        drawnow;
        pause(.05);
    end
end
end
