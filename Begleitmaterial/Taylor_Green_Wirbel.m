%% Visualisierung des Taylor-Green-Wirbels

function Taylor_Green_Wirbel(options)
arguments
    options.nu    double = 0   % Kinematische Viskosität
    options.t_end double = 5   % Endzeit
    options.t_nr  uint16 = 100 % Anzahl der Zeitschritte
    options.x_nr  uint16 = 100 % Auflösung in x-Richtung
    options.y_nr  uint16 = 20  % Auflösung in y-Richtung
end

% Gitter
x     = linspace(0,2*pi,options.x_nr);
y     = linspace(0,2*pi,options.y_nr);
dx    = x(2)-x(1);
dy    = y(2)-y(1);
[X,Y] = meshgrid(x,y);

% Lösung
U     = @(t) sin(X).*cos(Y).*exp(-2*options.nu*t);
V     = @(t) -cos(X).*sin(Y).*exp(-2*options.nu*t);
Psi   = @(t) sin(X).*sin(Y).*exp(-2*options.nu*t);

% Lagrange Partikel
subplot(1,2,1);
X_pos = X;
Y_pos = Y;
sctr  = scatter(X_pos',Y_pos','filled');
xlim([0 2*pi]);
ylim([0 2*pi]);
title('Lagrange Partikel');
xlabel('$x$','Interpreter','latex');
ylabel('$y$','Interpreter','latex');

% Stromfunktion
subplot(1,2,2);
cont = surfc(X,Y,Psi(0));
xlim([0 2*pi]);
ylim([0 2*pi]);
zlim([-1 1]);
title('Stromfunktion');
xlabel('$x$','Interpreter','latex');
ylabel('$y$','Interpreter','latex');
zlabel('$\Psi$','Interpreter','latex');

% Zeitschleife
time  = linspace(0,options.t_end,options.t_nr);
t_stp = time(2)-time(1);
for t = time

    % Aktuelle Strömungsgrößen
    U_now   = U(t);
    V_now   = V(t);
    Psi_now = Psi(t);
    
    % Plot
    for row = 1:options.y_nr
        sctr(row).XData = X_pos(row,:);
        sctr(row).YData = Y_pos(row,:);
    end
    cont(1).ZData = Psi_now;
    cont(2).ZData = Psi_now;
    sgtitle([ ...
        '$\nu=' num2str(options.nu) ...
        ',~t=' num2str(t,'%.2f') '$'],'Interpreter','latex');
    drawnow;
    
    % Nächster Zeitschritt
    X_pos = X_pos + t_stp*interp2(X,Y,U_now,X_pos,Y_pos);
    Y_pos = Y_pos + t_stp*interp2(X,Y,V_now,X_pos,Y_pos);

end
end
