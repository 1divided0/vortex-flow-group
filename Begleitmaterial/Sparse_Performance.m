%% Laufzeitvergleich für die Multiplikation dünnbesetzter Matrizen

function Sparse_Performance(N)
arguments
    N uint8 % Maximale Problemgröße
end

% Initialisierung
dims  = 1:N;
times = zeros(N,3);

% Iterationsschleife
for i=dims

    % Definition der Funktionen
    with_sparse    = @() speye(i)*ones(i,1);
    without_sparse = @() eye(i)*ones(i,1);
    no_zeros       = @() ones(i,1).*ones(i,1);
    
    % Laufzeitmessung
    times(i,1) = timeit(with_sparse);
    times(i,2) = timeit(without_sparse);
    times(i,3) = timeit(no_zeros);
end

% Plot
figure();
axes('XScale','log','YScale','log');
hold on;
plot(dims,times(:,1));
plot(dims,times(:,2));
plot(dims,times(:,3));
hold off;
xlabel('$N$','Interpreter','latex');
ylabel('Laufzeit in Sekunden');
xlim([-inf,inf]);
ylim([-inf,inf]);
axis padded;
grid on;
legend( ...
    '$\mathbf{I}_{N \times N}\cdot\mathbf{1}_{N \times 1}$ (mit CSC)', ...
    '$\mathbf{I}_{N \times N}\cdot\mathbf{1}_{N \times 1}$ (ohne Kompression)', ...
    '$\mathbf{1}_{N \times 1}\odot\mathbf{1}_{N \times 1}$', ...
    'Interpreter','latex','Location','best');
end
