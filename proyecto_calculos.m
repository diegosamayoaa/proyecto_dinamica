% Análisis dinámico 
% Jose Méndez, Christopher Padilla, DIego Samayoa

clc; clear;

%% Datos del mecanismo

A_punto = [0.0, 0.0];
D_punto = [0.0936, 0.08499];

AB = 0.050;
BC = 0.100;
CD = 0.100;
DP = 0.115;
DG = 0.02128;

omega_AB = -2.701947;
alpha_AB = 0.0;

g = 9.81;

% Masas de las barras. No se toman en cuenta tornillos.
m_AB  = 0.1  * 0.0283495;
m_BC  = 0.15 * 0.0283495;
m_CDP = 0.4  * 0.0283495;

W_AB  = m_AB  * g;
W_BC  = m_BC  * g;
W_CDP = m_CDP * g;

% Momentos de inercia respecto al centro de masa de cada barra
I_AB  = (1/12) * m_AB * AB^2;
I_BC  = (1/12) * m_BC * BC^2;
I_CDP = 57.823e-6;

thetas_deg = [-145.3, -185.3, -226.8, -246.9, -270.4, -316.4];
tiempos    = [0.900,  1.266,  1.633,  1.766,  1.899,  2.099];

npos = length(thetas_deg);

%% Matrices para guardar resultados

RES_CIN  = zeros(npos, 10);
RES_RXNS = zeros(npos, 9);
RES_VEL  = zeros(npos, 6);

fprintf('============================================================\n');
fprintf(' ANALISIS DINAMICO DE MECANISMO PLANO\n');
fprintf('============================================================\n\n');

%% Cálculo posición por posición

for i = 1:npos

    th = deg2rad(thetas_deg(i));

    %% Posición de los puntos

    B = AB * [cos(th), sin(th)];

    % Punto C por intersección de círculos:
    % C debe cumplir BC = 0.100 m y CD = 0.100 m.
    dBD = norm(D_punto - B);
    M   = (B + D_punto) / 2;
    h   = sqrt(BC^2 - (dBD/2)^2);

    C(1) = M(1) - h*(D_punto(2) - B(2)) / dBD;
    C(2) = M(2) + h*(D_punto(1) - B(1)) / dBD;

    r_CD = C - D_punto;
    P = D_punto + (-DP/CD) * r_CD;

    r_CB     = C - B;
    r_CD_vec = C - D_punto;
    r_PD     = P - D_punto;

    %% Velocidades

    % El punto B pertenece a AB, que gira alrededor del punto fijo A.
    vB = omega_AB * [-B(2), B(1)];

    % Se iguala la velocidad de C vista desde BC y desde CDP.
    M_vel = [-r_CB(2),  r_CD_vec(2);
              r_CB(1), -r_CD_vec(1)];

    b_vel = -vB';

    sol_vel = M_vel \ b_vel;

    w_BC = sol_vel(1);
    w_CD = sol_vel(2);

    vC = w_CD * [-r_CD_vec(2), r_CD_vec(1)];
    vP = w_CD * [-r_PD(2),     r_PD(1)];

    %% Aceleraciones

    % Como omega_AB es constante, alpha_AB = 0.
    % Por eso B solo tiene aceleración normal.
    aB = -omega_AB^2 * B;

    % Se iguala la aceleración de C vista desde BC y desde CDP.
    b_acc = (-aB + w_BC^2 * r_CB - w_CD^2 * r_CD_vec)';

    sol_acc = M_vel \ b_acc;

    alpha_BC = sol_acc(1);
    alpha_CD = sol_acc(2);

    aC = alpha_CD * [-r_CD_vec(2), r_CD_vec(1)] - w_CD^2 * r_CD_vec;
    aP = alpha_CD * [-r_PD(2),     r_PD(1)]     - w_CD^2 * r_PD;

    %% Centros de masa

    % AB es uniforme, por eso su centro de masa está a la mitad.
    G_AB  = B / 2;
    aG_AB = -omega_AB^2 * G_AB;

    % BC también se toma uniforme, sin tornillos.
    r_GB  = 0.5 * r_CB;
    aG_BC = aB + alpha_BC * [-r_GB(2), r_GB(1)] - w_BC^2 * r_GB;

    % El centro de masa de CDP está del lado de P respecto a D.
    r_GD   = (-DG/CD) * r_CD_vec;
    aG_CDP = alpha_CD * [-r_GD(2), r_GD(1)] - w_CD^2 * r_GD;

    %% Sistema dinámico

    % Incógnitas:
    % [Ax, Ay, Bx, By, Cx, Cy, Dx, Dy, MA]
    %
    % Se aplica:
    % sum Fx = m*aGx
    % sum Fy = m*aGy
    % sum MG = IG*alpha

    rA_gAB = A_punto - G_AB;
    rB_gAB = B       - G_AB;

    G_BC_abs = B + r_GB;
    rB_gBC   = B - G_BC_abs;
    rC_gBC   = C - G_BC_abs;

    G_CDP_abs = D_punto + r_GD;
    rC_gCDP   = C       - G_CDP_abs;
    rD_gCDP   = D_punto - G_CDP_abs;

    Amat = zeros(9,9);
    bvec = zeros(9,1);

    % Barra AB
    Amat(1,1) = 1;
    Amat(1,3) = 1;
    bvec(1) = m_AB * aG_AB(1);

    Amat(2,2) = 1;
    Amat(2,4) = 1;
    bvec(2) = m_AB * aG_AB(2) + W_AB;

    Amat(3,9) = 1;
    Amat(3,1) = -rA_gAB(2);
    Amat(3,2) =  rA_gAB(1);
    Amat(3,3) = -rB_gAB(2);
    Amat(3,4) =  rB_gAB(1);
    bvec(3) = I_AB * alpha_AB;

    % Barra BC
    Amat(4,3) = -1;
    Amat(4,5) =  1;
    bvec(4) = m_BC * aG_BC(1);

    Amat(5,4) = -1;
    Amat(5,6) =  1;
    bvec(5) = m_BC * aG_BC(2) + W_BC;

    Amat(6,3) =  rB_gBC(2);
    Amat(6,4) = -rB_gBC(1);
    Amat(6,5) = -rC_gBC(2);
    Amat(6,6) =  rC_gBC(1);
    bvec(6) = I_BC * alpha_BC;

    % Barra CDP
    Amat(7,5) = -1;
    Amat(7,7) =  1;
    bvec(7) = m_CDP * aG_CDP(1);

    Amat(8,6) = -1;
    Amat(8,8) =  1;
    bvec(8) = m_CDP * aG_CDP(2) + W_CDP;

    Amat(9,5) =  rC_gCDP(2);
    Amat(9,6) = -rC_gCDP(1);
    Amat(9,7) = -rD_gCDP(2);
    Amat(9,8) =  rD_gCDP(1);
    bvec(9) = I_CDP * alpha_CD;

    sol = Amat \ bvec;

    Ax = sol(1); Ay = sol(2);
    Bx = sol(3); By = sol(4);
    Cx = sol(5); Cy = sol(6);
    Dx = sol(7); Dy = sol(8);
    MA = sol(9);

    residuo = norm(Amat*sol - bvec);

    %% Guardar resultados

    RES_CIN(i,:) = [w_BC, w_CD, alpha_BC, alpha_CD, ...
                    norm(vB), norm(aB), norm(vC), norm(aC), ...
                    norm(vP), norm(aP)];

    RES_RXNS(i,:) = [Ax, Ay, Bx, By, Cx, Cy, Dx, Dy, MA];

    RES_VEL(i,:) = [vB(1), vB(2), vC(1), vC(2), vP(1), vP(2)];

    %% Imprimir resultados de cada posición

    fprintf('------ Posicion %d  (t = %.3f s,  theta_AB = %.1f deg) ------\n', ...
            i, tiempos(i), thetas_deg(i));

    fprintf('  Cinematica:\n');
    fprintf('    w_BC   = %10.4f rad/s    w_CD   = %10.4f rad/s\n', w_BC, w_CD);
    fprintf('    a_BC   = %10.4f rad/s^2  a_CD   = %10.4f rad/s^2\n', alpha_BC, alpha_CD);
    fprintf('    |vB|   = %10.4f m/s      |aB|   = %10.4f m/s^2\n', norm(vB), norm(aB));
    fprintf('    |vC|   = %10.4f m/s      |aC|   = %10.4f m/s^2\n', norm(vC), norm(aC));
    fprintf('    |vP|   = %10.4f m/s      |aP|   = %10.4f m/s^2\n', norm(vP), norm(aP));

    fprintf('  Velocidades por componente:\n');
    fprintf('    vBx = %8.4f m/s   vBy = %8.4f m/s\n', vB(1), vB(2));
    fprintf('    vCx = %8.4f m/s   vCy = %8.4f m/s\n', vC(1), vC(2));
    fprintf('    vPx = %8.4f m/s   vPy = %8.4f m/s\n', vP(1), vP(2));

    fprintf('  Reacciones:\n');
    fprintf('    Ax = %10.6f N    Ay = %10.6f N\n', Ax, Ay);
    fprintf('    Bx = %10.6f N    By = %10.6f N\n', Bx, By);
    fprintf('    Cx = %10.6f N    Cy = %10.6f N\n', Cx, Cy);
    fprintf('    Dx = %10.6f N    Dy = %10.6f N\n', Dx, Dy);
    fprintf('    MA = %10.6f N*m\n', MA);
    fprintf('    Residuo = %.2e\n\n', residuo);

end

%% Tabla resumen de cinemática

fprintf('============================================================\n');
fprintf(' TABLA RESUMEN - CINEMATICA\n');
fprintf('============================================================\n');
fprintf('%6s %8s %8s %8s %8s %7s %7s %7s %7s %7s %7s\n', ...
    't(s)','w_BC','w_CD','a_BC','a_CD','|vB|','|aB|','|vC|','|aC|','|vP|','|aP|');

for i = 1:npos
    fprintf('%6.3f %8.4f %8.4f %8.4f %8.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n', ...
        tiempos(i), RES_CIN(i,:));
end

%% Tabla resumen de velocidades por componente

fprintf('\n============================================================\n');
fprintf(' TABLA RESUMEN - VELOCIDADES POR COMPONENTE (m/s)\n');
fprintf('============================================================\n');
fprintf('%6s %9s %9s %9s %9s %9s %9s\n', ...
    't(s)','vBx','vBy','vCx','vCy','vPx','vPy');

for i = 1:npos
    fprintf('%6.3f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f\n', ...
        tiempos(i), RES_VEL(i,:));
end

%% Tabla resumen de reacciones

fprintf('\n============================================================\n');
fprintf(' TABLA RESUMEN - REACCIONES\n');
fprintf('============================================================\n');
fprintf('%6s %10s %10s %10s %10s %10s %10s %10s %10s %10s\n', ...
    't(s)','Ax(N)','Ay(N)','Bx(N)','By(N)','Cx(N)','Cy(N)','Dx(N)','Dy(N)','MA(N*m)');

for i = 1:npos
    fprintf('%6.3f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n', ...
        tiempos(i), RES_RXNS(i,:));
end