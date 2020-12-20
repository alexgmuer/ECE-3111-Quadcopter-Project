m = 0.5;                                %Constants used in the state space 
L = 0.25;                               %representation matrices
k = 3e-6;
b = 1e-7;
g = 9.81;
kd = 0.25;
Ixx = 5e-3;
Iyy = 5e-3;
Izz = 1e-2;
cm = 1e4;

%%%%%%%%%%%%%%%%%%%%%%
% STATE SPACE MATRICES 
%%%%%%%%%%%%%%%%%%%%%%

A = [ 0 0 0 1 0 0 0 0 0 0 0 0;          %Matrix A
    0 0 0 0 1 0 0 0 0 0 0 0;            %
    0 0 0 0 0 1 0 0 0 0 0 0;            %Dimensions: 12x12
    0 0 0 (-1)*kd/m 0 0 0 0 0 0 0 0;
    0 0 0 0 (-1)*kd/m 0 0 0 0 0 0 0;
    0 0 0 0 0 (-1)*kd/m 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 1 0 0;
    0 0 0 0 0 0 0 0 0 0 1 0;
    0 0 0 0 0 0 0 0 0 0 0 1;
    0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0];

B = [0 0 0 0;                           %Matrix B, each column corresponds to inputs u1, u2, u3, u4 
    0 0 0 0;                            %in order from left to right. 
    0 0 0 0;                            %
    0 0 0 0;                            %Dimensions: 12x4
    0 0 0 0;
    (k*cm)/m (k*cm)/m (k*cm)/m (k*cm)/m;
    0 0 0 0;
    0 0 0 0;
    0 0 0 0;
    (L*k*cm)/Ixx 0 (-1)*(L*k*cm)/Ixx 0;
    0 (L*k*cm)/Iyy 0 (-1)*(L*k*cm)/Iyy;
    (L*b*cm)/Izz (-1)*(L*b*cm)/Izz (L*b*cm)/Izz (-1)*(L*b*cm)/Izz];

C = [0 0 1 0 0 0 0 0 0 0 0 0;           %Matrix C, each column/entry corresponds to an output:
    0 0 0 0 0 0 1 0 0 0 0 0;            %from left to right:
    0 0 0 0 0 0 0 1 0 0 0 0;            %x, y, z, vx, vy, vz, phi, theta, psi, wx, wy, wz.
    0 0 0 0 0 0 0 0 1 0 0 0];           %The quadcopter system is only
                                        %concerned with the z, theta, phi,
                                        %and psi outputs, thus it is 4x12
                                        %in order to create a 4x4 T output
                                        %matrix.
                                        %
                                        %C matrix Dimensions: 6x12

D = [0 0 0 0;
    0 0 0 0;
    0 0 0 0;
    0 0 0 0];                            %Zero value for matrix D

%%%%%%%%%%%%%%%%%%
%TRANSFER FUNCTION
%%%%%%%%%%%%%%%%%%

 syms s                                  %Creates a T matrix using the state space to transfer function
                                         %equation. Each entry corresponds
%                                         %to the transfer functions as so: 
%                                         % [z(s)/u1(s)     z(s)/u2(s)     z(s)/u3(s)     z(s)/u4(s)    ]
 T = C/(s*eye(12)-A)*B                  % [phi(s)/u1(s)   phi(s)/u2(s)   phi(s)/u3(s)   phi(s)/u4(s)  ]
                                        % [theta(s)/u1(s) theta(s)/u2(s) theta(s)/u3(s) theta(s)/u4(s)]
                                        % [psi(s)/u1(s)   psi(s)/u2(s)   psi(s)/u3(s)   psi(s)/u4(s)  ]

%  system = ss(A,B,C,D);
%  system;
% T = tf(system);

%Z_tf = tf([0.24],[1, 0.5, 0])

%%%%%%%%%%%%%%%%%%%
%CONTROL CONVERSION
%%%%%%%%%%%%%%%%%%%

 M = [0.5 0.5 0.5 0.5;                  % Mapping matrix which can map the square voltages u1, u2, u3, u4
   0 0 0.5 -0.5;                       % to the control variables u_z, u_phi, u_theta, u_psi
   1 0 -0.5 -0.5;
   0.5 -0.5 0 0];
% %     
 H = T / M *(1/cm)                          % Using the mapping matrix M, and the relation between the square voltages
% %                                              and w1, w2, w3, w4, we can make a
% %                                             transfer function matrix H where
% %                                             each entry is:
% %                                             [z(s)/u_z(s)     z(s)/u_phi(s)     z(s)/u_theta(s)     z(s)/u_psi(s)    ]
% %                                             [phi(s)/u_z(s)   phi(s)/u_phi(s)   phi(s)/u_theta(s)   phi(s)/u_psi(s)  ]
% %                                             [theta(s)/u_z(s) theta(s)/u_phi(s) theta(s)/u_theta(s) theta(s)/u_psi(s)]
% %                                             [psi(s)/u_z(s)   psi(s)/u_phi(s)   psi(s)/u_theta(s)   psi(s)/u_psi(s)  ]
% %                                             The diagonal entries of H
% %                                             correspond to the transfer
% %                                             functions we are interested in: 
% %                                             Z(s)/Uz(s)
% %                                             Phi(s)/U_phi(s)
% %                                             Theta(s)/U_theta(s)
% %                                             Psi(s)/U_psi(s)
% 
 Z_tf_sym = H(1,1);                           % Extracting the specific transfer functions 
 Phi_tf_sym = H(2,2);                         % Z(s)/Uz(s), Phi(s)/U_phi(s), Theta(s)/U_theta(s), 
 Theta_tf_sym = H(3,3);                       % and Psi(s)/U_psi(s)
 Psi_tf_sym = H(4,4);

% [num_x, den_x] = numden(X_tf_sym);
% [num_y, den_y] = numden(Y_tf_sym);
[num_z, den_z] = numden(Z_tf_sym);
[num_phi, den_phi] = numden(Phi_tf_sym);
[num_theta, den_theta] = numden(Theta_tf_sym);
[num_psi, den_psi] = numden(Psi_tf_sym);
% 
% 
% % x_n = sym2poly(num_x);
% % x_d = sym2poly(den_x);
% % 
% % y_n = sym2poly(num_y);
% % y_d = sym2poly(den_y);
% 
z_n = sym2poly(num_z);
z_d = sym2poly(den_z);
% 
phi_n = sym2poly(num_phi);
phi_d = sym2poly(den_phi);

theta_n = sym2poly(num_theta);
theta_d = sym2poly(den_theta);

psi_n = sym2poly(num_psi);
psi_d = sym2poly(den_psi);
% 
% % X_tf = tf(x_n, x_d)
% % Y_tf = tf(y_n, y_d)
 Z_tf = tf(z_n, z_d);
Phi_tf = tf(phi_n, phi_d);
Theta_tf = tf(theta_n, theta_d);
Psi_tf = tf(psi_n, psi_d);
% 
% 
% figure('visible', 'off');
% rlocusplot(Z_tf);
% title('Z Translation Root Locus');
% saveas(gcf, 'Z_rootlocus.png');
% 
% figure('visible', 'off');
% rlocusplot(Phi_tf);
% title('Pitch Root Locus');
% saveas(gcf, 'Pitch_rootlocus.png');
% 
% figure('visible', 'off');
% rlocusplot(Theta_tf);
% title('Roll Root Locus');
% saveas(gcf, 'Roll_rootlocus.png');
% 
% figure('visible', 'off');
% rlocusplot(Psi_tf);
% title('Yaw Root Locus');
% saveas(gcf, 'Yaw_rootlocus.png');
% 
% 
% figure('visible', 'off');
% bode(Z_tf);
% title('Z Translation Bode Plots');
% saveas(gcf, 'Z_bode.png');
% 
% figure('visible', 'off');
% bode(Phi_tf);
% title('Pitch Bode Plots');
% saveas(gcf, 'Pitch_bode.png');
% 
% figure('visible', 'off');
% bode(Theta_tf);
% title('Roll Bode Plots');
% saveas(gcf, 'Roll_bode.png');
% 
% figure('visible', 'off');
% bode(Psi_tf);
% title('Yaw Bode Plots');
% saveas(gcf, 'Yaw_bode.png');

% [C_pi,info] = pidtune(system,'PI')

% figure('visible', 'on');
% opt = stepDataOptions('InputOffset',-1,'StepAmplitude',5);
% t = 0:0.1:15;
% [y t] = step(system,opt)
% saveas(gcf, 'System_response.png');

% system.InputName = {'u1','u2','u3','u4'};
% system.OutputName = {'x Position','y Position','z Position','Pitch','Roll','Yaw'};

