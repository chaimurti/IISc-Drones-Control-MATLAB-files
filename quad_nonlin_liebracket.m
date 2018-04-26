syms x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12
% syms Ix Iy Iz 
Ix = 2.5
Iy = 2.5
Iz = 2.5


X = [x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12];
% Cart_states = X(1:3);
% Ang_states = X(4:6);
% Cart_vels = X(7:9);
% Ang_vels = X(10:12);




f_n(1,1) = x7; %*(sin(x4)*sin(x6) + cos(x4)*cos(x6)*sin(x5)) - x8*(cos(x4)*sin(x6) -cos(x6)*sin(x4)*sin(x5)) + x7*(cos(x6)*cos(x5));
f_n(2,1) = x8; %*(cos(x4)*cos(x6)+sin(x4)*sin(x5)*sin(x6)) - x9*(cos(x6)*sin(x4) - cos(x4)*sin(x5)*sin(x6)) + x7*(cos(x5)*sin(x6));
f_n(3,1) = x9 % *(cos(x4)*cos(x5)) - x7*(sin(x5)) +x8*(cos(x5)*sin(x4));
f_n(4,1) = x10 + x12*(cos(x4)*tan(x5)) +x11*(sin(x4)*tan(x5));
f_n(5,1) = x11*(cos(x4)) - x12*sin(x4);
f_n(6,1) = x12*(cos(x4)/cos(x5)) +x11*(sin(x4)/cos(x5));
f_n(7,1) = 0;
f_n(8,1) = 0;
f_n(9,1) = 9.81;
f_n(10,1) = ((Iy - Iz)/Ix)*x11*x12 ;
f_n(11,1) = ((Iz-Ix)/Iy)*x10*x12;
f_n(12,1) = ((Ix-Iy)/Iz)*x10*x11;

g_n1 = sym(zeros(12,1));
g_n2 = sym(zeros(12,1));
g_n3 = sym(zeros(12,1));
g_n4 = sym(zeros(12,1));

g_n1(7,1) =-sin(x4)*sin(x6) - cos(x4)*cos(x6)*sin(x5);
g_n1(8,1) = -cos(x4)*sin(x6)-cos(x4)*sin(x5)*sin(x6);
g_n1(9,1) = -cos(x4)*cos(x5);

g_n2(10,1) = 1/Ix;
g_n3(11,1) = 1/Iy;
g_n4(12,1) = 1/Iz;

LB{1} = g_n1;
LB{2} = g_n2;
LB{3} = g_n3;
LB{4} = g_n4;

% F{1} = LB{1};
% F{2} = LB{2};
% F{3} = LB{3};
% F{4} = LB{4};

for i=1:5
    LB{4*i+1} = jacobian(f_n,X)*LB{4*i-3} - jacobian(LB{4*i-3},X)*f_n;
    LB{4*i+2} = jacobian(f_n,X)*LB{4*i-2} - jacobian(LB{4*i-2},X)*f_n;
    LB{4*i+3} = jacobian(f_n,X)*LB{4*i-1} - jacobian(LB{4*i-1},X)*f_n;
    LB{4*i+4} = jacobian(f_n,X)*LB{4*i} - jacobian(LB{4*i},X)*f_n;
end


T = zeros(12,32);

for i=1:22
    t = zeros(12,1);
    t = subs(LB{i},{x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12},{rand(1,1),rand(1,1),rand(1,1),pi/4,pi/4,pi/4,rand(1,1),rand(1,1),rand(1,1),rand(1,1),rand(1,1),rand(1,1)});
    T(:,i) = t;
end

