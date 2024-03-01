% Notes on theoretical part of the Project
% Part 1: D-H
syms l0 l1 l2 l3 a positive
syms q1 q2 q3 q4 q5 
syms c1 s1 c2 s2 c3 s3 c4 s4 c5 s5
syms s12 s13 s23 s14 s15 s24 s25 s34 s35 s45
syms c12 c13 c23 c14 c15 c24 c25 c34 c35 c45

% Assuming symbolic variables and their cosine/sine replacements are defined
cs2cossin = @(M) subs(M, {c1, s1, c2, s2, c3, s3, c4, s4, c5, s5, ...
                         s12, s13, s23, s14, s15, s24, s25, s34, s35, s45, ...
                         c12, c13, c23, c14, c15, c24, c25, c34, c35, c45}, ...
                      {cos(q1), sin(q1), cos(q2), sin(q2), cos(q3), sin(q3), ...
                       cos(q4), sin(q4), cos(q5), sin(q5), ...
                       sin(q1 + q2), sin(q1 + q3), sin(q2 + q3), sin(q1 + q4), sin(q1 + q5), ...
                       sin(q2 + q4), sin(q2 + q5), sin(q3 + q4), sin(q3 + q5), sin(q4 + q5), ...
                       cos(q1 + q2), cos(q1 + q3), cos(q2 + q3), cos(q1 + q4), cos(q1 + q5), ...
                       cos(q2 + q4), cos(q2 + q5), cos(q3 + q4), cos(q3 + q5), cos(q4 + q5)});

cossin2cs = @(M) subs(M, {cos(q1), sin(q1), cos(q2), sin(q2), cos(q3), sin(q3), ...
                         cos(q4), sin(q4), cos(q5), sin(q5), ...
                         sin(q1 + q2), sin(q1 + q3), sin(q2 + q3), sin(q1 + q4), sin(q1 + q5), ...
                         sin(q2 + q4), sin(q2 + q5), sin(q3 + q4), sin(q3 + q5), sin(q4 + q5), ...
                         cos(q1 + q2), cos(q1 + q3), cos(q2 + q3), cos(q1 + q4), cos(q1 + q5), ...
                         cos(q2 + q4), cos(q2 + q5), cos(q3 + q4), cos(q3 + q5), cos(q4 + q5)}, ...
                      {c1, s1, c2, s2, c3, s3, c4, s4, c5, s5, ...
                       s12, s13, s23, s14, s15, s24, s25, s34, s35, s45, ...
                       c12, c13, c23, c14, c15, c24, c25, c34, c35, c45});


getT = @(vec) [1 0 0 vec(1); 0 1 0 vec(2); 0 0 1 vec(3); 0 0 0 1];

getR = @(orientation, theta) ...
    (orientation == 'x') * [1 0 0 0; 0 cos(theta) -sin(theta) 0; 0 sin(theta) cos(theta) 0; 0 0 0 1] + ...
    (orientation == 'y') * [cos(theta) 0 sin(theta) 0; 0 1 0 0; -sin(theta) 0 cos(theta) 0; 0 0 0 1] + ...
    (orientation == 'z') * [cos(theta) -sin(theta) 0 0; sin(theta) cos(theta) 0 0; 0 0 1 0; 0 0 0 1];

getA = @(theta, d, alpha, a) vpa(cossin2cs(simplify(getR('z', theta) * getT([0, 0, d]) * getR('x', alpha) * getT([a, 0, 0]))),2);

vpa_cs = @(expr) vpa(cossin2cs(cs2cossin(expr)), 2);

vpa_simplify_cs = @(expr) vpa(cossin2cs(simplify(cs2cossin(expr))), 2);

getZ = @(M) vpa_simplify_cs(M(1:3, 3));

getP = @(M) vpa_simplify_cs(M(1:3, 4));

S = @(x) [
    0 -x(3) x(2);
    x(3) 0 -x(1);
    -x(2) x(1) 0];

% D-H MATRICES
a01 = getA(q1 + pi/2, l0, -pi/2, 0);
a12 = getA(q2-pi/2, -l1, 0 , l2);
a23 = getA(q3, 0, pi/2, 0);
a3E = getA(-pi/2, l3, -pi/2, 0);
% finding p2 z2
tmpa01 = getA(q1+pi/2, l0, -pi/2, 0);
tmpa12 = getA(-pi/2, -l1, 0, 0);
Oq2 = vpa(tmpa01*tmpa12,2) 
% finding p3 z3
Oq3 = vpa(a01 * a12,2)
% finding pE zE
OE = vpa_cs(a01*a12*a23*a3E)

z0 = [0;0;1];
z1 = [-c1; -s1; 0];
z2 = [-c1; -s1; 0];
p0 = [0;0;l0];
p1 = [l1*c1; l1*s1; l0];
p2 = [l1*c1 - l2*s1*s2; l1*s1 + l2*c1*s2; l0 + l2*c2];
pE = [
    l1*c1 + l3*s1*c23 - l2*s1*s2;
    l1*s1 - l3*c1*c23 + l2*c1*s2;
    l0 + l2*c2 + l3*s23;];

J = vpa_simplify_cs([
    S(z0)*(pE - p0) S(z1)*(pE - p1) S(z2)*(pE - p2);
    z0 z1 z2])

JL = J(1:3, :);

J_inv = vpa_simplify_cs(inv(JL))
J_det_exp = vpa_cs(det(JL))
JL_det = vpa_simplify_cs(det(JL))
assume(l2, 'positive')
assume(l3, 'positive')
JL_det_factors = solve(cs2cossin(l3*c23 - l2*s2) == 0)