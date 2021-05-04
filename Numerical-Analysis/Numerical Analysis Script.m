clear; clc; close all

width = 1.5; % Line thickness
lsize = 12; tsize = 12; legsize = 9; msize = 8; % Font size in label, title, legend, markers
set(0, 'defaultAxesTickLabelInterpreter', 'latex');
set(0, 'defaultTextInterpreter', 'latex'); set(0, 'defaultLegendInterpreter', 'latex');

%% Question 1 (d)
x = [1 1 1 1 1; -2 -1 0 1 2; 4/2 1/2 0 1/2 4/2; -8/6 -1/6 0 1/6 8/6; 16/24 1/24 0 1/24 16/24];
b = [0 0 0 0 1];
a = (inv(x)).*b;

%% Question 2 (a)
x0 = 0;                % Initial x value
g = sqrt((exp(x0))/3); % Initial g(x) value
i = 1;                 % Counter
while abs(g - x0) > 0.0001 % Error threshold
    z(i,1) = x0;           % Store data points
    z(i,2) = g;
    x0 = g;
    g = sqrt((exp(g))/3);  % g(x) points
    i = i + 1; % Counter advancement
end
x = [-0.2:0.01:1.2];
y = x;                 % y = x line
g2 = sqrt((exp(x))/3); % g(x) line
figure (1)
hold on; grid on
xlim([-0.2 1.2]); ylim([-0.2 1.2]);
plot1 = plot(x, y, 'LineWidth', width);
spl_1 = spline(x, g2);
plot2 = plot(x, ppval(spl_1, x), 'g-', 'LineWidth', width); % x, g2
for j=1:2:length(z)
    plot([z(j,1),z(j,1)],[z(j,2),z(j,1)],'r--', 'LineWidth', width)
    plot([z(j,2),z(j,1)],[z(j,2),z(j,2)],'r--', 'LineWidth', width)
end
for j=2:2:length(z)
    plot([z(j,1),z(j,1)],[z(j,1),z(j,2)],'r--', 'LineWidth', width)
    plot([z(j,1),z(j,2)],[z(j,2),z(j,2)],'r--', 'LineWidth', width)
end
hold off
title('g(x) = $\sqrt{\frac{exp(x)}{3}}$', 'FontSize', tsize);
xlabel('x', 'FontSize', lsize); ylabel('y', 'FontSize', lsize);
L = legend([plot1 plot2], 'x', 'g(x)', 'Location', 'best');
L.FontSize = legsize;
% labelx = 0.1; labely = 0.2; txt = 'from\ $x_0$';
% text(labelx, labely, txt, 'Interpreter', 'latex')
annotation('textarrow',[0.24,0.24],[0.23,0.45]);
hold off

%% Question 2 (a)
x0 = 0;                 % Initial x value
g = -sqrt((exp(x0))/3); % Initial g(x) value
i = 1;                  % Counter
while abs(g - x0) > 0.0001 % Error threshold
    z(i,1) = x0;           % Store data points
    z(i,2) = g;
    x0 = g;
    g = -sqrt((exp(g))/3);  % g(x) points
    i = i + 1; % Counter advancement
end
x = [-0.7:0.01:0.1];
y = x;                  % y = x line
g2 = -sqrt((exp(x))/3); % g(x) line
figure (2)
hold on; grid on
xlim([-0.7 0.1]); ylim([-0.9 0.2]);
plot1 = plot(x, y, 'LineWidth', width);
spl_1 = spline(x, g2);
plot2 = plot(x, ppval(spl_1, x), 'g-', 'LineWidth', width); % x, g2
for j=1:2:length(z)
    plot([z(j,1),z(j,1)],[z(j,2),z(j,1)],'r--', 'LineWidth', width)
    plot([z(j,2),z(j,1)],[z(j,2),z(j,2)],'r--', 'LineWidth', width)
end
for j=2:2:length(z)
    plot([z(j,1),z(j,1)],[z(j,1),z(j,2)],'r--', 'LineWidth', width)
    plot([z(j,1),z(j,2)],[z(j,2),z(j,2)],'r--', 'LineWidth', width)
end
hold off
title('g(x) = $-\sqrt{\frac{exp(x)}{3}}$', 'FontSize', tsize);
xlabel('x', 'FontSize', lsize); ylabel('y', 'FontSize', lsize);
L = legend([plot1 plot2], 'x', 'g(x)', 'Location', 'best');
L.FontSize = legsize;
% labelx = 0.1; labely = 0.2; txt = 'from\ $x_0$';
% text(labelx, labely, txt, 'Interpreter', 'latex')
annotation('textarrow',[0.81,0.81],[0.77,0.5]);
hold off

%% Question 2 (a)
x0 = 3.3;              % Initial x value
g = log(3 * x0^2);     % Initial g(x) value
i = 1;                 % Counter
while abs(g - x0) > 0.0001 % Error threshold
    z(i,1) = x0;           % Store data points
    z(i,2) = g;
    x0 = g;
    g = log(3 * g^2);  % g(x) points
    i = i + 1; % Counter advancement
end
x = [3.2:0.01:4];
y = x;              % y = x line
g2 = log(3 * x.^2); % g(x) line
figure (3)
hold on; grid on
xlim([3.2 4]); ylim([3.2 4]);
plot1 = plot(x, y, 'LineWidth', width);
spl_1 = spline(x, g2);
plot2 = plot(x, ppval(spl_1, x), 'g-', 'LineWidth', width); % x, g2
for j=1:2:length(z)
    plot([z(j,1),z(j,1)],[z(j,2),z(j,1)],'r--', 'LineWidth', width)
    plot([z(j,2),z(j,1)],[z(j,2),z(j,2)],'r--', 'LineWidth', width)
end
for j=2:2:length(z)
    plot([z(j,1),z(j,1)],[z(j,1),z(j,2)],'r--', 'LineWidth', width)
    plot([z(j,1),z(j,2)],[z(j,2),z(j,2)],'r--', 'LineWidth', width)
end
hold off
title('g(x) = $ln(3 x^2)$', 'FontSize', tsize);
xlabel('x', 'FontSize', lsize); ylabel('y', 'FontSize', lsize);
L = legend([plot1 plot2], 'x', 'g(x)', 'Location', 'best');
L.FontSize = legsize;
% labelx = 0.1; labely = 0.2; txt = 'from\ $x_0$';
% text(labelx, labely, txt, 'Interpreter', 'latex')
annotation('textarrow',[0.225,0.225],[0.21,0.35]);
hold off

%% Question 2 (a)
x0 = 4;                % Initial x value
g = sqrt((exp(x0))/3); % Initial g(x) value
i = 1;                 % Counter
while abs(g - x0) > 0.0001 % Error threshold
    z(i,1) = x0;           % Store data points
    z(i,2) = g;
    x0 = g;
    g = sqrt((exp(g))/3);  % g(x) points
    i = i + 1; % Counter advancement
end
x = [3:0.01:8];
y = x;                 % y = x line
g2 = sqrt((exp(x))/3); % g(x) line
figure (4)
hold on; grid on
xlim([3 8]); ylim([3 20]);
plot1 = plot(x, y, 'LineWidth', width);
spl_1 = spline(x, g2);
plot2 = plot(x, ppval(spl_1, x), 'g-', 'LineWidth', width); % x, g2
for j=1:2:length(z)
    plot([z(j,1),z(j,1)],[z(j,2),z(j,1)],'r--', 'LineWidth', width)
    plot([z(j,2),z(j,1)],[z(j,2),z(j,2)],'r--', 'LineWidth', width)
end
for j=2:2:length(z)
    plot([z(j,1),z(j,1)],[z(j,1),z(j,2)],'r--', 'LineWidth', width)
    plot([z(j,1),z(j,2)],[z(j,2),z(j,2)],'r--', 'LineWidth', width)
end
hold off
title('g(x) = $\sqrt{\frac{exp(x)}{3}}$', 'FontSize', tsize);
xlabel('x', 'FontSize', lsize); ylabel('y', 'FontSize', lsize);
L = legend([plot1 plot2], 'x', 'g(x)', 'Location', 'best');
L.FontSize = legsize;
% labelx = 0.1; labely = 0.2; txt = 'from\ $x_0$';
% text(labelx, labely, txt, 'Interpreter', 'latex')
annotation('textarrow',[0.42,0.55],[0.28,0.28]);
hold off

%% Question 2 (b)
clear;
width = 1.5; % Line thickness
lsize = 12; tsize = 12; legsize = 9; msize = 8; % Font size in label, title, legend, markers
set(0, 'defaultAxesTickLabelInterpreter', 'latex');
set(0, 'defaultTextInterpreter', 'latex'); set(0, 'defaultLegendInterpreter', 'latex');

q = 0;
p = (-3*pi)/2;
n = 2;          % number of intervals 
exact = (1/50)*(7*(exp((3*pi)/2)) - 1);
error = 1;
rel_err = 0.1; % relative error (0.1 or 0.01)

while error >= rel_err
    h = (q - p)/n;
    
    for j = 1:n+1
        x(j) = p + h*(j-1);
        f(j) = (cos(7*x(j)))/(exp(x(j)));
    end
    
    simpson = (h/3)*(2*sum(f) - f(1) - f(n+1) + 2*sum(f(2:2:n)));
    z(n) = simpson;
    error = 100 * abs((simpson - exact)/exact);
    n = n + 2;
end

fprintf("For %.2f %% error, ", rel_err)
fprintf("%.0f intervals ", n - 2)
fprintf("& integral value = %f \n", simpson)

figure (5)
hold on; grid on
plot1 = plot(2:2:n-2, z(2:2:end), 'kx', 'MarkerSize', msize, 'LineWidth', width);
plot2 = plot([0, n-2], [exact, exact], 'b-', 'LineWidth', width);
title('Convergence of the Simpson $\frac{1}{3}$ Scheme', 'FontSize', tsize);
xlabel('Number of sub-intervals', 'FontSize', lsize);
ylabel('Integral Approximation', 'FontSize', lsize);
L = legend([plot1 plot2], 'Simpson $\frac{1}{3}$', 'Exact Integral', 'Location', 'best');
L.FontSize = legsize;

%% Question 2 (c)
clear;
width = 1.5; % Line thickness
lsize = 12; tsize = 12; legsize = 9; msize = 8; % Font size in label, title, legend, markers
set(0, 'defaultAxesTickLabelInterpreter', 'latex');
set(0, 'defaultTextInterpreter', 'latex'); set(0, 'defaultLegendInterpreter', 'latex');

field = 'f';
value = {[1 2 0];
[1	1.0000000000000000	-0.5773502691896257
2	1.0000000000000000	0.5773502691896257];
[1	0.8888888888888888	0.0000000000000000
2	0.5555555555555556	-0.7745966692414834
3	0.5555555555555556	0.7745966692414834];
[1	0.6521451548625461	-0.3399810435848563
2	0.6521451548625461	0.3399810435848563
3	0.3478548451374538	-0.8611363115940526
4	0.3478548451374538	0.8611363115940526];
[1	0.5688888888888889	0.0000000000000000
2	0.4786286704993665	-0.5384693101056831
3	0.4786286704993665	0.5384693101056831
4	0.2369268850561891	-0.9061798459386640
5	0.2369268850561891	0.9061798459386640];
[1	0.3607615730481386	0.6612093864662645
2	0.3607615730481386	-0.6612093864662645
3	0.4679139345726910	-0.2386191860831969
4	0.4679139345726910	0.2386191860831969
5	0.1713244923791704	-0.9324695142031521
6	0.1713244923791704	0.9324695142031521];
[1	0.4179591836734694	0.0000000000000000
2	0.3818300505051189	0.4058451513773972
3	0.3818300505051189	-0.4058451513773972
4	0.2797053914892766	-0.7415311855993945
5	0.2797053914892766	0.7415311855993945
6	0.1294849661688697	-0.9491079123427585
7	0.1294849661688697	0.9491079123427585];
[1	0.3626837833783620	-0.1834346424956498
2	0.3626837833783620	0.1834346424956498
3	0.3137066458778873	-0.5255324099163290
4	0.3137066458778873	0.5255324099163290
5	0.2223810344533745	-0.7966664774136267
6	0.2223810344533745	0.7966664774136267
7	0.1012285362903763	-0.9602898564975363
8	0.1012285362903763	0.9602898564975363];
[1	0.3302393550012598	0.0000000000000000
2	0.1806481606948574	-0.8360311073266358
3	0.1806481606948574	0.8360311073266358
4	0.0812743883615744	-0.9681602395076261
5	0.0812743883615744	0.9681602395076261
6	0.3123470770400029	-0.3242534234038089
7	0.3123470770400029	0.3242534234038089
8	0.2606106964029354	-0.6133714327005904
9	0.2606106964029354	0.6133714327005904];
[1	0.2955242247147529	-0.1488743389816312
2	0.2955242247147529	0.1488743389816312
3	0.2692667193099963	-0.4333953941292472
4	0.2692667193099963	0.4333953941292472
5	0.2190863625159820	-0.6794095682990244
6	0.2190863625159820	0.6794095682990244
7	0.1494513491505806	-0.8650633666889845
8	0.1494513491505806	0.8650633666889845
9	0.0666713443086881	-0.9739065285171717
10	0.0666713443086881	0.9739065285171717];
[1	0.2729250867779006	0.0000000000000000
2	0.2628045445102467	-0.2695431559523450
3	0.2628045445102467	0.2695431559523450
4	0.2331937645919905	-0.5190961292068118
5	0.2331937645919905	0.5190961292068118
6	0.1862902109277343	-0.7301520055740494
7	0.1862902109277343	0.7301520055740494
8	0.1255803694649046	-0.8870625997680953
9	0.1255803694649046	0.8870625997680953
10	0.0556685671161737	-0.9782286581460570
11	0.0556685671161737	0.9782286581460570];
[1	0.2491470458134028	-0.1252334085114689
2	0.2491470458134028	0.1252334085114689
3	0.2334925365383548	-0.3678314989981802
4	0.2334925365383548	0.3678314989981802
5	0.2031674267230659	-0.5873179542866175
6	0.2031674267230659	0.5873179542866175
7	0.1600783285433462	-0.7699026741943047
8	0.1600783285433462	0.7699026741943047
9	0.1069393259953184	-0.9041172563704749
10	0.1069393259953184	0.9041172563704749
11	0.0471753363865118	-0.9815606342467192
12	0.0471753363865118	0.9815606342467192];
[1	0.2325515532308739	0.0000000000000000
2	0.2262831802628972	-0.2304583159551348
3	0.2262831802628972	0.2304583159551348
4	0.2078160475368885	-0.4484927510364469
5	0.2078160475368885	0.4484927510364469
6	0.1781459807619457	-0.6423493394403402
7	0.1781459807619457	0.6423493394403402
8	0.1388735102197872	-0.8015780907333099
9	0.1388735102197872	0.8015780907333099
10	0.0921214998377285	-0.9175983992229779
11	0.0921214998377285	0.9175983992229779
12	0.0404840047653159	-0.9841830547185881
13	0.0404840047653159	0.9841830547185881];
[1	0.2152638534631578	-0.1080549487073437
2	0.2152638534631578	0.1080549487073437
3	0.2051984637212956	-0.3191123689278897
4	0.2051984637212956	0.3191123689278897
5	0.1855383974779378	-0.5152486363581541
6	0.1855383974779378	0.5152486363581541
7	0.1572031671581935	-0.6872929048116855
8	0.1572031671581935	0.6872929048116855
9	0.1215185706879032	-0.8272013150697650
10	0.1215185706879032	0.8272013150697650
11	0.0801580871597602	-0.9284348836635735
12	0.0801580871597602	0.9284348836635735
13	0.0351194603317519	-0.9862838086968123
14	0.0351194603317519	0.9862838086968123];
[1	0.2025782419255613	0.0000000000000000
2	0.1984314853271116	-0.2011940939974345
3	0.1984314853271116	0.2011940939974345
4	0.1861610000155622	-0.3941513470775634
5	0.1861610000155622	0.3941513470775634
6	0.1662692058169939	-0.5709721726085388
7	0.1662692058169939	0.5709721726085388
8	0.1395706779261543	-0.7244177313601701
9	0.1395706779261543	0.7244177313601701
10	0.1071592204671719	-0.8482065834104272
11	0.1071592204671719	0.8482065834104272
12	0.0703660474881081	-0.9372733924007060
13	0.0703660474881081	0.9372733924007060
14	0.0307532419961173	-0.9879925180204854
15	0.0307532419961173	0.9879925180204854];};

s = struct(field, value);
exact = (1/50)*(7*(exp((3*pi)/2)) - 1);
error = 1;
i = 0;
rel_err = 0.01; % relative error

while error >= rel_err
    sum = 0;
    i = i+1;
    n = s(i).f;
    
    for j = 1:i
        sum = sum + n(j,2)*((cos(7*(((-3*pi)/4)-((3*pi*n(j,3))/4))))/(exp(((-3*pi)/4)-((3*pi*n(j,3))/4))));
    end
    
    z(i) = ((3*pi)/4) * sum;
    answ = ((3*pi)/4) * sum;
    error = 100 * abs((answ - exact)/exact);
end

% % Iterate 2 more times, change i to i + 2 in plots
% n14 = s(14).f ;
% sum14 = 0;
% 
% for j = 1:14
%     sum14 = sum14 + n14(j,2)*((cos(7*(((-3*pi)/4)-((3*pi*n14(j,3))/4))))/(exp(((-3*pi)/4)-((3*pi*n14(j,3))/4))));
% end
% 
% z(14) = ((3*pi)/4) * sum14;
% n15 = s(15).f ;
% sum15 = 0;
% 
% for j = 1:15
%     sum15 = sum15 + n15(j,2)*((cos(7*(((-3*pi)/4)-((3*pi*n15(j,3))/4))))/(exp(((-3*pi)/4)-((3*pi*n15(j,3))/4))));
% end
% z(15) = ((3*pi)/4) * sum15;

fprintf("For %.2f %% error, ", rel_err)
fprintf("%.0f gauss points ", i)
fprintf("& integral value = %f \n", answ)

figure (6)
hold on; grid on
plot1 = plot(1:i, z, 'kx', 'MarkerSize', msize, 'LineWidth', width); % i + 2
plot2 = plot([0, i], [answ, answ], 'b-', 'LineWidth', width); % i + 2
% plot([13,13],[-80,100],'g-')
% plot([14,14],[-80,100],'g-')
title('Convergence of the Gauss-Legendre Scheme', 'FontSize', tsize);
xlabel('Quadrature Points, i', 'FontSize', lsize);
ylabel('Integral Approximation', 'FontSize', lsize);
L = legend([plot1 plot2], 'Gauss-Legendre', 'Exact Integral', 'Location', 'best');
L.FontSize = legsize;

%% Question 3 (a) & (b)
clear;
width = 1.5; % Line thickness
lsize = 12; tsize = 12; legsize = 9; msize = 8; % Font size in label, title, legend, markers
set(0, 'defaultAxesTickLabelInterpreter', 'latex');
set(0, 'defaultTextInterpreter', 'latex'); set(0, 'defaultLegendInterpreter', 'latex');

sampling = [10 20 40 80 160]; % sampling at 10, 20, 40, 80, 160
l2_arr = zeros(length(sampling), 1);
linf_arr = zeros(length(sampling), 1);

for sample = 1:length(sampling)
    
    x = linspace(-1, 1, 1001);
    x2 = linspace(-1, 1, sampling(sample));
    y = 1./(1 + 25*(x2.^2));
    y1 = 1./(1 + 25*(x.^2));
    sum = 0;

    for i = 1:1001
        sum = 0;
        for j = 1:sampling(sample)
            L = 1;
            for k = 1:sampling(sample)
                if j ~= k
                    L = L*((x(i)-x2(k))/(x2(j)-x2(k)));
                end
            end
            sum = sum + (y(j)*L);
        end
        y2(i) = sum;
    end

    error = y2 - y1;
    l2 = norm(error, 2);
    linf = norm(error, inf);

    fprintf("Using %.0f equispaced points; ", sampling(sample))
    fprintf("L_2_norm = %f ", l2)
    fprintf("& L_inf_norm = %f \n", linf)
    
    l2_arr(sample, 1) = l2;     % [3.4658 56.2386 6.2097*(10^4) 2.1183*(10^11) 1.01128*(10^29)]
    linf_arr(sample, 1) = linf; % [0.3003 8.5756 1.4428*(10^4) 7.3032*(10^10) 8.4497*(10^28)]

    figure (6 + sample)
    hold on; grid on
    plot(x, y2, 'LineWidth', width)
    plot(x, y1, 'LineWidth', width)
    title(['Interpolatory Polynomial Approximation sampling at ', num2str(sampling(sample)), ' equispaced points'], 'FontSize', tsize);
    xlabel('x', 'FontSize', lsize);
    ylabel('y', 'FontSize', lsize);
    hold off

end

figure (12)
hold on; grid on
xlim([10 160]); % ylim([-0.5*(10^28) 11*(10^28)])
plot1 = plot(sampling, log(l2_arr), 'LineWidth', width);
plot2 = plot(sampling, log(linf_arr), 'LineWidth', width);
title('Norm (log scale) vs sampling points', 'FontSize', tsize);
xlabel('Number of equispaced points used', 'FontSize', lsize);
ylabel('Norm (log scale) of the error', 'FontSize', lsize);
L = legend([plot1 plot2], '$L_2$ norm', '$L_{\infty}$ norm', 'Location', 'best');
L.FontSize = legsize;

%% Question 3 (c)
clear;
width = 1.5; % Line thickness
lsize = 12; tsize = 12; legsize = 9; msize = 8; % Font size in label, title, legend, markers
set(0, 'defaultAxesTickLabelInterpreter', 'latex');
set(0, 'defaultTextInterpreter', 'latex'); set(0, 'defaultLegendInterpreter', 'latex');

sampling = [10 20 40 80 160]; % sampling at 10, 20, 40, 80, 160
l2_arr = zeros(length(sampling), 1);
linf_arr = zeros(length(sampling), 1);

for sample = 1:length(sampling)
    
    x = linspace(-1, 1, 1001);
    x2 = -cos(((2*(0:sampling(sample) - 1) + 1)/(2*(sampling(sample) - 1) + 2))*pi);
    y = 1./(1 + 25*(x2.^2));
    y1 = 1./(1 + 25*(x.^2));
    sum = 0;

    for i = 1:1001
        sum = 0;
        for j = 1:sampling(sample)
            L = 1;
            for k = 1:sampling(sample)
                if j ~= k
                    L = L*((x(i)-x2(k))/(x2(j)-x2(k)));
                end
            end
            sum = sum + (y(j)*L);
        end
        y2(i) = sum;
    end

    error = y2 - y1;
    l2 = norm(error, 2);
    linf = norm(error, inf);

    fprintf("Using %.0f Chebyshev points; ", sampling(sample))
    fprintf("L_2_norm = %f ", l2)
    fprintf("& L_inf_norm = %f \n", linf)
    
    l2_arr(sample, 1) = l2;     % [0.333119882806597 0.006255723108271 2.211423697165300e-06 2.765427572135907e-13]
    linf_arr(sample, 1) = linf; % [0.037590328892906 7.070159314979119e-04 2.499358247032291e-07 3.330669073875470e-14]

    figure (12 + sample)
    hold on; grid on
    plot(x, y2, 'LineWidth', width)
    plot(x, y1, 'LineWidth', width)
    title(['Interpolatory Polynomial Approximation sampling at ', num2str(sampling(sample)), ' Chebyshev points'], 'FontSize', tsize);
    xlabel('x', 'FontSize', lsize);
    ylabel('y', 'FontSize', lsize);
    hold off

end

figure (18)
hold on; grid on
xlim([10 160]); % ylim([-0.5*(10^28) 11*(10^28)])
plot1 = plot(sampling, log(l2_arr), 'LineWidth', width);
plot2 = plot(sampling, log(linf_arr), 'LineWidth', width);
title('Norm (log scale) vs sampling points', 'FontSize', tsize);
xlabel('Number of Chebyshev points used', 'FontSize', lsize);
ylabel('Norm (log scale) of the error', 'FontSize', lsize);
L = legend([plot1 plot2], '$L_2$ norm', '$L_{\infty}$ norm', 'Location', 'best');
L.FontSize = legsize;
