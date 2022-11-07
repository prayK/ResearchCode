clear;
disp('请输入介质折射率n1和n2')
n1 = input('n1 =' );
n2 = input('n2 =' );
theta = 0:0.1:90;
a = theta * pi / 180;
rp = (n2 * cos(a) - n1 * sqrt(1 - (n1 / n2 * sin(a)).^2))./...
    (n2 * cos(a) + n1 * sqrt(1 - (n1 / n2 * sin(a)).^2 ));
rs = (n1 * cos(a) - n2 * sqrt(1 - (n1 / n2 * sin(a)).^2))./...
    (n1 * cos(a) + n2 * sqrt(1 - (n1 / n2 * sin(a)).^2 ));

tp = 2 * n1 * cos(a)./(n2 * cos(a) + n1 * sqrt(1 - (n1 / n2 * sin(a).^2)));
ts = 2 * n1 * cos(a)./(n1 * cos(a) + n2 * sqrt(1 - (n1 / n2 * sin(a).^2)));

figure(1);
subplot(1,2,1);
plot(theta, rp, '-', theta, rs, '- -', theta, abs(rp), ':', theta, abs(rs), '-.', 'LineWidth', 2);
legend('rp', 'rs', '|rp|', '|rs|');
xlable('入射角\theta_i');
ylable('振幅');
title(['n_1 = ', num2str(n1), 'n_2 = ', num2str(n2), '时反射系数随入射角的变换曲线']);
axis([0 90 - 1 1]);
grid on;

subplot(1,2,2);
plot(theta, tp, '-', theta, ts, '- -', theta, abs(tp), ':',theta, abs(ts), '-.', 'LineWidth', 2);
legend('tp', 'ts', '|tp|', '|ts|');
xlable('入射角\theta_i');
ylable('振幅');
title(['n_1 = ', num2str(n1), 'n_2 = ', num2str(n2), '时反射系数随入射角的变换曲线']);
if n1 < n2
    axis([0 90 0 1]);
else
    axis([0 90 0 0.35]);
end

grid on;
