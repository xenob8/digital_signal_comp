%Фильтр Колмогорова-Винера
 
A=1; %амплитуда сигнала
Q=0.5; %СКО шума
N=1024;%количество точек расчета
kp1=5;%количество переиодов сигнала
 
clc;
q=Q*randn(1,N);%генерация одномерного массива нормально распределенного %шума с СКО=Q
for k=1:N %цикл вычисления зашумленного сигнала
  %s(k)=A*exp(-0.0003*(k-200)^2.0); %колоколообразный сигнал
  s(k)=A*sin(2*pi*kp1*k/N);%гармонический сигнал
  % s(k)=0; % сигнал прямоугольной формы
  % if (k>100)&(k<300) % сигнал прямоугольной формы
  % s(k)=A;
    x(k)=s(k)+q(k); % суммирование сигнала и шума
end
figure
plot(x(1:N));
title('Зашумленный сигнал до фильтра');
axis tight; 
 
Y=fft(x,N)/N; %БПФ  сигнала с шумом
SS1=Y.*conj(Y)/N; %спектр мощности
i=1:200;
figure
%plot(i,SS1(1:200)); 
semilogy(i,SS1(1:200)); %вывод спектра мощности сигнала с шумом
title('Частотный спектр сигнала с шумом');
 
Y=fft(s,N)/N; %БПФ сигнала без шума
SS1=Y.*conj(Y)/N; %спектр мощности сигнала без шума
 
Y1=fft(q,N)/N; %БПФ  шума
SS2=Y1.*conj(Y1)/N; %спектр мощности  шума
 
for i=1:N    
    H(i)=SS1(i)/(SS1(i)+SS2(i));%частотная характеристика оптимального фильтра
end
i=1:200;
figure
%plot(i,abs(H(1:200)));
semilogx(i,abs(H(1:200)));
%hold on
title('Частотная характеристика оптимального фильтра');
 
i=1:N;
XX1=fft(x,N); %частотный спектр сигнала с шумом
Z=ifft(XX1.*H);%свертка зашумленного сигнала с частотной хар-кой фильтра
axis tight;
 
figure
plot(i,s(1:N)); %вывод незашумленного сигнала до фильтра сигнала
title('Незашумленный сигнал до фильтра');
axis tight;    
figure
plot(i,Z(1:N)); %вывод отфильтрованного сигнала
title('Сигнал после свертки с част. хар-кой опт. фильтра');
axis tight;       
i=1:N;
DZ(i)=Z(i)-s(i);
DZ1=DZ*100/(max(s)-min(s));
SKO_total=std(DZ1)
 
i=1:N;
figure
plot(i,DZ1(1:N)); %вывод  погрешности отфильтрованного сигнала
title('Погрешность отфильтрованного сигнала');
ylabel('Полная погрешность, %'); % подпись по оси Y
axis tight;
 
pause;
close all; %закрытие всех окон графического вывода
clear;%очистка Workspace