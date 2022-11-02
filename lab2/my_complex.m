%Разложение функции t^p в комплексный ряд Фурье
%в дискретизированном виде на интервале [0,T]

clc;
%T=0.9*pi; 
T=5; 
N=64; %количество значений функции на интервале [0,T]
% M=64; %количество членов ряда Фурье
p=4; %показатель степени функции x^p 
kp=2.4;%количество периодов гармонического сигнала
C0=0;
for i=1:N+1 %генерация модельной функции
   x(i)=(2*T*(((i-1-N/2))/N)); % -T до T
   %x(i)= T*(i-1)/N;%для интервала от 0 до Т 
   f(i)=sin(2*pi*kp*(i-1)/N); % гармоническая функция 
   %f(i)=sin(x(i));
   % f(i)=(x(i)*cos(x(i)));
   % f(i)=abs(x(i));

   % x(i)=T*(((i-1))/N); %для интервала от 0 до Т
   %f(i)= (x(i))^p; %функция t^p  
   % f(i)=x(i)*exp(x(i)); 
   % f(i)=sinh(x(i));
    C0=C0+f(i);
end
C0=C0*(2/N);
max_freq = N/4;%N/4;
for M=1:max_freq

    for k=1:M
       C(k)=0; 
    end   
    for i=1:N+1
        for k=1:M
        C(k)=C(k)+f(i)*exp(-j*2*pi*k*(i-1)/N);    
        end
    end
    for k=1:M
    C(k)=C(k)*(2/N);
    end 
    %Вычисление и отображение спектра амплитуд (начало)
    for k=1:M 
    Cab(k)=abs(C(k));%коэффициенты Cab(k)- комплексные числа вида a+jb, 
    %функция abs вычисляет sqrt(a^2+b^2 )
    end
    k=1:M;
    %figure 
    %plot(k,Cab);
    %stem(Cab(1:M)); %вывод графика  дискретной последовательности данных
    %axis([1 8 -0.2 1.2]);%задание осей: [xmin xmax ymin ymax]
    %title('Амплитуды частотных составляющих спектра');
    %xlabel('Количество периодов')
    %axis tight;
    %Вычисление и отображение спектра амплитуд (конец) 
    for i=1:N+1
        y(i)=0;   
        for k=1:M    
        y(i)=y(i)+C(k)*exp(j*2*pi*k*(i-1)/N);  
        end
        y(i)=C0/2+y(i); 
    end 
    i=1:N+1;
   % figure
    %plot(i,f);
    %axis tight;
    %title('Исходная и восстановленная функция')
    %xlabel('Номер элемента массива')
    %hold on;
    %plot(i,real(y),'r-');
    %axis tight;
    %hold off;
    
    
    
    for i=2:N
      dy(i)=real(y(i))-f(i);%абсолютная погрешность восстановления
    end
    dy_proc=dy/(max(f)-min(f))*100;
    CKO=std(dy);
    CKO_proc=std(dy_proc)%СКО в процентах
    dispers(M) = CKO_proc
%     disp("for freq")
%     disp(M)
end
t = 1:max_freq
% figure
plot(t, dispers)
axis tight;
xlabel('Количество членов разложения')
ylabel("СКО, %")
T = table(t(:), dispers(:));
disp(T)


pause
close all;
clear;