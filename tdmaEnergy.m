V=3;
Il=18.8*0.001;
It=17.4*0.001;
Ir=19.7*0.001;
TCCA = 128*0.000001;
TCC = 128*0.00000001;
Ia=19.7*0.001;
TB=32*0.000001;
SE=827*0.000000001;
Is= 0.2* 0.001;
P = 12;
Ack = 11;
m = 10;
ST = 0.000016;
TB = 32*0.000001;
BP = 0.00032;
TCCA = 128*0.000001;
macMinBE=4;
macMaxBE=5;
run = 10;
  N = 11;
  counter = 1;
E4 = zeros(1,N);
 for i = 1 : 20
     while counter < run
         for j = 1 : N
             if j == 1
                  E4(j) = 0;
             else
                 E4(j) = E4(j) + m*P*TB*It*V;
             end
         end
         counter = counter + 1;
     end
 end
                        plot(cumsum(E4),'g')
    hold on
 grid on
 legend('Energy Consumption "TDMA Scheme"')
 xlabel('Node Index');
 ylabel('Energy Consumption [milliJoule]');

