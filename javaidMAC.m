%************************************************************************
%MATLAB SOURCE CODE OF AR-MAC. 
% PROPOSED AND IMPLEMENTED BY Aziz ur Rahim
% STUDENT OF MS(EE), AND MEMBER OF COMSENSE(COMMUNICATION OVER SENSOR)
% RESEARCH GROUP.
% HTTP:\\WWW.NJAVAID.COM
%************************************************************************

clc
clear all
% Supply Voltage (V)
V=3;

% Listening Mode Current (Il)
Il=18.8*0.001;

% Transmit Mode Current (It)
It=17.4*0.001;

% Receive Mode Current (Ir)
Ir=19.7*0.001;

% Clear Channel Assessment (TCCA)
TCCA = 128*0.000001;
TCC = 128*0.00000001;

% Active Radio Current Iactive (Idle Listeing) 
Ia=19.7*0.001;

% Byte Transmission Time (TB)
TB=32*0.000001;

% Sleep Current Is
Is= 0.2* 0.001;

%Data Packet Length in bytes
P = 120;

% ACK Packet Length in Bytes
ACK=11;


% The maximum time to wait for an acknowledgment frame to
% arrive following a transmitted data frame

macACKWaitDuration = 120*16e-6;

%Switching Energy (Eidleswitch)

SE=827*0.000000001;

%The minimum value of the backoff exponent (BE) in the CSMA-CA algorithm
macMinBE=3;

%The maximum value of the backoff exponent (BE) in the CSMA-CA algorithm
macMaxBE=5;

%The maximum number of backoffs the CSMA-CA algorithm will attempt before
%declaring a channel access failure

macMaxCSMABackoffs=4;

% The maximum time to wait for an acknowledgment frame to
% arrive following a transmitted data frame

macACKWaitDuration = 120*16e-6;

% RX-to-TX or TX-to-RX maximum turnaround time

aTurnaroudTime = 12*16e-6;

% number of nodes
numn=12;

% Number of Cycles
numc=1000;

% %Probability of Packet drops
PPD = [1 2 3 4 5 6 7 8 9 10 11 12]/100;

E1=zeros(1,12);
E2=zeros(1,12);
PD=zeros(1,12);

D1=zeros(1,20);
D2=zeros(1,20);

Transmit=0;
att=0;

%Time Frame 
T=1;

for i=1:1:20;
  for d=1:1:numn
    for n=1:1:numc
      
      E1(i)=E1(i)+(2*SE);
      NB=0;
      BE=macMinBE;
        while (NB<=macMaxCSMABackoffs)
          att=att+1 ;
          D1(i)=D1(i)+TCC;
          E1(i)=E1(i)+TCCA*Il*V;
          Tactive=TCCA;
            m=rand();
            if (m >(i/100))
                  Transmit=Transmit+1;
                  E1(i)=E1(i)+((P+(i)*m)*TB*It*V);
                  E1(i)=E1(i)+ (macACKWaitDuration*Il*V);
                  E1(i)=E1(i)+ (ACK*TB*Ir*V);
                  Tactive=Tactive+(P*TB)+(macACKWaitDuration)+(TB*ACK);
                  Tsleep=1-Tactive;
                  E1(i)=E1(i)+ Tsleep*Is*V;
                  NB=macMaxCSMABackoffs+1;
            else
                NB=NB+1; 
                if (NB> macMaxCSMABackoffs)
                    PD(i)=PD(i)+1;
                end
            end
        end
    end
  end
     
end


% for i=1:1:20;
%   for d=1:1:numn
%     for n=1:1:numc
%       
%       E2(i)=E2(i)+(2*SE);
%       Tactive=0;
%       
%      
%          g=rand();
%             if (g <(i/100))
%                 E2(i)=E2(i)+TCCA*Il*V;
%                 D2(i)=D2(i)+TCC;
%             end
%      
%                   E2(i)=E2(i)+(P*TB*It*V);
%                   E2(i)=E2(i)+ (macACKWaitDuration*Il*V);
%                   E2(i)=E2(i)+ (ACK*TB*Ir*V);
%                   Tactive=Tactive+(P*TB)+(macACKWaitDuration)+(TB*ACK);
%                   Tsleep=1-Tactive;
%                   E2(i)=E2(i)+ Tsleep*Is*V;
%                   
%     end
%   end
%      
% end

plot( E1,'r:+')
hold
plot( E2, 'b:*')
legend('802.15.4','AR-MAC')
xlabel('Packet Error Rate[%]');
ylabel('Energy Consumption [milliJoule]');
grid on

% figure,
% plot(D1,'r:+')
% hold
% plot( D2, 'b:*')
% legend('802.15.4','AR-MAC')
% grid on 

%**************************************************************************
% %ZIGBEE MAC
% %**************************************************************************
% %node = [0.5 1.2; 0.4 1.9; 0.2 0.5; 0.4 0.6; 0.3 1.4; 0.2 1.9; 0.6 0.3; 0.7 0.8; 0.55 0.9; 0.5 1.6; 0.6 1.9]; %10 nodes
% % plot(node(:,1),node(:,2), 'b*')
% % hold
% % grid on
% for j=1:N
%     D(1,j) =  topo_dist(1, j);
% end
% srcPktsT = zeros(1,N);
% E2 = zeros(1,N);
% counter = 1;
% macMaxCSMABackoffs = 4;
% numberMessagesTr = 0;
% pktsTr = zeros(1,N);
% pktsLo = zeros(1,N);
% while counter < 1000
%     Isources = [];
%     NB = zeros(1,N);
%     csmaTime = zeros(1,N);
%     pktStartTXTime = zeros(1,N);
%     pktEndTXTime = zeros(1,N);
%     CW = [2 2 2 2 2 2 2 2 2 2];
%     sourceState = floor(3*rand(1,N));
%     %     sourceState(1) = 20;
%     %     sourceState(2) = 1;
%     %     sourceState(3) = 1;
%     %     sourceState(4) = 0;
%     %     sourceState(5) = 0;%relay
%     %     sourceState(6) = 1;%source
%     %     sourceState(7) = 1;
%     %     sourceState(8) = 1;
%     %     sourceState(9) = 0;
%     %     sourceState(10) = 0;
%     %     sourceState(11) = 1;
%     %     sourceState(1) = 20;
%     sources = find(sourceState>=1 & sourceState < 20);
%     S = length(sources); %sources size
%     baseChannel = zeros(1,S);
%     backoff = floor(2^minBE*rand(1,S));
%     for s=1:S
%         Isources(s).srcID = sources(s);
%         i = 1;
%         for t =1:S
%             if backoff(s) == backoff(t)
%                 Isources(s).b(i) = sources(t);
%                 i = i + 1;
%             end
%         end
%     end
%     %**************************************************************************
%     %create an array structure of each source and its corresponding
%     %interferers
%     tr = zeros(1,S);
%     for s=1:S
%         if length(Isources(s).b) == 1
%             tr(s) = sources(s);
%         end
%     end
%     %scheduling the event by ascending sorting the backoff array
%     sortMatrix = zeros(S,3);
%     for i = 1:S
%         sortMatrix(i,1) = backoff(i);
%         sortMatrix(i,2) = sources(i);
%         sortMatrix(i,3) = D(sources(i))*100;
%     end
%     sortMatrix = sortrows(sortMatrix,1);
%     for i = 1:S
%         backoff(i)= sortMatrix(i,1);
%         sources(i)= sortMatrix(i,2);
%     end
%     fprintf('\nZigbee Case: Sources = ');
%     disp(sources);
%     fprintf('\nZigbee Case: Backoff = ');
%     disp(backoff);
%     src2Cdist = sortMatrix(:,3)';
%     %**************************************************************************
%     %number of transmitters on basechannel;
%     trans = find(tr>0);
%     sizeT = length(trans);
%     transmitters = zeros(1,sizeT);
%     sortedTransmitters = zeros(1,sizeT);
%     sortedBackoff = zeros(1,sizeT);
%     p = 1;
%     for s=1:S
%         if tr(s) > 0
%             transmitters(p) = tr(s);
%             p = p + 1;
%         end
%     end
%     [c,ia]=intersect(sources,transmitters);
%     trans = sort(ia);
%     sortedTransmitters = sources(trans);
%
%     fprintf('Zigbee Case: Sorted backoff: ');
%     sortedBackoff = backoff(trans);
%     disp(sortedBackoff);
%     D = D * 100;
%     sizeSB = length(sortedBackoff);
%     %for all transmitters set the channel as clear
%
%     intSrcs = setxor(sources,transmitters); %complement of transmitters in sources
%
%     baseChannel(sortedTransmitters)= 1;
%     baseChannel(intSrcs ) = 0;
%
%     fprintf('\nZigbee Case channel : ');
%     disp(baseChannel);
%     minBE = 4;
%     %initialize the sources backoff array
%     for g =1:S
%         csmaTime(g) = backoff(g);
%     end
%     %Initialize sensors energy
%     for i = 1:N
%         if r == 1 % coordinator
%             E2(i)= 0;
%         else
%             E2(i)=E2(i)+(2*SE);
%         end
%     end
%     for i = 1:S % for all sources
%         src = sources(i);
%         NB(src) = 0;
%         CW(src) = 2;
%         csmaTime(src) = backoff(i)*BP + TCCA;
%         E2(src)=E2(src)+ TCCA*Il*V; % stays listening for TCCA time
%         E2(src)=E2(src)+ SE;
%         fprintf('\nZigbee Case: Node %d performs CCA', src);
%         while NB(src) <= macMaxCSMABackoffs
%             if baseChannel(src) == 1
%                 CW(src) = CW(src) - 1;
%                 if CW(src) == 0
%                     %Source transmits m data packets
%                     pktStartTXTime(src) = csmaTime(src) + ST + SIFS; % pkts start transmission time
%                     pktEndTXTime(src) = pktStartTXTime(src) + m*P*TB; % pkts end transmission time
%                     fprintf('\nZigbee Case: Node %d transmits data packet', src);
%                     E2(src)=E2(src) + 2*m*SIFS*Il*V + m*P*TB*It*V + m*Ack*TB*Ir*V + 2*(m+1)*SE + m*(macACKWaitDuration*Il*V);
%                     pktsTr(src) = pktsTr(src) + 1;
%                     numberMessagesTr = numberMessagesTr + 1;
%                     NB(src) = 0;
%                     CW(src) = 2;
%                     srcPktsT(src) = srcPktsT(src) + 1;
%                     break;
%                 else
%                     csmaTime(src) = csmaTime(src) + TCCA;
%                     fprintf('\nZigbee Case: Node %d performs CCA,%d', src,CW(src));
%                     E2(src)=E2(src)+ TCCA*Il*V; % stays listening for TCCA time
%                 end
%             else
%                 CW(src) = 2;
%                 NB(src) = NB(src) + 1;
%                 backoff(i) = min(backoff(i)*BP + 1, macMaxBE*BP);
%                 if NB(src) <= macMaxCSMABackoffs
%                     backoff(i) = floor(2^backoff(i)*BP*rand);
%                     csmaTime(src) = backoff(i)*BP + TCCA;
%                     E2(src)=E2(src)+ TCCA*Il*V; % stays listening for TCCA time
%                     fprintf('\nZigbee Busy Channel Case: NB[%d] = %d', src,NB(src));
%                 else
%                     fprintf('\nZigbee Case: Node %d terminates; NB[ %d ] exceeds macMaxCSMABackoffs = 5', src, src);
%                     pktsLo(src) = pktsLo(src) + 1;
%                     break;
%                 end
%             end
%         end
%     end%end for
%     counter = counter + 1;
%     if counter == 1000
%         break;
%     end
%     if numberMessagesTr == 1000
%         break;
%     end
% end%end while
% plr = zeros(1,N);
% for n = 2 : N
%     plr(n) = pktsLo(n) / ( pktsLo(n) + srcPktsT(n))*100;
% end




