close all
clear all 
clc
% identifying Simulation parameters:
n=1e6;                   % Number of bits/SNR=1 million bits
SNR=0:2:30;              % Signal to noise ratio range in dB with 2Db steps
snr=10.^(SNR/10);        % Signal to noise ratio
%initializing empty matices to be filled with the sequences 
%OOK=3;
OOK=zeros(1,16);
PRK=zeros(1,16);
FSK=zeros(1,16);
OOK_mat=zeros(1,16);
PRK_mat=zeros(1,16);
FSK_mat=zeros(1,16);
QAM_mat=zeros(1,16);
% Transmitter side variables  
bits=randi([0,1],1,n);      %Generate random binary data vector(1 x n)
%we will work on ASK PSK FSK modulation according to the specified cases in the manual

%ASK special case: OOK_Modulation(0,1)
signal_OOK=bits; %No change in the bits will be required

%PSK special case: PRK modulation(-1,1)
signal_PRK=(bits.*2)-1;  %represent the 1 by 1 and the 0 bit by -1 using provided formula 

%FSK special case: Orthogonal-FSK modulation
ones=find(bits);     %find all elements that are equal to one in the bits sequence
zeros=find(~bits);   %find all elements that are equal to zero in the bits sequence
signal_FSK(ones)=i;  %replace all the ones with i which is an imaginary number
signal_FSK(zeros)=1; %replace all the zeros with 1

% modulating using matlab built in function 
OOK_mat_mod=genqammod(bits,[0  1]);   %OOK modulation using built in functions  
PRK_mat_mod=pskmod(bits,2);           %PRK modulation using built in functions  
FSK_mat_mod=genqammod(bits,[1  1i]);  %FSK modulation using built in functions 
rbdv=randi([0 16-1],1,n);              %Generate random data vector(1 x n)
QAM_mat_mod=qammod(rbdv,16);        %QAM modulation using built in functions    

%Applying noise to the bit stream with different values of SNR,receiving
%the bits stream,recovering the original bits by demodulating it 
%and calculating the BER and Pe for each type of modulated signal

for i=1:16  
    Ptx_OOK=mean(signal_OOK.^2)
    noise=(sqrt(mean(abs(signal_OOK).^2)/(2*snr(i))))*(randn(1,n)+(1j*randn(1,n)));%applying the noise on OOK sequence
    OOK_sequence=signal_OOK+noise;
  
    Ptx_PRK=mean(signal_PRK.^2)
    noise=(sqrt(mean(abs(signal_PRK).^2)/(2*snr(i))))*(randn(1,n)+(1j*randn(1,n)));%applying the noise on PRK sequence
    PRK_sequence=signal_PRK+noise;
    
    Ptx_FSK=mean(signal_FSK.^2)
    noise=(sqrt(mean(abs(signal_FSK).^2)/(2*snr(i))))*(randn(1,n)+(1j*randn(1,n)));%applying the noise on FSK sequence
    FSK_sequence=signal_FSK+noise;
    
    %aplying noise on signals after modulation using functions
    noise=(sqrt(mean(abs(OOK_mat_mod).^2)/(2*snr(i))))*(randn(1,n)+(1j*randn(1,n)));%applying the noise on OOK_tb sequence
    OOK_sequence_mat=OOK_mat_mod+noise;
    
    noise=(sqrt(mean(abs(PRK_mat_mod).^2)/(2*snr(i))))*(randn(1,n)+(1j*randn(1,n)));%applying the noise on PRK_tb sequence
    PRK_sequence_mat=PRK_mat_mod+noise;
    
    noise=(sqrt(mean(abs(FSK_mat_mod).^2)/(2*snr(i))))*(randn(1,n)+(1j*randn(1,n)));%applying the noise on FSK_tb sequence
    FSK_sequence_mat=FSK_mat_mod+noise;
    
    noise=(sqrt(mean(abs(QAM_mat_mod).^2)/(2*snr(i))))*(randn(1,n)+(1j*randn(1,n)));%applying the noise on QAM_tb sequence
    QAM_sequence_mat=QAM_mat_mod+noise;
    
    %Decide whether the Rx_sequence is ‘1’ or ‘0’ by comparing each bit with a threshold 
    OOK_demod=(real(OOK_sequence) >=0.5);%0.5 for OOK
    PRK_demod=(real(PRK_sequence) >=0);%0 for PRK
    FSK_demod=(real(FSK_sequence)<imag(FSK_sequence));%in FSK if real > imaginary component then it equal zero else equals one
  
    %demodulation using built in functions 
    OOK_mat_demod = genqamdemod(OOK_sequence_mat,[0  1]);
    PRK_mat_demod = pskdemod(PRK_sequence_mat,2);
    FSK_mat_demod = genqamdemod(FSK_sequence_mat,[1  1i]);
    rbdv_mat_demod = qamdemod(QAM_sequence_mat,16);
    
    %Calculate #of errors (BER) %[NUMBER,RATIO] = BITERR(X,Y)
    [~,OOK(i)]=biterr(bits,OOK_demod);
    [~,PRK(i)]=biterr(bits,PRK_demod);  
    [~,FSK(i)]=biterr(bits,FSK_demod);
    %for modulated using matlab 
    [~,OOK_mat(i)] = biterr(bits,OOK_mat_demod);
    [~,PRK_mat(i)] = biterr(bits,PRK_mat_demod);
    [~,FSK_mat(i)] = biterr(bits,FSK_mat_demod);
    [~,QAM_mat(i)] = symerr(rbdv,rbdv_mat_demod);
end

% constellation
n=1e6;
snrq=15;
refC = qammod(0:15,16);
constDiagram = comm.ConstellationDiagram('ReferenceConstellation',refC , 'XLimits',[-4 4],'YLimits',[-4 4]);
data = randi([0 15],n,1);
sym = qammod(data,16);
rcv = awgn(sym,snrq);
 step(constDiagram,rcv)
% Ploting all in one curve
figure
semilogy(SNR,OOK,'r-*','LineWidth',2)
hold on;
semilogy(SNR,OOK_mat,'m--o','LineWidth',2)
hold on;
semilogy(SNR,PRK,'g-*','LineWidth',2)
hold on;
semilogy(SNR,PRK_mat,'k--+','LineWidth',2)
hold on;
semilogy(SNR,FSK,'b-*','LineWidth',2)
hold on;
semilogy(SNR,FSK_mat,'y--s','LineWidth',2)
hold on;
semilogy(SNR,QAM_mat,'c-p','LineWidth',2)
title('BER vs. SNR ')
ylabel('BER')
xlabel('SNR')
legend('OOK','OOK-Function','PRK','PRK-Function','FSK','FSK-Function','QAM-Function')
grid on;
