clear; close all
rng default
M = 2;                 % Modulation order
k = log2(M);            % Bits per symbol
EbNoVec = (-8:0)';       % Eb/No values (dB)
numSymPerFrame = 1e4;   % Number of QAM symbols per frame
berEstSoft = zeros(size(EbNoVec)); 
berEstHard = zeros(size(EbNoVec));
berEstUncoded = zeros(size(EbNoVec));
trellis = poly2trellis(9,[557 663 711]);
tbl = 32;
rate = 1/3;
for n = 1:length(EbNoVec)
    % Convert Eb/No to SNR
    snrdB = EbNoVec(n);%+ 10*log10(k*rate);
    % Noise variance calculation for unity average signal power.
    noiseVar = 10.^(-snrdB/10);
    % Reset the error and bit counters
    [numErrsSoft,numErrsHard,numErrsUncoded,numBits] = deal(0);
    for packet = 1:1e2
     fprintf('SNR=%d, packet=%d\n',EbNoVec(n),packet);
   % while numErrsSoft < 100 && numBits < 1e7
        % Generate binary data and convert to symbols
        dataIn = randi([0 1],numSymPerFrame*k,1);        
        % Convolutionally encode the data
        dataEnc = convenc(dataIn,trellis);        
        % MPSK modulate
           dataTx = reshape(dataEnc,length(dataEnc)/k,k);
           dataMtx = bi2de(dataTx);
           
           dataUncoded = reshape(dataIn,length(dataIn)/k,k);
           UncodedMtx = bi2de(dataUncoded);
           
           if M ==4
            txSig = pskmod(dataMtx,M,pi/4,'gray');
            rxSig = awgn(txSig,snrdB,'measured');       
            rxDemod= pskdemod(rxSig,M,pi/4,'gray'); %Ó²ÅÐ¾ö
            
            txUncoded = pskmod(UncodedMtx,M,pi/4,'gray');
            rxUncoded = awgn(txUncoded,snrdB,'measured');
            UncodedDemod = pskdemod(rxUncoded,M,pi/4,'gray');% ÎÞ±àÂë
           elseif M==2
            txSig = pskmod(dataMtx,M);  
            rxSig = awgn(txSig,snrdB,'measured');
            rxDemod= pskdemod(rxSig,M); %Ó²ÅÐ¾ö
            
            txUncoded = pskmod(UncodedMtx,M);
            rxUncoded = awgn(txUncoded,snrdB,'measured');
            UncodedDemod = pskdemod(rxUncoded,M); % ÎÞ±àÂë
           end
%%        % Hard decision   
            rxMtx = de2bi(rxDemod);
            rxDataHard = reshape(rxMtx,length(dataEnc),1);
            % Viterbi decode the demodulated data
            dataHard = vitdec(rxDataHard,trellis,tbl,'cont','hard');
            numErrsInFrameHard = biterr(dataIn(1:end-tbl),dataHard(tbl+1:end));    
            
%%       % soft decision
         if M == 4
            [~,qcodeR]=quantiz(real(rxSig),[-0.75,-0.5,-0.25,0,0.25,0.5,0.75],7:-1:0);
            [~,qcodeI]=quantiz(imag(rxSig),[-0.75,-0.5,-0.25,0,0.25,0.5,0.75],7:-1:0);
            rxDataSoft =[qcodeR qcodeI]';
            dataSoft= vitdec(rxDataSoft,trellis,tbl,'cont','soft',3);

         elseif M==2
            [~,qcodeR]=quantiz(real(rxSig),[-0.75,-0.5,-0.25,0,0.25,0.5,0.75],7:-1:0);
            rxDataSoft = qcodeR';
            dataSoft= vitdec(rxDataSoft,trellis,tbl,'cont','soft',3);
         end
        %dataSoft = vitdec(real(rxSig),trellis,tbl,'cont','unquant');      
        numErrsInFrameSoft = biterr(dataIn(1:end-tbl),dataSoft(tbl+1:end));
%%      % uncoded decision
        rxDataUncoded=de2bi(UncodedDemod);
        dataUncoded = reshape(rxDataUncoded,length(dataIn),1);
        numErrsInFrameUncoded=biterr(dataUncoded,dataIn);
        
%%        % Increment the error and bit counters
        numErrsHard = numErrsHard + numErrsInFrameHard;
        numErrsSoft = numErrsSoft + numErrsInFrameSoft;
        numErrsUncoded = numErrsUncoded + numErrsInFrameUncoded;
        numBits = numBits + numSymPerFrame*k;

    end
    
    % Estimate the BER for both methods
    berEstSoft(n) = numErrsSoft/numBits;
    berEstHard(n) = numErrsHard/numBits;
    berEstUncoded(n) = numErrsUncoded/numBits;
end
semilogy(EbNoVec,[berEstSoft berEstHard berEstUncoded],'-*')
hold on
if M == 4 
% semilogy(EbNoVec,berawgn(EbNoVec,'psk',M,'diff'))
% semilogy(EbNoVec,berawgn(EbNoVec,'psk',M,'nondiff'))
% semilogy(EbNoVec,berawgn(EbNoVec,'oqpsk','diff'))
semilogy(EbNoVec,berawgn(EbNoVec,'dpsk',M))
elseif M == 2
semilogy(EbNoVec,berawgn(EbNoVec,'psk',M,'nondiff'))
end
legend('Soft','Hard','Uncoded','theoretic')
grid
xlabel('Eb/No (dB)')
ylabel('Bit Error Rate')