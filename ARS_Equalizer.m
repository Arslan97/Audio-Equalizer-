%%
clc
clear all
close all
%%
% [voice_t,Fs]= audioread('Arslan.wav');
file_name = input('Enter the File Name..\n','s')
[voice_t,Fs]= audioread(file_name);
sound (voice_t,Fs);
T = length(voice_t);

%original signal time domain plot
timeAxis = linspace(0, T/Fs , T) ;
figure;
subplot(2,1,1); plot(timeAxis,voice_t) ; grid on ; xlabel('time(s)') ; title('Original Signal Time Domain');
 
%Take F.T to analyze original signal in frequency domain
voice_f = abs(fft(voice_t));
voice_f = voice_f(1:T/2+1);
f = 0:Fs/length(voice_t):Fs/2;
subplot(2,1,2); plot(f,voice_f); grid on ; xlabel('freq(Hz)'); title('Original Signal Freq Domain');

%normalize sample freq
Fn = Fs/2 ;
%%  
go=true; %flag
while go
    FilterOption=input('Choose Filter Type \n (1)IIR Butterworth \n (2)FIR  \n (3)Exit \n ');
    switch FilterOption
        case (1)
            N1=4;                                %IIR filter order
            
            %low pass IIR (0-170)

            wc = 170/(Fn);                        %normalize digital cutoff freq
            [n, d] = butter(N1, wc,'low');        %evaluate coef
            [h1, f1] = freqz(n, d, 512 , Fs);     %freq response
            mag1 = abs(h1);                       %magnitude response
            phase1 = angle(h1) * 180 / pi;        %phase response
            [h11,t11] = impz(n,d);                %impulse response
            [h12,t12] = stepz(n,d);               %step response
            transf1 = tf(n, d);                   %evaluate transfer function
            
            %Enter Gain Input
            promptString = sprintf('Please enter the gain of (%d) - (%d) band : ', 0, 170);
            gain=input(promptString);
            
            %plot LPF IIR
            figure;
            suptitle('0-170 Hz Butterworth IIR');
            subplot(3,2,1);plot(f1,mag1); grid on ; title('Magnitude');xlabel('Freq Hz') ;
            subplot(3,2,2);plot(f1,phase1); grid on ;title('Phase');xlabel('Freq Hz');
            subplot(3,2,3);stem(t11,h11); grid on ;title('Impulse Response');xlabel('samples'); ylabel('Amplitude');
            subplot(3,2,4);stem(t12,h12); grid on ;title('Step Response');xlabel('samples'); ylabel('Amplitude');
            subplot(3,2,[5 6]); zplane(n,d) ; grid on ;title('Zeros and Poles');

            %composite signal            
            y= db2mag(gain)*filter(n, d, voice_t);
 
            %IIR (170Hz--16kHz)
            %storing cutoff freq vector
            cutOff =[170 310 600 1000 3000 6000 12000 14000 16000];          
            for i=1:8

                wc2 = [cutOff(i) cutOff(i+1)]/(Fn);   %normalized digital cutoff freq
                [n2,d2] = butter(N1,wc2,'bandpass');  %evaluate coef
                [h2,f2] = freqz(n2,d2,512,Fs);        %freq response
                mag2= abs(h2);                        %magnitude response
                phase2= angle(h2)*180/pi;             %phase response
                [h21,t21] = impz(n2,d2);              %impulse response
                [h22,t22] = stepz(n2,d2);             %step response
                transf2=tf(n2,d2);                    %evaluate transfer function
                
                %Enter Gain Input
                promptString = sprintf('Please enter the gain of (%d) - (%d) band : ', cutOff(i), cutOff(i+1));
                gain=input(promptString);
                
                figure;   
                suptitle([num2str(cutOff(i)),'~', num2str(cutOff(i+1)),'Hz Butterworth IIR']);
                subplot(3,2,1);plot(f2,mag2); grid on ; title('Magnitude');xlabel('Freq(Hz)');
                subplot(3,2,2);plot(f2,phase2); grid on ;title('Phase');xlabel('Freq(Hz)');
                subplot(3,2,3);stem(t21,h21) ; grid on ;title('Impulse Response');xlabel('samples'); ylabel('Amplitude');
                subplot(3,2,4);stem(t22,h22); grid on ;title('Step Response');xlabel('samples'); ylabel('Amplitude');
                subplot(3,2,[5 6]);zplane(n2,d2); grid on ;title('Zeros and Poles');

                %composite signal                
                x= db2mag(gain)*filter(n2, d2, voice_t);
                y= y+x ;
            end
            
            %resampling of the output
            promptString2 = sprintf('Please enter the output sample rate multiplicator :');
            Fs_mul=input(promptString2);
                
            if (Fs_mul==2) %if output sample rate is doubled
                Y = downsample(voice_t,2);
                sound(Y,Fs)
                y = Y;
                 Ymag = abs(fft(y));
                Ymagh = Ymag(1:length(Ymag)/2+1);
                %play and save
                sound (y,2*Fs);
                filename1 = 'ericIIR.wav';
                audiowrite(filename1,y,Fs);
                figure;
               timeAxis = downsample(timeAxis/2,2);
                subplot(2,1,1);plot(timeAxis,y);grid on ; title('filtered signal (IIR) time domain');xlabel('time(s)');
                %Ymagh = upsample(Ymagh,2);
                YO = upsample(Ymagh,2);
                YO = YO(1:length(YO)-1);
               subplot(2,1,2);plot(f,YO);grid on ;title('filtered signal (IIR) freq domain');xlabel('freq(Hz)');
            
            elseif (Fs_mul==1/2) %if output sample rate is halved                              
                Y = upsample(voice_t,2);
                sound(Y,Fs)
                y = Y;
                Ymag = abs(fft(y));
                Ymagh = Ymag(1:length(Ymag)/2+1);
                %play and save
                %sound (y,2*Fs);
                filename1 = 'ericIIR.wav';
                audiowrite(filename1,y,Fs);
                figure;
                timeAxis = upsample(timeAxis*2,2);
                
                subplot(2,1,1);plot(timeAxis,y);grid on ; title('filtered signal (IIR) time domain');xlabel('time(s)');
                %Ymagh = downsample(Ymagh,2);
                YO = downsample(Ymagh,2);
                subplot(2,1,2);plot(f,YO);grid on ;title('filtered signal (IIR) freq domain');xlabel('freq(Hz)');
            else   %if output sample rate is normal
                y = resample(y , 1 ,1);
                Ymag = abs(fft(y));
                Ymagh = Ymag(1:length(Ymag)/2+1);
                %play and save
                sound (y,Fs);
                filename1 = 'ericIIR.wav';
                audiowrite(filename1,y,Fs);
                figure;
                subplot(2,1,1);plot(timeAxis,y);grid on ; title('filtered signal (IIR) time domain');xlabel('time(s)');
                subplot(2,1,2);plot(f,Ymagh);grid on ;title('filtered signal (IIR) freq domain');xlabel('freq(Hz)');
            end
            
            
            %%%%%%%%%% FIR Filter %%%%%%%%%
        case(2)
            N = 64;   %FIR filter order
            
            %choose window method 
            Window_Option={'hamming','kaiser','rectwin'};
            disp('Here is the Windowing option:');
            disp('1: hamming');
            disp('2: kaiser');
            disp('3: rectwin');
            WindowOption=input('WindowOption=');
            
            %Warning Error
            while WindowOption ~=1 && WindowOption ~=2 && WindowOption ~=3 
                disp('Choose only 1,2,or 3');
                WindowOption=input('WindowOption=');
            end
            
            %defining the user window option
            switch Window_Option{WindowOption}
                case 'hamming'        
                    userWindow=hamming(N+1);
                case 'kaiser'
                    userWindow=kaiser(N+1);
                case 'rectwin'
                    userWindow=rectwin(N+1);
                otherwise
                    disp('Wrong option. EXIT');
                    exit 
            end
            
            
            %low pass FIR (0-170)
            wc = 170/(Fn);                      %normalize digital cutoff freq            
            n = fir1(N,wc,'low',userWindow);    %Find coefficient
            [h1, f1] = freqz(n,1,512,Fs);       %frequency response
            [h11,t11] = impz(n,1);              %impulse response
            [h12,t12] = stepz(n,1);             %step response
            mag1 = abs(h1);                     %magnitude response 
            phase1 = angle(h1)*180/pi;          %phase response
            transf1 = tf(n,1);                  %transfer function
            
            %Enter gain input
            promptString = sprintf('Please enter the gain of (%d) - (%d) band : ', 0, 170);
            gain=input(promptString);
            
            figure;
            suptitle('0-170 Hz FIR');
            subplot(3,2,1);plot(f1,mag1); grid on ; title('Magnitude'); xlabel('Freq(Hz)');
            subplot(3,2,2);plot(f1,phase1); grid on ;title('Phase');xlabel('Freq(Hz)');
            subplot(3,2,3);stem(t11,h11); grid on ;title('Impulse Response');xlabel('samples'); ylabel('Amplitude');
            subplot(3,2,4);stem(t12,h12); grid on ;title('Step Response');xlabel('samples'); ylabel('Amplitude');
            subplot(3,2,[5 6]);zplane(n,1); grid on ;title('Zeros and Poles');

            %composite signal
            y= db2mag(gain)*filter(n,1,voice_t);
            
            %BPF FIR (170--16kHz)
            %storing cutoff freq vector
            cutOff =[170 310 600 1000 3000 6000 12000 14000 16000]; 
            for i=1:8                                                            
                wc2 = [cutOff(i) cutOff(i+1)]/(Fn);               %normalized digital cutoff freq
                n2 =fir1(N ,wc2,'bandpass',userWindow) ;          %normalize coef
                [h2, f2] = freqz(n2,1,512,Fs);                    %frequency response
                [h21,t21] = impz(n2,1);                           %impulse response
                [h22,t22] = stepz(n2,1);                          %step response
                mag2 = abs(h2);                                   %magnitude response
                phase2 = angle(h2)*180/pi;                        %phase response
                transf2 = tf(n2,1);                               %transfer function           
                
                %Enter Gain Input
                promptString = sprintf('Please enter the gain of (%d) - (%d) band : ', cutOff(i), cutOff(i+1));
                gain=input(promptString);
                
                figure;
                suptitle([num2str(cutOff(i)),'~', num2str(cutOff(i+1)),'Hz FIR']);
                subplot(3,2,1);plot(f2,mag2); grid on ; title('Magnitude');xlabel(' Freq Hz');
                subplot(3,2,2);plot(f2,phase2); grid on ;title('Phase'); xlabel(' Freq Hz');
                subplot(3,2,3);stem(t21,h21); grid on ;title('Impulse Response');xlabel('samples'); ylabel('Amplitude');
                subplot(3,2,4);stem(t22,h22); grid on ;title('Step Response');xlabel('samples'); ylabel('Amplitude');
                subplot(3,2,[5 6]);zplane(n2,1); grid on ;title('Zeros and Poles');

                %composite signal               
                x= db2mag(gain) * filter(n2, 1, voice_t);
                y= y+x ;
            end
            y = resample(y , 1 ,1);            
            Ymag = abs(fft(y));
            Ymagh = Ymag(1:length(Ymag)/2+1);
            %play and save
            sound (y,Fs)
            filename2 = 'ericFIR.wav';
            audiowrite(filename2,y,Fs);
            figure;
            subplot(2,1,1);plot(timeAxis,y); grid on ;title('filtered signal (FIR) time domain'); xlabel('time(s)');
            subplot(2,1,2);plot(f,Ymagh); grid on ;title('filtered signal (FIR) freq domain'); xlabel('freq(Hz)');
      
        
        case(3) 
            go = false;            
    end
end