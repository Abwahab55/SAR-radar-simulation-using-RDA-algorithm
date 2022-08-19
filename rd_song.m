
clear all

thetaT=0;%ƽ̨����б�ӽ�
thetaT=thetaT*pi/180;%rad
c=3e8;%����
fc=1.5e9;%��Ƶ
lambda=c/fc;%����

%%��������
X0=400;%��λ��[-X0,X0]
H=4000;%�߶�
Y=10000;
Rb=sqrt(H^2+Y^2);
Rc=Rb/cos(thetaT);
R0=500;%������2*R0


%%��λ��(Azimuth,Cross-Range),x/u domain
v=100;%SAR ƽ̨�ٶ�
D=4;
Lsar=lambda*Rc/D;%�ϳɿ׾�����
Na=1024;
x=linspace(-X0,X0,Na);%u������
u=x/v;%��λ��
du=2*X0/v/Na;
fu=linspace(-1/2/du,1/2/du,Na);%fu������
fdc=2*v*sin(thetaT)/lambda;%Doppler��Ƶ����Ƶ��
fdr=-2*(v*cos(thetaT))^2/lambda/Rc;%Doppler��Ƶб��

%%������(Range),r/t domain
Tr=5e-6;%LFM�ź����� 5us 
Br=50e6; %LFM�źŴ��� 50MHz
Kr=Br/Tr; %��Ƶб��
Nr=512;
Rmin=sqrt(H^2+(Y-R0/2)^2);
Rmax=sqrt(H^2+(Y+R0/2)^2+(Lsar/2)^2);
dt=(2*(Rmax-Rmin)/c+Tr)/Nr; %sampling frequency
t=linspace(2*Rmin/c-Tr/2,2*Rmax/c+Tr/2,Nr);%��ʱ��
r=t*c/2;
f=linspace(-1/2/dt,1/2/dt,Nr);%f������


%%Ŀ��λ��
Ptar=[Y,0;Y-50,0;Y,50];%����������,��λ������
[Ntar,n]=size(Ptar);

%%�����ز�
s_ut=zeros(Nr,Na);
U=ones(Nr,1)*u;%����Ϊ���� %��λ��
T=t'*ones(1,Na);%��ʱ��
for k=1:Ntar
    rn=Ptar(k,1);
    xn=Ptar(k,2);
    R=sqrt(H^2+rn^2+(xn-v*U).^2);
    DT=T-2*R/c;
    phase=pi*Kr*DT.^2-4*pi/lambda*R;
    s_ut=s_ut+exp(j*phase).*(abs(DT)<Tr/2).*(abs(v*U-xn)<Lsar/2);
end

%�Ӻ�����
%wr=hamming(Nr);
Nwr=floor(Tr/dt);%��ʱ��
wr=hamming(Nwr);
N=ceil((Nr-Nwr)/2);
wr=[zeros(N,1);wr;zeros(Nr-Nwr-N,1)];
% s_ut=s_ut%.*(wr*ones(1,Na));%����Ӵ�
 s_ut=s_ut.*(wr*ones(1,Na));%����Ӵ�
wa=hamming(Na);
% s_ut=s_ut%.*(ones(Nr,1)*wa');%��λ�Ӵ�
 s_ut=s_ut.*(ones(Nr,1)*wa');%��λ�Ӵ�


%%����ѹ��
p0_t=exp(j*pi*Kr*(t-2*Rc/c).^2).*(abs(t-2*Rc/c)<Tr/2);%������LFM�ź�
p0_f=fftshift(fft(fftshift(p0_t)));
s_uf=fftshift(fft(fftshift(s_ut)));%������FFT
src_uf=s_uf.*(conj(p0_f).'*ones(1,Na));%����ѹ��
src_ut=fftshift(ifft(fftshift(src_uf)));%����ѹ������ź�?
src_fut=fftshift(fft(fftshift(src_ut).')).';%�����������

%%�����㶯У��ԭ��
F=f'*ones(1,Na);%����Ϊ���󣬾�����
FU=ones(Nr,1)*fu; %��λ��
Hrcc=exp(-j*pi/fc/fdr*FU.^2.*F);%�����㶯У������
src_fuf=fftshift(fft(fftshift(src_fut)));%����ѹ����Ķ�άƵ��?
s2rc_fuf=src_fuf.*Hrcc;
s2rc_fut=fftshift(ifft(fftshift(s2rc_fuf)));%�����������

%%��λѹ��
p0_2fu=exp(j*pi/fdr*(FU-fdc).^2);%��λ��ѹ������%%%%%��
s2rcac_fut=s2rc_fut.*p0_2fu;%��λѹ��
s2rcac_ut=fftshift(ifft(fftshift(s2rcac_fut).')).';%��λ��IFFT%%%%%ת�ã�


%
figure(1)
Sr=abs(src_ut(Nr/4:3*Nr/4,Na/2));
Srmax=max(Sr);
Sr=20*log10(Sr/Srmax);
rs=((Nr/4:3*Nr/4)*dt+2*Rmin/c-Tr/2)*c/2; %%%%%��
ys=sqrt(rs.^2-H^2);
plot(ys,Sr)
axis tight
xlabel('����/m')
ylabel('���ȹ�һ��/dB')
grid on

figure(2)
G=(abs(s2rcac_ut));
gm=max(max(G));
gn=min(min(G));
G=255/(gm-gn)*(G-gn); %%%%%��
yr=sqrt((r*cos(thetaT)).^2-H^2);
imagesc(x,yr,255-G),%%%%%��
colormap(gray)
axis tight,
xlabel('Azimuth/m')
ylabel('Range/m')
title('Ŀ��ͼ��')

figure(3)
subplot(121)
G=(abs(src_fut(Nr/2-30:Nr/2+20,:)));
gm=max(max(G));
gn=min(min(G));
G=255/(gm-gn)*(G-gn);
yr=sqrt((r*cos(thetaT)).^2-H^2);
imagesc(x,yr,255-G),colormap(gray)
grid on,axis tight,
xlabel('Doppler/Hz')
ylabel('Range/m')
title('(a)�����������')

subplot(122)
G=(abs(s2rc_fut(Nr/2-30:Nr/2+20,:)));
gm=max(max(G));
gn=min(min(G));
G=255/(gm-gn)*(G-gn);
yr=sqrt((r*cos(thetaT)).^2-H^2);
imagesc(x,yr,255-G),colormap(gray)
grid on,axis tight,
xlabel('Doppler/Hz')
ylabel('Range/m')
title('(b)У�����RD��')

figure(4)
Src_ut=abs(s2rcac_ut).';
yr=((Nr/2-50:Nr/2+50)*dt+2*Rmin/c-Tr/2)*c/2;
ys=sqrt(yr.^2-H^2);
xa=(Na/2-50:Na/2+100)*du*v-X0;
waterfall(ys,xa,Src_ut(Na/2-50:Na/2+100,Nr/2-50:Nr/2+50));
axis tight
xlabel('Range/m')
ylabel('Azimuth/m')
title('Ŀ��3Dͼ��')
