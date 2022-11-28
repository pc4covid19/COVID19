%% supply values High and Low
PicL=0; %pro-inflam
debL=50; %debris
AicL=0;%anti-inflammatory cytokine
IgL=0;%Antibody
CD8L=0;%CD8
CD4L=0;%CD4
DCL=28;%DC
DML=0;
TCL=10;
TH1L=1;
TH2L=1;
BcL=20;
PsL=0;
DLL=0;
L=[PicL debL AicL IgL CD8L CD4L DCL DML TCL TH1L TH2L BcL PsL DLL];
PicH=0.07; %pro-inflam
debH=84; %debris
AicH=1.2;%anti-inflammatory cytokine
IgH=0.19;%Antibody
CD8H=250;%CD8
CD4H=9;%CD4
DCH=48;%DC
DMH=0.1;
TCH=1.1E5;
TH1H=13;
TH2H=3.5E3;
BcH=5.2e03;
PsH=2.5e04;
DLH=0.084;
H=[PicH debH AicH IgH CD8H CD4H DCH DMH TCH TH1H TH2H BcH PsH DLH];
%% supply x value
x=0.1;
%% solve for new values
F=(1-x).*L+x.*H;

Fciel=[F(1:4) ceil(F(5:7)) F(8:end)];
disp(['pro ',num2str(Fciel(1))])
disp(['mac ',num2str(Fciel(2))])
disp(['anti-im ',num2str(Fciel(3))])
disp(['Ig ',num2str(Fciel(4))])
disp(['CD8 ',num2str(Fciel(5))])
disp(['CD4 ',num2str(Fciel(6))])
disp(['DC ',num2str(Fciel(7))])
disp(['DClymph ',num2str(Fciel(8))])
disp(['TC ',num2str(Fciel(9))])
disp(['TH1 ',num2str(Fciel(10))])
disp(['TH2 ',num2str(Fciel(11))])
disp(['BC ',num2str(Fciel(12))])
disp(['PS ',num2str(Fciel(13))])
disp(['DL ',num2str(Fciel(14))])