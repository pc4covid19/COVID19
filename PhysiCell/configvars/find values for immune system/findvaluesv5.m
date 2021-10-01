%% supply values High and Low
PicL=0; %pro-inflam
debL=0; %debris
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
L=[PicL debL AicL IgL CD8L CD4L DCL DML TCL TH1L TH2L BcL PsL];
PicH=0.07; %pro-inflam
debH=0.01; %debris
AicH=1.2;%anti-inflammatory cytokine
IgH=0.65;%Antibody
CD8H=140;%CD8
CD4H=20;%CD4
DCH=2;%DC
DMH=0.13;
TCH=52000;
TH1H=19;
TH2H=617;
BcH=3e08;
PsH=6e07;
H=[PicH debH AicH IgH CD8H CD4H DCH DMH TCH TH1H TH2H BcH PsH];
%% supply x value
x=0.01;
%% solve for new values
F=(1-x).*L+x.*H;

Fciel=[F(1:4) ceil(F(5:7)) F(8:end)];
disp(['pro ',num2str(Fciel(1))])
disp(['debris ',num2str(Fciel(2))])
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