%Comparison of the number of stable poles used by smiaaa and FastAAA to
%achieve target different accuracies
function [smiaaa_poles]=polecompare(f,s);
smiaaa_poles=[];
Fastaaa_poles=[];
tol2=[];
for ii=3:8
    tol=10^(-ii);
    tol2= [tol tol2];
    [faaa_poles]=FastAAAcompare(f,s,tol,false,0);
    if(isempty(faaa_poles))
        Fastaaa_poles=[Fastaaa_poles nan];
    else
        Fastaaa_poles=[Fastaaa_poles length(faaa_poles)];
    end
    
    [symaaal,pwj,symaaa,pzj,~,pfj] = symmetricsmiaaa(f,s,tol,false,1,1);
    nn=length(pwj)/2;
    [ppoles_aaa,~,ppfaaaf,~,~]=properrational(pzj.',pwj(nn+1:end),pwj(1:nn),pfj.',f,s);
    smiaaa_poles=[smiaaa_poles length(ppoles_aaa)];
end
figure()
semilogx(tol2,flip(smiaaa_poles),'b',tol2,flip(Fastaaa_poles),'r','Linewidth',1.5);
set(gca, 'XDir','reverse')
legend('Number of smiaaa poles','Number of FastAAA poles');
ylabel('Number of Stable Poles');
xlabel('Accuracy');
title("Accuracy vs Poles ISS example")

    
end