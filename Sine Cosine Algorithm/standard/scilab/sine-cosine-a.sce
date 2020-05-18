
clc
clear

exec('objective.sce');
disp('RUNNING...')
rand('seed',getdate('s'))

POPSIZE = 50
DIM = 2
UPPER = 100
LOWER = -100
MAXITER = 50

SA = rand(POPSIZE,DIM).*(UPPER-LOWER)+LOWER
//  EVALUATE FITNESS
FIT = F1(SA)
[BESTFIT IND] = min(FIT)
GBESTPOS = SA(IND,:)
GBESTFIT = BESTFIT

xtitle("INITIALIZATION")
square(LOWER,LOWER,UPPER,UPPER)
plot(SA(:,1),SA(:,2),'b.')
plot(GBESTPOS(:,1),GBESTPOS(:,2),'g.')
xs2png(gcf(),'gif/sca/ITER0.png')
clf()

for ITER=1:MAXITER
    // UPDATE r1, r2, r3, r4
    a=2
    // r1 linearly decreases from a to 0
    r1 = a - (a-0)*((ITER-1)/MAXITER)
    r2 = rand(POPSIZE,DIM).*(2*%pi)
    r3 = rand(POPSIZE,DIM).*2
    r4 = rand(POPSIZE,1)
    IND1 = find(r4<0.5)
    IND2 = find(r4>=0.5)
    ABSVAL = abs( r3.*repmat(GBESTPOS,[POPSIZE 1]) - SA )
    if length(IND1)>0
        SA(IND1,:) = SA(IND1,:) + r1.*sin(r2(IND1,:)).*ABSVAL(IND1,:)
    end
    if length(IND2)>0
        SA(IND2,:) = SA(IND2,:) + r1.*cos(r2(IND2,:)).*ABSVAL(IND2,:)
    end
    SA = ammend(SA,UPPER,LOWER)
    //  EVALUATE FITNESS
    FIT = F1(SA)
    [BESTFIT IND] = min(FIT)
    if BESTFIT < GBESTFIT
        GBESTFIT = BESTFIT
        GBESTPOS = SA(IND,:)
    end
    FITRUN(ITER) = GBESTFIT
    
    xtitle("ITER "+string(ITER))
    square(LOWER,LOWER,UPPER,UPPER)
    plot(SA(:,1),SA(:,2),'b.')
    plot(GBESTPOS(:,1),GBESTPOS(:,2),'g.')
    xs2png(gcf(),'gif/sca/ITER'+string(ITER)+'.png')
    clf()
end
//plot((1:MAXITER)',GBESTPLOT,'g-')
