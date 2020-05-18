
clc
clear

exec('objective.sce');
disp('RUNNING...')
rand('seed',getdate('s'))

POPSIZE = 300
DIM = 100
UPPER = 100
LOWER = -100
MAXITER = 200

b=1.5

WHALE = rand(POPSIZE,DIM).*(UPPER-LOWER)+LOWER
FIT = F1(WHALE)
[BESTFIT IND] = min(FIT, 'r')
GBESTPOS = WHALE(IND,:)
GBESTFIT = BESTFIT
/*
xtitle("INITIALIZATION")
square(LOWER,LOWER,UPPER,UPPER)
plot(WHALE(:,1),WHALE(:,2),'b.')
plot(GBESTPOS(:,1),GBESTPOS(:,2),'g.')
//xs2png(gcf(),'gif/woa/ITER0.png')
clf()
*/
disp(GBESTFIT)

for ITER=1:MAXITER
    //a decreases linearly fron 2 to 0 in Eq. (2.3)
    a=2-(ITER-1)*((2-0)/MAXITER)
    //a2 linearly dicreases from -1 to -2 to calculate t in Eq. (3.12)
    a2=-1-(ITER-1)*((-1-(-2))/MAXITER)
    A = rand(POPSIZE,1).*(2*a)-a
    C = rand(POPSIZE,1).*2
    PROBABILITY = rand(POPSIZE,1)
    IND1 = find(PROBABILITY<0.5 & abs(A)<1)
    IND2 = find(PROBABILITY<0.5 & abs(A)>=1)
    IND3 = find(PROBABILITY>=0.5)
    
    //  BUBBLE-NET ATTACKING METHOD (EXPLOITATION)
    //  SHRINKING CIRCLING MECHANISM
    if length(IND1)>0
        GBESTREP = repmat(GBESTPOS,[length(IND1) 1])
        AREP = repmat(A(IND1),[1 DIM])
        CREP = repmat(C(IND1),[1 DIM])
        D = abs(GBESTREP.*CREP - WHALE(IND1,:))
        WHALE(IND1,:)=GBESTREP-AREP.*D
    end
    //  SEARCH FOR PREY (EXPLORATION)
    if length(IND2)>0
        RW = grand(1,'prm',1:POPSIZE)
        AREP = repmat(A(IND2),[1 DIM])
        CREP = repmat(C(IND2),[1 DIM])
        RANDWHALE = WHALE(RW(1:length(IND2)),:)
        D = abs(RANDWHALE.*CREP - WHALE(IND2,:))
        WHALE(IND2,:) = RANDWHALE - AREP.*D
    end
    //  SPIRAL UPDATING POSITION
    if length(IND3)>0
        L=grand(length(IND3),1,"unf",-1,1)
        //L = (a2-1).*rand(length(IND3),1) + 1
        GBESTREP = repmat(GBESTPOS,[length(IND3) 1])
        EXPO = repmat(exp(b.*L),[1 DIM])
        COSF = repmat(cos(L.*(2*%pi)),[1 DIM])
        D = abs(GBESTREP - WHALE(IND3,:))
        WHALE(IND3,:) = D.*EXPO.*COSF + GBESTREP
    end
    
    WHALE = ammend(WHALE,UPPER,LOWER)
    FIT = F1(WHALE)
    [BESTFIT IND] = min(FIT,'r')
    if BESTFIT < GBESTFIT
        GBESTFIT = BESTFIT
        GBESTPOS = WHALE(IND,:)
    end
    FITRUN(ITER) = GBESTFIT
    /*
    xtitle("ITER "+string(ITER))
    square(LOWER,LOWER,UPPER,UPPER)
    plot(WHALE(:,1),WHALE(:,2),'b.')
    plot(GBESTPOS(:,1),GBESTPOS(:,2),'g.')
    //xs2png(gcf(),'gif/woa/ITER'+string(ITER)+'.png')
    clf()
    */
end

disp(GBESTFIT)
//plot((1:MAXITER)',FITRUN,'g-')
