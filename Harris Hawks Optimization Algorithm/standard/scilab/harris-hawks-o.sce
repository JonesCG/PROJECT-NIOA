
clear
clc

exec('objective.sce');

POPSIZE = 50
DIM = 2
UPPER = 100
LOWER = -100
MAXITER = 50

disp("RUNNING...")
rand('seed',getdate('s'))
HAWK = rand(POPSIZE,DIM).*(UPPER-LOWER) + LOWER
FIT = F1(HAWK)
[BESTFIT IND] = min(FIT)
GBESTFIT = BESTFIT
GBESTPOS = HAWK(IND,:)

mul = 0.01
// FOR LEVY FLIGHT
bet = 1   // beta constant
num = gamma(1+bet)*sin(%pi*bet*0.5)
den = gamma(0.5*(1+bet))*bet*2^(0.5*(bet-1))
sigma = (num/den)^(1/bet)

xtitle("INITILIZATION")
square(LOWER,LOWER,UPPER,UPPER)
plot(HAWK(:,1),HAWK(:,2),'b.')
plot(GBESTPOS(:,1),GBESTPOS(:,2),'g.')
xs2png(gcf(),'gif/hho/ITER0.png')
clf()

for ITER=1:MAXITER
    E0 = grand(POPSIZE,1,"unf",-1,1)
    J = (1 - rand(POPSIZE,1)).*2
    E = E0.*(2*(1-((ITER-1)/MAXITER)))  // UPDATE E
    R = rand(POPSIZE,1)
    
    IND1 = find(abs(E)>=1)
    IND2 = find(abs(E)<1 & R>=0.5 & abs(E)>=0.5 )
    IND3 = find(abs(E)<1 & R>=0.5 & abs(E)<0.5 )
    IND4 = find(abs(E)<1 & R<0.5 & abs(E)>=0.5 )
    IND5 = find(abs(E)<1 & R<0.5 & abs(E)<0.5 )
    
    //  EXPLORATION
    if length(IND1)>0
        Q = rand(length(IND1),1)
        Q1 = find(Q>=0.5)
        Q2 = find(Q<0.5)
        
        if length(Q1)>0
            R1 = rand(length(Q1),DIM)
            R2 = rand(length(Q1),DIM)
            PERM = grand(1,'prm',1:POPSIZE)
            RAND = PERM(1:length(Q1))
            HAWK(IND1(Q1),:) = HAWK(RAND,:) - abs( HAWK(RAND,:) - HAWK(IND1(Q1),:).*R2.*2 ).*R1
        end
        if length(Q2)>0
            M = mean(HAWK,'r') // Mean Position of all Hawks
            R3 = rand(length(Q2),DIM)
            R4 = rand(length(Q2),DIM)
            RABMEAN = repmat(GBESTPOS-M,[length(Q2) 1])
            HAWK(IND1(Q2),:) = RABMEAN - (LOWER + R4.*(UPPER-LOWER)).*R3
        end
        HAWK(IND1,:) = ammend(HAWK(IND1,:),UPPER,LOWER)
        FIT(IND1,:) = F1(HAWK(IND1,:))
    end
    
    //  EXPLOITATION
    //  SOFT BESIEGE
    if length(IND2)>0
        RABREP = repmat(GBESTPOS,[length(IND2) 1]) // repeat rabbit
        DELTA_POS = RABREP - HAWK(IND2,:)
        E_VAL = repmat(E(IND2),[1 DIM])
        J_VAL = repmat(J(IND2),[1 DIM])
        HAWK(IND2,:,:) = DELTA_POS - E_VAL.*abs( J_VAL.*RABREP - HAWK(IND2,:) )
        
        HAWK(IND2,:) = ammend(HAWK(IND2,:),UPPER,LOWER)
        FIT(IND2,:) = F1(HAWK(IND2,:))
    end
    //  HARD BESIEGE
    if length(IND3)>0
        RABREP = repmat(GBESTPOS,[length(IND3) 1])
        DELTA_POS = RABREP - HAWK(IND3,:)
        E_VAL = repmat(E(IND3),[1 DIM])
        HAWK(IND3,:) = RABREP - E_VAL.*abs(DELTA_POS)
        
        HAWK(IND3,:) = ammend(HAWK(IND3,:),UPPER,LOWER)
        FIT(IND3,:) = F1(HAWK(IND3,:))
    end
    //  SOFT BESIEGE WITH RAPID PROGRESSIVE DIVE
    if length(IND4)>0
        RABREP = repmat(GBESTPOS,[length(IND4) 1])
        E_VAL = repmat(E(IND4),[1 DIM])
        J_VAL = repmat(J(IND4),[1 DIM])
        HY = RABREP - E_VAL.*abs(J_VAL.*RABREP - HAWK(IND4,:))
        // LEVY FLIGHT 
        u = rand(length(IND4),DIM)
        v = rand(length(IND4),DIM)
        LF = ( (u.*sigma)./(abs(v).^(1/bet)) ).*mul
        HZ = HY + rand(length(IND4),DIM).*LF
        
        HY = ammend(HY,UPPER,LOWER)
        HZ = ammend(HZ,UPPER,LOWER)
        FY = F1(HY)
        FZ = F1(HZ)
        IMPY = find(FY<FIT(IND4))
        if length(IMPY)>0
            HAWK(IND4(IMPY),:) = HY(IMPY,:)
            FIT(IND4(IMPY)) = FY(IMPY)
        end
        IMPZ = find(FZ<FIT(IND4))
        if length(IMPZ)>0
            HAWK(IND4(IMPZ),:) = HZ(IMPZ,:)
            FIT(IND4(IMPZ)) = FZ(IMPZ)
        end
    end
    //  HARD BESIEGE WITH RAPID PROGRESSIVE DIVE
    if length(IND5)>0
        M = mean(HAWK,'r')
        MREP = repmat(M,[length(IND5) 1])
        RABREP = repmat(GBESTPOS,[length(IND5) 1 1])
        E_VAL = repmat(E(IND5),[1 DIM])
        J_VAL = repmat(J(IND5),[1 DIM])
        HY = RABREP - E_VAL.*abs(J_VAL.*RABREP - MREP)
        u = rand(length(IND5),DIM)
        v = rand(length(IND5),DIM)
        LF = ( (u.*sigma)./(abs(v).^(1/bet)) ).*mul
        HZ = HY + rand(length(IND5),DIM).*LF
        HY = ammend(HY,UPPER,LOWER)
        HZ = ammend(HZ,UPPER,LOWER)
        FY = F1(HY)
        FZ = F1(HZ)
        IMPY = find(FY<FIT(IND5))
        if length(IMPY)>0
            HAWK(IND5(IMPY),:) = HY(IMPY,:)
            FIT(IND5(IMPY)) = FY(IMPY)
        end
        IMPZ = find(FZ<FIT(IND5))
        if length(IMPZ)>0
            HAWK(IND5(IMPZ),:) = HZ(IMPZ,:)
            FIT(IND5(IMPZ)) = FZ(IMPZ)
        end
    end
    
    [BESTFIT IND] = min(FIT,'r')
    if BESTFIT<GBESTFIT
        GBESTFIT = BESTFIT
        GBESTPOS = HAWK(IND,:)
    end
    FITRUN(ITER) = GBESTFIT
    
    //gcf().axes_size = [500 500]
    xtitle("ITER "+string(ITER))
    square(LOWER,LOWER,UPPER,UPPER)
    plot(HAWK(:,1),HAWK(:,2),'b.')
    plot(GBESTPOS(:,1),GBESTPOS(:,2),'g.')
    xs2png(gcf(),'gif/hho/ITER'+string(ITER)+'.png')
    clf()
    
    
end
