
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

a = 0.1 // power exponent
c = 0.01 // sensory modality
prob = 0.5 // or 0.8 // switch probability

BUTTERFLY = rand(POPSIZE,DIM).*(UPPER-LOWER)+LOWER
FIT = F1(BUTTERFLY)
[BESTFIT IND] = min(FIT, 'r')
GBESTPOS = BUTTERFLY(IND,:)
GBESTFIT = BESTFIT

xtitle("INITIALIZATION")
square(LOWER,LOWER,UPPER,UPPER)
plot(BUTTERFLY(:,1),BUTTERFLY(:,2),'b.')
plot(GBESTPOS(:,1),GBESTPOS(:,2),'g.')
xs2png(gcf(),'gif/boa/ITER0.png')
clf()

for ITER=1:MAXITER
    //CALCULATE FRAGRANCE
    FRAG = (FIT.^a).*c
    RANDOM = rand(POPSIZE,1)
    IND1 = find(RANDOM<prob)
    IND2 = find(RANDOM>=prob)
    if length(IND1)>0
        R1 = rand(length(IND1),DIM)
        DIST1 = (R1.^2).*repmat(GBESTPOS,[length(IND1) 1]) - BUTTERFLY(IND1,:)
        STEP1 = DIST1.*repmat(FRAG(IND1),[1 DIM])
        BUTTERFLY(IND1,:) = BUTTERFLY(IND1,:) + STEP1
    end
    if length(IND2)>0
        INDRAND = grand(1,'prm',(1:POPSIZE)')
        R2 = rand(length(IND2),DIM)
        DIST2 = (R2.^2).*BUTTERFLY(INDRAND(1:length(IND2)),:) -  BUTTERFLY(INDRAND(1:length(IND2)),:)
        STEP2 = DIST2.*repmat(FRAG(IND2),[1 DIM])
        BUTTERFLY(IND2,:) = BUTTERFLY(IND2,:) + STEP2
    end
    
    BUTTERFLY = ammend(BUTTERFLY,UPPER,LOWER)
    FIT = F1(BUTTERFLY)
    [BESTFIT IND] = min(FIT,'r')
    if BESTFIT < GBESTFIT
        GBESTFIT = BESTFIT
        GBESTPOS = BUTTERFLY(IND,:)
    end
    
    a = 0.1 - (0.1-0.3)*(ITER/MAXITER) // UPDATE POWER EXPONENT
    c = c + 0.025/(c*MAXITER)
    
    FITRUN(ITER) = GBESTFIT
    
    //gcf().axes_size = [500 500]
    xtitle("ITER "+string(ITER))
    square(LOWER,LOWER,UPPER,UPPER)
    plot(BUTTERFLY(:,1),BUTTERFLY(:,2),'b.')
    plot(GBESTPOS(:,1),GBESTPOS(:,2),'g.')
    xs2png(gcf(),'gif/boa/ITER'+string(ITER)+'.png')
    clf()
    
end

//plot((1:MAXITER)',FITRUN,'g-')
