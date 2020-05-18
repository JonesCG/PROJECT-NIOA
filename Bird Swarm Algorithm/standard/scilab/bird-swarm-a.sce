
clear
clc

exec('objective.sce');
disp('RUNNING...')
rand('seed',getdate('s'))

POPSIZE = 50
DIM = 2
UPPER = 1000
LOWER = -1000
MAXITER = 300

FQ = 3
PB = 0.9 // 0.8 - 1.0
C = 1.5
S = 1.5
A1 = 1
A2 = 1
FL = 0.7 // 0.5 - 0.9

BIRD = rand(POPSIZE,DIM).*(UPPER-LOWER)+LOWER
FIT = F1(BIRD)
[BESTFIT IND] = min(FIT, 'r')
PBESTPOS = BIRD
PBESTFIT = FIT
GBESTPOS = BIRD(IND,:)
GBESTFIT = BESTFIT

for ITER=1:MAXITER
    disp(ITER)
    if modulo(ITER,FQ)~=0
        R = rand(POPSIZE,1)
        IND1 = find(R<PB)
        IND2 = find(R>=PB)
        // EQUATION 1
        if length(IND1) > 0
            PBESTCALC = PBESTPOS(IND1,:) - BIRD(IND1,:)
            GBESTCALC = repmat(GBESTPOS(1,:),[length(IND1) 1]) - BIRD(IND1,:)
            BIRD(IND1,:) = BIRD(IND1,:) + (PBESTCALC.*C).*rand(length(IND1),DIM) + (GBESTCALC.*S).*rand(length(IND1),DIM)
        end
        
        // EQUATION 2
        if length(IND2) > 0
            A1CALC = A1.*exp( (-1).*(PBESTFIT(IND2)./(sum(PBESTFIT)+%eps).*POPSIZE) )
            ch = 1
            while ch>0
                INDK = grand(1,'prm',1:POPSIZE)(1:length(IND2))
                ch = sum(IND2==INDK)
            end
            A2CALC = A2.*exp( ( (PBESTFIT(IND2)-PBESTFIT(INDK))./(abs(PBESTFIT(INDK)-PBESTFIT(IND2))+%eps) ).*( (POPSIZE.*PBESTFIT(INDK))./(sum(PBESTFIT)+%eps) ) )
            MEANREP = repmat(mean(BIRD,'r'),[length(IND2) 1])
            BIRD(IND2,:) = BIRD(IND2,:) + repmat(A1CALC,[1 DIM]).*(MEANREP-BIRD(IND2,:)).*rand(length(IND2),DIM) + repmat(A2CALC,[1 DIM]).*(PBESTPOS(INDK,:)-BIRD(IND2,:)).*grand(length(IND2),DIM,'unf',-1,1)
        end
    else
        [FIT SORTIND] = gsort(FIT,'lr','i')
        BIRD = BIRD(SORTIND,:)
        PBESTPOS = PBESTPOS(SORTIND,:)
        PBESTFIT = PBESTFIT(SORTIND)
        // RANDOM NUMBER OF PRODUCERS
        NUMPROD = floor(rand()*POPSIZE)
        // INDEX OF PRODUCERS
        INDPROD = [1 grand(1,'prm',(2:POPSIZE-1)) POPSIZE](1:NUMPROD)
        // INDEX OF SCROUNGERS
        INDSCR = setdiff(2:POPSIZE,INDPROD)
        
        // UPDATE PRODUCER - EQUATION 5
        BIRD(INDPROD,:) = BIRD(INDPROD,:) + BIRD(INDPROD,:).*grand(length(INDPROD),DIM,"nor",0,1)
        // UPDATE SCROUNGER - EQUATION 6
        ch = 1
        while ch>0
            INDK = grand(1,'prm',1:POPSIZE)(1:length(INDSCR))
            ch = sum(INDSCR==INDK)
        end
        BIRD(INDSCR,:) = BIRD(INDSCR,:) + (BIRD(INDK,:)-BIRD(INDSCR,:)).*rand(length(INDSCR),DIM).*FL
    end
    BIRD = ammend(BIRD,UPPER,LOWER)
    FIT = F1(BIRD)
    I = find(FIT<PBESTFIT)
    PBESTPOS(I,:) = BIRD(I,:)
    PBESTFIT(I) = FIT(I)
    [BESTFIT IND] = min(PBESTFIT)
    if BESTFIT < GBESTFIT
        GBESTFIT = BESTFIT
        GBESTPOS = PBESTPOS(IND,:)
    end
    FITRUN(ITER) = GBESTFIT
    
    xtitle("ITER "+string(ITER))
    square(LOWER,LOWER,UPPER,UPPER)
    plot(BIRD(:,1),BIRD(:,2),'b.')
    plot(GBESTPOS(:,1),GBESTPOS(:,2),'g.')
    clf()
    
end
