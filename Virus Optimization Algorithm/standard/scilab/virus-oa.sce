
clear
clc

exec('objective.sce');

// POPSIZE = 100
INSIZE = 50
DIM = 2
MAXITER = 50
UPPER = 100
LOWER = -100
NSTRONGV = 10
GRATESV = 8
GRATECV = 2
MAXPOP = 1000
IN = 1
disp("RUNNING...")
rand('seed',getdate('s'))

VIRUS = rand(INSIZE,DIM).*(UPPER-LOWER)+LOWER
FIT = F1(VIRUS)
[BESTFIT IND] = min(FIT)
GBESTPOS = VIRUS(IND,:)
GBESTFIT = BESTFIT
AV = mean(FIT)

xtitle("INITIALIZATION")
square(LOWER,LOWER,UPPER,UPPER)
plot(VIRUS(:,1),VIRUS(:,2),'b.')
plot(GBESTPOS(:,1),GBESTPOS(:,2),'g.')
xs2png(gcf(),'gif/voa/ITER0.png')
clf()

for ITER=1:MAXITER
    //disp(ITER)
	POPSIZE = size(VIRUS)(1)
	[FIT SORTIND] = gsort(FIT,'lr','i')
	VIRUS = VIRUS(SORTIND,:)
	STRONGV = VIRUS(1:NSTRONGV,:)
	COMV = VIRUS(NSTRONGV+1:POPSIZE,:)
	RANDS = grand(NSTRONGV*GRATESV,DIM,'unf',-1,1)
	RANDC = grand((POPSIZE-NSTRONGV)*GRATECV,DIM,'unf',-1,1)
	REPSV = repmat(STRONGV,[GRATESV 1])
	REPCV = repmat(COMV,[GRATECV 1])
	VIRS = REPSV + (RANDS./IN).*REPSV
	VIRC = REPCV + (RANDC).*REPCV
	NEWVIRUS = [VIRS; VIRC]
	NEWVIRUS = ammend(NEWVIRUS,UPPER,LOWER)
	NEWFIT = F1(NEWVIRUS)
	
	VIRUS = [VIRUS; NEWVIRUS]
	FIT = [FIT; NEWFIT]

    // ===================
    xtitle("ITER "+string(ITER))
    square(LOWER,LOWER,UPPER,UPPER)
    plot(VIRUS(:,1),VIRUS(:,2),'b.')
    plot(GBESTPOS(:,1),GBESTPOS(:,2),'g.')
    xs2png(gcf(),'gif/voa/ITER'+string(ITER)+'_B.png')
    clf()
    // ========================


	NAV = mean(FIT)
	
	if NAV > AV
		IN = IN + 1
	end
	AV = NAV
	
	POPSIZE = size(VIRUS)(1)
	AMT = round(rand()*(POPSIZE-NSTRONGV))
	REMOVE_LOW = length(find(FIT>AV))
	[FIT SORTIND] = gsort(FIT,'lr','i')
	VIRUS = VIRUS(SORTIND,:)
	if AMT <= REMOVE_LOW then
        VIRUS = VIRUS(1:(POPSIZE-AMT),:)
        FIT = FIT(1:(POPSIZE-AMT))
	else
		REMOVE_HI = AMT-REMOVE_LOW
		INDPERM = grand(1,'prm',(NSTRONGV+1:(POPSIZE-REMOVE_LOW)))
		INDREM = INDPERM(1:REMOVE_HI)
		STORE = setdiff((1:(POPSIZE - REMOVE_LOW)),INDREM)
		VIRUS = VIRUS(STORE,:,:)
        FIT = FIT(STORE)
	end
	
	if size(VIRUS)(1) > MAXPOP
        VIRUS = VIRUS(1:INSIZE,:)
        FIT = FIT(1:INSIZE)
    end
	
	if FIT(1) < GBESTFIT
        GBESTFIT = FIT(1)
        GBESTPOS = VIRUS(1,:)
    end
    FITRUN(ITER) = GBESTFIT
    
    //gcf().axes_size = [500 500]
    xtitle("ITER "+string(ITER))
    square(LOWER,LOWER,UPPER,UPPER)
    plot(VIRUS(:,1),VIRUS(:,2),'b.')
    plot(GBESTPOS(:,1),GBESTPOS(:,2),'g.')
    xs2png(gcf(),'gif/voa/ITER'+string(ITER)+'_C.png')
    clf()
    
end
//plot((1:MAXITER)',FITRUN,'c-')
disp('DONE')
