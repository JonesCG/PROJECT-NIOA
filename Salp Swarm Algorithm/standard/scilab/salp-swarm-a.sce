
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

SALP = rand(POPSIZE,DIM).*(UPPER-LOWER)+LOWER
FIT = F1(SALP)
[FIT SORTIND] = gsort(FIT,'lr','i')
SALP = SALP(SORTIND,:)
GBESTFIT = FIT(1)
GBESTPOS = SALP(1,:)
BESTSALP = SALP(1,:)

xtitle("INITIALIZATION")
square(LOWER,LOWER,UPPER,UPPER)
plot(SALP(:,1),SALP(:,2),'b.')
plot(GBESTPOS(:,1),GBESTPOS(:,2),'g.')
xs2png(gcf(),'gif/ssa/ITER0.png')
clf()

for ITER=1:MAXITER
    C1 = 2*exp(-1*(4*ITER/MAXITER)^2)
    //C1 = 2*%e^((-1)*(4*i/MAXITER)^2)
    C2 = rand(1,DIM)
    //C2 = grand(1,p,2,'unf',-1,1)
    C3 = rand()
    //C3 = grand(1,1,"unf",-1,1)
    
    SALP(1,:) = BESTSALP(1,:) + C1.*(C2(1,:).*(UPPER-LOWER) + LOWER)
    if C3 < 0.5
        SALP(1,:) = BESTSALP(1,:) - C1.*(C2(1,:).*(UPPER-LOWER) + LOWER)
    end
    SALP((2:POPSIZE),:) = (SALP((2:POPSIZE),:) + SALP((1:POPSIZE-1),:)).*0.5
    SALP = ammend(SALP,UPPER,LOWER)
    FIT = F1(SALP)
    
    [FIT SORTIND] = gsort(FIT,'lr','i')
    SALP = SALP(SORTIND,:)
    BESTSALP = SALP(1,:)
    if FIT(1) < GBESTFIT
        GBESTFIT = FIT(1)
        GBESTPOS = SALP(1,:)
    end
    FITRUN(ITER) = GBESTFIT
    xtitle("ITER "+string(ITER))
    square(LOWER,LOWER,UPPER,UPPER)
    plot(SALP(:,1),SALP(:,2),'b.')
    plot(GBESTPOS(:,1),GBESTPOS(:,2),'g.')
    xs2png(gcf(),'gif/ssa/ITER'+string(ITER)+'.png')
    clf()
end
