
clc
clear

exec('objective.sce');
disp('RUNNING...')

POPSIZE = 50
DIM = 2
UPPER = 100
LOWER = -100
CSIZE = POPSIZE/2
ASIZE = POPSIZE/2
MAXITER = 50
rand('seed',getdate('s'))
C = rand(CSIZE,DIM).*(UPPER-LOWER) + LOWER
CFIT = F1(C)
A = rand(ASIZE,DIM).*(UPPER-LOWER) + LOWER
AFIT = F1(A)
[BESTFIT IND] = min([CFIT;AFIT])
GBESTFIT = BESTFIT
GBESTPOS = [C;A](IND,:,:)

xtitle("INITIALIZATION")
square(LOWER,LOWER,UPPER,UPPER)
plot([A;C](:,1),[A;C](:,2),'b.')
plot(GBESTPOS(:,1),GBESTPOS(:,2),'g.')
xs2png(gcf(),'gif/imo/ITER0.png')
clf()

for ITER=1:MAXITER
    //disp(ITER)
    // sort cations and anions
    [CFIT CSORTIND] = gsort(CFIT,'lr','i')
    C = C(CSORTIND,:)
    [AFIT ASORTIND] = gsort(AFIT,'lr','i')
    A = A(ASORTIND,:)
    // get best and worst cation
    CBESTFIT = CFIT(1)
    CBEST = C(1,:)
    CWORSTFIT = CFIT(CSIZE)
    CWORST = C(CSIZE,:)
    // get best and worst anion
    ABESTFIT = AFIT(1)
    ABEST = A(1,:)
    AWORSTFIT = AFIT(ASIZE)
    AWORST = A(ASIZE,:)
    //calculate force of cation and anion
    AD = abs(A - repmat(CBEST,[ASIZE 1]))
    AF = 1./(1 + exp((-0.1)./(AD)))
    CD = abs(C - repmat(ABEST,[CSIZE 1]))
    CF = 1./(1 + exp((-0.1)./(CD)))
    // update cation and anion
    A = A + AF.*(repmat(CBEST,[ASIZE 1]) - A)
    C = C + CF.*(repmat(ABEST,[CSIZE 1]) - C)
    
    if CBESTFIT>=CWORSTFIT/2 & ABESTFIT>=AWORSTFIT/2
        Q1 = grand(ASIZE,DIM,"unf",-1,1)
        Q2 = grand(CSIZE,DIM,"unf",-1,1)
        AT = A + Q1.*repmat(CBEST,[ASIZE 1])
        CT = C + Q2.*repmat(ABEST,[CSIZE 1])
        if rand()>0.5
            AT = A + Q1.*(repmat(CBEST,[ASIZE 1])-1)
        end
        if rand()>0.5
            CT = C + Q2.*(repmat(ABEST,[CSIZE 1])-1)
        end
        A = AT
        C = CT
        if rand()<0.1
            C = rand(CSIZE,DIM).*(UPPER-LOWER) + LOWER
            A = rand(ASIZE,DIM).*(UPPER-LOWER) + LOWER
        end
    end
    
    A = ammend(A,UPPER,LOWER)
    C = ammend(C,UPPER,LOWER)
    //A = penalty(A,ubx,lbx,lby,uby)
    //C = penalty(C,ubx,lbx,lby,uby)
    AFIT = F1(A)
    CFIT = F1(C)
    // get GBEST
    [BESTFIT IND] = min([AFIT;CFIT])
    if BESTFIT < GBESTFIT
        GBESTFIT = BESTFIT
        GBESTPOS = [A;C](IND,:,:)
    end
    FITRUN(ITER) = GBESTFIT
    ALL = [A;C]
    
    xtitle("ITER "+string(ITER))
    square(LOWER,LOWER,UPPER,UPPER)
    //plot(A(:,1), A(:,2), 'bx')
    //plot(C(:,1),C(:,2),'bx')
    plot(ALL(:,1),ALL(:,2),'b.')
    plot(GBESTPOS(:,1),GBESTPOS(:,2),'g.')
    xs2png(gcf(),'gif/imo/ITER'+string(ITER)+'.png')
    clf()
    
    
end
//plot((1:MAXITER)',FITRUN,'g-')
