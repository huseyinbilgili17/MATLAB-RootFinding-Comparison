
%% Newton's Method

rootsf=[];
tstart1 = cputime %for cpu time
syms x
f1(x) = x^2-4*exp(-x) -3
df(x) = diff(f1);
double(df(3));
x0 = 1.6;
i1=0;
errors= [];
if abs(double(f1(x0))) > 10^(-5)
    while abs(double(f1(x0)))> 10^-5
        error = abs(f1(x0));
        errors = [errors, error]; % I wrote this because we need errors to calculate rate of convergence

        x0 = x0- double(f1(x0))/double(df(x0));
        i1= i1+1;
    end
    fprintf('%f is the root obtained by Newtons Method',x0)
    rootsf(end+1)= x0;
    q = double(log(errors(end) / errors(end-1)) / log(errors(end-1) / errors(end-2)));
else
    fprintf('%f is the root',x0)
    
end
tEnd1 = cputime -tstart1

%Secant Method

tstart2 = cputime

f =@(x) x^2-4*exp(-x) -3
x0 = 1.6;
x1 = 1.7;
i2=0;
errors2= [];
if abs(f(x1)) > 10^(-5)
    while abs(f(x1)) > 10^-5
       x1 = x1- f(x1)/((f(x1)- f(x0))/(x1-x0));
       i2=i2+1;
       error = abs(f(x1));
       errors2(end +1) = error;
    end
   
    fprintf('%f is the root obtained by Secant Method\n',x1)
    rootsf(end+1)=x1;
    q2 = log(errors2(end) / errors2(end-1)) / log(errors2(end-1) / errors2(end-2));
    
else
    fprintf('%f is the root',x1)

end
tEnd2 = cputime -tstart2

% Regula Falsı

tstart3 = cputime
f =@(x) x^2 -4*exp(-x)-3;
a = 1;
b= 4;
wfunc = @(z,y) (y*f(z)-z*f(y))/(f(z)-f(y))
w= wfunc(a,b); %weighted point
i3=0;
errors3 = [];
if f(a)*f(b) <0
    while abs(f(w)) > 10^-5
        if f(a)*f(w)<0
            b=w;
        else
            a=w;
        end
        w= wfunc(a,b);
        f(w);
        i3=i3+1;
        error = abs(f(w));
        errors3(end +1) =  error;
    end
    fprintf('%f is the root find by Regula Falsı Method\n',w)
    q3 = log(errors3(end) / errors3(end-1)) / log(errors3(end-1) / errors3(end-2));
    rootsf(end+1)=w;
else
    fprintf('There is no root')
end
tEnd3 = cputime -tstart3
q3

% Fix Point Method

tstart4=cputime
syms x
f=@(x)  x.^2-4.*exp(-x) -3;
g(x) = sqrt(4*exp(-x)+3);
x1=1.6;
dfG(x)= diff(g,x);
i4=0;
errors4 = [];
if double(dfG(x1)) <1 % we checked wether it is convergent or not
    while abs(f(x1)) >10^-5
       x1=double(g(x1));
       error = abs(f(x1));
       errors4(end +1) =  error;
       i4=i4+1;
    end
    q4 = log(errors4(end) / errors4(end-1)) / log(errors4(end-1) / errors4(end-2));
    fprintf('%f is the root obtained by fixed poit method.\n',x1)
    rootsf(end+1)=x1;
else
    fprintf('This funciton diverges');
end
tEnd4= cputime - tstart4;

% Muller's Method

tstart5 = cputime;
f= @(x) x.^2-4.*exp(-x) -3;
xi= [1.6 1.7 1.8];
z= polyfit(xi,f(xi),2);
a=[];
l1=@(x) (-x(2)+sqrt(x(2)^2-4*x(1)*x(3)))/(2*x(1));
l2=@(x) (-x(2)-sqrt(x(2)^2-4*x(1)*x(3)))/(2*x(1)); % Found the possible roots
a(1)= l1(z);
a(2)= l2(z);
i5=0;

if abs(f(a(1))) < abs(f(a(2)))% choose the one that is close to 0
            root = a(1);
else
            root= a(2);
end

if ~all(f(xi))
    fprintf('%f is the root',f(xi).find(f(xi)==0));
else
    
    while abs(f(root))> 10^-5
        z=polyfit(xi,f(xi),2);
        a(1)= l1(z);
        a(2)= l2(z);
        if abs(f(a(1))) < abs(f(a(2)))
            root = a(1)
        else
            root= a(2)
        end
        xi(1)=xi(2); 
        xi(2)=xi(3); 
        xi(3)=root;
        i5=i5+1;
        
    end
   
    fprintf('%f is the root obtained by Mullers method\n',root)
    rootsf(end+1)= root;
end
tEnd5= cputime - tstart5;


% Table

Methods =["Newton Raphson";"Secant";"Regula Falsı";"Fix Point";"Müller"];

cputimes = [tEnd1;tEnd2;tEnd3;tEnd4;tEnd5];

iteration_num=[i1;i2;i3;i4;i5];
ConvergenceRates= [2;1.618;1;1;1.84];
RateOfConvergenceCalculated= [q;q2;q3;q4;NaN];
RootsCalculated = rootsf';
table1 = table(Methods,cputimes,iteration_num,ConvergenceRates,RateOfConvergenceCalculated,RootsCalculated)

% https://en.wikipedia.org/wiki/Rate_of_convergence (Reference)