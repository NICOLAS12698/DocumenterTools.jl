using Plots
###### primal simplex method
println("#######################################################")
println("Prof. Vanderbei's experiment")
println("#######################################################")

### Tolerance
eps = 1e-7;
### Parameter to create data
sigma = 10;
### Number of iterations (= number of LPs to be solved)
niters = 100
### Variables initialization
global k=0; # Main iteration counter

ms = zeros(niters,1); # array of m's
ns = zeros(niters,1); # array of n's
iters = zeros(niters,1); # array of number of Simplex method iterations
opts = zeros(niters,1); # array of optimization status

##################################
######### Beginning of main loop
#################################
while k != niters
   m = Int(round(10*exp(log(100)*rand())));
   n = Int(round(10*exp(log(100)*rand())));
   global A = Int.(round.(sigma*(randn(m,n))));
   global b = Int.(round.(sigma*abs.(randn(m))));
   global c = Int.(round.(sigma*randn(n)));

   A = -A;    # to reflect subtraction from right-hand side


   nonbasics = collect(1:n);
   basics = collect(n+1:n+m);

   global iter = 0; # Simplex method iteration counter

 # Optimization status: opt=1 (optimal),-1 (unbounded), 0 (too hard)
   global opt = 1;
   ##################################
   ######### Beginning of Simplex method loop
   #################################
   while  maximum(c) > eps
    ### Largest positive reduced costs
      col = argmax(c); # position of the entering nonbasic variable
      Acol = A[:,col];  # column associated to the entering nonbasic variable

    ### Condition for unboundeness of the objective function
    ### (The ratio test cannot be performed!)
     if minimum(Acol)>=-eps
      #if sum(1 for i in 1:length(Acol) if Acol[i]<-eps) == 0
	    opt = -1;  ### unbounded
	    break;
      end

    ### Computing numerators and denominators to perform the ratio test
    ### (In order to obtain an mx1 vector I fill the components with Acol[i]>=-eps with a large number)
      nums = [Acol[i]<-eps ? b[i] : 2*sum(b[i]/(-Acol[i]) for i in 1:m if Acol[i]<-eps) for i in 1:m];
      dens = [Acol[i]<-eps ? -Acol[i] : 1 for i in 1:m];

     ### Performing the ratio test

      row = argmin(nums./dens); # position of the leaving basic variable

      j = nonbasics[col]; # index of entering nonbasic variable
      i = basics[row]; # index of leaving basic variable

      ### Updating matrix A
      Arow = A[row,:]; # row of the leaving basic variable
      a = A[row,col]; # pivot
      global A = A - Acol*(Arow/a)';
      A[row,:] = -Arow/a;
      A[:,col] = Acol/a;
      A[row,col] = 1/a;
      ### Updating r.h.s. b
      brow = b[row];
      global b = b - brow*Acol/a;
      b[row] = -brow/a;
      ### Updating objective function c
      ccol = c[col];
      global c = c - ccol*Arow/a;
      c[col] = ccol/a;
      ### Updating basic and nonbasic variables
      basics[row] = j;
      nonbasics[col] = i;

      ### Going to next Simplex method iteration
      global iter = iter+1;

      ### Checking if we go above the iteration limit (problem is too hard)
      if iter > 100*(m+n)
      	  opt = 0;  # iter limit
          break;
      end
   end
   ##################################
   ######### End of Simplex method loop
   #################################

   ### Saving the data if optimization status is optimal (opt=1) or unbounded (opt=-1)
   if iter>0 && opt!=0
      println("k = ", k)
      global k = k+1; # going to next main iteration
      ms[k] = m; # saving m
      ns[k] = n; # saving n
      iters[k] = iter; # saving the number of iterations
      opts[k] = opt; # saving the optimization status
  end

end
##################################
######### End of main loop
#################################

### Separating the data into problems with optimal value and unbounded problems
opt_ms=[ms[i] for i in 1:niters if opts[i]==1]
opt_ns=[ns[i] for i in 1:niters if opts[i]==1]
opt_iters=[iters[i] for i in 1:niters if opts[i]==1]

unb_ms=[ms[i] for i in 1:niters if opts[i]==-1]
unb_ns=[ns[i] for i in 1:niters if opts[i]==-1]
unb_iters=[iters[i] for i in 1:niters if opts[i]==-1]

### Plots environment
gr()
### Plotting
display(scatter(min.(opt_ms,opt_ns), opt_iters, label="Problems with an opt. sol",legend=:bottomright))
display(scatter!(min.(unb_ms,unb_ns), unb_iters,label="Problems that are unb.",legend=:bottomright))
