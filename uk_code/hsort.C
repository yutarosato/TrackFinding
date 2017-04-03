void hsort(int n,int ra[]){
  int l,j,ir,i;
  int rra;

  l=(n >> 1)+1;
  ir=n-1;
  for(;;){
    if(l>0){
      rra=ra[--l];
    }
    else{
      rra=ra[ir];
      ra[ir]=ra[0];
      if(--ir == 0){
	ra[0]=rra;
	return;
      }
    }
    i=l;
    j=(l<< 1);

    while(j <=ir){
      if(j < ir && ra[j] < ra[j+1]) ++j;
      if(rra < ra[j]){
	ra[i]=ra[j];
	if(j==0){
	i=j;
	j++;
	}
	else{
	j += (i=j);
	}	
      }
      else{
	j=ir+1;
      }
    }
    ra[i]=rra;
  }
}
