double fExact(double x, double *xP, double *yP, int Npoint)
{
  double f, lagrange = 0.0;
  f = 0;
  for(int i = 0; i <Npoint ; i++)
    {
      lagrange = 1.0;	  
      for(int j = 0; j< Npoint ;  j++)
	{
	  if(j != i) lagrange *= ( x - xP[j] )/( xP[i] - xP[j] );
	}
      f += yP[i]*lagrange;
    }  
  return f;
}
