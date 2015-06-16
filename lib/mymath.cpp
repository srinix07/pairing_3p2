#include "mymath.h"

using namespace std;




double interpol(double *x,double *y,int N,double x_des)
{

	gsl_interp_accel *acc = gsl_interp_accel_alloc ();
    gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, N);

    gsl_spline_init (spline, x, y, N);

	double yi;

    if(x_des<x[0])
    	yi = y[0];
    else if(x_des>x[N-1])
    	yi = y[N-1];
    else
        yi = gsl_spline_eval (spline, x_des, acc);

    gsl_spline_free (spline);
    gsl_interp_accel_free (acc);

    return yi;
}

int readv(char *fi, double* &K, double** &V, int &ksize)
{


	fstream F;

	string line;
	F.open(fi,ios::in);
	int i;
	int size = 0;

	while(getline(F,line))
	{
		if(line.length()!=0 && line.at(0)!='#')
			size++;
	}

F.clear();
F.seekg(0, ios::beg);


double *X=new double [size];
double *Y=new double [size];
double *Z=new double [size];

	while(getline(F,line))
	{
		if(line.length()!=0 && line.at(0)!='#')
	  {
		string st=line;
		istringstream ss(st);
		int count=1;
		while (!ss.eof())
		    {
				  string t;
				  getline( ss, t , ' ' );
				 if(t.length()!=0)
				 {
				  if(count==1)
				  	X[i]=atof(t.c_str());
				  if(count==2)
				  	Y[i]=atof(t.c_str());
				  if(count==3)
				  	Z[i]=atof(t.c_str());
				  	count++;
				  }

			}

			i++;
	   }

	}

	F.close();
	ksize=(int)sqrt(size);



			K = new double [ksize];


			for(int l=0;l<ksize;l++)
				{
				K[l]=Y[l];
				}

				V=new double* [ksize];
				for(int h=0;h<ksize;h++)
				{
					V[h]=new double [ksize];
				}


			for(i=0;i<ksize;i++)
				for(int j=0;j<ksize;j++)
				{
					//V[i][j]=0.0378769364844850*Z[i*ksize+j];
					V[i][j]=Z[i*ksize+j];
				}
delete [] X;
delete [] Z;
delete [] Y;
return 0;
}



/*
int gauss_leg(double* x, double* w,int N,double a,double b)
{
	//x=new double [N];
	//w=new double [N];
if(N)
{
	char* file_x;


		if(N==20)
			{
				file_x="gauss/l20";

			}
		if(N==64)
			{
				file_x="gauss/l64";

			}
		if(N==120)
			{
				file_x="gauss/l120";

			}
		if(N==54)
			{
				file_x="gauss/l54";

			}
			if(N==240)
			{
				file_x="gauss/l240";

			}
			if(N==320)
			{
				file_x="gauss/l320";

			}
			if(N==400)
			{
				file_x="gauss/l400";

			}

		vector<double> D;

		readV(file_x,D);
		N =  D.size();
		for (int i = 0; i <N/2 ; i++)
			{
				x[i]=D[2*i];
				x[i]=((b+a)*0.5)+((b-a)*0.5)*x[i];

				w[i]=(b-a)*0.5*D[1 + 2*i];

			}

		//delete D;
			return 0;
}
else
		return 0;

}

	double int_gauss(double (*f)(double,double*), double* params,double *x,double *w, double a,double b,int N)
	{
		//double *x,*w;
		//x=params[3];
		//w=params[4];
		double I=0;
		//gauss_leg(x,w,N);
		for(int i=0;i<N;i++)
			{
				x[i]=((b+a)/2)+((b-a)/2)*x[i];
				I=I+w[i]*f(x[i],params);
			}

			I=I*((b-a)/2);


			delete [] x;
			delete [] w;

		return I;
	}
*/
	double interpol2(double* x, double *V, int N,double xd,double yd)
	{
		if(xd<=x[N-1] && yd<=x[N-1])
		{
			double *Y=new double [N];
			int i=0;
		#pragma omp parallel shared(Y,V,x) private(i)
			{

   			#pragma omp for schedule(dynamic)
			for(i=0;i<N;i++)
				{
					Y[i]=interpol(x,&V[N*i],N,yd);
				}
			}
				double r= interpol(x,Y,N,xd);
				delete [] Y;

			return r;
		}
		else
			return 0;
	}



int readV(string fi, vector<double> &D)
{

fstream F;
	string l,w;

	//vector<double> v;

	F.open(fi,ios::in);

	while(getline(F,l))
	{
		if(l.length()!=0 && l.at(0)!='#')
		{
			istringstream ss(l);

			while(getline(ss,w,' '))
			{
				if(w.length()!=0 && w.at(0)!='\t')
					{

						D.push_back(atof(w.c_str()));
					}
			}


		}


	}

	F.close();

return 0;
}

int getV(vector<double> &D,vector<double> &K,vector<double> &V,int col)
{
	int i=0;
	int size = D.size();

	do
	{
			K.push_back(D[1+col*i]);
			//std::cout<<D[1+i+5*i]<<endl;

		i++;
	}while(D[1+col*i]!=D[1] || i==0);



	//vector<double> VV;

	for (int j = 0; j < col-2; j++)
		{
			i=0;
			while((col-1+col*i)<size)
				{

					V.push_back(D[2 + j + col*i]);

					i++;
				}
		}


return 0;
}

/*
int grid_split(double kf,double kstart,double kend,double offset,double offset2, int N1,int N2,int N3, double* &K, double* &W,int &N)
{

	if(kf - offset < 0)
		N1=0;
	if(kf + offset > kend)
		N3=0;

	N = N1 + N2 + N3;
//cerr<<N;
	K = new double[N];
	W = new double[N];

	gauss_leg(K,W,N1,kstart,kf-offset);

	gauss_leg(&K[N1],&W[N1],N2,kf-offset,kf+offset2);

	gauss_leg(&K[N1+N2],&W[N1+N2],N3,kf+offset2,kend);

//cerr<<N;



	return 0;
}
*/

std::string shell_execute(string cmd)
{
    FILE* pipe = popen(cmd.c_str(), "r");
    if (!pipe) return "ERROR";
    char buffer[128];
    std::string result = "";
    while(!feof(pipe)) {
    	if(fgets(buffer, 128, pipe) != NULL)
    		result += buffer;
    }
    pclose(pipe);
    return result;
}
