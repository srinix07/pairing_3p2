#include "gauss.h"


gauss::gauss():size(0),x(),w(),xr(),wr()
{
	

}

gauss::gauss(std::string gsize):size(0),x(),w(),xr(),wr()
{
	std::string file_x = "../lib/gauss/gauss_";
	file_x += gsize;
	file_x += ".dat";
	
				
		std::vector<double> D;

		readV(file_x,D);
		
		for (unsigned int i = 0; i <D.size()*0.5 ; i++)
			{
				x.push_back(D[2*i]);
				//x[i]=((b+a)*0.5)+((b-a)*0.5)*x[i];
		
				w.push_back(D[1 + 2*i]);
			
			}
		
		
		size = D.size()*0.5;	
	
		
}



void gauss::rescale(double aa, double bb)
{
	xr.clear();
	wr.clear();
	for (int i = 0; i < size; i++)
	{
		xr.push_back(((bb+aa)*0.5)+((bb-aa)*0.5)*x[i]);
		wr.push_back((bb-aa)*0.5*w[i]);
	}
	a=aa;
	b=bb;
}

std::vector<double> gauss::get(int j)
{
	if(j)
		return wr;
	else
		return xr;

}

void gauss::app(gauss g)
{
	
	for (int i = 0; i < g.size ; ++i)
	{
		xr.push_back(g.xr[i]);
		wr.push_back(g.wr[i]);
	}

	size += g.size;
	b = g.b;
}

void gauss::check()
{
	std::cerr.precision(15);
	std::cerr.setf(std::ios::scientific | std::ios::showpos);
	
	double I0=0,I1=0,I2=0;
	for (int i = 0; i < size; ++i)
	{
		I0 += wr[i];
		I1 += wr[i]*xr[i];
		I2 += wr[i]*xr[i]*xr[i];
	}
	std::cerr<<"size"<<"  "<<"a"<<"  "<<"b"<<"  "<<"I0"<<"  "<<"b-a"<<"  "<<"I1"<<"  "<<"(b*b-a*a)*0.5"<<"  "<<"I2"<<"  "<<"(b*b*b-a*a*a)/3"<<std::endl;
	std::cerr<<size<<"  "<<a<<"  "<<b<<"  "<<I0<<"  "<<b-a<<"  "<<I1<<"  "<<(b*b-a*a)*0.5<<"  "<<I2<<"  "<<(b*b*b-a*a*a)/3<<std::endl;
}

void gauss::print()
{	
	std::cerr.precision(15);
	std::cerr.setf(std::ios::scientific | std::ios::showpos);
	
	for (int i = 0; i < size; ++i)
	{
		std::cerr<<xr[i]<<"  "<<wr[i]<<std::endl;	
	}
	
}

void gauss::stitch(gauss g1,gauss g2)
{
	xr.clear();
	wr.clear();
	for (int i = 0; i < g1.size ; ++i)
	{
		xr.push_back(g1.xr[i]);
		wr.push_back(g1.wr[i]);
	}
	for (int i = 0; i < g2.size; ++i)
	{
		xr.push_back(g2.xr[i]);
		wr.push_back(g2.wr[i]);
	}
	
	size = g1.size + g2.size;
	a = g1.a;
	b = g2.b;
}

gauss::~gauss()
{

}
