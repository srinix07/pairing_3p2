#include "../lib/mymath.h"
#include "../lib/gauss.h"
#include "../lib/const.h"

using namespace std;
using namespace arma;
using namespace TCLAP;


vec getf(vec &DD, double* &K, int ksize,double *VKF, double *we, double kf, vec &SHP11, vec &SHP12, vec &SHP21, vec &SHP22, bool relativistic)
{
  double gap1 =0;
  double gap2 =0;
  double E    =0,ek=0,ef=0;
  double dkf  =0;
  
  mat vkf(2,2);
  vec F(4);
  
  F[0]=0;F[1]=0;F[2]=0;F[3]=0;
  
  for (int i = 0; i < ksize; ++i)
    {
      vkf(0,0) = VKF[ksize*0 + i];
      vkf(1,0) = VKF[ksize*1 + i];
      vkf(0,1) = VKF[ksize*2 + i];
      vkf(1,1) = VKF[ksize*3 + i];
      
      gap1 = DD[0]*SHP11(2*i) + DD[1]*SHP12(2*i) + DD[2]*SHP21(2*i) + DD[3]*SHP22(2*i);
      gap2 = DD[0]*SHP11(2*i+1) + DD[1]*SHP12(2*i+1) + DD[2]*SHP21(2*i+1) + DD[3]*SHP22(2*i+1);
      
      dkf = gap1*gap1 + gap2*gap2;
      
      if(relativistic)
	{
	  ek = sqrt(Mn*K[i]*K[i] + Mn*Mn);
	  ef = sqrt(Mn*kf*kf + Mn*Mn);
	}
      else
	{
	  ek = K[i]*K[i]*0.5;
	  ef = kf*kf*0.5;
	}
      
      E      = -1*we[i]*K[i]*K[i]*(1/PI)*(1/(sqrt((ek - ef)*(ek - ef) + dkf)));
      
      F[0] += E*vkf(0,0)*gap1;
      F[1] += E*vkf(0,1)*gap2;
      F[2] += E*vkf(1,0)*gap1;
      F[3] += E*vkf(1,1)*gap2;
      
    }
  
  return F;
}

int main(int argc, char** argv)
{
  string potential;
  bool relativistic;
  bool np;
  bool mev;
  double del;
    
  //Command line Parameter parsing.
  try
    {
    
      CmdLine cmd("NN-Pairing By Khodel Method.", ' ', "jun-15");
      
      ValueArg<string> Potential("p","potential","Provide the NN Potential",true,"potentia","string");
      cmd.add(Potential);
      
      ValueArg<double> Focus("d","delta","Sets the Grid width around kf.(default:0.005)",false,0.005,"double");
      cmd.add(Focus);
      
      SwitchArg Relativistic("r","relativistic","Toggle Relativistic Calculation.(default: false)", cmd, false);
      SwitchArg Np("P","npsinglet","Toggle NP Calculation.(default: NN)", cmd, false);
      SwitchArg MeV("M","MeV","Print Gap in Mev.(default: false)", cmd, false);
      
      cmd.parse( argc, argv );
      
      potential    = Potential.getValue();
      relativistic = Relativistic.getValue();
      np           = Np.getValue();
      mev          = MeV.getValue();
      del          = Focus.getValue();
      
    }
  catch(ArgException &e)
    {
      std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    }
  
  
  cout.precision(15);
  cout.setf(ios::scientific | ios::showpos);
  
  
  std::vector<double> D;
  std::vector<double> K;
  std::vector<double> V;
  
  readV(potential,D);
  getV(D,K,V,6);
  
  double *VV = &V[0];
  
  int ksize = K.size();
  cerr << ksize<<endl;
  
  double *VKF;
  double su    = (np)?NP:NN;
  double units = (mev)?((np)?HFM_NP:HFM_NN):1;
  double E     = 0,ek=0,ef=0;
  
  double *x,*we,*VG;
  double *KK = &K[0];
    
  int ii =0;
  int Ng;
  
  gauss g64("64"),g20("30");     // #gauss-points Change grid size here. allowed sizes 20,30,54,64,120,240,320,400,2000.
  
  
  
  for (double kf = 0.1; kf <3.6; kf+=0.1)
    {
      
      gauss gt;
      
      g20.rescale(0,kf-del);       //grid splitting 0---kf-del---kf---kf+del----kmax.
      gt.app(g20);
      g64.rescale(kf-del,kf);
      gt.app(g64);
      g64.rescale(kf,kf+del);
      gt.app(g64);
      g20.rescale(kf+del,K[ksize-1]);
      gt.app(g20);
      Ng = gt.size;
      cerr<<"Grid Size:"<<Ng<<endl;
      x = &(gt.xr[0]);
      we = &(gt.wr[0]);
      VG=new double [Ng*Ng*4];

      //Interpolation
#pragma omp parallel shared(VG,VV,KK,x) private(ii)  //OpenMp #parallel 
      {
#pragma omp for schedule(dynamic)                    //Openmp #parallel
	for (ii = 0; ii < Ng*Ng; ii++)
	  {
	    VG[ii+Ng*Ng*0] = su*interpol2(KK,&VV[ksize*ksize*0],ksize,x[ii/Ng],x[ii%Ng]);
	    VG[ii+Ng*Ng*1] = -su*interpol2(KK,&VV[ksize*ksize*1],ksize,x[ii/Ng],x[ii%Ng]);
	    VG[ii+Ng*Ng*2] = -su*interpol2(KK,&VV[ksize*ksize*2],ksize,x[ii/Ng],x[ii%Ng]);
	    VG[ii+Ng*Ng*3] = su*interpol2(KK,&VV[ksize*ksize*3],ksize,x[ii/Ng],x[ii%Ng]);
	    
	  }
      }
      
      VKF = new double[4*Ng];
      
      for (int i = 0; i < Ng; i++)
	{
	  
	  VKF[Ng*0 + i] = 	interpol(x,&VG[Ng*Ng*0 + Ng*i],Ng,kf);
	  VKF[Ng*1 + i] = 	interpol(x,&VG[Ng*Ng*1 + Ng*i],Ng,kf);
	  VKF[Ng*2 + i] = 	interpol(x,&VG[Ng*Ng*2 + Ng*i],Ng,kf);
	  VKF[Ng*3 + i] = 	interpol(x,&VG[Ng*Ng*3 + Ng*i],Ng,kf);
	  
	}
      
      mat vff(2,2);
      
      vff(0,0)=interpol(x,&VKF[Ng*0],Ng,kf);  //V(kf,kf)
      vff(0,1)=interpol(x,&VKF[Ng*1],Ng,kf);
      vff(1,0)=interpol(x,&VKF[Ng*2],Ng,kf);
      vff(1,1)=interpol(x,&VKF[Ng*3],Ng,kf);
      
      double *W;
      W = new double[4*Ng*Ng];
      
      double dkf  = 1e-8;                     //initial guess
      double dkf2 = 1;
      
      while(fabs(dkf-dkf2)>1e-12)             // iteration loop
	{
	  std::vector<double> B11,B12,B21,B22;
	  
	  for (int l = 0; l <Ng*Ng; l++)
	    {
	      //matrix loop
	      int i=l/Ng;
	      int j=l%Ng;
	      
	      mat vv(2,2);
	      mat vkf(2,2);
	      mat vfk(2,2);
	      
	      vv(0,0) = VG[j + i*Ng + 0*Ng*Ng];  //V(i,j)
	      vv(0,1) = VG[j + i*Ng + 1*Ng*Ng];
	      vv(1,0) = VG[j + i*Ng + 2*Ng*Ng];
	      vv(1,1) = VG[j + i*Ng + 3*Ng*Ng];
	      
	      vkf(0,0) = VKF[Ng*0 + i];          //V(i,kf)
	      vkf(0,1) = VKF[Ng*1 + i];
	      vkf(1,0) = VKF[Ng*2 + i];
	      vkf(1,1) = VKF[Ng*3 + i];
	      
	      vfk(0,0) = VKF[Ng*0 + j];          //V(kf,j);
	      vfk(1,0) = VKF[Ng*1 + j];
	      vfk(0,1) = VKF[Ng*2 + j];
	      vfk(1,1) = VKF[Ng*3 + j];
	      
	      mat B_temp;
	      B_temp= vkf/vff;
	      
	      if(j==0)
		{
		  B11.push_back(B_temp(0,0));
		  B11.push_back(0);
		  B12.push_back(B_temp(0,1));
		  B12.push_back(0);
		  B21.push_back(0);
		  B21.push_back(B_temp(1,0));
		  B22.push_back(0);
		  B22.push_back(B_temp(1,1));
		}
	      
	      mat w_temp;
	      
	      ek     = (relativistic)?sqrt(Mn*x[j]*x[j] + Mn*Mn):(x[j]*x[j]*0.5);   //starting energies
	      ef     = (relativistic)?sqrt(Mn*kf*kf + Mn*Mn)    :(kf*kf*0.5);
	      
	      E      = we[j]*x[j]*x[j]*(1/PI)*(1/(sqrt((ek - ef)*(ek - ef) + dkf)));
	      
	      w_temp = E*(vv - (vkf/vff)%vfk);      //Khodel Equation Note: '/' does element wise division, '%' does element wise multiplication. 
	      
	      W[2*i+2*Ng*(2*j)]       = w_temp(0,0);  //Making the Fredholm matrix.
	      W[1 + 2*i+2*Ng*(2*j)]   = w_temp(1,0);
	      W[2*i+2*Ng*(2*j+1)]     = w_temp(0,1);
	      W[1 + 2*i+2*Ng*(2*j+1)] = w_temp(1,1);
	      
	    }  //matrix loop end
	  
	  mat WW(W,2*Ng,2*Ng,false);
	  mat	WWt;
	  WWt = WW + eye<mat>(2*Ng,2*Ng);  //Adding I matrix.
	  
	  vec b11(B11);
	  vec b12(B12);
	  vec b21(B21);
	  vec b22(B22);
	  
	  vec SHP11 = solve(WWt,b11);      //Solving for 4 shape factors.
	  vec SHP12 = solve(WWt,b12);
	  vec SHP21 = solve(WWt,b21);
	  vec SHP22 = solve(WWt,b22);
	  
	  vec DD(4);
	  DD.fill(1e-3);
	  vec F(4);
	  F.fill(10);
	  vec ff(2);
	  
	  int Ndd=0;
	  double norm = 1;
	  while(sqrt(norm)  > 1e-12)       //iteraratively solving for  D Co-efficients.
	    {
	      
	      F    = getf(DD,x,Ng,VKF,we,kf,SHP11,SHP12,SHP21,SHP22,relativistic);
	      
	      norm = (DD(0)-F(0))*(DD(0)-F(0)) + (DD(1)-F(1))*(DD(1)-F(1)) + (DD(2)-F(2))*(DD(2)-F(2)) + (DD(3)-F(3))*(DD(3)-F(3));
	      
	      DD   = F;
	      
	      Ndd++;
	      
	    }
	  
	  double gapf1=DD[0]+DD[1];
	  double gapf2=DD[2]+DD[3];
	  
	  dkf2 = dkf;
	  dkf  = gapf1*gapf1 + gapf2*gapf2;
	  
	}//gap iteration
      
      cout<<kf<< "  " << sqrt(dkf)*units<< endl;
      cerr<< del << "   "<< kf <<"  "<< sqrt(dkf)*units << endl;
      
      delete [] W;
      delete [] VKF;
      
    } //kf loop
  delete [] VG;
  
  return 0;
  
}
