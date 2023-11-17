#include <stdio.h>
#include <math.h>
#include <time.h>
//#include <curses.h>

#define SG 1555168.
#define PI 3.14159265

float n1, n2, de, d, r1, r2, d1, d2, z1, z2, c1, c2, zz[3], v01, v02, gele, surf;
float hi(void), n11(float z), n22(float z), n3(float z), dn1(float z), dn2(float z),
      dn3(float z), fi1(float z), fi2(float z), fi3(float z), w(float nx, float fi, float dn, float dnn),
      w1(float x), w2(float x), w3(float x), sei(float dr,float rc), ei(void),
      ii(void), me(void), wku(float nx,float fi), wki(float nx), wob(float nx),
      wcor(float nx), wvk(float nx,float dn), wku1(float x), wku2(float x),
      wku3(float x), wki1(float x), wki2(float x), wki3(float x), wob1(float x),
      wob2(float x), wob3(float x), wcor1(float x), wcor2(float x), wcor3(float x),
      wvk1(float x), wvk2(float x), wvk3(float x), sum1(float fn(float x)), sum2(float fn(float x)),
      sum3(float fn(float x)), t(float x0,float x1, float f(float z)), ww(float nx,float dn),
      ww1(float nx), ww2(float nx), ww3(float nx), eid(void), dn21(float z), dn22(float z),
      dn23(float z), svk41(float x), svk42(float x), svk43(float x), sww41(float x),
      sww42(float x), sww43(float x), svk4(float nx, float dn, float dn22), sww4(float nx,float dn22),
      a_exit(void);

int m;
FILE *vh, *vv, *fopen( );

int main(void)
{
    extern float n1, n2, r1, r2, d, d1, d2, z1, z2, c1, c2, de, zz[3], v01, v02;
    extern int m;
    int i, j, mm, ps, bs, N, jj, i1;
    float step, y[3], b[3], p[3], d_max;
    float s7[42], be2[42], be3[42], s9[42], s10[42], s11[42];
    float ep, d0, h, sku1, sku2, sku3, ski1, ski2, ski3, beta, f1, f2, f3, v1, g, s8[42], sww1, sww2, sww3;
    float r[42], be[42], s[42], s1[42], s2[42], s3[42], s4[42], s5[42], s6[42], ea[42];
    float sob1, sob2, sob3, scor1, scor2, scor3, svk1, svk2, svk3, se1, se2;
    float ssvk41, ssvk42, ssvk43, ssww41, ssww42, ssww43;
    float av;

    float rs, t, t1, t2, t3, t4, t5;

    int k, l, M = 0, obrazec = 0, stop = 0, net = 0;
    float hbd[] = {1.,1.}, betadelta[] = {1.4,0}, beta1, Pbd[2], zz_temp[2], bd1, bd2, f0;

    //////////////////////////////////////////////////////////
    time_t start;  time_t finish;  float et;  unsigned let, min, sec;
    time( &start );
    //////////////////////////////////////////////////////////

    N = 1;
    n1 = 0.0038;
    n2 = 0;
    d1 = 5.71;
    d2 = 3.85;
    r1 = 1.74;
    r2 = 0.92;
    c1 = 6.99;
    c2 = 4.70;
    z1 = 1;
    z2 = 0;
    ep = 1;
    d0 = 0.0;
    h = 0.2;
    se1 = 1000;
    se2 = 1000;
    mm = 0;
    m = 16;
    gele = 0.0;

    printf(" Number variational parameters? \n");
    scanf("%d",&N);
    printf(" Bulk electron density 1 material n1 \n");
    scanf("%f",&n1);
    printf(" Uniform ion density? 0-Yes  1-No \n");
    scanf("%d",&gele);
    if (gele > 0.0) { gele=1.0;
        printf(" interplaner spacing 1 material d1 \n");
        scanf("%f",&d1);
        printf("plane nearest neighbors distance  1 material c1 \n");
        scanf("%f",&c1);
        printf("ion charge  1 material Z1 \n");
        scanf("%f",&z1);
        printf(" cutoff radius of the Ashkroft pseudopotential 1 material rc1 \n");
        scanf("%f",&r1);
        if (r1==0) {
            rs=pow((double)(3./(4.*PI*n1)),(double)(1./3.));
            t=rs*rs;
            t=t/(3*r1*r1*r1);
            t1=2.21/rs;
            t2=0.9*pow((double)(z1),(double)(2./3.));
            t3=(4.5*r1*r1)/(rs*rs);
            t4=0.458;
            t5=0.00711*rs*rs;
            t5=t5/(1+0.127*rs);
            t5=t5/(1+0.127*rs);
            r1=pow((double)(t4-t1+t2+t5),(double)(0.5))*(sqrt(2.0)/3.0)*rs; //this is rc1
        }
        printf(" Surface or Adgesion? 0-surface  1-Adgesion \n");
        scanf("%d",&surf);
        if (surf>0) { mm=100;
        printf(" Valume electron density 2 material n2 \n");
        scanf("%f",&n2);
        printf(" interplaner spacing 2 material d2 \n");
        scanf("%f",&d2);
        printf("plane nearest neighbors distance  2 material c2 \n");
        scanf("%f",&c2);
        printf(" cutoff radius of the Ashkroft pseudopotential 2 material rc2 \n");
        scanf("%f",&r2);
        printf("ion charge  2 material Z2 \n");
        scanf("%f",&z2);
        printf("dielectric const epsilon? \n");
        scanf("%f",&ep);
        printf("Distance step h? \n");
        scanf("%f",&h);
        printf("Max distance D? \n");
        scanf("%f",&d_max);
        printf("surface Energy 1 material se1? \n");
        scanf("%f",&se1);
        printf("surface Energy 2 material se2? \n");
        scanf("%f",&se2);
        }
    }

    printf("n1=%f,d1=%f,r1=%f,z1=%f,mm=%f\n",n1,d1,r1,z1,mm);

    vv = fopen("result.txt","w");

    fprintf(vv, "n1=%f,d1=%f,r1=%f,z1=%f\n", n1, d1, r1, z1);
    fprintf(vv, "epsilon=%f\n", ep);
    fprintf(vv, "n2=%f,d2=%f,r2=%f,z2=%f\n", n2, d2, r2, z2);
    fprintf(vv, "c1=%f,c2=%f,se1=%f,se2=%f\n", c1, c2, se1, se2);
    fprintf(vv, "d0=%f,h=%f\n", d0, h);
    //for(jj=0;jj<=0;jj++){
    de = sqrt((double)ep);
    d = d0;
    for(i1 = 0; i1 <= mm; i1++){
      do { hbd[0]=1.; hbd[1]=1.;
                    do { net=0;
    for (j = 0; j < N; j++) {
     k = 0; l = 0;

     zz[j]=betadelta[j];
     f1=me();
     zz[j]=betadelta[j]+hbd[j]; f2=me();
     zz[j]=betadelta[j]-hbd[j]; f3=me();
     zz[j]=betadelta[j];

     if (f2<f1) { bd1=betadelta[j]+hbd[j]; k=1;}
     if (f3<f1) { bd2=betadelta[j]-hbd[j]; l=1;}

     if (k==1) { if (l==1) {if (f3<f2) zz[j]=bd2;
                            else zz[j]=bd1;
                           }
                 else zz[j]=bd1;
               }
     else { if (l==1) zz[j]=bd2;}
     if ( (k==0) && (l==0) ) net++;
                      }

    if ( net==N ) {
      if (M==0) {hbd[0]/=10; M=1;}
      else      {hbd[1]/=10; M=0;}
                  }

    if ( (hbd[0]<1.e-4) && (hbd[1]<1.e-4) ) { if (obrazec==0) stop=1;break;}

    			      } while(net==N);
                           M=0; net=0;
    	if (stop==1) break;

    if (obrazec==0) {
    f0=me(); obrazec=1;
     for (i=0;i<N;i++) {
      Pbd[i]=betadelta[i]+2*(zz[i]-betadelta[i]);
      betadelta[i]=Pbd[i];
      zz_temp[i]=zz[i];
    	               }
    				}
     else { if (me()<f0)  {

    f0=me();
    for (i=0;i<N;i++) {
      Pbd[i]=zz_temp[i]+2*(zz[i]-zz_temp[i]);
      betadelta[i]=Pbd[i];
      zz_temp[i]=zz[i];

    	              }   }
    		else { obrazec=0; betadelta[0]=zz_temp[0]; betadelta[1]=zz_temp[1];}
    	   }


    				   } while(stop==0);
    				   stop=0;


    /*
    step=1;
    for(j=1;j<=N-1;j++)y[j]=p[j]=b[j]=zz[j]=0.0;y[0]=p[0]=b[0]=zz[0]=1.4;
    f1=f2=me( );ps=0;bs=1;
    aaa:for(j=0;j<=N-1;j++){zz[j]=y[j]-step;f3=me( );
    			if(f3<f1){y[j]=zz[j];f1=f3;}
    			else{ zz[j]=y[j]+step;f3=me( );
    				   if(f3<f1){y[j]=zz[j];f1=f3;}
    				   else{zz[j]=y[j];f1=me( );}}}
    				   f1=me( );
    printf("f1=%f  f2=%f  %f  %f  %f\n",f1,f2,zz[0],zz[1],zz[2]);
    g=f2-0.00001;
    if(f1<g) goto ddd;
    if((ps==1)&&(bs==0)){for(j=0;j<=N-1;j++) zz[j]=y[j]=p[j]=b[j];
    		     f1=f2=me( );bs=1;ps=0;printf("vvv\n");goto aaa;}
    step/=10;
    printf("ШАГ УМЕНЬШЕН В 10 РАЗ\n");
    if(step>1.e-4)goto aaa;
    	      else goto min;
    ddd:for(j=0;j<=N-1;j++){p[j]=2*y[j]-b[j];b[j]=y[j];zz[j]=y[j]=p[j]; }
    f2=f1=me( );ps=1;bs=0;printf("sss\n");goto aaa;
    min:
    /**///putch(0x07);//putch(0x07);
    printf("          РЕЗУЛЬТАТ МИНИМИЗАЦИИ     \n ");

    printf("==============================================================\n");
    //fprintf(vv,"          РЕЗУЛЬТАТ МИНИМИЗАЦИИ     \n ");
    fprintf(vv,"==============================================================\n");

    for(j=0;j<=N-1;j++)printf("%f\n",zz[j]);

         v1 = f0;
         sku1 = sum1(wku1);
         sku2 = sum2(wku2);
         sku3 = sum3(wku3);
         ski1 = sum1(wki1);
         ski2 = sum2(wki2);
         ski3 = sum3(wki3);
         sob1 = sum1(wob1);
         sob2 = sum2(wob2);
         sob3 = sum3(wob3);
         scor1 = sum1(wcor1);
         scor2 = sum2(wcor2);
         scor3 = sum3(wcor3);
         svk1 = sum1(wvk1);
         svk2 = sum2(wvk2);
         svk3 = sum3(wvk3);
         sww1 = sum1(ww1);
         sww2 = sum2(ww2);
         sww3 = sum3(ww3);
         ssvk41 = sum1(svk41);
         ssvk42 = sum2(svk42);
         ssvk43 = sum3(svk43);
         ssww41 = sum1(sww41);
         ssww42 = sum2(sww42);
         ssww43 = sum3(sww43);
         s1[i1] = sku1+sku2+sku3;
         s2[i1] = ski1+ski2+ski3;
         s3[i1] = sob1+sob2+sob3;
         s4[i1] = scor1+scor2+scor3;
         s5[i1] = svk1+svk2+svk3;
         s6[i1] = ei( );
         s7[i1] = ii( );
         s8[i1] = sww1+sww2+sww3;
     //    s10[i]=ssvk41+ssvk42+ssvk43;
      //  s11[i]=ssww41+ssww42+ssww43;
         s9[i1] = eid( );
         ea[i1] = se1+se2-v1;
         s[i1] = s1[i1]+s2[i1]+s3[i1]+s4[i1]+s5[i1]+gele*s6[i1]+gele*s7[i1]+s8[i1]+s9[i1];//+s10[i]+s11[i];

         r[i1] = d;
         be[i1] = beta = zz[0];
         be2[i1] = zz[1];
         be3[i1] = zz[2];
         av = a_exit( );
         fprintf(vv,"d=%f,beta=%f,m.s.n.=%f,en.adg.=%f\n",d,beta,v1,se1+se2-v1);
    //     fprintf(vv,"sku=%f,ski=%f,sob=%f,scor=%f\n",s1[i],s2[i],s3[i],s4[i]);
    //     fprintf(vv,"svk=%f,ei=%f,ii=%f,seid=%f\n",s5[i],s6[i],s7[i],s9[i]);
         d=d+h;
         printf("d=%f i=%f\n",d,i1);
         if (d>d_max) break;}
    fprintf(vv," d !beta!mf.e.!en.ad!kulon!kinet!obmen!corel!neodn!elion!eid  !ii   !SV   ! fa  \n");
    for(i=0;i<=mm;i++)
    fprintf(vv,"%3.1f!%4.2f!%5.0f!%5.0f!%5.0f!%5.0f!%5.0f!%5.0f!%5.0f!%5.0f!%5.0f!%5.0f!%5.0f!\n",
    r[i],be[i],s[i],ea[i],s1[i],s2[i],s3[i],s4[i],s5[i],s6[i],s9[i],s7[i],s8[i]);
    //for(i=0;i<=mm;i++)fprintf(vv," %f  %f  %f %f %f\n",r[i],be2[i],be3[i],v01,r1);
    //for(i=0;i<=mm;i++)fprintf(vv," W4kin=%f  \n",s10[i]);
    fprintf(vv,"beta=%f     delta=%f    \n",zz[0],zz[1]);
    printf("beta=%f     delta=%f    \n",zz[0],zz[1]);
    //fprintf(vv,"Rm=%f     V0=%f    \n",r1,v01);
    //printf("Rm=%f     V0=%f    \n",r1,v01);
    fprintf(vv,"W(ЁрсюЄр т√їюфр)а =%f\n",av);
    printf("работа выхода =%f\n",av);
    //r1=r1+0.001;
    //}
    fclose(vv);
    ////////////////////////////////////////
    time( &finish ); et = difftime( finish, start );
    let=(unsigned)et; min=let/60; sec=let%60;
    printf("\n\n   *** Рассчет занял %4u мин. и %4u сек. ***\n",min,sec);
    //putch(0x07);// putch(0x07); putch(0x07); putch(0x07);
    ////////////////////////////////////////
}

float a_exit(void)
{//extern float zz[3];
float w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,w10_, y,b,alpha,beta,gama,dt;
float u,u0,u1,u2,u3,y_,ys;
b=zz[0];
dt=zz[1];
alpha=1.0;
beta=1.0;
gama=1.0;

       w1=4.0*PI*n1/(b*b);
       w2=4.0*PI*n1*d1*dt;
       w3=0.5*pow(3.0*PI*PI*n1, 2.0/3.0);
       w4=pow(3.0*n1/PI,1./3.);
       w5=0.056*pow(n1,2./3.)+0.0059*pow(n1,1./3.);
       w5=w5/((0.079+pow(n1,1./3.))*(0.079+pow(n1,1./3.)));
       w6=beta*0.4*pow(z1,2./3.);
       w6=w6*pow(4.*PI*n1/3.0,1./3.);
       w7=alpha*4*PI*n1*r1*r1;
       w8=gama*8.0*PI*v01*r1*r1*r1*n1/3.;
       w9=2.*PI*n1*(d1*d1*d1/12.+(2./3.)*v01*d1*r1*r1*r1-r1*r1*d1)/d1;
u=4.*PI*n1/(b*b);
u0=2-exp(-b*d1);
u1=-2*sinh(0.5*b*d1);
u2=v01*d1*(sinh(b*r1)-b*r1*cosh(b*r1));
u3=b*d1*cosh(b*r1);
w10=-u*(u1+u2+u3)*exp(-0.5*b*d1)/u0;
       y=w1/*-w2*/-w3+w4+w5+gele*(w6/*-w7+w8+w10*/);
       y=y*27.2;
u=4.*PI*n1*exp(b*dt)/(b*b);
u1=-2*sinh(0.5*b*d1)*exp(-b*dt)-b*b*dt*d1;
w10_=-u*(u1+u2+u3)*exp(-0.5*b*d1)/u0;
       y_=w1-w2-w3+w4+w5+w6-w7+w8+w10_;
       y_=y_*27.2;
       ys=w1-w2-w3+w4+w5+w6-w7+w8+w9;
       ys=ys*27.2;
/*fprintf(vv,"w1=%f\n",w1*27.2);
fprintf(vv,"w2=%f\n",w2*27.2);
fprintf(vv,"w3=%f\n",w3*27.2);
fprintf(vv,"w4=%f\n",w4*27.2);
fprintf(vv,"w5=%f\n",w5*27.2);
fprintf(vv,"w6=%f\n",w6*27.2);
fprintf(vv,"w7=%f\n",w7*27.2);
fprintf(vv,"w8=%f\n",w8*27.2);
fprintf(vv,"w9=%f\n",w9*27.2);
fprintf(vv,"w10=%f\n",w10*27.2);*/
//fprintf(vv,"w10_=%f\n",w10_*27.2);
//fprintf(vv,"раб.вых.с учетом релакс.=%f\n",y_);
//fprintf(vv,"раб.вых.ср.=%f\n",ys);
//printf("w1=%f\n",w1*27.2);
//printf("w2=%f\n",w2*27.2);
//printf("w3=%f\n",w3*27.2);
//printf("w4=%f\n",w4*27.2);
//printf("w5=%f\n",w5*27.2);
//printf("w6=%f\n",w6*27.2);
//printf("w7=%f\n",w7*27.2);
//printf("w8=%f\n",w8*27.2);
//printf("w9=%f\n",w9*27.2);
//printf("w10=%f\n",w10*27.2);
//printf("w10_=%f\n",w10_*27.2);
//printf("раб.вых.с учетом релакс.=%f\n",y_);
printf("раб.вых.ср.=%f\n",y);
return(y);
}

float hi(void)
{extern float d,de,zz[3];
float u,y,beta;
beta=zz[0];
u=(2.*beta*d)/de;
y=(1.+de*de)*sinh((double)u)+2.*de*cosh((double)u);
return(y);
}

float n11(float z)
{extern float d,n1,n2,zz[3],de;
float u,u1,u2,y,beta;
beta=zz[0];
u=(2.*beta*d)/de;
u1=hi( );
u2=(n1*exp((double)(-u))*(de-1.)+n2*(de+1.))/u1;
u2=n1-u2;
y=n1-(de/(de+1.))*z*u2;
return(y);
}

float n22(float z)
{extern float d,n1,n2,zz[3],de;
float u,u1,u2,y,beta;
beta=zz[0];
u=beta/de;
u1=n1*(de-1.)*exp((double)(u*(z-d)))+n2*(de+1.)*exp((double)(u*(z+d)));
u2=n1*(de+1.)*exp((double)(-u*(z-d)))+n2*(de-1.)*exp((double)(-u*(z+d)));
y=0.5*(u1+u2)/hi( );
return(y);
}

float n3(float z)
{extern float d,n1,n2,zz[3],de;
float u,u1,u2,y,beta;
beta=zz[0];
u=(2.*beta*d)/de;
u1=(de/(1.+de))*z;
u2=n2*(de-1.)*exp((double)(-u))+n1*(de+1.);
u2=u2/hi( );
y=n2-u1*(n2-u2);
return(y);
}

float dn1(float z)
{extern float d,n1,n2,de,zz[3];
float u1,u2,y,beta;
beta=zz[0];
u1=n1*(de-1.)*exp((double)(-2.*beta*d/de))+n2*(1.+de);
u1=n1-u1/hi( );
u2=-beta*(de/(de+1.))*z;
y=u2*u1;
return(y);
}

float dn2(float z)
{extern float d,n1,n2,de,zz[3];
float u,u1,u2,y,beta;
beta=zz[0];
u=beta/de;
u1=n1*(de-1.)*exp((double)(u*(z-d)))+n2*(de+1.)*exp((double)(u*(z+d)));
u2=n1*(de+1.)*exp((double)(-u*(z-d)))+n2*(de-1.)*exp((double)(-u*(z+d)));
y=0.5*u*(u1-u2)/hi( );
return(y);
}

float dn3(float z)
{extern float d,n1,n2,de,zz[3];
float u,u1,u2,y,beta;
beta=zz[0];
u=((beta*de)/(1.+de))*z;
u1=n2*(de-1.)*exp((double)(-2.*beta*d/de))+n1*(1.+de);
u2=n2-u1/hi( );
y=u*u2;
return(y);
}

float dn21(float z)
{extern float d,n1,n2,de,zz[3];
float u1,u2,y,beta;
beta=zz[0];
u1=n1*(de-1.)*exp((double)(-2.*beta*d/de))+n2*(1.+de);
u1=n1-u1/hi( );
u2=-beta*(de/(de+1.))*z;
y=beta*u2*u1;
return(y);
}

float dn22(float z)
{extern float d,n1,n2,de,zz[3];
float u,u1,u2,y,beta;
beta=zz[0];
u=beta/de;
u1=n1*(de-1.)*exp((double)(u*(z-d)))+n2*(de+1.)*exp((double)(u*(z+d)));
u2=n1*(de+1.)*exp((double)(-u*(z-d)))+n2*(de-1.)*exp((double)(-u*(z+d)));
y=0.5*u*u*(u1+u2)/hi( );
return(y);
}

float dn23(float z)
{extern float d,n1,n2,de,zz[3];
float u,u1,u2,y,beta;
beta=zz[0];
u=((beta*de)/(1.+de))*z;
u1=n2*(de-1.)*exp((double)(-2.*beta*d/de))+n1*(1.+de);
u2=n2-u1/hi( );
y=-beta*u*u2;
return(y);
}

float fi1(float z)
{extern float zz[3];
float y,beta;
beta=zz[0];
y=-(4*PI*n11(z))/(beta*beta);
return(y);
}

float fi2(float z)
{extern float zz[3];
float beta,y;
beta=zz[0];
y=-(4*PI*n22(z))/(beta*beta);
return(y);
}

float fi3(float z)
{extern float zz[3];
float beta,y;
beta=zz[0];
y=-(4*PI*n3(z))/(beta*beta);
return(y);
}

float w(float nx,float fi,float dn,float dnn)
{float y;
y=wku(nx,fi)+wki(nx)+wob(nx)+wcor(nx)+wvk(nx,dn)+ww(nx,dn) + svk4(nx,dn,dnn);
y+=sww4(nx,dnn);
return(y);
}

float w1(float x)
{extern float n1;
float nx,dn,fi,y,dnn;
nx=n11(x);
dn=dn1(x);
dnn=dn21(x);
fi=fi1(x);
y=(w(nx,fi,dn,dnn)-w(n1,fi,0.,0.))/x;
return(y);
}

float w2(float x)
{float fi,nx,dn,dnn,y;
nx=n22(x);
dn=dn2(x);
dnn=dn22(x);
fi=fi2(x);
y=w(nx,fi,dn,dnn);
return(y);
}

float w3(float x)
{extern float n2;
float fi,nx,dn,dnn,y;
nx=n3(x);
dn=dn3(x);
dnn=dn23(x);
fi=fi3(x);
y=(w(nx,fi,dn,dnn)-w(n2,fi,0.,0.))/x;
return(y);
}

float ei(void)
{extern float z1,n1,n2,de,d,d1,d2,r1,r2,zz[3],v01,v02;
float u0,u,u1,u2,u3,u4,u5,u6,y1,y,t1,t2,t3,t4,t5,t,beta,rs;
rs=pow((double)(3./(4.*PI*n1)),(double)(1./3.));
t=rs*rs;
t=t/(3*r1*r1*r1);
t1=2.21/rs;
t2=0.9*pow((double)(z1),(double)(2./3.));
t3=(4.5*r1*r1)/(rs*rs);
t4=0.458;
t5=0.00711*rs*rs;
t5=t5/(1+0.127*rs);
t5=t5/(1+0.127*rs);
v01=0;//t*(t1-t2+t3-t4-t5);
beta=zz[0];
u=(2.*PI*n1*n1)/(beta*beta*beta);
u0=beta*d1*exp(-(double)(0.5*beta*d1));
u1=1-exp(-(double)(beta*d1));
u2=u*(1.-u0*cosh((double)(beta*r1))/u1);
u3=r1*cosh((double)(beta*r1))-sinh((double)(beta*r1))/beta;
u4=v01*u*u3*u0/u1;
y=u2+u4;
return(SG*y);
}

float eid(void)
{extern float n1,n2,de,d,d1,d2,r1,r2,zz[3],v01,v02;
float u0,u,u1,u2,u3,u4,u5,u6,u7,u8,u9,y,beta;
beta=zz[0];
u=(4.*PI*de)/((de+1.)*beta*beta*beta);
u0=exp((double)(-2.*beta*d1/de));
u1=hi( );
u2=n1*(n1-(u0*n1*(de-1.)+n2*(de+1.))/u1);
u3=n2*(n2-(u0*n2*(de-1.)+n1*(de+1.))/u1);
y=u2*d1*(1.-exp((double)(beta*zz[1])))*exp((double)(-0.5*beta*d1))*cosh((double)(beta*r1));
y+=u3*d2*(1.-exp((double)(beta*zz[2])))*exp((double)(-0.5*beta*d2))*cosh((double)(beta*r2));
u4=2*PI*(n1*n1*d1*zz[1]*zz[1]+n2*n2*d2*zz[2]*zz[2]);
/*u5=2.*PI*(n1*n1*d1*zz[1]+n2*n2*d2*zz[2]);*/
/*u5/=beta;*/
y*=u*beta;
y+=u4;
/*y+=2*PI*n1*n1*d1*zz[1]*exp(-beta*(0.5*d1-zz[1]))/beta;*/
u5=2*PI*n1*n1*d1;
u5/=beta*beta;
u6=(1-exp(beta*zz[1]))*exp(-0.5*beta*d1);
u7=sinh(beta*r1);
u8=cosh(beta*r1);
u9=u5*u6*((v01/beta)*u7-v01*r1*u8);
y+=u9;
return(SG*y);
}

float ii(void)
{extern float z1,z2,c1,c2,d1,d2,d,de,zz[3];
float u1,u2,u3,u4,y;
u1=1.73205*z1*z1/(c1*c1*c1);
u1=u1*exp((double)(-(4.*PI*(d1-2*zz[1]))/(1.73205*c1)));
u2=1.73205*z2*z2/(c2*c2*c2);
u2=u2*exp((double)(-(4.*PI*(d2-2*zz[2]))/(1.73205*c2)));
u3=-2.*1.73205*z1*z2;
u3=u3/(pow((double)(c1*c2),3./2.));
u4=-(2.*PI/(1.73205))*((d1+d/de-2*zz[1]/de)/c1+(d2+d/de-2*zz[2]/de)/c2);
y=(u1+u2+u3*exp((double)(u4)))*SG;
return(y);
}

float me(void)
{extern float zz[3],d;
float me1,me2,me3,y,beta;
if(zz[0]<=0.) y=1.e25;else{
beta=zz[0];
me1=(1./beta)*t(0.,1.,w1);
me2=t(-d,d,w2);
me3=(1./beta)*t(0.,1.,w3);
y=SG*(me1+me2+me3)+gele*ei( )+gele*ii( )+gele*eid( );}
return(y);
}

float wku(float nx,float fi)
{
float y;
y=0.5*fi*nx;
return(y);
}

float wki(float nx)
{
float y;
y=0.3*pow((double)(3.*PI*PI),(double)(2./3.))*pow((double)nx,(double)(5./3.));
return(y);
}

float wob(float nx)
{
float y;
y=-0.75*pow((double)(3./PI),(double)(1./3.))*pow((double)nx,(double)(4./3.));
return(y);
}

float wcor(float nx)
{
float y;
y=-0.056*pow((double)nx,(double)(4./3.));
y=y/(0.079+pow((double)nx,(double)(1./3.)));
return(y);
}

float wvk(float nx,float dn)
{float y;
y=0.;
if(fabs((double)nx)>0.00000000001)
y=(1./72.)*(dn*dn)/nx;
return(y);
}

float wku1(float x)
{extern float n1;
float nx,fi,y;
nx=n11(x);
fi=fi1(x);
y=(wku(nx,fi)-wku(n1,fi))/x;
return(y);
}

float wku2(float x)
{float fi,nx,y;
nx=n22(x);
fi=fi2(x);
y=wku(nx,fi);
return(y);
}

float wku3(float x)
{extern float n2;
float fi,nx,y;
nx=n3(x);
fi=fi3(x);
y=(wku(nx,fi)-wku(n2,fi))/x;
return(y);
}

float wki1(float x)
{extern float n1;
float nx,y;
nx=n11(x);
y=(wki(nx)-wki(n1))/x;
return(y);
}

float wki2(float x)
{float nx,y;
nx=n22(x);
y=wki(nx);
return(y);
}

float wki3(float x)
{extern float n2;
float nx,y;
nx=n3(x);
y=(wki(nx)-wki(n2))/x;
return(y);
}

float wob1(float x)
{extern float n1;
float nx,y;
nx=n11(x);
y=(wob(nx)-wob(n1))/x;
return(y);
}

float wob2(float x)
{float nx,y;
nx=n22(x);
y=wob(nx);
return(y);
}

float wob3(float x)
{extern float n2;
float nx,y;
nx=n3(x);
y=(wob(nx)-wob(n2))/x;
return(y);
}

float wcor1(float x)
{extern float n1;
float nx,y;
nx=n11(x);
y=(wcor(nx)-wcor(n1))/x;
return(y);
}

float wcor2(float x)
{float nx,y;
nx=n22(x);
y=wcor(nx);
return(y);
}

float wcor3(float x)
{extern float n2;
float nx,y;
nx=n3(x);
y=(wcor(nx)-wcor(n2))/x;
return(y);
}

float wvk1(float x)
{extern float n1;
float nx,dn,y;
nx=n11(x);
dn=dn1(x);
y=(wvk(nx,dn)-wvk(n1,0.))/x;
return(y);
}

float wvk2(float x)
{float nx,dn,y;
nx=n22(x);
dn=dn2(x);
y=wvk(nx,dn);
return(y);
}

float wvk3(float x)
{extern float n2;
float nx,dn,y;
nx=n3(x);
dn=dn3(x);
y=(wvk(nx,dn)-wvk(n2,0.))/x;
return(y);
}

float sum1(float fn(float x))
{extern float zz[3];
float y,beta;
beta=zz[0];
y=SG*t(0.,1.,fn)/beta;
return(y);
}

float sum2(float fn(float x))
{extern float d;
float y;
y=SG*t(-d,d,fn);
return(y);
}

float sum3(float fn(float x))
{extern float zz[3];
float y,beta;
beta=zz[0];
y=SG*t(0.,1.,fn)/beta;
return(y);
}

float t(float x0,float x1,float f(float z))
{extern int m;
float h,a,b,a1,a2,g,z,y;
int i,j;
static float ag[8]={0.10122854,0.22238104,
		     0.31370664,0.36268378,
		     0.36268378,0.31370664,
		     0.22238104,0.10122854},
	     xg[8]={-0.96028986,-0.79666648,
		    -0.52553242,-0.18343464,
		     0.18343464,0.52553242,
		     0.79666648,0.96028986};
h=(x1-x0)/m;
a=x0;
b=x0+h;
y=0.;
 for(j=1;j<=m;j++){
 g=0.;
a1=0.5*(b+a);
a2=0.5*(b-a);
for(i=0;i<=7;i++){
     z=a1+a2*xg[i];
     g=g+ag[i]*f(z);
     }
y=y+g*a2;
a=b;
b=a+h;}
return(y);
}

float ww(float nx,float dn)
{float a,b,u,y;
if(fabs((double)nx)<0.000000000001) y=0;             else{
a=0.4666+0.3735*pow((double)(3.*PI*PI*nx),(double)(-2./9.));
b=-0.0085+0.3318*pow((double)(3.*PI*PI*nx),(double)(1./15.));
u=pow((double)PI,(double)(-5./3.));
y=(u*a*b*b*dn*dn);
y/=pow((double)(3.*nx),(double)(4./3.));}
return(y);
}

float ww1(float x)
{extern float n1;
float nx,dn,y;
nx=n11(x);
dn=dn1(x);
y=(ww(nx,dn)-ww(n1,0.))/x;
return(y);
}

float ww2(float x)
{float nx,dn,y;
nx=n22(x);
dn=dn2(x);
y=ww(nx,dn);
return(y);
}

float ww3(float x)
{extern float n2;
float nx,dn,y;
nx=n3(x);
dn=dn3(x);
y=(ww(nx,dn)-ww(n2,0.))/x;
return(y);
}

float svk4(float nx,float dn,float dnn)
{float a,b,b2,c,d,e,y;
if(fabs((double)nx)<0.0000001) y=0;else {
e=(dn/nx)*(dn/nx);
a=nx/(540*pow((double)(3*PI*PI*nx),(double)(2./3.)));
b=(dnn/nx);
b2=b*b;
c=1.125*b*e;
d=0.3333333333*e*e;
y=1.336*a*(b2-c+d);
}
return(y);
}

float sww4(float nx,float dnn)
{float a,b,y;
if(fabs((double)nx)<0.0000001) y=0;else {
a=0.001766;
b=-0.2986*pow((double)nx,-0.26);
y=dnn*dnn*a*exp((double)b);
y/=(nx*nx);
}
return(0.1666666667*y);
}

float svk41(float x)
{extern float n1;
float nx,dn,y,dnn;
nx=n11(x);
dn=dn1(x);
dnn=dn21(x);
y=(svk4(nx,dn,dnn)-svk4(n1,0.,0.))/x;
return(y);
}

float svk42(float x)
{float nx,dn,y,dnn;
nx=n22(x);
dn=dn2(x);
dnn=dn22(x);
y=svk4(nx,dn,dnn);
return(y);
}

float svk43(float x)
{extern float n2;
float nx,dn,y,dnn;
nx=n3(x);
dn=dn3(x);
dnn=dn23(x);
y=(svk4(nx,dn,dnn)-svk4(n2,0.,0.))/x;
return(y);
}

float sww41(float x)
{extern float n1;
float nx,dnn,y;
nx=n11(x);
dnn=dn21(x);
y=(sww4(nx,dnn)-sww4(n1,0.))/x;
return(y);
}

float sww42(float x)
{float nx,dnn,y;
nx=n22(x);
dnn=dn22(x);
y=sww4(nx,dnn);
return(y);
}

float sww43(float x)
{extern float n2;
float nx,dnn,y;
nx=n3(x);
dnn=dn23(x);
y=(sww4(nx,dnn)-sww4(n2,0.))/x;
return(y);
}