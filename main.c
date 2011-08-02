#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

/* get the orintation of the focal plane and position of the optical center using 3 fudical points*/

/* You will need to initialize Objectx, Objecty, Objectz, focallengthm (meters),
   Imagex, Imagey, Testpixelx, Testpixely (pixels)
   fieldofviewd (degrees)
   image coordinates are assumed to be 0, 0 for the image center */

double Objectx[3],Objecty[3],Objectz[3],Imagex[3], Imagey[3],focallengthp,Object_Image[8][3][3],Pinhole[8][3],ray[3];
int realcount;

int main()
{
    int resection(void);

    double focallengthm,fieldofviewr,fieldofviewd,pixels_wide,pixels_high,pixels_diagonal;
    double Testpixelx,Testpixely,Testobjectx,Testobjecty,Testobjectz,Xpoint,Ypoint,Zpoint,testimage[3][3],testobject[3][3];
    double OriginalImagex[3],OriginalImagey[3],ii,io,oo;
    int index,row,rvalue;

pixels_wide=4288;
pixels_high=2848;

OriginalImagex[0] =  -419;   OriginalImagey[0] = -714;
OriginalImagex[1] = -1741;   OriginalImagey[1] =  -87;
OriginalImagex[2] =   452;   OriginalImagey[2] =  601;

Objectx[0] =  0.0;    Objecty[0] =  0.00;   Objectz[0] = 0.0;
Objectx[1] =  0.0;    Objecty[1] = 0.908;   Objectz[1] = 0.0;
Objectx[2] = 1.83;    Objecty[2] = 0.908;   Objectz[2] = 0.0;

/*
  Camera focal length was set to 28 mm
*/

focallengthm = 0.028; // 28 millimeters
fieldofviewd = 40.0;  // degrees

Testpixelx =   124.0;
Testpixely =   -84.0;
Testobjectx=    .823;
Testobjecty=    .288;
Testobjectz=    .022;

//------------------ end Mick's changes -------------------------------


/* put the focal length and image coordinates into the same dimensions
   for proper scaling; if you know the physical dimension of imaging array this could be simplified */

    fieldofviewr=fieldofviewd*M_PI/180;
    focallengthp=pixels_wide/tan(fieldofviewr/2)/2;
    printf("focal length in pixels being taken as %f\n",focallengthp);

/* field of view could refer to the horizontal or vertical angle or diagonal
this choice could be part of the initialization

    pixels_diagonal=sqrt(pixels_high*pixels_high+pixels_wide*pixels_wide);
    focallengthp=pixels_diagonal/tan(fieldofviewr/2)/2;

    focallengtp=pixels_high/tan(fieldofviewr/2)/2;

    focallengtp=vpixels/tan(fieldofviewr/2)/2; */


    for(index=0;index<3;++index){
        Imagex[index]=OriginalImagex[index]/focallengthp;
        Imagey[index]=OriginalImagey[index]/focallengthp;
    }

/* The results that depend on both the transform and the pinhole will be the same but either one separately will look wrong.
Think of it like this: if you lay the picture on the plane of the object triangle then the
sense of of going from positive x to positive y should be the same for image and object
coordinates. */

    rvalue=resection();

    printf("There are %i real solution out 8 possible.\n",realcount);
    if(realcount==0)
        return(-1);

    for(index=0;index<realcount;++index)
        for(row=0;row<3;++row){
            testimage[row][0]=Object_Image[index][0][0]*Imagex[row]+Object_Image[index][0][1]*Imagey[row]+Object_Image[index][0][2];
            testimage[row][1]=Object_Image[index][1][0]*Imagex[row]+Object_Image[index][1][1]*Imagey[row]+Object_Image[index][1][2];
            testimage[row][2]=Object_Image[index][2][0]*Imagex[row]+Object_Image[index][2][1]*Imagey[row]+Object_Image[index][2][2];
            testobject[row][0]=Objectx[row]-Pinhole[index][0];
            testobject[row][1]=Objecty[row]-Pinhole[index][1];
            testobject[row][2]=Objectz[row]-Pinhole[index][2];
            ii=testimage[row][0]*testimage[row][0]+testimage[row][1]*testimage[row][1]+testimage[row][2]*testimage[row][2];
            io=testimage[row][0]*testobject[row][0]+testimage[row][1]*testobject[row][1]+testimage[row][2]*testobject[row][2];
            oo=testobject[row][0]*testobject[row][0]+testobject[row][1]*testobject[row][1]+testobject[row][2]*testobject[row][2];
            printf("for scale %i point %i the least square error is %e\n",index,row,oo-io*io/ii);
        }
    printf("\n");

    printf("For image point is %f %f\n",Testpixelx,Testpixely);
    printf("The object point was %f %f %f\n",Testobjectx,Testobjecty,Testobjectz);
    printf("using the object z as %f to set the scale\n",Testobjectz);
    for(index=0;index<realcount;++index){
        testimage[0][0]=Xpoint=(Object_Image[index][0][0]*Testpixelx/focallengthp+Object_Image[index][0][1]*Testpixely/focallengthp+focallengthp*Object_Image[index][0][2]);
        testimage[0][1]=Ypoint=(Object_Image[index][1][0]*Testpixelx/focallengthp+Object_Image[index][1][1]*Testpixely/focallengthp+focallengthp*Object_Image[index][1][2]);
        testimage[0][2]=Zpoint=(Object_Image[index][2][0]*Testpixelx/focallengthp+Object_Image[index][2][1]*Testpixely/focallengthp+focallengthp*Object_Image[index][2][2]);
        Xpoint=Xpoint*(Testobjectz-Pinhole[index][2])/Zpoint+Pinhole[index][0];
        Ypoint=Ypoint*(Testobjectz-Pinhole[index][2])/Zpoint+Pinhole[index][1];
        printf("for scale %i the object position is %f %f 0.022\n",index,Xpoint,Ypoint);
        testobject[0][0]=Testobjectx-Pinhole[index][0];
        testobject[0][1]=Testobjecty-Pinhole[index][1];
        testobject[0][2]=Testobjectz-Pinhole[index][2];
        ii=testimage[0][0]*testimage[0][0]+testimage[0][1]*testimage[0][1]+testimage[0][2]*testimage[0][2];
        io=testimage[0][0]*testobject[0][0]+testimage[0][1]*testobject[0][1]+testimage[0][2]*testobject[0][2];
        oo=testobject[0][0]*testobject[0][0]+testobject[0][1]*testobject[0][1]+testobject[0][2]*testobject[0][2];
        printf("for the test point the least error is %e\n",sqrt(oo-io*io/ii));
    }
    printf("\n");

    return 0;
}

int resection()
{
    double cf[11][10],CubicCoef[3],adjustor,zerotest,generator[3];
    double DiffObject[3][3],con1,con2,con3,offset,scale;
    double mx1[3][3],mx2[3][3],mx3[3][3],mx4[3][3],mx5[3][3],mx6[3][3];
    double cubic1,scalelist[8][3],ScaledImage[3][3],DifScaImg[3][3],twist[2];
    complex double eigenvalue1,eigenvalue2,eigenvalue3,eigenvector1[3],eigenvector2[3],eigenvector3[3];
    complex double cmx1[3][3],cmx2[3][3],cmx3[3][3],ccf[3][3],cadj,complexscales[8][3],testscales,cubic2;
    int index,row,column;

    zerotest=1.e-14;

/*The base vector equation is;

scale[point_index]*rotation_matrix.image_vector[point_index]+pinhole_position-object_vector[point_index]

*/

    DiffObject[0][0]=Objectx[1]-Objectx[2];
    DiffObject[0][1]=Objecty[1]-Objecty[2];
    DiffObject[0][2]=Objectz[1]-Objectz[2];
    DiffObject[1][0]=Objectx[2]-Objectx[0];
    DiffObject[1][1]=Objecty[2]-Objecty[0];
    DiffObject[1][2]=Objectz[2]-Objectz[0];
    DiffObject[2][0]=Objectx[0]-Objectx[1];
    DiffObject[2][1]=Objecty[0]-Objecty[1];
    DiffObject[2][2]=Objectz[0]-Objectz[1];

/* initialize the coefficients
if we subtract the third base equation from the second, the first from the third and the second from the first
then take the vector magnitude then mx's are coefficients for the quadratic scale terms and the con's for objects */

    mx1[0][0]=0;
    mx1[1][1]=Imagex[1]*Imagex[1]+Imagey[1]*Imagey[1]+1;
    mx1[2][2]=Imagex[2]*Imagex[2]+Imagey[2]*Imagey[2]+1;
    mx1[0][1]=mx1[1][0]=0;
    mx1[0][2]=mx1[2][0]=0;
    mx1[1][2]=mx1[2][1]=-Imagex[1]*Imagex[2]-Imagey[1]*Imagey[2]-1;
    con1=-(DiffObject[0][0]*DiffObject[0][0]+DiffObject[0][1]*DiffObject[0][1]+DiffObject[0][2]*DiffObject[0][2]);

    mx2[0][0]=Imagex[0]*Imagex[0]+Imagey[0]*Imagey[0]+1;
    mx2[1][1]=0;
    mx2[2][2]=Imagex[2]*Imagex[2]+Imagey[2]*Imagey[2]+1;
    mx2[0][1]=mx2[1][0]=0;
    mx2[0][2]=mx2[2][0]=-Imagex[0]*Imagex[2]-Imagey[0]*Imagey[2]-1;
    mx2[1][2]=mx2[2][1]=0;
    con2=-(DiffObject[1][0]*DiffObject[1][0]+DiffObject[1][1]*DiffObject[1][1]+DiffObject[1][2]*DiffObject[1][2]);

    mx3[0][0]=Imagex[0]*Imagex[0]+Imagey[0]*Imagey[0]+1;
    mx3[1][1]=Imagex[1]*Imagex[1]+Imagey[1]*Imagey[1]+1;
    mx3[2][2]=0;
    mx3[0][1]=mx3[1][0]=-Imagex[0]*Imagex[1]-Imagey[0]*Imagey[1]-1;
    mx3[0][2]=mx3[2][0]=0;
    mx3[1][2]=mx3[2][1]=0;
    con3=-(DiffObject[2][0]*DiffObject[2][0]+DiffObject[2][1]*DiffObject[2][1]+DiffObject[2][2]*DiffObject[2][2]);

/* to solve these we need a llinear combination of the three equations that depends on just two variables
first we find two equations that elimate the con's */

    for(row=0;row<3;++row)
        for(column=0;column<3;++column){
            mx4[row][column]=mx1[row][column]-con1*mx3[row][column]/con3;
            mx5[row][column]=mx2[row][column]-con2*mx3[row][column]/con3;
        }

/* It is easier to find the right combination of two equations if we first convert one to a diagonal form */

    cf[0][0]=mx5[1][1]*mx5[2][2]-mx5[1][2]*mx5[2][1];
    cf[1][1]=mx5[0][0]*mx5[2][2]-mx5[0][2]*mx5[2][0];
    cf[2][2]=mx5[1][1]*mx5[0][0]-mx5[1][0]*mx5[0][1];
    cf[0][1]=cf[1][0]=mx5[0][2]*mx5[2][1]-mx5[0][1]*mx5[2][2];
    cf[0][2]=cf[2][0]=mx5[0][1]*mx5[1][2]-mx5[0][2]*mx5[1][1];
    cf[1][2]=cf[2][1]=mx5[1][0]*mx5[0][2]-mx5[1][2]*mx5[0][0];
    adjustor=mx5[0][0]*cf[0][0]+mx5[0][1]*cf[1][0]+mx5[0][2]*cf[2][0];

    for(row=0;row<3;++row)
        for(column=0;column<3;++column)
            mx6[row][column]=(cf[row][0]*mx4[0][column]+cf[row][1]*mx4[1][column]+cf[row][2]*mx4[2][column])/adjustor;

    cf[0][0]=mx6[1][1]*mx6[2][2]-mx6[1][2]*mx6[2][1];
    cf[1][1]=mx6[0][0]*mx6[2][2]-mx6[0][2]*mx6[2][0];
    cf[2][2]=mx6[1][1]*mx6[0][0]-mx6[1][0]*mx6[0][1];
    cf[0][1]=mx6[0][2]*mx6[2][1]-mx6[0][1]*mx6[2][2];
    cf[1][0]=mx6[1][2]*mx6[2][0]-mx6[1][0]*mx6[2][2];
    cf[0][2]=mx6[0][1]*mx6[1][2]-mx6[0][2]*mx6[1][1];
    cf[2][0]=mx6[2][1]*mx6[1][0]-mx6[2][0]*mx6[1][1];
    cf[1][2]=mx6[1][0]*mx6[0][2]-mx6[1][2]*mx6[0][0];
    cf[2][1]=mx6[2][0]*mx6[0][1]-mx6[2][1]*mx6[0][0];

/* build a cubic equation Determinent(matrix-scalar) the coefficients of which are */

    CubicCoef[0]=mx6[0][0]*cf[0][0]+mx6[0][1]*cf[1][0]+mx6[0][2]*cf[2][0];
    CubicCoef[1]=cf[0][0]+cf[1][1]+cf[2][2];
    CubicCoef[2]=mx6[0][0]+mx6[1][1]+mx6[2][2];

    offset=CubicCoef[2]/3;
    CubicCoef[0]+=2*offset*offset*offset-CubicCoef[1]*offset;
    CubicCoef[1]-=offset*offset*3;
    scale=CubicCoef[1]*4/3;
    if((scale<zerotest)&&(scale>-zerotest))
        scale=1;
    if(scale<-zerotest)
        scale*=-1;
    scale=sqrt(scale);
    if(CubicCoef[0]>0)
        scale*=-1;
    CubicCoef[0]/=scale*scale*scale;

/* the solution has four forms depending the sign of CubicCoef[1] and for the negative case the magnitude of CubicCoef[0] */

    if((CubicCoef[1]<zerotest)&&(CubicCoef[1]>-zerotest)){
        cubic1=pow(CubicCoef[0],1./3);
        cubic2=cubic1*csqrt(-3)/2;
    }
    if(CubicCoef[1]>zerotest){
        cubic1=asinh(-4*CubicCoef[0])/3;
        cubic2=cosh(cubic1)*csqrt(-3)/2;
        cubic1=sinh(cubic1);
    }
    if(CubicCoef[1]<zerotest){
        if(CubicCoef[0]<-1./4){
            cubic1=acosh(-4*CubicCoef[0])/3;
            cubic2=sinh(cubic1)*csqrt(-3)/2;
            cubic1=cosh(cubic1);
        }
        else{
            cubic1=acos(-4*CubicCoef[0])/3;
            cubic2=sin(cubic1)*sqrt(3)/2;
            cubic1=cos(cubic1);
        }
    }

    eigenvalue1=scale*cubic1-offset;
    eigenvalue2=scale*(-cubic1/2+cubic2)-offset;
    eigenvalue3=scale*(-cubic1/2-cubic2)-offset;

/* now we combine mx4 and mx5 using the eigenvalues the cadj is a transformed scale item from the
base equations the results to eliminate of the two of the three scale unknowns
understanding how this actually works is a linear algebra thing */

    for(row=0;row<3;++row)
        for(column=0;column<3;++column){
            cmx1[row][column]=mx4[row][column]+eigenvalue1*mx5[row][column];
            cmx2[row][column]=mx4[row][column]+eigenvalue2*mx5[row][column];
            cmx3[row][column]=mx4[row][column]+eigenvalue3*mx5[row][column];
        }

    eigenvector1[0]=cmx1[1][1]*cmx1[2][2]-cmx1[1][2]*cmx1[2][1];
    eigenvector1[1]=cmx1[0][2]*cmx1[2][1]-cmx1[0][1]*cmx1[2][2];
    eigenvector1[2]=cmx1[0][1]*cmx1[1][2]-cmx1[0][2]*cmx1[1][1];

    eigenvector2[0]=cmx2[1][1]*cmx2[2][2]-cmx2[1][2]*cmx2[2][1];
    eigenvector2[1]=cmx2[0][2]*cmx2[2][1]-cmx2[0][1]*cmx2[2][2];
    eigenvector2[2]=cmx2[0][1]*cmx2[1][2]-cmx2[0][2]*cmx2[1][1];

    eigenvector3[0]=cmx3[1][1]*cmx3[2][2]-cmx3[1][2]*cmx3[2][1];
    eigenvector3[1]=cmx3[0][2]*cmx3[2][1]-cmx3[0][1]*cmx3[2][2];
    eigenvector3[2]=cmx3[0][1]*cmx3[1][2]-cmx3[0][2]*cmx3[1][1];

    ccf[1][1]=ccf[2][2]=0;
    for(row=0;row<3;++row)
        for(column=0;column<3;++column){
            ccf[1][1]+=eigenvector2[row]*cmx1[row][column]*eigenvector2[column];
            ccf[2][2]+=eigenvector3[row]*cmx1[row][column]*eigenvector3[column];
        }
    cadj=csqrt(-ccf[2][2]/ccf[1][1]);
    for(row=0;row<3;++row)
        eigenvector2[row]*=cadj;

    ccf[0][0]=ccf[2][2]=0;
    for(row=0;row<3;++row)
        for(column=0;column<3;++column){
            ccf[0][0]+=eigenvector1[row]*cmx2[row][column]*eigenvector1[column];
            ccf[2][2]+=eigenvector3[row]*cmx2[row][column]*eigenvector3[column];
        }
    cadj=csqrt(-ccf[2][2]/ccf[0][0]);
    for(row=0;row<3;++row)
        eigenvector1[row]*=cadj;

/* now use the weighted eigenvectors to build a new set of vectors leaving on unknown
using the even cvector slots to allow for sign changes in the last scale veriable */

    for(row=0;row<3;++row){
        complexscales[0][row]=eigenvector1[row]+eigenvector2[row]+eigenvector3[row];
        complexscales[2][row]=-eigenvector1[row]+eigenvector2[row]+eigenvector3[row];
        complexscales[4][row]=eigenvector1[row]-eigenvector2[row]+eigenvector3[row];
        complexscales[6][row]=-eigenvector1[row]-eigenvector2[row]+eigenvector3[row];
    }

/* use the complexscales to reduce the mx1 to get a single equation on the remaining scale variable
the result would be the same for mx2 or mx3, another niffty linear algebra thing */

    for(index=0;index<4;++index){
        cadj=0;
        for(row=0;row<3;++row)
            for(column=0;column<3;++column){
                cadj+=complexscales[2*index][row]*mx1[row][column]*complexscales[2*index][column];
            }
        cadj=csqrt(-con1/cadj);
        for(row=0;row<3;++row){
            complexscales[2*index][row]*=cadj;
            complexscales[2*index+1][row]=-complexscales[2*index][row];
        }
    }

/* for(index=0;index<4;++index)
printf("(%f,%f) (%f,%f) (%f,%f)\n",creal(complexscales[2*index][0]),cimag(complexscales[2*index][0]),
       creal(complexscales[2*index][1]),cimag(complexscales[2*index][1]),creal(complexscales[2*index][2]),cimag(complexscales[2*index][2])); */

/* we only want the real scales */

    realcount=0;
    for(index=0;index<8;++index){
        testscales=complexscales[index][0]*complexscales[index][0]+complexscales[index][1]*complexscales[index][1]+
            complexscales[index][2]*complexscales[index][2];
        adjustor=cimag(testscales)/creal(testscales);
        if((adjustor<zerotest)&&(adjustor>-zerotest)){
            for(row=0;row<3;++row)
                scalelist[realcount][row]=creal(complexscales[index][row]);
            ++realcount;
        }
    }

/* got the scales now calcultate the transform
for(index=0;index<realcount;++index)
printf("%f %f %f\n",scalelist[index][0],scalelist[index][1],scalelist[index][2]); */

    for(index=0;index<realcount;++index){
        for(row=0;row<3;++row){
            ScaledImage[row][0]=scalelist[index][row]*Imagex[row];
            ScaledImage[row][1]=scalelist[index][row]*Imagey[row];
            ScaledImage[row][2]=scalelist[index][row];
        }
        for(row=0;row<3;++row){
            DifScaImg[0][row]=ScaledImage[1][row]-ScaledImage[2][row];
            DifScaImg[1][row]=ScaledImage[2][row]-ScaledImage[0][row];
            DifScaImg[2][row]=ScaledImage[0][row]-ScaledImage[1][row];
        }
        for(row=0;row<3;++row)
            for(column=0;column<3;++column)
                mx1[row][column]=DifScaImg[row][column]-DiffObject[row][column];

        generator[0]=mx1[1][1]*mx1[2][2]-mx1[1][2]*mx1[2][1];
        generator[1]=mx1[1][2]*mx1[2][0]-mx1[1][0]*mx1[2][2];
        generator[2]=mx1[1][0]*mx1[2][1]-mx1[1][1]*mx1[2][0];
        adjustor=sqrt(generator[0]*generator[0]+generator[1]*generator[1]+generator[2]*generator[2]);
        for(row=0;row<3;++row)
            generator[row]/=adjustor;

/* could probably safely limit this to six (2?) equations but three (7) extra calculation are cheap */

        cf[0][0]=generator[2]*DifScaImg[0][1]-generator[1]*DifScaImg[0][2];
        cf[1][0]=generator[0]*DifScaImg[0][2]-generator[2]*DifScaImg[0][0];
        cf[2][0]=generator[1]*DifScaImg[0][0]-generator[0]*DifScaImg[0][1];
        cf[3][0]=generator[2]*DifScaImg[1][1]-generator[1]*DifScaImg[1][2];
        cf[4][0]=generator[0]*DifScaImg[1][2]-generator[2]*DifScaImg[1][0];
        cf[5][0]=generator[1]*DifScaImg[1][0]-generator[0]*DifScaImg[1][1];
        cf[6][0]=generator[2]*DifScaImg[2][1]-generator[1]*DifScaImg[2][2];
        cf[7][0]=generator[0]*DifScaImg[2][2]-generator[2]*DifScaImg[2][0];
        cf[8][0]=generator[1]*DifScaImg[2][0]-generator[0]*DifScaImg[2][1];

        cf[0][1]=generator[0]*generator[0]*DifScaImg[0][0] + generator[0]*generator[1]*DifScaImg[0][1] + generator[0]*generator[2]*DifScaImg[0][2]-DifScaImg[0][0];
        cf[1][1]=generator[0]*generator[1]*DifScaImg[0][0] + generator[1]*generator[1]*DifScaImg[0][1] + generator[1]*generator[2]*DifScaImg[0][2]-DifScaImg[0][1];
        cf[2][1]=generator[0]*generator[2]*DifScaImg[0][0] + generator[1]*generator[2]*DifScaImg[0][1] + generator[2]*generator[2]*DifScaImg[0][2]-DifScaImg[0][2];
        cf[3][1]=generator[0]*generator[0]*DifScaImg[1][0] + generator[0]*generator[1]*DifScaImg[1][1] + generator[0]*generator[2]*DifScaImg[1][2]-DifScaImg[1][0];
        cf[4][1]=generator[0]*generator[1]*DifScaImg[1][0] + generator[1]*generator[1]*DifScaImg[1][1] + generator[1]*generator[2]*DifScaImg[1][2]-DifScaImg[1][1];
        cf[5][1]=generator[0]*generator[2]*DifScaImg[1][0] + generator[1]*generator[2]*DifScaImg[1][1] + generator[2]*generator[2]*DifScaImg[1][2]-DifScaImg[1][2];
        cf[6][1]=generator[0]*generator[0]*DifScaImg[2][0] + generator[0]*generator[1]*DifScaImg[2][1] + generator[0]*generator[2]*DifScaImg[2][2]-DifScaImg[2][0];
        cf[7][1]=generator[0]*generator[1]*DifScaImg[2][0] + generator[1]*generator[1]*DifScaImg[2][1] + generator[1]*generator[2]*DifScaImg[2][2]-DifScaImg[2][1];
        cf[8][1]=generator[0]*generator[2]*DifScaImg[2][0] + generator[1]*generator[2]*DifScaImg[2][1] + generator[2]*generator[2]*DifScaImg[2][2]-DifScaImg[2][2];

        cf[0][2]=mx1[0][0];
        cf[1][2]=mx1[0][1];
        cf[2][2]=mx1[0][2];
        cf[3][2]=mx1[1][0];
        cf[4][2]=mx1[1][1];
        cf[5][2]=mx1[1][2];
        cf[6][2]=mx1[2][0];
        cf[7][2]=mx1[2][1];
        cf[8][2]=mx1[2][2];

        for(row=0;row<9;++row)
            if((cf[row][0]>zerotest)||(cf[row][0]<-zerotest)){
                cf[9][0]=1;
                cf[9][1]=cf[row][1]/cf[row][0];
                cf[9][2]=cf[row][2]/cf[row][0];
                break;
            }
        for(row=0;row<9;++row){
            adjustor=cf[row][0];
            for(column=0;column<3;++column)
                cf[row][column]-=cf[9][column]*adjustor;
        }

        for(row=0;row<9;++row)
            if((cf[row][1]>zerotest)||(cf[row][1]<-zerotest)){
                cf[10][0]=cf[row][0]/cf[row][1];
                cf[10][1]=1;
                cf[10][2]=cf[row][2]/cf[row][1];
                break;
            }
        for(row=0;row<10;++row){
            adjustor=cf[row][1];
            for(column=0;column<3;++column)
                cf[row][column]-=cf[10][column]*adjustor;
        }

/* the coefficients represented by twists could be derived analyticly but this is easier */

        twist[0]=-cf[9][2];
        twist[1]=-cf[10][2];
        Object_Image[index][0][0]=1-twist[1]+generator[0]*generator[0]*twist[1];
        Object_Image[index][0][1]=generator[0]*generator[1]*twist[1]+generator[2]*twist[0];
        Object_Image[index][0][2]=generator[2]*generator[0]*twist[1]-generator[1]*twist[0];
        Object_Image[index][1][0]=generator[0]*generator[1]*twist[1]-generator[2]*twist[0];
        Object_Image[index][1][1]=1-twist[1]+generator[1]*generator[1]*twist[1];
        Object_Image[index][1][2]=generator[1]*generator[2]*twist[1]+generator[0]*twist[0];
        Object_Image[index][2][0]=generator[2]*generator[0]*twist[1]+generator[1]*twist[0];
        Object_Image[index][2][1]=generator[1]*generator[2]*twist[1]-generator[0]*twist[0];
        Object_Image[index][2][2]=1-twist[1]+generator[2]*generator[2]*twist[1];

/* got the scale and rotation matrix how get the pinhole position */

        Pinhole[index][0]=Objectx[0]-(Object_Image[index][0][0]*ScaledImage[0][0]+Object_Image[index][0][1]*ScaledImage[0][1]+
            Object_Image[index][0][2]*ScaledImage[0][2]);
        Pinhole[index][1]=Objecty[0]-(Object_Image[index][1][0]*ScaledImage[0][0]+Object_Image[index][1][1]*ScaledImage[0][1]+
            Object_Image[index][1][2]*ScaledImage[0][2]);
        Pinhole[index][2]=Objectz[0]-(Object_Image[index][2][0]*ScaledImage[0][0]+Object_Image[index][2][1]*ScaledImage[0][1]+
            Object_Image[index][2][2]*ScaledImage[0][2]);

    }

    return(0);
}
