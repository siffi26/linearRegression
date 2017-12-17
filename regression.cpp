/*Program to do regularized linear model (polynomial basis) regression*/
/*Author: Siffi Singh */
/*Dated: 27/9/2017 */

/* standard Headers */
#include<iostream>
#include<cmath>
#include<cstdio>
#include<cstring>
#include<algorithm>
#include<fstream>
using namespace std;
/* The function that calculates the LU decomposed matrix.*/
void LU(float(*D)[3][3],int n)
{
    int i,j,k,m,an,am;
    float x;
    printf("The matrix: \n");
    for(j=0; j<=2; j++)
    {
        printf(" %f  %f  %f \n",(*D)[j][0],(*D)[j][1],(*D)[j][2]);
    }
    for(k=0; k<=n-1; k++)
    {
        for(j=k+1; j<=n; j++)
        {
            x=(*D)[j][k]/(*D)[k][k];
            for(i=k; i<=n; i++)
            {
                (*D)[j][i]=(*D)[j][i]-x*(*D)[k][i];
            }
            (*D)[j][k]=x;
        }
    }
}
int main ()
{
    /* Variable declarations */
    int n,m,an,am;
    float d[3],C[3][3]= {0}, B[3][1]= {0};
    float x,s[3][3]= {0},y[3], X[3][1]= {0};
    FILE *FP,*fp1;
    an=3;
    am=3;
    n=2;
    string line;
    int i=0,j=0,k=0,len,last=0;
    string num = "";
    float A[10][10]= {0};
    /*Storing the data from csv file to a matrix*/
    ifstream infile ("data.csv");
    if(infile.is_open())
    {
        while ( getline(infile,line) )
        {
            // Takes complete row
            // cout << line<< '\n';
            k=0,last=0,j=0;
            len=line.length();

            while(k!=len) {
                if(line[k]==','||k==len-1) {
                    //for converting string into number
                    A[i][j]=atof(num.append(line,last,k).c_str());
                    //Emtying string for getting next data
                    num="";
                    //increasing column number after saving data
                    j++;
                    if(k!=len-1)
                        last=k+1;
                }
                k++;
            }
            //increase row number for array
            i++;
        }
        //close the file
        infile.close();
    }
    else cout << "Unable to open file";
    cout<<"Matrix of CSV file below:\n";
    int col=0;

    //For finding the number of columns
    for(int x=0; x<len; x++) {
        if(line[x]==',')
            col++;
    }
    col++;
    // i= number of rows
    // col = number of columns
    for(int l=0; l<i; l++) {
        for(int m=0; m<col; m++) {
            if(m==col-1)
            {
                B[l][0]=A[l][col-1];
                A[l][col-1]=1;
            }
            cout<<A[l][m]<<"\t";
        }
        cout<<"\n";
    }
    /*Inputing no. of bases*/
    int bases;
    cout<<"No. of bases: ";
    cin>>bases;
    /*Inputing value of lambda*/
    int lambda;
    cout<<"Lambda: ";
    cin>>lambda;
    /*Creating Transpose of 'A' matrix in 'AT'*/
    int AT[i][col]= {0};
    for(int l=0; l<i; l++)
    {
        for(int m=0; m<col; m++)
        {
            AT[l][m]=A[m][l];
        }
    }

    /*'P' matrix stores the product of AT('A' transpose) and A*/
    float P[i][i]= {0};
    for(int d=0; d<i; d++)
    {
        for(int e=0; e<col; e++)
        {
            for(int f=0; f<col; f++)
            {
                P[d][e] += AT[d][f] * A[f][e];
            }
        }
    }
    /*'L' matrix stores lambda*I(Identity) matrix*/
    float L[i][i]= {0};
    for(int l=0; l<i; l++)
    {
        for(int m=0; m<i; m++)
        {
            if(l==m)
                L[l][m]=lambda;
        }
    }
    /*'Q' matrix stores AT*A + lambda*I. Hence, Q[][] = AT*A + lambda*I*/
    float Q[3][3]= {0};
    for(int l=0; l<i; l++)
    {
        for(int m=0; m<i; m++)
        {
            Q[l][m] = P[l][m] + L[l][m];
        }
    }
    /*Finding Inverse of 'Q' matrix using LU decmposition*/
    float (*b)[3][3];
    b=&Q;
    for(m=0; m<=i; m++)
    {
        for(j=0; j<=i; j++)
        {
            C[m][j]=Q[m][j];
        }
    }
    /* Call a sub-function to calculate the LU decomposed matrix.*/
    LU(b,n);
    printf(" \n");
    printf("The matrix after LU decomposition is:\n");
    for(m=0; m<=2; m++)
    {
        for(j=0; j<=2; j++)
        {
            cout<<Q[m][j]<<" ";
        }
        cout<<endl;
    }
    /*Finding the inverse with the above computed LU decomposd matrix.*/
    for(m=0; m<=2; m++)
    {
        d[0]=0.0;
        d[1]=0.0;
        d[2]=0.0;
        d[m]=1.0;
        for(i=0; i<=n; i++)
        {
            x=0.0;
            for(j=0; j<=i-1; j++)
            {
                x=x+Q[i][j]*y[j];
            }
            y[i]=(d[i]-x);
        }
        for(i=n; i>=0; i--)
        {
            x=0.0;
            for(j=i+1; j<=n; j++)
            {
                x=x+Q[i][j]*s[j][m];
            }
            s[i][m]=(y[i]-x)/Q[i][i];
        }
    }
    /* Print the inverse matrix */
    printf("\nThe Inverse Matrix of AT*A + lambda*I:\n");
    for(m=0; m<=2; m++)
    {
        for(j=0; j<=2; j++)
        {
            cout<<s[m][j]<<" ";
        }
        cout<<endl;
    }
    /* check that the product of the matrix with its inverse results is indeed a unit matrix
    printf("\nThe product of Q and Q^-1 is:\n");
    for(m=0;m<=2;m++)
    {
    	for(j=0;j<=2;j++)
    	{
    	x=0.0;
    		for(i=0;i<=2;i++)
    		{
    			x=x+C[m][i]*s[i][j];
    		}
    	printf(" %d    %d    %f \n", m,j,x );
    	}
    } */


    /*After computing inverse of (AT*A + lambda*I), we need to find (AT*A + lambda*I).AT */
    float R[3][3]= {0};
    for(int d=0; d<3; d++)
    {
        for(int e=0; e<3; e++)
        {
            for(int f=0; f<3; f++)
            {
                R[d][e] += s[d][f] * AT[f][e];
            }
        }
    }
    /*After computing inverse of (AT*A + lambda*I).AT, we need to find (AT*A + lambda*I).AT.b */
    for(int d=0; d<3; d++)
    {
        for(int e=0; e<3; e++)
        {
            for(int f=0; f<1; f++)
            {
                X[d][e] += R[d][f] * B[f][e];
            }
        }
    }
    /*Printing (AT*A + lambda*I).AT.b */
    cout<<"\nThe final parameters are: \n";
    for(int l=0; l<3; l++)
    {
        for(int m=0; m<1; m++)
        {
            cout<<(char)(97+l)<<": "<<X[l][m]<<" ";
        }
        cout<<endl;
    }
    /*Calculating ||AX-b||^2 */
    float E[3][1]= {0};
    for(int d=0; d<3; d++)
    {
        for(int e=0; e<3; e++)
        {
            for(int f=0; f<1; f++)
            {
                E[d][e] += A[d][f] * X[f][e];
            }
        }
    }
    for(int l=0; l<3; l++)
    {
        E[l][0]-=B[l][0];
    }
    /*Printing ||AX-b||^2 */
    cout<<"\nError is: \n";
    for(int l=0; l<3; l++)
    {
        for(int m=0; m<1; m++)
        {
            cout<<pow(E[l][m],2)<<" ";
        }
        cout<<endl;
    }
    return 0;
}


