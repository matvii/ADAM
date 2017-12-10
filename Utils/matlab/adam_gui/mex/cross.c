#include"prepare.h"
/*Calculate the cross product of two vectors A and B, place the result
 * in the vector C
 */
void cross(double *A,double *B,double *C)
{
  C[0]=A[1]*B[2]-A[2]*B[1];
  C[1]=A[2]*B[0]-A[0]*B[2];
  C[2]=A[0]*B[1]-A[1]*B[0];
}