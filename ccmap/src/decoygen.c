#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>



//--------------------------------------------------
// Function rotateAtom
//--------------------------------------------------
void rotateAtom (float oldX, float oldY, float oldZ,
                 float *newX, float *newY, float *newZ,
                 float psi, float theta, float phi )
{
    float r11, r21, r31, r12, r22, r32, r13, r23, r33;
    r11 = cos(psi)*cos(phi)  -  sin(psi)*cos(theta)*sin(phi);
    r21 = sin(psi)*cos(phi)  +  cos(psi)*cos(theta)*sin(phi);
    r31 = sin(theta)*sin(phi);

    r12 = -cos(psi)*sin(phi)  -  sin(psi)*cos(theta)*cos(phi);
    r22 = -sin(psi)*sin(phi)  +  cos(psi)*cos(theta)*cos(phi);
    r32 = sin(theta)*cos(phi);

    r13 = sin(psi)*sin(theta);
    r23 = -cos(psi)*sin(theta);
    r33 = cos(theta);

    *newX = r11 * oldX + r12 * oldY + r13 * oldZ;
    *newY = r21 * oldX + r22 * oldY + r23 * oldZ;
    *newZ = r31 * oldX + r32 * oldY + r33 * oldZ;

}

void translateAtom (float oldX, float oldY, float oldZ,
                 float *newX, float *newY, float *newZ,
                 float xOffset, float yOffset, float zOffset )
{

    *newX = oldX + xOffset;
    *newY = oldY + yOffset;
    *newZ = oldZ + zOffset;
}

