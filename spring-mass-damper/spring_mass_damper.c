#include <stdio.h>

void BuildSystem(int *dim, float sys[][102]);
void PrintSystem(int dim, float sys[][102]);
void SolveGaussJordan(int dim, float sys[][102]);

int main(void)
{
    int c, dimension;
    float system[101][102];

    BuildSystem(&dimension, system);
    printf("\n\nAssembled system:\n\n");
    PrintSystem(dimension, system);
    SolveGaussJordan(dimension, system);
    printf("\n\nSolution:\n");
    for (c = 1; c <= dimension; c++)
        printf("\n U%d = %f\n", c, system[c][dimension + 1]);

    return 0;
}

void BuildSystem(int *dim, float sys[][102])
{
    int i, j;

    printf("\n|| FEM solution: spring-mass-damper system ||");
    printf("\n\nEnter the number of nodes: ");
    scanf("%d", dim);

    for (i = 0; i <= *dim; i++)
        for (j = 0; j <= *dim; j++)
            sys[i][j] = 0;

    /* Master element coefficients derived from Rayleigh-Ritz:
       x'' + 0.8x' + x = 0,  h = 6.4 / (n-1) */
    sys[1][1] = 1.5796;
    sys[1][2] = 1.6239;
    sys[2][1] = 0.8239;
    sys[2][2] = 2.3796;

    /* Assemble global stiffness matrix */
    for (i = 1; i < *dim; i++)
        for (j = 1; j < *dim; j++) {
            if (j == i)
                sys[i+1][j+1] = sys[1][1] + sys[i+1][j+1];
            break;
        }

    for (i = 1; i < *dim; i++)
        for (j = 1; j < *dim; j++) {
            if (j == i) {
                sys[i+1][i+1] = sys[2][2];
                if (i < *dim - 1) {
                    sys[i+1][j+2] = sys[1][2];
                    sys[i+2][j+1] = sys[2][1];
                }
            }
        }

    sys[*dim][*dim] = sys[2][2] - sys[1][1];

    /* Apply boundary conditions: x(0) = 1, dx/dt at last node = 0 */
    sys[1][1] = 1;
    sys[2][1] = 0;
    for (i = 0; i <= *dim; i++) {
        if (i == 0)
            sys[i+1][*dim+1] = -1.5756;
        else if (i == 1)
            sys[i+1][*dim+1] = -0.8239;
        else
            sys[i+1][*dim+1] = 0;
    }
}

void PrintSystem(int dim, float sys[][102])
{
    int i, j;
    for (i = 1; i <= dim; i++) {
        for (j = 1; j <= dim + 1; j++) {
            printf("%.4f\t", sys[i][j]);
            if (j == dim) printf("   |");
        }
        printf("\n");
    }
}

void SolveGaussJordan(int dim, float sys[][102])
{
    int nonzero, col, c1, c2, row;
    float pivot, temp;

    for (col = 1; col <= dim; col++) {
        nonzero = 0;
        row = col;
        while (nonzero == 0) {
            if (sys[row][col] > 0.0000001 || sys[row][col] < -0.0000001)
                nonzero = 1;
            else
                row++;
        }
        pivot = sys[row][col];
        for (c1 = 1; c1 <= dim + 1; c1++) {
            temp = sys[row][c1];
            sys[row][c1] = sys[col][c1];
            sys[col][c1] = temp / pivot;
        }
        for (c2 = col + 1; c2 <= dim; c2++) {
            temp = sys[c2][col];
            for (c1 = col; c1 <= dim + 1; c1++)
                sys[c2][c1] = sys[c2][c1] - temp * sys[col][c1];
        }
    }

    for (col = dim; col >= 1; col--)
        for (c1 = col - 1; c1 >= 1; c1--) {
            sys[c1][dim+1] = sys[c1][dim+1] - sys[c1][col] * sys[col][dim+1];
            sys[c1][col] = 0;
        }
}
