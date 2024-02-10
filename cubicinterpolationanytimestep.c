#include <stdio.h>
#include <stdlib.h>
#include <math.h>


typedef struct PointCoords {
    double xval;
    double yval;
} DataPoint;

void inputCoordinatesManual(DataPoint *points, int n) { //This function entails manual input of coords
    for (int i = 0; i < n; ++i) {
        printf("\nEnter coordinates for Point %d:\n", i + 1);
        printf("X Value: ");
        scanf("%lf", &points[i].xval);
        printf("Y Value: ");
        scanf("%lf", &points[i].yval);
    }
}

void inputCoordinatesFromCSV(DataPoint *points, int n, const char *filename) { //This function entails input from csv
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        fprintf(stderr, "Error opening file: %s\n", filename);
        return;
    }

    for (int i = 0; i < n; ++i) {
        if (fscanf(file, "%lf,%lf", &points[i].xval, &points[i].yval) != 2) {
            fprintf(stderr, "Error reading data from CSV file.\n");
            break;
        }
    }

    fclose(file);
}

void calculate_sdys(int m, double sdy[], DataPoint points[]) { //This function calculates second derivates of y
    int n = m - 2;

    double a[n], b[n], c[n], d[n], cdash[n], ddash[n];

    printf("Enter second derivative at first point: ");
    scanf("%lf", &sdy[0]);
    printf("Enter second derivative at last point: ");
    scanf("%lf", &sdy[m - 1]);

    for (int i = 0; i < n; i++) {
        a[i] = (points[i + 1].xval - points[i].xval) / 6;
        b[i] = (points[i + 2].xval - points[i].xval) / 3;
        c[i] = (points[i + 2].xval - points[i + 1].xval) / 6;
        d[i] = ((points[i + 2].yval - points[i + 1].yval) / (points[i + 2].xval - points[i + 1].xval)) -
               ((points[i + 1].yval - points[i].yval) / (points[i + 1].xval - points[i].xval));

        if (i == 0) {
            cdash[i] = c[i] / b[i];
            ddash[i] = (d[i] - a[i] * sdy[0]) / b[i];
        } else if (i != 0 && i != n - 1) {
            cdash[i] = c[i] / (b[i] - cdash[i - 1] * a[i]);
            ddash[i] = (d[i] - ddash[i - 1] * a[i]) / (b[i] - cdash[i - 1] * a[i]);
        } else if (i == n - 1) {
            cdash[i] = c[i] / (b[i] - cdash[i - 1] * a[i]);
            ddash[i] = (d[i] - c[i] * sdy[m - 1] - ddash[i - 1] * a[i]) / (b[i] - cdash[i - 1] * a[i]);
        }
    }

    for (int i = m - 2; i > 0; i--) {
        if (i == m - 2) {
            sdy[m - 2] = ddash[n - 1];
        } else {
            sdy[i] = ddash[n - m + 1 + i] - cdash[n - m + 1 + i] * sdy[i + 1];
        }
    }

    printf("\nInterpolated Second Derivatives:\n");
    for (int i = 0; i < m; i++) {
        printf("Sdy of point %d is: %lf\n", i + 1, sdy[i]);
    }
}

void interpolateCubicSpline(int m, double sdy[], DataPoint points[], DataPoint vels_all[], DataPoint acc_all[], double timestep) {
    // This function will calculate values of y through interpolation and also plot the curve

    FILE *gnuplotPipe = popen("gnuplot -persistent", "w");
    if (!gnuplotPipe) {
        fprintf(stderr, "Error opening gnuplot pipe.\n");
        exit(EXIT_FAILURE);
    }

    fprintf(gnuplotPipe, "plot '-' title 'Input Points' with points pointtype 7 pointsize 1.5, '-' with lines title 'Interpolation', '-' title 'Interpolated Points' with points\n");
    for (int i = 0; i < m; i++) {
        fprintf(gnuplotPipe, "%lf %lf\n", points[i].xval, points[i].yval);
    }
    fprintf(gnuplotPipe, "e\n");

    int temp=0;
    for (int i = 0; i < m - 1; i++) {
        for (double t = 0.0; t <= 1.0; t += timestep) {
            // Calculate the x and y coordinates of the interpolated point
            double A,B,C,D;  // Coefficents of y(j), y(j+1), sdy(j), sdy(j+1)
            double x = points[i].xval + t;
            vels_all[temp].xval = x;
            acc_all[temp].xval = x;

            A = (points[i+1].xval - x) / (points[i+1].xval - points[i].xval);
            B = (x - points[i].xval) / (points[i+1].xval - points[i].xval) ;
            C = 0.166667 * (A*A*A - A) * (points[i+1].xval - points[i].xval)* (points[i+1].xval - points[i].xval);
            D = 0.166667 * (B*B*B - B) * (points[i+1].xval - points[i].xval) * (points[i+1].xval - points[i].xval);

            double y = A * points[i].yval + B * points[i+1].yval + C * sdy[i] + D * sdy[i+1];
            vels_all[temp].yval = y;
            
            if (temp != 0) {
                acc_all[temp].yval = (vels_all[temp].yval - vels_all[temp-1].yval) / (vels_all[temp].xval - vels_all[temp-1].xval);   // acc is first der of vel. Therefore 
            }
            else {acc_all[temp].yval = 0;}
            fprintf(gnuplotPipe, "%lf %lf\n", x, y);
            temp += 1;
        }
    }
    fprintf(gnuplotPipe, "e\n");

    for (int i=0; i<(m-1)*((1/timestep) +1);i++){
        fprintf(gnuplotPipe, "%lf %lf\n", vels_all[i].xval, vels_all[i].yval);
    }

    fclose(gnuplotPipe);
}

double lagrange_term(int n,DataPoint points[], double x_val) {
    double result = 0.0;
    for (int i = 0; i < n; i++) {
        double term = points[i].yval;
        for (int j = 0; j < n; j++) {
            if (j != i) {
                term *= (x_val - points[j].xval) / (points[i].xval - points[j].xval);
            }
        }
        result += term;
    }
    return result;
}

void interpolateLagrange(int m, DataPoint points[], DataPoint vels_all[], DataPoint acc_all[], double timestep){

    FILE *gnuplotPipe = popen("gnuplot -persistent", "w");
    if (!gnuplotPipe) {
        fprintf(stderr, "Error opening gnuplot pipe.\n");
        exit(EXIT_FAILURE);
    }

    fprintf(gnuplotPipe, "plot '-' title 'Input Points' with points pointtype 7 pointsize 1.5, '-' with lines title 'Interpolation', '-' title 'Interpolated Points' with points\n");
    for (int i = 0; i < m; i++) {
        fprintf(gnuplotPipe, "%lf %lf\n", points[i].xval, points[i].yval);
    }
    fprintf(gnuplotPipe, "e\n");


    int temp = 0;
    for (int i = 0; i < m - 1; i++) {
        for (double t = 0.0; t <= 1.0; t += timestep){
            double x = points[i].xval + t;
            double y = lagrange_term(m, points, x);

            vels_all[temp].xval = x;
            acc_all[temp].xval = x;
            vels_all[temp].yval = y;
            if (temp != 0) {
                acc_all[temp].yval = (vels_all[temp].yval - vels_all[temp-1].yval) / (vels_all[temp].xval - vels_all[temp-1].xval) ;   // acc is first der of vel. Therefore (x2-x1)/(y2-y1)
            }
            else {acc_all[temp].yval = 0;}
            fprintf(gnuplotPipe, "%lf %lf\n", x, acc_all[temp].yval);
            temp += 1;
       }
    }   
    fprintf(gnuplotPipe, "e\n");

    for (int i=0; i<(m-1)*((1/timestep) +1);i++){
       fprintf(gnuplotPipe, "%lf %lf\n", vels_all[i].xval, vels_all[i].yval);
    }

    fclose(gnuplotPipe);
}

void inputCoordinates(int n, DataPoint *points) {
    int choice;

    printf("Type 1 to enter coordinates manually or Type 2 to import coordinates from a CSV file: ");
    scanf("%d", &choice);

    do {
        if (choice == 1) {
            inputCoordinatesManual(points, n);
        } else if (choice == 2) {
            char filename[100];
            printf("Enter the path of your file: ");
            scanf("%s", filename);
            inputCoordinatesFromCSV(points, n, filename);
        } else {
            printf("Wrong Value entered, Please Try again\n");
        }

    } while (choice != 1 && choice != 2);
}

int main() {
    int m;
    double timestep;
    printf("Enter the number of points: ");
    scanf("%d", &m);
    printf("Enter the timestep: ");
    scanf("%lf", &timestep);
    int numOfInterpolatedPoints;
    numOfInterpolatedPoints = (m-1)*((1/timestep) +1);

    DataPoint points[m];
    DataPoint vels_all[numOfInterpolatedPoints];
    DataPoint acc_all[numOfInterpolatedPoints];
    double sdy_nodes[m]; // Array of second derivatives at each point
    double fdy[m]; // Array of first derivatives at each point

    inputCoordinates(m, points);

    int choice;

    printf("Type 1 for Lagrange Interpolation or Type 2 to for Cubic Interpolation: ");
    scanf("%d", &choice);

    do {
        if (choice == 1) {
            interpolateLagrange(m, points, vels_all, acc_all, timestep);
        } else if (choice == 2) {
            calculate_sdys(m, sdy_nodes, points);
            interpolateCubicSpline(m, sdy_nodes, points, vels_all, acc_all,timestep);
        } else {
            printf("Wrong Value entered, Please Try again\n");
        }

    } while (choice != 1 && choice != 2);

    char filename[100];
    printf("Enter the path of your output file: ");
    scanf("%s", filename);
    
    FILE *csvFile = fopen(filename, "w");
    if (csvFile == NULL) {
        fprintf(stderr, "Error opening file for writing.\n");
        return 1;
    }

    fprintf(csvFile, "Timestep,Second Derivative at First Point,Second Derivative at Last Point\n");

    fprintf(csvFile, "%lf,%lf,%lf \n\n", timestep, sdy_nodes[0], sdy_nodes[m-1]);

    fprintf(csvFile, "X-value, Velocity, Acceleration\n");

    for (int i = 0; i < numOfInterpolatedPoints; i++) {
        fprintf(csvFile, "%lf,%lf,%lf\n",
                vels_all[i].xval, vels_all[i].yval,
                acc_all[i].yval);
    }

    printf("File Saved Successfully!!");
    fclose(csvFile);


    return 0;
}

