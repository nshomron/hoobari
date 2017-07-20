#include <stdio.h>
#include <math.h>

double calculatePhred(double* probs);

int main(int argc, char** argv)
{
    double probs[] = {0.5,1.45,3.23521};
    double res = calculatePhred(probs);
    printf("%f", res);
    return 0;
}

double calculatePhred(double* probs)
{
    // exponent
    for(int i=0; i<3; i++){
        probs[i]=exp(probs[i]);
    }

    // Sum exponents
    double sum_exp=0.0;
    for(int i=0; i<3; i++){
        sum_exp+=probs[i];
    }

    // Divide by sum
    double posteriors[3];
    for(int i=0; i<3; i++){
        posteriors[i] = probs[i]/sum_exp;
    }

    // Find maximal posterior value
    double max = 0.0;
    for(int i=0; i<3; i++){
        if(max < posteriors[i]){
            max = posteriors[i];
        }
    }

    double phred = (log10(1-max))*(-10);

    return phred;
}
