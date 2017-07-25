#include <stdio.h>
#include <math.h>

long double calculatePhred(long double* probs);

long double calculatePhred(long double* probs)
{
    // Find minimal value
    long double min = probs[0];
    for(int i=1; i<3; i++){
        if (probs[i] < min){
            min=probs[i];
        }
    }

    // Substract minimal value from all values
    for(int i=0; i<3; i++){
        probs[i]-=min;
    }

    // exponent
    for(int i=0; i<3; i++){
        probs[i]=exp(probs[i]);
    }

    // Sum exponents
    long double sum_exp=0.0;
    for(int i=0; i<3; i++){
        sum_exp+=probs[i];
    }

    // Divide by sum
    long double posteriors[3];
    for(int i=0; i<3; i++){
        posteriors[i] = probs[i]/sum_exp;
    }

    // Find maximal posterior value
    long double max = posteriors[0];
    for(int i=1; i<3; i++){
        if(max < posteriors[i]){
            max = posteriors[i];
        }
    }

    long double phred = (log10(1-max))*(-10);

    return phred;
}
