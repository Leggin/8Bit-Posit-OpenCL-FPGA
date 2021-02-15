#include <math.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdlib.h>

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

typedef unsigned char posit8; // posit with 8 bits and exponent size 0

typedef struct posit_values
{
    bool sign;
    int k;
    unsigned char exp;
    unsigned char frac;
    unsigned char fracLength;

    bool inf;
    bool zero;
} posit_values;


int clz(char c)
{
    unsigned int n = 0;
    if (c == 0)
        return 8;

    if ((c & 0xF0) == 0)
    {
        n = n + 4;
        c = c << 4;
    }
    if ((c & 0xC0) == 0)
    {
        n = n + 2;
        c = c << 2;
    }
    if ((c & 0x80) == 0)
    {
        n = n + 1;
    }
    return n;
}

void printBits(size_t const size, void const *const ptr)
{
    unsigned char *b = (unsigned char *)ptr;
    unsigned char byte;
    int i, j;

    for (i = size - 1; i >= 0; i--)
    {
        for (j = 7; j >= 0; j--)
        {
            byte = (b[i] >> j) & 1;
            printf("%u", byte);
        }
    }
    puts("");
}

posit8 twosComplement(posit8 input)
{
    input = ~input;
    input = input + 1;
    return input;
}

void extractPositValues(posit8 input, posit_values *output)
{
    output->exp = 0x0; //allways zero by definition

    if (input == 0x0)
    {
        output->zero = true;
        output->inf = false;
        output->sign = false;
        output->k = 0;
        output->frac = 0x0;
        output->fracLength = 0;
    }
    else if (input == 0x80)
    {
        output->zero = false;
        output->inf = true;
        output->sign = true;
        output->k = 0;
        output->frac = 0x0;
        output->fracLength = 0;
    }
    else
    {
        output->sign = input >> 7;
        // if negative make two's complement
        if (output->sign)
        {
            input = twosComplement(input);
        }

        // regime bits consisting of zeros
        int regime_bits = clz(input << 1);
        output->k = -1 * regime_bits;

        if (regime_bits == 0)
        {
            // regime bits consisting of ones
            unsigned char temp = ~(input << 1); // make bitshift before inverting for all-1-case
            regime_bits = clz(temp);
            output->k = regime_bits - 1;
        }

        // flush sign and regime
        if (regime_bits < 7)
        {
            output->fracLength = 8 - (regime_bits + 2); //mask all the regime bits (regime_bits + 1) + the sign bit (1)
        }
        else
        {
            output->fracLength = 0;
        }
        if (output->fracLength > 0)
        {
            unsigned char mask = ((1 << output->fracLength) - 1);
            output->frac = input & mask;
        }
        else
        {
            output->frac = 0x0;
        }
    }
}

void posit8ToDouble(posit8 posit, double *output)
{
    posit_values pval;
    extractPositValues(posit, &pval);
    double fraction_max = pow(2, pval.fracLength);
    double fracVal = 1 + ((double)pval.frac / fraction_max);
    if (posit == 0x0)
    {
        *output = 0;
    }
    else
    {
        *output = pow(-1, pval.sign) * pow(2, pval.k) * fracVal;
    }
}

char convertFraction(double doubleFraction, unsigned int fracLength)
{
    char positFraction = 0x0;
    unsigned int fracCounter = fracLength;
    double temp = 1;
    while (fracCounter > 0)
    {
        temp /= 2;

        if (doubleFraction >= temp)
        {
            // shift in 1
            positFraction <<= 1;
            positFraction += 1;
            doubleFraction -= temp;
        }
        else
        {
            // shift in 0
            positFraction <<= 1;
        }
        fracCounter--;
    }
    return positFraction;
}

void doubleToPosit8(double doubleInput, posit8 *out)
{
    char output = 0x0;
    unsigned int sign;
    unsigned int regimeLength;
    double tempDoubleInput;
    unsigned int regimeIsOnes;
    // extract sign
    (doubleInput >= 0) ? (sign = 0) : (sign = 1);

    if (doubleInput == 0)
    {
        // check for zero
        output = 0;
    }
    else if (isinf(doubleInput) || isnan(doubleInput))
    {
        // check for infinity
        output = 0x80;
    }
    else if (doubleInput >= 64)
    {
        // check for maxpos
        output = 0x7F;
    }
    else if (doubleInput <= -64)
    {
        // check for minpos
        output = 0x81;
    }
    else if (doubleInput <= 0.015625 && !sign)
    {
        // check for +minpos
        output = 0x1;
    }
    else if (doubleInput >= -0.015625 && sign)
    {
        // check for -minpos
        output = 0xFF;
    }
    else
    {
        if (sign)
        {
            // Make negative numbers positive for easier computation
            doubleInput = -doubleInput;
        }

        tempDoubleInput = doubleInput;

        // check if regime is composed of ones
        if (doubleInput > 1 || doubleInput < -1)
        {
            regimeIsOnes = 1;
            regimeLength = 1; // because k = m - 1 we need to add one
            while (tempDoubleInput >= 2)
            {
                tempDoubleInput /= 2;
                regimeLength++;
            }
        }
        // regime is composed of zeros
        else
        {
            regimeIsOnes = 0;
            regimeLength = 0; // because k = m
            while (tempDoubleInput < 1)
            {
                tempDoubleInput *= 2;
                regimeLength++;
            }
        }

        double doubleFraction = tempDoubleInput - 1; // remove hidden bit (1.0101010)
        unsigned int fracLength = 6 - regimeLength;
        char positFraction = convertFraction(doubleFraction, fracLength);

        output = 0x0;

        for (int i = 0; i < regimeLength; i++)
        {
            output = (output << 1) + regimeIsOnes;
        }
        output = (output << 1) + !regimeIsOnes;
        output = (output << fracLength);
        output |= positFraction;

        if (sign)
        {
            output = twosComplement(output);
        }
    }
    *out = output;
}

char kToRegime(int k)
{
    if (k >= 6)
    {
        // regime is only 1s
        return (1 << 7) - 1;
    }
    if (k <= -7)
    {
        // regime is only 0s
        return 0x0;
    }
    if (k >= 0)
    {
        // k is positive, regime consists of 1s ending in 0
        return (1 << (k + 2)) - 2;
    }
    // k is negative and regime ends therefore with 1
    return 0x1;
}

int regimeLengthFromK(int k, int positSize)
{
    int result = 0;

    if (k >= 0)
    {
        result = k + 2;
    }
    else if (k < 0)
    {
        result = -k + 1;
    }

    if (result > positSize - 1)
    {
        result = positSize - 1;
    }
    return result;
}

void addHiddenBitToFraction(posit_values *a)
{
    char hiddenBit = 0x1;
    hiddenBit <<= a->fracLength;
    a->frac |= hiddenBit;
}

void multPosit8(posit8 a, posit8 b, posit8 *result)
{
    *result = 0x0;
    if (a == 0x0)
    {
        *result = 0x0;
    }
    else if (b == 0x0)
    {
        *result = 0x0;
    }
    else if (a == 0x80 || b == 0x80)
    {
        *result = 0x80;
    }
    else
    {
        // extracting posit values
        posit_values valuesA;
        extractPositValues(a, &valuesA);
        posit_values valuesB;
        extractPositValues(b, &valuesB);
        bool sign = valuesA.sign ^ valuesB.sign;

        // adding scales of a and b
        int newK = valuesA.k + valuesB.k;

        // adding hidden bit to fraction
        addHiddenBitToFraction(&valuesA);
        addHiddenBitToFraction(&valuesB);

        // keeping track of shift value for position of hidden bit
        int fracLength;

        // shifting the mantissa of the smaller factor so both are of equal scaling
        if (valuesA.fracLength > valuesB.fracLength)
        {
            valuesB.frac <<= valuesA.fracLength - valuesB.fracLength;
            //fracLength = 2 * valuesA.fracLength;
            fracLength = valuesA.fracLength << 1;
        }
        else
        {
            valuesA.frac <<= valuesB.fracLength - valuesA.fracLength;
            //fracLength = 2 * valuesB.fracLength;
            fracLength = valuesB.fracLength << 1;
        }

        // calculating product of aligned mantissas
        unsigned int intResult = valuesA.frac * valuesB.frac;

        {
            if (intResult >> fracLength + 1 > 0)
            {
                newK++;
                fracLength++;
            }

            if (fracLength > 0)
            {
                // mask hidden bit
                intResult &= 0xffffffff >> 32 - fracLength;
            }
            else
            {
                intResult = 0x0;
            }
        }

        int newRegimeLength = regimeLengthFromK(newK, 8);
        int resultFracLength = 8 - newRegimeLength - 1;

        unsigned int carryCondition = intResult & 0xffffffff >> 32 - (fracLength - resultFracLength);
        bool carryBit = 0;

        // check if rest of mantissa that is not used in new posit is bigger than (0.)100
        // if so carry bit is one and will be added in the end
        if ((fracLength - resultFracLength) > 1 && carryCondition > (0x1 << (fracLength - resultFracLength) - 1))
        {
            carryBit = 1;
        }

        // determine new regime length
        newRegimeLength = regimeLengthFromK(newK, 8);

        // transform k value to regime bits
        char newRegime = kToRegime(newK);

        // calculate length of mantissa bits
        resultFracLength = 8 - newRegimeLength - 1;
        int fracShift = fracLength - resultFracLength;
        unsigned char newFrac;
        if (fracShift > 0)
        {
            newFrac = intResult >> fracShift;
        }
        else
        {
            newFrac = intResult << -fracShift;
        }

        *result |= newRegime << resultFracLength;

        if (resultFracLength > 0)
        {
            *result |= newFrac;
        }

        if (*result != 0x7f && (carryBit || ((*result & 0x1) == 1 && carryCondition == (0x1 << (fracLength - resultFracLength) - 1))))
        {
            *result += 1;
        }

        // per definition posits never drop to zero
        if (*result == 0x0)
        {
            *result = 0x1;
        }
        if (sign)
        {
            *result = twosComplement(*result);
        }
    }
}

void divPosit8(posit8 a, posit8 b, posit8 *result)
{
    *result = 0x0;
    if (a == 0x0)
    {
        *result = 0x0;
    }
    else if (b == 0x0)
    {
        *result = NAN;
    }
    else if (a == 0x80 || b == 0x80)
    {
        *result = 0x80;
    }
    else
    {
        // extracting posit values
        posit_values valuesA;
        extractPositValues(a, &valuesA);
        posit_values valuesB;
        extractPositValues(b, &valuesB);
        bool sign = valuesA.sign ^ valuesB.sign;

        // subtracting scales of a and b
        int newK = valuesA.k - valuesB.k;

        // adding hidden bit to fraction
        addHiddenBitToFraction(&valuesA);
        addHiddenBitToFraction(&valuesB);

        valuesB.frac <<= 7 - valuesB.fracLength;
        valuesA.frac <<= 7 - valuesA.fracLength;

        unsigned int dividend = valuesA.frac << 7;

        unsigned int quot = dividend / valuesB.frac;
        unsigned int rem = dividend % valuesB.frac;

        if (quot != 0)
        {
            bool rcarry = quot >> 7; // this is the hidden bit (7th bit) , extreme right bit is bit 0
            if (!rcarry)
            {
                newK--;
                quot <<= 1;
            }
        }

        // determine new regime length
        unsigned int newRegimeLength = regimeLengthFromK(newK, 8);

        // transform k value to regime bits
        char newRegime = kToRegime(newK);

        int resultFracLength = 8 - newRegimeLength - 1;

        newRegime <<= resultFracLength;

        *result = newRegime;

        int scale;
        if (newK < 0)
        {
            scale = -newK;
        }
        else
        {
            scale = newK + 1;
        }

        quot &= 0x7F;
        unsigned char newFrac = quot;
        newFrac >>= scale + 1;

        *result |= newFrac;

        bool guard = (bool)(0x1 & quot >> scale);
        if (guard)
        {
            bool roundSticky = (((1 << scale) - 1) & quot) ? 1 : 0;
            if (rem > 0 || roundSticky || ((*result & 0x1) == 1 && !roundSticky))
            {
                *result += 1;
            }
        }

        if (*result == 0x0)
        {
            *result = 0x1;
        }

        if (sign)
        {
            *result = twosComplement(*result);
        }
    }
}

void addPosit8(posit8 a, posit8 b, posit8 *result)
{
    *result = 0x0;
    if (a == 0x0)
    {
        *result = b;
    }
    else if (b == 0x0)
    {
        *result = a;
    }
    else if (a == 0x80 || b == 0x80)
    {
        *result = 0x80;
    }
    else if (twosComplement(a) == b)
    {
        *result = 0x0;
    }
    else
    {

        posit_values valuesA;
        extractPositValues(a, &valuesA);
        posit_values valuesB;
        extractPositValues(b, &valuesB);

        bool newSign = 0;

        posit8 compA = a;
        if (valuesA.sign == 1)
        {
            compA = twosComplement(a);
        }

        posit8 compB = b;
        if (valuesB.sign == 1)
        {
            compB = twosComplement(b);
        }
        newSign = (compA > compB) ? valuesA.sign : valuesB.sign;

        // compute scale factor
        int newK = valuesA.k;
        addHiddenBitToFraction(&valuesA);
        addHiddenBitToFraction(&valuesB);

        if (valuesA.k > valuesB.k)
        {
            valuesB.fracLength += valuesA.k - valuesB.k;
            newK = valuesA.k;
        }
        else if (valuesA.k < valuesB.k)
        {
            valuesA.fracLength += valuesB.k - valuesA.k;
            newK = valuesB.k;
        }

        int fracLength;
        unsigned int fracA, fracB;
        if (valuesA.fracLength > valuesB.fracLength)
        {
            fracB = valuesB.frac << valuesA.fracLength - valuesB.fracLength;
            fracA = valuesA.frac;
            fracLength = valuesA.fracLength;
        }
        else
        {
            fracA = valuesA.frac << valuesB.fracLength - valuesA.fracLength;
            fracB = valuesB.frac;
            fracLength = valuesB.fracLength;
        }
        unsigned int tempFrac;

        if (valuesA.sign == valuesB.sign)
        {
            tempFrac = fracA + fracB;
        }
        else
        {
            if (fracA > fracB)
            {
                tempFrac = fracA - fracB;
            }

            else
            {
                tempFrac = fracB - fracA;
            }
        }

        {
            if (tempFrac >> fracLength + 1 > 0)
            {
                newK++;
                fracLength++;
            }
            else
            {
                while (tempFrac >> fracLength == 0 && fracLength > 0)
                {
                    newK--;
                    fracLength--;
                }
            }

            if (fracLength > 0)
            {
                // mask hidden bit
                tempFrac &= 0xffffffff >> 32 - fracLength;
            }
            else
            {
                tempFrac = 0x0;
            }
        }

        int newRegimeLength = regimeLengthFromK(newK, 8);
        int resultFracLength = 8 - newRegimeLength - 1;
        int fractionShiftValue = fracLength - resultFracLength;

        unsigned int cond = tempFrac & (0xffffffff >> 32 - (fracLength - resultFracLength));

        bool carryBit = 0;
        if ((fracLength - resultFracLength) > 0 && cond > (0x1 << (fracLength - resultFracLength) - 1))
        {
            carryBit = 1;
        }
        char newRegime = kToRegime(newK);
        unsigned char newFrac = tempFrac;

        *result = 0x0 << 7;
        *result |= newRegime << resultFracLength;

        if (resultFracLength > 0)
        {
            if (fractionShiftValue > 0)
            {
                *result |= (tempFrac >> fractionShiftValue);
            }
            else
            {
                *result |= (tempFrac << -fractionShiftValue);
            }
        }

        if (carryBit || (((*result & 0x1) == 1 && cond == (0x1 << (fracLength - resultFracLength) - 1)) && *result != 0x7f))
        {
            *result += 1;
        }

        if (newSign == 1)
        {
            *result = twosComplement(*result);
        }
    }
}

void subPosit8(posit8 a, posit8 b, posit8 *result)
{
    *result = 0x0;
    posit8 negative = twosComplement(b);
    addPosit8(a, negative, result);
}

void sigmoidPosit8(posit8 *x)
{
    *x = *x ^ (1 << 7);
    *x >>= 2;
}


void main(int argc, char *argv[])
{
    double doubleInputA;
    double doubleInputB;
    int mode;
    char a = 0x0;
    char b = 0x0;

    sscanf(argv[1], "%d", &mode);


    if (argc > 2)
    {
        sscanf(argv[2], "%lf", &doubleInputA);
        doubleToPosit8(doubleInputA, &a);
    }

    if (argc > 3)
    {
        sscanf(argv[3], "%lf", &doubleInputB);
        doubleToPosit8(doubleInputB, &b);
    }

    char out;
    if (mode == -2)
    {
        // from double to posit binary string
        printBits(sizeof(a), &a);
        return;
    }
    else if (mode == -1)
    {
        // from binary string to posit

        char input[8];
        sscanf(argv[2], "%s", input);
        for (int i = 0; i < 8; i++)
        {
            if (input[i] == '1')
            {
                out |= 1;
            }
            if (i < 7)
            {
                out <<= 1;
            }
        }
    }
    else if (mode == 0)
    {
        addPosit8(a, b, &out);
    }
    else if (mode == 1)
    {
        subPosit8(a, b, &out);
    }
    else if (mode == 2)
    {
        multPosit8(a, b, &out);
    }
    else if (mode == 3)
    {
        divPosit8(a, b, &out);
    }
    else if (mode == 4)
    {
        sigmoidPosit8(&a);
        out = a;
    }

    double output;
    posit8ToDouble(out, &output);
    printf("%lf", output);
}
