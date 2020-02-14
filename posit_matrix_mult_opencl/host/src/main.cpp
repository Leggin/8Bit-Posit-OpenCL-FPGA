#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cstring>
#include "CL/opencl.h"
#include "AOCLUtils/aocl_utils.h"

using namespace aocl_utils;

#define STRING_BUFFER_LEN 1024

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

double doubleInput = 1.0;

// OpenCL runtime configuration
static cl_platform_id platform = NULL;
static cl_device_id device = NULL;
static cl_context context = NULL;
static cl_command_queue queue = NULL;
static cl_kernel doubleToPositKernel = NULL;
static cl_kernel positToDoubleKernel = NULL;
static cl_kernel computationKernel = NULL;

static cl_program program = NULL;
static cl_mem output_buf1;
static cl_mem output_buf2;

static cl_mem input_buf;  // num_devices elements
static cl_mem output_buf; // num_devices elements

typedef char _posit8;

unsigned int N = 4; // problem size
unsigned NUM_ITERATIONS = 1;
static scoped_aligned_ptr<_posit8> input;  // num_devices elements
static scoped_aligned_ptr<_posit8> output; // num_devices elements
_posit8 INITIAL_TEMPERATURE = 0x50;        // 01010000: 1,5

// Function prototypes
bool init();
void cleanup();
void init_problem();
void posit8ToDouble(_posit8 p, double *d);
int _clz(char c);
void extractPositValues(_posit8 input, posit_values *output);
void doubleToPosit8(double doubleInput, _posit8 *out);

// Entry point.
int main(int argc, char **argv)
{
    cl_int status;

    if (argc > 1)
    {
        sscanf(argv[1], "%d", &N);
    }

    if (!init())
    {
        return -1;
    }
    init_problem();

    cl_event write_event;

    status = clEnqueueWriteBuffer(queue, input_buf, CL_FALSE, 0, N * N * sizeof(_posit8), input, 0, NULL, &write_event);
    checkError(status, "Failed to transfer input A");

    scoped_aligned_ptr<cl_event> kernel_event;
    kernel_event.reset(NUM_ITERATIONS);
    cl_event finish_event;

    for (unsigned i = 0; i < NUM_ITERATIONS; i++)
    {
        // Set kernel arguments.
        unsigned argi = 0;

        size_t global_work_size[2];
        global_work_size[0] = N;
        global_work_size[1] = N;

        status = clSetKernelArg(computationKernel, argi++, sizeof(cl_mem), &input_buf);
        checkError(status, "Failed to set argument %d", argi - 1);

        status = clSetKernelArg(computationKernel, argi++, sizeof(cl_mem), &input_buf);
        checkError(status, "Failed to set argument %d", argi - 1);

        status = clSetKernelArg(computationKernel, argi++, sizeof(cl_int), (void *)&N);
        checkError(status, "Failed to set argument %d", argi - 1);

        status = clSetKernelArg(computationKernel, argi++, sizeof(cl_mem), &output_buf);
        checkError(status, "Failed to set argument %d", argi - 1);

        status = clEnqueueNDRangeKernel(queue, computationKernel, 2, NULL, global_work_size, NULL, 1, &write_event, &kernel_event[i]);
        checkError(status, "Failed to launch kernel");

        status = clEnqueueReadBuffer(queue, output_buf, CL_FALSE, 0, N * N * sizeof(_posit8), output, 1, &kernel_event[i], &finish_event);

        //status = clEnqueueWriteBuffer(queue, input_buf, CL_FALSE, 0, N * N * sizeof(_posit8), output, 0, NULL, &write_event);
    }

    // Release local events.
    clReleaseEvent(write_event);

    // Wait for all devices to finish.
    clWaitForEvents(1, &finish_event);

    cl_ulong time_ns = 0;
    for (unsigned i = 0; i < NUM_ITERATIONS; i++)
    { 
        time_ns += getStartEndTime(kernel_event[i]);
    } 
    double seconds = double(time_ns) * 1e-9;
    //printf("Kernel time: %0.3f ms\n", double(time_ns) * 1e-6);

    double operations = 2 * pow(N, 3);
    double gflops = (operations / seconds) * 1.0e-9;
    printf("%lf,%lf\n", seconds, gflops);  

    // Release all events.
    for (unsigned i = 0; i < NUM_ITERATIONS; ++i)
    {
        clReleaseEvent(kernel_event[i]);
    }
    clReleaseEvent(finish_event);

  

    cleanup();

    return 0;
}

/////// HELPER FUNCTIONS ///////

double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

bool init()
{
    cl_int status;

    if (!setCwdToExeDir())
    {
        return false;
    }

    // Get the OpenCL platform.
    platform = findPlatform("Intel");
    if (platform == NULL)
    {
        printf("ERROR: Unable to find Intel FPGA OpenCL platform.\n");
        return false;
    }
    // Query the available OpenCL devices.
    scoped_array<cl_device_id> devices;
    cl_uint num_devices;

    devices.reset(getDevices(platform, CL_DEVICE_TYPE_ALL, &num_devices));

    // We'll just use the first device.
    device = devices[0];

    // Create the context.
    context = clCreateContext(NULL, 1, &device, &oclContextCallback, NULL, &status);
    checkError(status, "Failed to create context");

    // Create the command queue.
    queue = clCreateCommandQueue(context, device, CL_QUEUE_PROFILING_ENABLE, &status);
    checkError(status, "Failed to create command queue");

    // Create the program.
    std::string binary_file = getBoardBinaryFile("device", device);
    program = createProgramFromBinary(context, binary_file.c_str(), &device, 1);

    // Build the program that was just created.
    status = clBuildProgram(program, 0, NULL, "", NULL, NULL);
    checkError(status, "Failed to build program");

    // Create the kernel - name passed in here must match kernel name in the
    // original CL file, that was compiled into an AOCX file using the AOC tool

    const char *computationKernelName = "matrix_mult"; // Kernel name, as defined in the CL file
    computationKernel = clCreateKernel(program, computationKernelName, &status);
    checkError(status, "Failed to create computationKernel");

    // Input buffers.
    input_buf = clCreateBuffer(context, CL_MEM_READ_ONLY, N * N * sizeof(_posit8), NULL, &status);
    checkError(status, "Failed to create buffer for input A");

    // Output buffer.
    output_buf = clCreateBuffer(context, CL_MEM_WRITE_ONLY, N * N * sizeof(_posit8), NULL, &status);
    checkError(status, "Failed to create buffer for output");

    return true;
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

// Initialize the data for the problem. Requires num_devices to be known.
void init_problem()
{
    // Generate input vector A and the output consisting
    // of a total of N * N elements.
    input.reset(N * N);
    output.reset(N * N);
    for (unsigned j = 0; j < N * N; ++j)
    {
        _posit8 p;
        doubleToPosit8(fRand(0, 1.0), &input[j]);
        output[j] = 0;
    }
}

// Free the resources allocated during initialization
void cleanup()
{
    if (computationKernel)
    {
        clReleaseKernel(computationKernel);
    }
    if (input_buf)
    {
        clReleaseMemObject(input_buf);
    }
    if (output_buf)
    {
        clReleaseMemObject(output_buf);
    }
    if (program)
    {
        clReleaseProgram(program);
    }
    if (queue)
    {
        clReleaseCommandQueue(queue);
    }
    if (context)
    {
        clReleaseContext(context);
    }
}

_posit8 _twosComplement(_posit8 input)
{
    input = ~input;
    input = input + 1;
    return input;
}

void extractPositValues(_posit8 input, posit_values *output)
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
            input = _twosComplement(input);
        }

        // regime bits consisting of zeros
        int regime_bits = _clz(input << 1);
        output->k = -1 * regime_bits;

        if (regime_bits == 0)
        {
            // regime bits consisting of ones
            unsigned char temp = ~(input << 1); // make bitshift before inverting for all-1-case
            regime_bits = _clz(temp);
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

void posit8ToDouble(_posit8 posit, double *output)
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

int _clz(char c)
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

char convertFraction(double doubleFraction, unsigned int fracLength)
{
    char _positFraction = 0x0;
    unsigned int _fracCounter = fracLength;
    double _temp = 1;
    while (_fracCounter > 0)
    {
        _temp /= 2;

        if (doubleFraction >= _temp)
        {
            // shift in 1
            _positFraction = (_positFraction << 1) + 1;
            doubleFraction -= _temp;
        }
        else
        {
            // shift in 0
            _positFraction = _positFraction << 1;
        }
        _fracCounter--;
    }
    return _positFraction;
}

void doubleToPosit8(double doubleInput, _posit8 *out)
{
    char output = 0x0;
    unsigned int sign;
    unsigned int _regimeLength;
    double _tempDoubleInput;
    unsigned int _regimeIsOnes;
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

        _tempDoubleInput = doubleInput;

        // check if regime is composed of ones
        if (doubleInput > 1 || doubleInput < -1)
        {
            _regimeIsOnes = 1;
            _regimeLength = 1; // because k = m - 1 we need to add one
            while (_tempDoubleInput >= 2)
            {
                _tempDoubleInput /= 2;
                _regimeLength++;
            }
            if (_regimeLength > 6)
            {
                output = 0x7F;
                *out = output;
                return;
            }
        }
        // regime is composed of zeros
        else
        {
            _regimeIsOnes = 0;
            _regimeLength = 0; // because k = m
            while (_tempDoubleInput < 1)
            {
                _tempDoubleInput *= 2;
                _regimeLength++;
            }

            if (_regimeLength > 6)
            {
                output = 0x1;
                *out = output;
                return;
            }
        }

        double _doubleFraction = _tempDoubleInput - 1; // remove hidden bit (1.0101010)
        unsigned int _fracLength = 6 - _regimeLength;
        char _positFraction = convertFraction(_doubleFraction, _fracLength);

        output = 0x0;

        for (int i = 0; i < _regimeLength; i++)
        {
            output = (output << 1) + _regimeIsOnes;
        }
        output = (output << 1) + !_regimeIsOnes;
        output = (output << _fracLength);
        output |= _positFraction;

        if (sign)
        {
            output = _twosComplement(output);
        }
    }
    *out = output;
}
