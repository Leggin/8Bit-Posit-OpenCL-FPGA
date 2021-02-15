import softposit as sp
import numpy as np
import csv
import subprocess
import parmap


def csv_list(csv_file_name):
    with open(csv_file_name) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        for row in csv_reader:
            yield row


def run_device_c(mode=0, a="0", b="0"):
    command = ['./device.out', str(mode), a, b]
    output = subprocess.check_output(command, encoding='UTF-8')
    return output


def addition_comparison(row_b, row_a):
    c_library_result = run_device_c(mode=0, a=row_a[0], b=row_b[0])
    python_posit_result = sp.posit8(float(row_a[0])) + sp.posit8(float(row_b[0]))
    if float(c_library_result) == float(python_posit_result):
        return 1
    else:
        print(f'{c_library_result} != {python_posit_result}')
        return 0


def division_comparison(row_b, row_a):
    c_library_result = run_device_c(mode=3, a=row_a[0], b=row_b[0])
    python_posit_result = sp.posit8(float(row_a[0])) / sp.posit8(float(row_b[0]))
    if float(c_library_result) == float(python_posit_result):
        return 1
    else:
        print(f'{c_library_result} != {python_posit_result}')  # 0.00 != NaR -> no standardized nan in c-code print possible
        return 0


def multiplication_comparison(row_b, row_a):
    c_library_result = run_device_c(mode=2, a=row_a[0], b=row_b[0])
    python_posit_result = sp.posit8(float(row_a[0])) * sp.posit8(float(row_b[0]))
    if float(c_library_result) == float(python_posit_result):
        return 1
    else:
        print(f'{c_library_result} != {python_posit_result}')
        return 0


def get_comparator(mode):
    if mode == 'addition':
        return addition_comparison
    elif mode == 'multiplication':
        return multiplication_comparison
    elif mode == 'division':
        return division_comparison


def unit_test(mode='addition'):
    comparison = get_comparator(mode)
    successcount = errcount = i = 0
    possible_values_a = iter(csv_list('8_bit.csv'))

    for row_a in possible_values_a:
        possible_values_b = iter(csv_list('8_bit.csv'))

        # skip rows that were checked in previous iterations
        for _ in range(i):
            next(possible_values_b)

        results = parmap.map(comparison, possible_values_b, row_a)
        successcount += np.count_nonzero(np.array(results) == 1)
        errcount += np.count_nonzero(np.array(results) == 0)
        i += 1

    print(f'we have {successcount} successes, and {errcount} errors in {mode} mode')


print(f'conversion from posit to real (10101111): {run_device_c(-1, a="10101111")}')
unit_test(mode='division')
unit_test(mode='addition')
unit_test(mode='multiplication')
