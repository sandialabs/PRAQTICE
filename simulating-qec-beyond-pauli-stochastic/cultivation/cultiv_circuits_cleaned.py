import stim 


#Gidney cultivation circuits w/ d=3 color code
cultiv_stimc_d3 = stim.Circuit('''
QUBIT_COORDS(0, 0) 0
    QUBIT_COORDS(0, 1) 1
    QUBIT_COORDS(1, 0) 2
    QUBIT_COORDS(1, 1) 3
    QUBIT_COORDS(1, 2) 4
    QUBIT_COORDS(1, 3) 5
    QUBIT_COORDS(2, 0) 6
    QUBIT_COORDS(2, 1) 7
    QUBIT_COORDS(2, 2) 8
    QUBIT_COORDS(2, 3) 9
    QUBIT_COORDS(3, 0) 10
    QUBIT_COORDS(3, 1) 11
    QUBIT_COORDS(3, 2) 12
    QUBIT_COORDS(4, 0) 13
    QUBIT_COORDS(4, 1) 14
    R 7 10 12 3 2 4 14 0 5 13 11 8 6
    TICK
    H 0 2 3 4 5 6 7 8 10 11 12 13 14
    TICK
    CZ 6 10 7 11 8 12
    TICK
    CZ 7 8 10 11
    TICK
    CZ 6 7 11 12
    TICK
    CZ 2 6 4 8 11 14
    TICK
    H 2 4 6 8 11 14
    TICK
    CZ 2 6 4 8 11 14
    TICK
    CZ 3 4 13 14
    TICK
    H 4 13 14
    TICK
    CZ 2 3 13 14
    TICK
    H 2 3
    TICK
    CZ 3 4
    TICK
    CZ 2 3
    TICK
    X 3
    TICK
    CZ 2 3
    TICK
    CZ 3 4
    TICK
    H 2 4
    TICK
    CZ 0 2 4 5
    TICK
    H 0 2 4 5
    TICK
    CZ 0 2 4 5
    TICK
    R 2 11 4 6 14 8 9
    TICK
    H 2 4 5 6 8 9 11 14
    TICK
    CZ 2 6 4 8 5 9 11 14
    TICK
    H 5 6 7 8 9 10 12 14
    TICK
    CZ 0 2 5 9 6 10 7 11 8 12
    TICK
    H 0 5 9
    TICK
    CZ 2 3 6 7 8 9 11 12
    TICK
    CZ 3 4 7 8 10 11 13 14
    TICK
    H 2 3 4 6 7 8 9 10 11 12 13 14
    TICK
    CZ 3 4 7 8 10 11 13 14
    TICK
    CZ 0 2 6 10 7 11 8 12
    TICK
    CZ 2 3 6 7 8 9 11 12
    TICK
    H 2 4 11
    TICK
    CZ 2 6 4 8 11 14
    TICK
    H 2 4 6 8 11 14
    TICK
    M 6 14 8 5 2 11 4
    DETECTOR(2, 0, 0) rec[-7]
    DETECTOR(4, 1, 0) rec[-6]
    DETECTOR(2, 2, 0) rec[-5]
    DETECTOR(1, 3, 0) rec[-4]
    DETECTOR(1, 0, 0) rec[-3]
    DETECTOR(3, 1, 0) rec[-2]
    DETECTOR(1, 2, 0) rec[-1]
    TICK
    R 14 11 6 2 8 1
    TICK
    H 1 2 4 6 8 11 14
    SQRT_Y 0 3 7 9 10 12 13
    TICK
    CZ 0 1 2 3 6 7 8 9 10 11 13 14
    TICK
    H 1 3 7 8 14
    TICK
    CZ 1 3 7 8 11 14
    TICK
    H 3
    TICK
    CZ 3 7 11 12
    TICK
    H 11
    TICK
    CZ 7 11
    TICK
    H 7
    TICK
    M 7
    TICK
    R 7
    TICK
    H 7
    TICK
    CZ 7 11
    TICK
    H 11
    TICK
    CZ 3 7 11 12
    TICK
    H 3
    SQRT_Y_DAG 12
    TICK
    CZ 1 3 7 8 11 14
    TICK
    H 1 3 7 8 14
    TICK
    CZ 0 1 2 3 6 7 8 9 10 11 13 14
    TICK
    H 1 2 6 8 11 14
    SQRT_Y_DAG 0 3 7 9 10 13
    TICK
    M 14 11 6 2 8 1
    OBSERVABLE_INCLUDE(0) rec[-6] rec[-5] rec[-2] rec[-1]
    DETECTOR(2.14286, 1, 1, -1, -9) rec[-9] rec[-8] rec[-7]
    DETECTOR(4, 1, 2) rec[-6]
    DETECTOR(3, 1, 2) rec[-5]
    DETECTOR(2, 1, 2) rec[-4] rec[-7]
    DETECTOR(1, 0, 2) rec[-3]
    DETECTOR(2, 2, 2) rec[-2]
    DETECTOR(0, 1, 2) rec[-1]
    H 1 2 6 8 11 14
    TICK
    MPP !X0*X3*X7*X9*X10*X12*X13 X0*X3*X7*X10 Z0*Z3*Z7*Z10 X3*X7*X9*X12 Z3*Z7*Z9*Z12 X7*X10*X12*X13 Z7*Z10*Z12*Z13
    OBSERVABLE_INCLUDE(0) rec[-7]
    DETECTOR(0.75, 0.25, 4, -1, -9) rec[-14] rec[-6]
    DETECTOR(0.75, 0.25, 5, -1, -9) rec[-5]
    DETECTOR(1.5, 1.375, 4, -1, -9) rec[-14] rec[-4]
    DETECTOR(1.5, 1.375, 5, -1, -9) rec[-3]
    DETECTOR(2.5, 0.875, 4, -1, -9) rec[-14] rec[-2]
    DETECTOR(2.5, 0.875, 5, -1, -9) rec[-1]
''')

injection_stimc_d3 = stim.Circuit('''QUBIT_COORDS(0, 0) 0
    QUBIT_COORDS(1, 0) 1
    QUBIT_COORDS(1, 1) 2
    QUBIT_COORDS(1, 2) 3
    QUBIT_COORDS(1, 3) 4
    QUBIT_COORDS(2, 0) 5
    QUBIT_COORDS(2, 1) 6
    QUBIT_COORDS(2, 2) 7
    QUBIT_COORDS(2, 3) 8
    QUBIT_COORDS(3, 0) 9
    QUBIT_COORDS(3, 1) 10
    QUBIT_COORDS(3, 2) 11
    QUBIT_COORDS(4, 0) 12
    QUBIT_COORDS(4, 1) 13
    R 10 7 5 6 9 11 2 1 3 13 0 4 12
    TICK
    H 0 1 2 3 4 5 6 7 9 10 11 12 13
    TICK
    CZ 5 9 6 10 7 11
    TICK
    CZ 6 7 9 10
    TICK
    CZ 5 6 10 11
    TICK
    CZ 1 5 3 7 10 13
    TICK
    H 1 3 5 7 10 13
    TICK
    CZ 1 5 3 7 10 13
    TICK
    CZ 2 3 12 13
    TICK
    H 3 12 13
    TICK
    CZ 1 2 12 13
    TICK
    H 1 2
    TICK
    CZ 2 3
    TICK
    CZ 1 2
    TICK
    SQRT_X_DAG 2
    TICK
    CZ 1 2
    TICK
    CZ 2 3
    TICK
    H 1 3
    TICK
    CZ 0 1 3 4
    TICK
    H 0 1 3 4
    TICK
    CZ 0 1 3 4
    TICK
    R 5 13 7 8 1 10 3
    TICK
    H 1 3 4 5 7 8 10 13
    TICK
    CZ 1 5 3 7 4 8 10 13
    TICK
    H 4 5 6 7 8 9 11 13
    TICK
    CZ 0 1 4 8 5 9 6 10 7 11
    TICK
    H 0 4 8
    TICK
    CZ 1 2 5 6 7 8 10 11
    TICK
    CZ 2 3 6 7 9 10 12 13
    TICK
    H 1 2 3 5 6 7 8 9 10 11 12 13
    TICK
    CZ 2 3 6 7 9 10 12 13
    TICK
    CZ 0 1 5 9 6 10 7 11
    TICK
    CZ 1 2 5 6 7 8 10 11
    TICK
    H 1 3 10
    TICK
    CZ 1 5 3 7 10 13
    TICK
    H 1 3 5 7 10 13
    TICK
    M 5 13 7 4 1 10 3
    DETECTOR(2, 0, 0) rec[-7]
    DETECTOR(4, 1, 0) rec[-6]
    DETECTOR(2, 2, 0) rec[-5]
    DETECTOR(1, 0, 0) rec[-3]
    DETECTOR(3, 1, 0) rec[-2]
    DETECTOR(1, 2, 0) rec[-1]''')

def create_d3_injection_circ():
    return injection_stimc_d3

def create_d3_cultiv_circ():
    return cultiv_stimc_d3