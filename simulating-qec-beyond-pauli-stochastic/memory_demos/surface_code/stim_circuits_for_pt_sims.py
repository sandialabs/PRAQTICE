import pygsti
import numpy as np
import stim 
from scipy.linalg import hadamard

init_perfect_circuit = {3:stim.Circuit('''R 1 3 5 8 10 12 15 17 19 2 9 11 13 14 16 18 25
    TICK
    H 2 11 16 25
    TICK
    CX 2 3 16 17 11 12 15 14 10 9 19 18
    TICK
    CX 2 1 16 15 11 10 8 14 3 9 12 18
    TICK
    CX 16 10 11 5 25 19 8 9 17 18 12 13
    TICK
    CX 16 8 11 3 25 17 1 9 10 18 5 13
    TICK
    H 2 11 16 25
    TICK
    MR 2 9 11 13 14 16 18 25'''),
    5: stim.Circuit('''R 1 3 5 7 9 12 14 16 18 20 23 25 27 29 31 34 36 38 40 42 45 47 49 51 53 2 6 13 15 17 19 21 22 24 26 28 30 35 37 39 41 43 44 46 48 50 52 59 63
    TICK
    H 2 6 15 19 24 28 37 41 46 50 59 63
    TICK
    CX 2 3 24 25 46 47 15 16 37 38 6 7 28 29 50 51 19 20 41 42 23 22 45 44 14 13 36 35 27 26 49 48 18 17 40 39 31 30 53 52
    TICK
    CX 2 1 24 23 46 45 15 14 37 36 6 5 28 27 50 49 19 18 41 40 12 22 34 44 3 13 25 35 16 26 38 48 7 17 29 39 20 30 42 52
    TICK
    CX 24 14 46 36 15 5 37 27 59 49 28 18 50 40 19 9 41 31 63 53 12 13 34 35 25 26 47 48 16 17 38 39 29 30 51 52 20 21 42 43
    TICK
    CX 24 12 46 34 15 3 37 25 59 47 28 16 50 38 19 7 41 29 63 51 1 13 23 35 14 26 36 48 5 17 27 39 18 30 40 52 9 21 31 43
    TICK
    H 2 6 15 19 24 28 37 41 46 50 59 63
    TICK
    MR 2 6 13 15 17 19 21 22 24 26 28 30 35 37 39 41 43 44 46 48 50 52 59 63'''),
    7: stim.Circuit('''R 1 3 5 7 9 11 13 16 18 20 22 24 26 28 31 33 35 37 39 41 43 46 48 50 52 54 56 58 61 63 65 67 69 71 73 76 78 80 82 84 86 88 91 93 95 97 99 101 103 2 6 10 17 19 21 23 25 27 29 30 32 34 36 38 40 42 47 49 51 53 55 57 59 60 62 64 66 68 70 72 77 79 81 83 85 87 89 90 92 94 96 98 100 102 109 113 117
    TICK
    H 2 6 10 19 23 27 32 36 40 49 53 57 62 66 70 79 83 87 92 96 100 109 113 117
    TICK
    CX 2 3 32 33 62 63 92 93 19 20 49 50 79 80 6 7 36 37 66 67 96 97 23 24 53 54 83 84 10 11 40 41 70 71 100 101 27 28 57 58 87 88 31 30 61 60 91 90 18 17 48 47 78 77 35 34 65 64 95 94 22 21 52 51 82 81 39 38 69 68 99 98 26 25 56 55 86 85 43 42 73 72 103 102
    TICK
    CX 2 1 32 31 62 61 92 91 19 18 49 48 79 78 6 5 36 35 66 65 96 95 23 22 53 52 83 82 10 9 40 39 70 69 100 99 27 26 57 56 87 86 16 30 46 60 76 90 3 17 33 47 63 77 20 34 50 64 80 94 7 21 37 51 67 81 24 38 54 68 84 98 11 25 41 55 71 85 28 42 58 72 88 102
    TICK
    CX 32 18 62 48 92 78 19 5 49 35 79 65 109 95 36 22 66 52 96 82 23 9 53 39 83 69 113 99 40 26 70 56 100 86 27 13 57 43 87 73 117 103 16 17 46 47 76 77 33 34 63 64 93 94 20 21 50 51 80 81 37 38 67 68 97 98 24 25 54 55 84 85 41 42 71 72 101 102 28 29 58 59 88 89
    TICK
    CX 32 16 62 46 92 76 19 3 49 33 79 63 109 93 36 20 66 50 96 80 23 7 53 37 83 67 113 97 40 24 70 54 100 84 27 11 57 41 87 71 117 101 1 17 31 47 61 77 18 34 48 64 78 94 5 21 35 51 65 81 22 38 52 68 82 98 9 25 39 55 69 85 26 42 56 72 86 102 13 29 43 59 73 89
    TICK
    H 2 6 10 19 23 27 32 36 40 49 53 57 62 66 70 79 83 87 92 96 100 109 113 117
    TICK
    MR 2 6 10 17 19 21 23 25 27 29 30 32 34 36 38 40 42 47 49 51 53 55 57 59 60 62 64 66 68 70 72 77 79 81 83 85 87 89 90 92 94 96 98 100 102 109 113 117''')}

circuit_to_add_noise_base = {3:stim.Circuit('''H 2 11 16 25
    TICK
    CX 2 3 16 17 11 12 15 14 10 9 19 18
    TICK
    CX 2 1 16 15 11 10 8 14 3 9 12 18
    TICK
    CX 16 10 11 5 25 19 8 9 17 18 12 13
    TICK
    CX 16 8 11 3 25 17 1 9 10 18 5 13
    TICK
    H 2 11 16 25
    TICK
    MR 2 9 11 13 14 16 18 25
    DETECTOR rec[-8] rec[-16]
    DETECTOR rec[-7] rec[-15]
    DETECTOR rec[-6] rec[-14]
    DETECTOR rec[-5] rec[-13]
    DETECTOR rec[-4] rec[-12]
    DETECTOR rec[-3] rec[-11]
    DETECTOR rec[-2] rec[-10]
    DETECTOR rec[-1] rec[-9]'''), 
                             5: stim.Circuit('''H 2 6 15 19 24 28 37 41 46 50 59 63
    TICK
    CX 2 3 24 25 46 47 15 16 37 38 6 7 28 29 50 51 19 20 41 42 23 22 45 44 14 13 36 35 27 26 49 48 18 17 40 39 31 30 53 52
    TICK
    CX 2 1 24 23 46 45 15 14 37 36 6 5 28 27 50 49 19 18 41 40 12 22 34 44 3 13 25 35 16 26 38 48 7 17 29 39 20 30 42 52
    TICK
    CX 24 14 46 36 15 5 37 27 59 49 28 18 50 40 19 9 41 31 63 53 12 13 34 35 25 26 47 48 16 17 38 39 29 30 51 52 20 21 42 43
    TICK
    CX 24 12 46 34 15 3 37 25 59 47 28 16 50 38 19 7 41 29 63 51 1 13 23 35 14 26 36 48 5 17 27 39 18 30 40 52 9 21 31 43
    TICK
    H 2 6 15 19 24 28 37 41 46 50 59 63
    TICK
    MR 2 6 13 15 17 19 21 22 24 26 28 30 35 37 39 41 43 44 46 48 50 52 59 63
    DETECTOR(2, 0, 0) rec[-24] rec[-48]
    DETECTOR(6, 0, 0) rec[-23] rec[-47]
    DETECTOR(2, 2, 0) rec[-22] rec[-46]
    DETECTOR(4, 2, 0) rec[-21] rec[-45]
    DETECTOR(6, 2, 0) rec[-20] rec[-44]
    DETECTOR(8, 2, 0) rec[-19] rec[-43]
    DETECTOR(10, 2, 0) rec[-18] rec[-42]
    DETECTOR(0, 4, 0) rec[-17] rec[-41]
    DETECTOR(2, 4, 0) rec[-16] rec[-40]
    DETECTOR(4, 4, 0) rec[-15] rec[-39]
    DETECTOR(6, 4, 0) rec[-14] rec[-38]
    DETECTOR(8, 4, 0) rec[-13] rec[-37]
    DETECTOR(2, 6, 0) rec[-12] rec[-36]
    DETECTOR(4, 6, 0) rec[-11] rec[-35]
    DETECTOR(6, 6, 0) rec[-10] rec[-34]
    DETECTOR(8, 6, 0) rec[-9] rec[-33]
    DETECTOR(10, 6, 0) rec[-8] rec[-32]
    DETECTOR(0, 8, 0) rec[-7] rec[-31]
    DETECTOR(2, 8, 0) rec[-6] rec[-30]
    DETECTOR(4, 8, 0) rec[-5] rec[-29]
    DETECTOR(6, 8, 0) rec[-4] rec[-28]
    DETECTOR(8, 8, 0) rec[-3] rec[-27]
    DETECTOR(4, 10, 0) rec[-2] rec[-26]
    DETECTOR(8, 10, 0) rec[-1] rec[-25]'''),
                                             7: stim.Circuit('''H 2 6 10 19 23 27 32 36 40 49 53 57 62 66 70 79 83 87 92 96 100 109 113 117
    TICK
    CX 2 3 32 33 62 63 92 93 19 20 49 50 79 80 6 7 36 37 66 67 96 97 23 24 53 54 83 84 10 11 40 41 70 71 100 101 27 28 57 58 87 88 31 30 61 60 91 90 18 17 48 47 78 77 35 34 65 64 95 94 22 21 52 51 82 81 39 38 69 68 99 98 26 25 56 55 86 85 43 42 73 72 103 102
    TICK
    CX 2 1 32 31 62 61 92 91 19 18 49 48 79 78 6 5 36 35 66 65 96 95 23 22 53 52 83 82 10 9 40 39 70 69 100 99 27 26 57 56 87 86 16 30 46 60 76 90 3 17 33 47 63 77 20 34 50 64 80 94 7 21 37 51 67 81 24 38 54 68 84 98 11 25 41 55 71 85 28 42 58 72 88 102
    TICK
    CX 32 18 62 48 92 78 19 5 49 35 79 65 109 95 36 22 66 52 96 82 23 9 53 39 83 69 113 99 40 26 70 56 100 86 27 13 57 43 87 73 117 103 16 17 46 47 76 77 33 34 63 64 93 94 20 21 50 51 80 81 37 38 67 68 97 98 24 25 54 55 84 85 41 42 71 72 101 102 28 29 58 59 88 89
    TICK
    CX 32 16 62 46 92 76 19 3 49 33 79 63 109 93 36 20 66 50 96 80 23 7 53 37 83 67 113 97 40 24 70 54 100 84 27 11 57 41 87 71 117 101 1 17 31 47 61 77 18 34 48 64 78 94 5 21 35 51 65 81 22 38 52 68 82 98 9 25 39 55 69 85 26 42 56 72 86 102 13 29 43 59 73 89
    TICK
    H 2 6 10 19 23 27 32 36 40 49 53 57 62 66 70 79 83 87 92 96 100 109 113 117
    TICK
    MR 2 6 10 17 19 21 23 25 27 29 30 32 34 36 38 40 42 47 49 51 53 55 57 59 60 62 64 66 68 70 72 77 79 81 83 85 87 89 90 92 94 96 98 100 102 109 113 117
    DETECTOR(2, 0, 0) rec[-48] rec[-96]
    DETECTOR(6, 0, 0) rec[-47] rec[-95]
    DETECTOR(10, 0, 0) rec[-46] rec[-94]
    DETECTOR(2, 2, 0) rec[-45] rec[-93]
    DETECTOR(4, 2, 0) rec[-44] rec[-92]
    DETECTOR(6, 2, 0) rec[-43] rec[-91]
    DETECTOR(8, 2, 0) rec[-42] rec[-90]
    DETECTOR(10, 2, 0) rec[-41] rec[-89]
    DETECTOR(12, 2, 0) rec[-40] rec[-88]
    DETECTOR(14, 2, 0) rec[-39] rec[-87]
    DETECTOR(0, 4, 0) rec[-38] rec[-86]
    DETECTOR(2, 4, 0) rec[-37] rec[-85]
    DETECTOR(4, 4, 0) rec[-36] rec[-84]
    DETECTOR(6, 4, 0) rec[-35] rec[-83]
    DETECTOR(8, 4, 0) rec[-34] rec[-82]
    DETECTOR(10, 4, 0) rec[-33] rec[-81]
    DETECTOR(12, 4, 0) rec[-32] rec[-80]
    DETECTOR(2, 6, 0) rec[-31] rec[-79]
    DETECTOR(4, 6, 0) rec[-30] rec[-78]
    DETECTOR(6, 6, 0) rec[-29] rec[-77]
    DETECTOR(8, 6, 0) rec[-28] rec[-76]
    DETECTOR(10, 6, 0) rec[-27] rec[-75]
    DETECTOR(12, 6, 0) rec[-26] rec[-74]
    DETECTOR(14, 6, 0) rec[-25] rec[-73]
    DETECTOR(0, 8, 0) rec[-24] rec[-72]
    DETECTOR(2, 8, 0) rec[-23] rec[-71]
    DETECTOR(4, 8, 0) rec[-22] rec[-70]
    DETECTOR(6, 8, 0) rec[-21] rec[-69]
    DETECTOR(8, 8, 0) rec[-20] rec[-68]
    DETECTOR(10, 8, 0) rec[-19] rec[-67]
    DETECTOR(12, 8, 0) rec[-18] rec[-66]
    DETECTOR(2, 10, 0) rec[-17] rec[-65]
    DETECTOR(4, 10, 0) rec[-16] rec[-64]
    DETECTOR(6, 10, 0) rec[-15] rec[-63]
    DETECTOR(8, 10, 0) rec[-14] rec[-62]
    DETECTOR(10, 10, 0) rec[-13] rec[-61]
    DETECTOR(12, 10, 0) rec[-12] rec[-60]
    DETECTOR(14, 10, 0) rec[-11] rec[-59]
    DETECTOR(0, 12, 0) rec[-10] rec[-58]
    DETECTOR(2, 12, 0) rec[-9] rec[-57]
    DETECTOR(4, 12, 0) rec[-8] rec[-56]
    DETECTOR(6, 12, 0) rec[-7] rec[-55]
    DETECTOR(8, 12, 0) rec[-6] rec[-54]
    DETECTOR(10, 12, 0) rec[-5] rec[-53]
    DETECTOR(12, 12, 0) rec[-4] rec[-52]
    DETECTOR(4, 14, 0) rec[-3] rec[-51]
    DETECTOR(8, 14, 0) rec[-2] rec[-50]
    DETECTOR(12, 14, 0) rec[-1] rec[-49]''') }

eoc_meas = {3:stim.Circuit('''M 1 3 5 8 10 12 15 17 19
    DETECTOR rec[-3] rec[-6] rec[-13]
    DETECTOR rec[-5] rec[-6] rec[-8] rec[-9] rec[-16]
    DETECTOR rec[-1] rec[-2] rec[-4] rec[-5] rec[-11]
    DETECTOR rec[-4] rec[-7] rec[-14]
    OBSERVABLE_INCLUDE(0) rec[-7] rec[-8] rec[-9]'''), 
           5: stim.Circuit('''M 1 3 5 7 9 12 14 16 18 20 23 25 27 29 31 34 36 38 40 42 45 47 49 51 53
    DETECTOR(0, 4, 1) rec[-15] rec[-20] rec[-42]
    DETECTOR(0, 8, 1) rec[-5] rec[-10] rec[-32]
    DETECTOR(2, 2, 1) rec[-19] rec[-20] rec[-24] rec[-25] rec[-47]
    DETECTOR(2, 6, 1) rec[-9] rec[-10] rec[-14] rec[-15] rec[-37]
    DETECTOR(4, 4, 1) rec[-13] rec[-14] rec[-18] rec[-19] rec[-40]
    DETECTOR(4, 8, 1) rec[-3] rec[-4] rec[-8] rec[-9] rec[-30]
    DETECTOR(6, 2, 1) rec[-17] rec[-18] rec[-22] rec[-23] rec[-45]
    DETECTOR(6, 6, 1) rec[-7] rec[-8] rec[-12] rec[-13] rec[-35]
    DETECTOR(8, 4, 1) rec[-11] rec[-12] rec[-16] rec[-17] rec[-38]
    DETECTOR(8, 8, 1) rec[-1] rec[-2] rec[-6] rec[-7] rec[-28]
    DETECTOR(10, 2, 1) rec[-16] rec[-21] rec[-43]
    DETECTOR(10, 6, 1) rec[-6] rec[-11] rec[-33]
    OBSERVABLE_INCLUDE(0) rec[-21] rec[-22] rec[-23] rec[-24] rec[-25]'''),
           7:stim.Circuit('''M 1 3 5 7 9 11 13 16 18 20 22 24 26 28 31 33 35 37 39 41 43 46 48 50 52 54 56 58 61 63 65 67 69 71 73 76 78 80 82 84 86 88 91 93 95 97 99 101 103
    DETECTOR(0, 4, 1) rec[-35] rec[-42] rec[-87]
    DETECTOR(0, 8, 1) rec[-21] rec[-28] rec[-73]
    DETECTOR(0, 12, 1) rec[-7] rec[-14] rec[-59]
    DETECTOR(2, 2, 1) rec[-41] rec[-42] rec[-48] rec[-49] rec[-94]
    DETECTOR(2, 6, 1) rec[-27] rec[-28] rec[-34] rec[-35] rec[-80]
    DETECTOR(2, 10, 1) rec[-13] rec[-14] rec[-20] rec[-21] rec[-66]
    DETECTOR(4, 4, 1) rec[-33] rec[-34] rec[-40] rec[-41] rec[-85]
    DETECTOR(4, 8, 1) rec[-19] rec[-20] rec[-26] rec[-27] rec[-71]
    DETECTOR(4, 12, 1) rec[-5] rec[-6] rec[-12] rec[-13] rec[-57]
    DETECTOR(6, 2, 1) rec[-39] rec[-40] rec[-46] rec[-47] rec[-92]
    DETECTOR(6, 6, 1) rec[-25] rec[-26] rec[-32] rec[-33] rec[-78]
    DETECTOR(6, 10, 1) rec[-11] rec[-12] rec[-18] rec[-19] rec[-64]
    DETECTOR(8, 4, 1) rec[-31] rec[-32] rec[-38] rec[-39] rec[-83]
    DETECTOR(8, 8, 1) rec[-17] rec[-18] rec[-24] rec[-25] rec[-69]
    DETECTOR(8, 12, 1) rec[-3] rec[-4] rec[-10] rec[-11] rec[-55]
    DETECTOR(10, 2, 1) rec[-37] rec[-38] rec[-44] rec[-45] rec[-90]
    DETECTOR(10, 6, 1) rec[-23] rec[-24] rec[-30] rec[-31] rec[-76]
    DETECTOR(10, 10, 1) rec[-9] rec[-10] rec[-16] rec[-17] rec[-62]
    DETECTOR(12, 4, 1) rec[-29] rec[-30] rec[-36] rec[-37] rec[-81]
    DETECTOR(12, 8, 1) rec[-15] rec[-16] rec[-22] rec[-23] rec[-67]
    DETECTOR(12, 12, 1) rec[-1] rec[-2] rec[-8] rec[-9] rec[-53]
    DETECTOR(14, 2, 1) rec[-36] rec[-43] rec[-88]
    DETECTOR(14, 6, 1) rec[-22] rec[-29] rec[-74]
    DETECTOR(14, 10, 1) rec[-8] rec[-15] rec[-60]
    OBSERVABLE_INCLUDE(0) rec[-43] rec[-44] rec[-45] rec[-46] rec[-47] rec[-48] rec[-49]''')}

