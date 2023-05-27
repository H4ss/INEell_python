from math import atan, cos, sqrt
import sys
import numpy as np

FluSor = np.zeros((0, 0))
MatAbb = np.zeros((0, 0))
MatCar = np.zeros((0, 0))
MatFleXY = np.zeros((0, 0))
MatFleYX = np.zeros((0, 0))
MatTor = np.zeros((0, 0))
MatLap = np.zeros((0, 0))
MatTagXZ = np.zeros((0, 0))
MatTagYZ = np.zeros((0, 0))
MatTop = np.zeros((0, 0))


def Command_click1():
    global FluSor, MatAbb, MatCar, MatFleXY, MatFleYX, MatTor, MatLap, MatTagXZ, MatTagYZ, MatTop
    NumRig = 76
    NumCol = 101
    FluSor = np.zeros((NumRig+1, NumCol+1))
    MatAbb = np.zeros((NumRig+1, NumCol+1))
    MatCar = np.zeros((NumRig+1, NumCol+1))
    MatFleXY = np.zeros((NumRig+1, NumCol+1))
    MatFleYX = np.zeros((NumRig+1, NumCol+1))
    MatTor = np.zeros((NumRig+1, NumCol+1))
    MatLap = np.zeros((NumRig+1, NumCol+1))
    MatTagXZ = np.zeros((NumRig+1, NumCol+1))
    MatTagYZ = np.zeros((NumRig+1, NumCol+1))
    MatTop = np.zeros((NumRig+1, NumCol+1), dtype=int)
    ParRig = np.zeros((NumRig+1, NumCol), dtype=float)
    ParCol = np.zeros((NumRig, NumCol+1), dtype=float)

    LunPia = 10.0     # length of the plate
    LarPia = 7.5      # width of the plate
    ModYou = 25000000000.0  # Young's modulus
    CoePoi = 0.2      # Poisson's ratio
    SpePia = 0.4      # thickness of the plate
    PesCal = 25000.0  # specific weight
    CarDis0 = 5000.0  # ???
    DifCarR = 2500.0  # ???
    DifCarC = 5000.0  # ???
    PesPro = PesCal * SpePia
    CarCon1 = 100000.0  # ???
    CarCon2 = 50000.0   # ???
    CarCon3 = 50000.0   # ???
    CarCon4 = 100000.0  # ???
    # flexural stiffness D (eq. 2)
    RigFle = ModYou * SpePia ** 3 / (12 * (1 - CoePoi ** 2))
    DelRig = LarPia / float(NumRig - 1)
    DelCol = LunPia / float(NumCol - 1)

    for IndRig in range(1, NumRig+1):
        for IndCol in range(1, NumCol+1):
            if IndRig == 1 or IndRig == NumRig or IndCol == 1 or IndCol == NumCol:
                MatTop[IndRig][IndCol] = -1
            else:
                MatTop[IndRig][IndCol] = 1

    for IndRig in range(1, NumRig+1):
        for IndCol in range(1, NumCol+1):
            if IndCol < NumCol:
                ParRig[IndRig][IndCol] = 1
            if IndRig < NumRig:
                ParCol[IndRig][IndCol] = 1

    if "TipCar1" in sys.argv:
        TipCar = 1
    elif "TipCar2" in sys.argv:
        TipCar = 2
    else:
        raise ValueError("TipCar argument not provided")

    # I don't actually understand the code about CarDis
    CarDis(TipCar, NumRig, NumCol, DelRig, DelCol, PesPro, CarDis0,
           DifCarR, DifCarC, CarCon1, CarCon2, CarCon3, CarCon4, MatTop)

    for IndRig in range(1, NumRig+1):
        for IndCol in range(1, NumCol+1):
            FluSor[IndRig, IndCol] = MatCar[IndRig, IndCol] / RigFle
            MatLap[IndRig, IndCol] = 0

            # anche sul bordo si mette zero, condizione di Dirichlet di valore nullo
            if IndRig == 1 or IndRig == NumRig or IndCol == 1 or IndCol == NumCol:
                MatLap[IndRig, IndCol] = 0

    FatSur = 0.0
    TolCon = 0.000001
    IndErr = 1
    MatLap, ParRig, ParCol = INEell(NumRig, NumCol, DelRig, DelCol, ParRig,
           ParCol, FatSur, TolCon, MatLap, IndErr)

    for IndRig in range(1, NumRig+1):
        for IndCol in range(1, NumCol+1):
            FluSor[IndRig, IndCol] = MatLap[IndRig, IndCol] * DelRig * DelCol
            MatAbb[IndRig, IndCol] = 0.0

    TolCon = 0.000001
    IndErr = 1
    MatAbb, ParRig, ParCol = INEell(NumRig, NumCol, DelRig, DelCol, ParRig,
           ParCol, FatSur, TolCon, MatAbb, IndErr)

    CoePoi = 0.3
    RigFle = 1.0
    TaFlTo(NumRig, NumCol, DelRig, DelCol, CoePoi, RigFle, MatTagXZ)

    MaxFleXY = 0.0
    MaxFleYX = 0.0
    MxTaPoXZ = 0.0
    MxTaNeXZ = 0.0
    MxTaPoYZ = 0.0
    MxTaNeYZ = 0.0
    MxToPo = 0.0
    MxToNe = 0.0
            
    AsMaFlXY = 0.0
    OrMaFlXY = 0.0
    AsMaFlYX = 0.0
    OrMaFlYX = 0.0
    AsMxTaPoXZ = 0.0
    OrMxTaPoXZ = 0.0
    AsMxTaNeXZ = 0.0
    OrMxTaNeXZ = 0.0
    AsMxTaPoYZ = 0.0
    OrMxTaPoYZ = 0.0
    AsMxTaNeYZ = 0.0
    OrMxTaNeYZ = 0.0
    AsMxToPo = 0.0
    OrMxToPo = 0.0
    AsMxToNe = 0.0
    OrMxToNe = 0.0

    for IndRig in range(1, NumRig+1):
        for IndCol in range(1, NumCol+1):
            if MatFleXY[IndRig, IndCol] > MaxFleXY:
                MaxFleXY = MatFleXY[IndRig, IndCol]
                AsMaFlXY = float(IndCol - 1) * DelCol
                OrMaFlXY = float(IndRig - 1) * DelRig
            if MatFleYX[IndRig, IndCol] > MaxFleYX:
                MaxFleYX = MatFleYX[IndRig, IndCol]
                AsMaFlYX = float(IndCol - 1) * DelCol
                OrMaFlYX = float(IndRig - 1) * DelRig
            if MatTagXZ[IndRig, IndCol] > MxTaPoXZ:
                MxTaPoXZ = MatTagXZ[IndRig, IndCol]
                AsMxTaPoXZ = float(IndCol - 1) * DelCol
                OrMxTaPoXZ = float(IndRig - 1) * DelRig
            if MatTagXZ[IndRig, IndCol] < MxTaNeXZ:
                MxTaNeXZ = MatTagXZ[IndRig, IndCol]
                AsMxTaNeXZ = float(IndCol - 1) * DelCol
                OrMxTaNeXZ = float(IndRig - 1) * DelRig
            if MatTagYZ[IndRig, IndCol] > MxTaPoYZ:
                MxTaPoYZ = MatTagYZ[IndRig, IndCol]
                AsMxTaPoYZ = float(IndCol - 1) * DelCol
                OrMxTaPoYZ = float(IndRig - 1) * DelRig
            if MatTagYZ[IndRig, IndCol] < MxTaNeYZ:
                MxTaNeYZ = MatTagYZ[IndRig, IndCol]
                AsMxTaNeYZ = float(IndCol - 1) * DelCol
                OrMxTaNeYZ = float(IndRig - 1) * DelRig
            if MatTor[IndRig, IndCol] > MxToPo:
                MxToPo = MatTor[IndRig, IndCol]
                AsMxToPo = float(IndCol - 1) * DelCol
                OrMxToPo = float(IndRig - 1) * DelRig
            if MatTor[IndRig, IndCol] < MxToNe:
                MxToNe = MatTor[IndRig, IndCol]
                AsMxToNe = float(IndCol - 1) * DelCol
                OrMxToNe = float(IndRig - 1) * DelRig

    # Create arrays of text boxes
    Text1 = np.empty(3, dtype=object)
    Text2 = np.empty(3, dtype=object)
    Text3 = np.empty(3, dtype=object)
    Text4 = np.empty(3, dtype=object)
    Text5 = np.empty(3, dtype=object)
    Text6 = np.empty(3, dtype=object)
    Text7 = np.empty(3, dtype=object)
    Text8 = np.empty(3, dtype=object)

    # Format values and store in arrays of text
    Text1[0] = "{:.1f}".format(MaxFleXY / 1000)
    Text1[1] = "{:.1f}".format(AsMaFlXY)
    Text1[2] = "{:.1f}".format(OrMaFlXY)

    Text2[0] = "{:.1f}".format(MaxFleYX / 1000)
    Text2[1] = "{:.1f}".format(AsMaFlYX)
    Text2[2] = "{:.1f}".format(OrMaFlYX)

    Text3[0] = "{:.1f}".format(MxTaPoXZ / 1000)
    Text3[1] = "{:.1f}".format(AsMxTaPoXZ)
    Text3[2] = "{:.1f}".format(OrMxTaPoXZ)

    Text4[0] = "{:.1f}".format(MxTaNeXZ / 1000)
    Text4[1] = "{:.1f}".format(AsMxTaNeXZ)
    Text4[2] = "{:.1f}".format(OrMxTaNeXZ)

    Text5[0] = "{:.1f}".format(MxTaPoYZ / 1000)
    Text5[1] = "{:.1f}".format(AsMxTaPoYZ)
    Text5[2] = "{:.1f}".format(OrMxTaPoYZ)

    Text6[0] = "{:.1f}".format(MxTaNeYZ / 1000)
    Text6[1] = "{:.1f}".format(AsMxTaNeYZ)
    Text6[2] = "{:.1f}".format(OrMxTaNeYZ)

    Text7[0] = "{:.1f}".format(MxToPo / 1000)
    Text7[1] = "{:.1f}".format(AsMxToPo)
    Text7[2] = "{:.1f}".format(OrMxToPo)

    Text8[0] = "{:.1f}".format(MxToNe / 1000)
    Text8[1] = "{:.1f}".format(AsMxToNe)
    Text8[2] = "{:.1f}".format(OrMxToNe)

    # Print the values stored in the arrays of text
    print("MaxFleXY: {}".format(Text1[0]))
    print("AsMaFlXY: {}".format(Text1[1]))
    print("OrMaFlXY: {}".format(Text1[2]))
    print("MaxFleYX: {}".format(Text2[0]))
    print("AsMaFlYX: {}".format(Text2[1]))
    print("OrMaFlYX: {}".format(Text2[2]))
    print("MxTaPoXZ: {}".format(Text3[0]))
    print("AsMxTaPoXZ: {}".format(Text3[1]))
    print("OrMxTaPoXZ: {}".format(Text3[2]))
    print("MxTaNeXZ: {}".format(Text4[0]))
    print("AsMxTaNeXZ: {}".format(Text4[1]))
    print("OrMxTaNeXZ: {}".format(Text4[2]))
    print("MxTaPoYZ: {}".format(Text5[0]))
    print("AsMxTaPoYZ: {}".format(Text5[1]))
    print("OrMxTaPoYZ: {}".format(Text5[2]))
    print("MxTaNeYZ: {}".format(Text6[0]))
    print("AsMxTaNeYZ: {}".format(Text6[1]))
    print("OrMxTaNeYZ: {}".format(Text6[2]))
    print("MxToPo: {}".format(Text7[0]))
    print("AsMxToPo: {}".format(Text7[1]))
    print("OrMxToPo: {}".format(Text7[2]))
    print("MxToNe: {}".format(Text8[0]))
    print("AsMxToNe: {}".format(Text8[1]))
    print("OrMxToNe: {}".format(Text8[2]))

    # Memorizzazione dei risultati
    PatFil = "./"
    if "TipCar1" in sys.argv:
        with open(PatFil + "MatAbb1.txt", "w") as f1, \
             open(PatFil + "MatFleXY1.txt", "w") as f2, \
             open(PatFil + "MatFleYX1.txt", "w") as f3, \
             open(PatFil + "MatTagXZ1.txt", "w") as f4, \
             open(PatFil + "MatTagYZ1.txt", "w") as f5, \
             open(PatFil + "MatTor1.txt", "w") as f6:
     
            for IndRig in range(1, NumRig + 1):
                for IndCol in range(1, NumCol + 1):
                    f1.write("{:.6f}  ".format(1000 * MatAbb[IndRig, IndCol]))
                    f1.write("\n")
                    f2.write("{:08.3f}  ".format(MatFleXY[IndRig, IndCol] / 1000))
                    f2.write("\n")
                    f3.write("{:08.3f}  ".format(MatFleYX[IndRig, IndCol] / 1000))
                    f3.write("\n")
                    if MatTagXZ[IndRig, IndCol] >= 0:
                        f4.write("{:08.3f}  ".format(MatTagXZ[IndRig, IndCol] / 1000))
                        f4.write("\n")
                    else:
                        f4.write("#{:07.3f}  ".format(MatTagXZ[IndRig, IndCol] / 1000))
                        f4.write("\n")
                    if MatTagYZ[IndRig, IndCol] >= 0:
                        f5.write("{:08.3f}  ".format(MatTagYZ[IndRig, IndCol] / 1000))
                        f5.write("\n")
                    else:
                        f5.write("#{:07.3f}  ".format(MatTagYZ[IndRig, IndCol] / 1000))
                        f5.write("\n")
                    if MatTor[IndRig, IndCol] >= 0:
                        f6.write("{:08.3f}  ".format(MatTor[IndRig, IndCol] / 1000))
                        f6.write("\n")
                    else:
                        f6.write("#{:07.3f}  ".format(MatTor[IndRig, IndCol] / 1000))
                        f6.write("\n")
                # for IndFil in range(1, 7):


    elif "TipCar2" in sys.argv:
        with open(PatFil + "MatAbb2.txt", "w") as f1, \
             open(PatFil + "MatFleXY2.txt", "w") as f2, \
             open(PatFil + "MatFleYX2.txt", "w") as f3, \
             open(PatFil + "MatTagXZ2.txt", "w") as f4, \
             open(PatFil + "MatTagYZ2.txt", "w") as f5, \
             open(PatFil + "MatTor2.txt", "w") as f6:
            
            for IndRig in range(1, NumRig + 1):
                for IndCol in range(1, NumCol + 1):
                    f1.write("{:.6f}  ".format(1000 * MatAbb[IndRig, IndCol]))
                    f1.write("\n")
                    f2.write("{:08.3f}  ".format(MatFleXY[IndRig, IndCol] / 1000))
                    f2.write("\n")
                    f3.write("{:08.3f}  ".format(MatFleYX[IndRig, IndCol] / 1000))
                    f3.write("\n")
                    if MatTagXZ[IndRig, IndCol] >= 0:
                        f4.write("{:08.3f}  ".format(MatTagXZ[IndRig, IndCol] / 1000))
                        f4.write("\n")
                    else:
                        f4.write("#{:07.3f}  ".format(MatTagXZ[IndRig, IndCol] / 1000))
                        f4.write("\n")
                    if MatTagYZ[IndRig, IndCol] >= 0:
                        f5.write("{:08.3f}  ".format(MatTagYZ[IndRig, IndCol] / 1000))
                        f5.write("\n")
                    else:
                        f5.write("#{:07.3f}  ".format(MatTagYZ[IndRig, IndCol] / 1000))
                        f5.write("\n")
                    if MatTor[IndRig, IndCol] >= 0:
                        f6.write("{:08.3f}  ".format(MatTor[IndRig, IndCol] / 1000))
                        f6.write("\n")
                    else:
                        f6.write("#{:07.3f}  ".format(MatTor[IndRig, IndCol] / 1000))
                        f6.write("\n")
                
    # Creazione file per Surfer /abbassamenti/
    PatFil = "./"
    with open(PatFil + "S-MatAbb{}.txt".format(str(TipCar)[1]), "w") as f:
        for IndRig in range(1, NumRig + 1):
            for IndCol in range(1, NumCol + 1):
                ValAsc = (IndCol - 1) * DelCol
                ValOrd = (IndRig - 1) * DelRig
                f.write("{:.6f}  {:.6f}  {:.6f}\n".format(ValAsc, ValOrd, 1000 * MatAbb[IndRig, IndCol]))

    

def TaFlTo(NumRig, NumCol, DelRig, DelCol, CoePoi, RigFle, MatTagXZ):
    global FluSor, MatAbb, MatCar, MatFleXY, MatFleYX, MatTor, MatLap, MatTagYZ, MatTop
    DerAbbXY = 0
    for IndRig in range(1, NumRig+1):
        for IndCol in range(1, NumCol+1):
            TopNod = MatTop[IndRig][IndCol]
            if TopNod == -1:
                MatFleXY[IndRig][IndCol] = 0
                MatFleYX[IndRig][IndCol] = 0
                if IndRig == 1:
                    DelLapY = MatLap[2][IndCol] - MatLap[1][IndCol]
                    DerLapY = DelLapY / DelRig
                    if IndCol == 1:
                        DerAbbXY = MatAbb[2][2] / (DelRig * DelCol)
                    elif IndCol == NumCol:
                        DerAbbXY = -MatAbb[2][NumCol-1] / (DelRig * DelCol)
                    else:
                        DerAbbXY = (MatAbb[2][IndCol+1] - MatAbb[2]
                                    [IndCol-1]) / (2 * DelRig * DelCol)
                elif IndRig == NumRig:
                    DelLapY = MatLap[NumRig][IndCol] - MatLap[NumRig-1][IndCol]
                    DerLapY = DelLapY / DelRig
                    if IndCol == 1:
                        DerAbbXY = -MatAbb[NumRig-1][2] / (DelRig * DelCol)
                    elif IndCol == NumCol:
                        DerAbbXY = MatAbb[NumRig -
                                          1][NumCol-1] / (DelRig * DelCol)
                    else:
                        DerAbbXY = - \
                            (MatAbb[NumRig-1][IndCol+1] - MatAbb[NumRig-1]
                             [IndCol-1]) / (2 * DelRig * DelCol)
                else:
                    DelLapY = MatLap[IndRig+1][IndCol] - \
                        MatLap[IndRig-1][IndCol]
                    DerLapY = DelLapY / (2 * DelRig)
                    if IndCol == 1:
                        DerAbbXY = (
                            MatAbb[IndRig+1][2] - MatAbb[IndRig-1][2]) / (2 * DelRig * DelCol)
                    elif IndCol == NumCol:
                        DerAbbXY = - \
                            (MatAbb[IndRig+1][NumCol-1] - MatAbb[IndRig-1]
                             [NumCol-1]) / (2 * DelRig * DelCol)
                MatTagYZ[IndRig][IndCol] = -RigFle * DerLapY
                MatTor[IndRig][IndCol] = -RigFle * (1 - CoePoi) * DerAbbXY
                if IndCol == 1:
                    DelLapX = MatLap[IndRig][2] - MatLap[IndRig][1]
                    DerLapX = DelLapX / DelCol
                elif IndCol == NumCol:
                    DelLapX = MatLap[IndRig][NumCol] - MatLap[IndRig][NumCol-1]
                    DerLapX = DelLapX / DelCol
                else:
                    DelLapX = MatLap[IndRig][IndCol+1] - \
                        MatLap[IndRig][IndCol-1]
                    DerLapX = DelLapX / (2 * DelCol)
                MatTagXZ[IndRig][IndCol] = -RigFle * DerLapX
            else:
                DelLapY = MatLap[IndRig+1][IndCol] - MatLap[IndRig-1][IndCol]
                MatTagYZ[IndRig][IndCol] = -RigFle * DelLapY / (2 * DelRig)
                DelLapX = MatLap[IndRig][IndCol+1] - MatLap[IndRig][IndCol-1]
                MatTagXZ[IndRig][IndCol] = -RigFle * DelLapX / (2 * DelCol)
                DifAbb1 = MatAbb[IndRig+1][IndCol+1] - \
                    MatAbb[IndRig+1][IndCol-1]
                DifAbb2 = MatAbb[IndRig-1][IndCol+1] - \
                    MatAbb[IndRig-1][IndCol-1]
                DerAbbXY = (DifAbb1 - DifAbb2) / (4 * DelRig * DelCol)
                MatTor[IndRig][IndCol] = -RigFle * (1 - CoePoi) * DerAbbXY
                SomAbbYP = MatAbb[IndRig+1][IndCol-1] + \
                    MatAbb[IndRig-1][IndCol-1]
                DerAbbY2 = 0.25 * \
                    (SomAbbYP - 2 * MatAbb[IndRig][IndCol-1]) / (DelRig ** 2)
                SomAbbYC = MatAbb[IndRig+1][IndCol] + MatAbb[IndRig-1][IndCol]
                DerAbbY2 = DerAbbY2 + 0.5 * \
                    (SomAbbYC - 2 * MatAbb[IndRig][IndCol]) / (DelRig ** 2)
                SomAbbYS = MatAbb[IndRig+1][IndCol+1] + \
                    MatAbb[IndRig-1][IndCol+1]
                DerAbbY2 = DerAbbY2 + 0.25 * \
                    (SomAbbYS - 2 * MatAbb[IndRig][IndCol+1]) / (DelRig ** 2)
                SomAbbXP = MatAbb[IndRig-1][IndCol+1] + \
                    MatAbb[IndRig-1][IndCol-1]
                DerAbbX2 = 0.25 * \
                    (SomAbbXP - 2 * MatAbb[IndRig-1][IndCol]) / (DelCol ** 2)
                SomAbbXC = MatAbb[IndRig][IndCol+1] + MatAbb[IndRig][IndCol-1]
                DerAbbX2 = DerAbbX2 + 0.5 * \
                    (SomAbbXC - 2 * MatAbb[IndRig][IndCol]) / (DelCol ** 2)
                SomAbbXS = MatAbb[IndRig+1][IndCol+1] + \
                    MatAbb[IndRig+1][IndCol-1]
                DerAbbX2 = DerAbbX2 + 0.25 * \
                    (SomAbbXS - 2 * MatAbb[IndRig+1][IndCol]) / (DelCol ** 2)
                MatFleXY[IndRig][IndCol] = -RigFle * \
                    (DerAbbX2 + CoePoi * DerAbbY2)
                MatFleYX[IndRig][IndCol] = -RigFle * \
                    (DerAbbY2 + CoePoi * DerAbbX2)


def CarDis(TipCar, NumRig, NumCol, DelRig, DelCol, PesPro, CarDis0, DifCarR, DifCarC, CarCon1, CarCon2, CarCon3, CarCon4, MatTop):
    global MatCar
    FasRig = (NumRig - 1) // 3
    FasCol = (NumCol - 1) // 4
    for IndRig in range(1, NumRig+1):
        for IndCol in range(1, NumCol+1):
            if IndRig <= FasRig:
                if IndCol <= FasCol:
                    MatCar[IndRig][IndCol] = (
                        PesPro + CarDis0) * DelRig * DelCol
                elif IndCol <= 2*FasCol:
                    VarCol = IndCol - FasCol - 1
                    CarLoc = PesPro + CarDis0 + DifCarC * VarCol / FasCol
                    MatCar[IndRig][IndCol] = CarLoc * DelRig * DelCol
                else:
                    MatCar[IndRig][IndCol] = (
                        PesPro + CarDis0 + DifCarC) * DelRig * DelCol
            elif IndRig <= 2*FasRig:
                if IndCol <= FasCol:
                    VarRig = IndRig - FasRig - 1
                    CarLoc = PesPro + CarDis0 + DifCarR * VarRig / FasRig
                    MatCar[IndRig][IndCol] = CarLoc * DelRig * DelCol
                elif IndCol <= 2*FasCol:
                    VarRig = IndRig - FasRig - 1
                    VarCol = IndCol - FasCol - 1
                    CarLoc = PesPro + CarDis0 + DifCarR * VarRig / FasRig
                    CarLoc += DifCarC * VarCol / FasCol
                    MatCar[IndRig][IndCol] = CarLoc * DelRig * DelCol
                else:
                    VarRig = IndRig - FasRig - 1
                    CarLoc = PesPro + CarDis0 + DifCarC + DifCarR * VarRig / FasRig
                    MatCar[IndRig][IndCol] = CarLoc * DelRig * DelCol
            else:
                if IndCol <= FasCol:
                    MatCar[IndRig][IndCol] = (
                        PesPro + CarDis0 + DifCarR) * DelRig * DelCol
                elif IndCol <= 2*FasCol:
                    VarCol = IndCol - FasCol - 1
                    CarLoc = PesPro + CarDis0 + DifCarR + DifCarC * VarCol / FasCol
                    MatCar[IndRig][IndCol] = CarLoc * DelRig * DelCol
                else:
                    CarLoc = PesPro + CarDis0 + DifCarC + DifCarR
                    MatCar[IndRig][IndCol] = CarLoc * DelRig * DelCol
            if MatTop[IndRig][IndCol] == -1:
                MatCar[IndRig][IndCol] /= 2
    MatCar[1][1] /= 2
    MatCar[1][NumCol] /= 2
    MatCar[NumRig][1] /= 2
    MatCar[NumRig][NumCol] /= 2
    if TipCar == 2:
        for IndRig in range(1, NumRig+1):
            for IndCol in range(1, NumCol+1):
                if 26 <= IndCol <= 27:
                    if 25 <= IndRig <= 27:
                        MatCar[IndRig][IndCol] += CarCon1 * DelCol * DelRig
                if 26 <= IndCol <= 27:
                    if 50 <= IndRig <= 52:
                        MatCar[IndRig][IndCol] += CarCon2 * DelCol * DelRig
                if 51 <= IndCol <= 52:
                    if 25 <= IndRig <= 27:
                        MatCar[IndRig][IndCol] += CarCon3 * DelCol * DelRig
                if 76 <= IndCol <= 77:
                    if 25 <= IndRig <= 27:
                        MatCar[IndRig][IndCol] += CarCon4 * DelCol * DelRig


def INEell(NumRig, NumCol, DelRig, DelCol, ParRig, ParCol, FatSur, TolCon, FunInc, IndErr):
    # MatTop, ParCol, ParRig, FluSor
    global FluSor, MatTop
    MasIte = 2500
    if NumRig < 1 or NumCol < 1 or NumRig + NumCol < 2:
        MesErr(IndErr, "INEell", -1)
        return
    if DelRig <= 0 or DelCol <= 0:
        MesErr(IndErr, "INEell", -2)
        return
    if TolCon <= 0:
        MesErr(IndErr, "INEell", -3)
        return
    ConDir = 1
    NodInt = 0
    for IndRig in range(1, NumRig+1):
        for IndCol in range(1, NumCol+1):
            TopNod = MatTop[IndRig-1][IndCol-1]
            # print(f'Top {MatTop[IndRig-1][IndCol-1]}')
            # print(f'Right {MatTop[IndRig-1][IndCol]}')
            # print(f'Bottom {MatTop[IndRig][IndCol-1]}')
            # print(f'Left {MatTop[IndRig-1][IndCol-2]}')
            print(f'Adjacent {MatTop[IndRig-2][IndCol-1]} rig {IndRig} col {IndCol}')

            # print(f'TopNod value: {TopNod}')
            # print(f'ConDir value: {ConDir}')

            if TopNod == -1:
                TopNod *= -1
            elif TopNod == 0 and IndRig != 1:
                TopNod = 1
            ConDir = (1+TopNod) / abs(1+TopNod) * ConDir
            if TopNod != 0 and TopNod != -1:
                NodInt += 1
                if IndRig == 1:
                    if TopNod != -2:
                        print(f'exit5 case1 {TopNod}')
                        MesErr(IndErr, "INEell", -5)
                        return
                elif TopNod == 1 and MatTop[IndRig-2][IndCol-1] == 0:
                    print(f'case1: {MatTop[IndRig-2][IndCol-1]}')
                    MesErr(IndErr, "INEell", -6)
                    return
                if IndCol == 1:
                    if TopNod != -2:
                        print(f'exit5 case2')
                        MesErr(IndErr, "INEell", -5)
                        return
                elif TopNod == 1 and MatTop[IndRig-1][IndCol-2] == 0:
                    print(f'case2')
                    MesErr(IndErr, "INEell", -6)
                    return
                if IndCol == NumCol:
                    if TopNod != -2:
                        print(f'exit5 case3')
                        MesErr(IndErr, "INEell", -5)
                        return
                elif TopNod == 1 and MatTop[IndRig-1][IndCol] == 0:
                    print(f'case3')
                    MesErr(IndErr, "INEell", -6)
                    return
                if IndRig == NumRig:
                    if TopNod != -2:
                        print(f'exit5 case4')
                        MesErr(IndErr, "INEell", -5)
                        return
                elif TopNod == 1 and MatTop[IndRig][IndCol-1] == 0:
                    print(f'case4')
                    MesErr(IndErr, "INEell", -6)
                    return
    if ConDir != 0:
        MesErr(IndErr, "INEell", -4)
        return
    RapDel = DelCol / DelRig
    if FatSur < 1 or FatSur >= 2:
        PiGre = 4 * atan(1)
        ColRig = NumCol / NumRig
        CosUno = cos(PiGre / sqrt(NodInt * ColRig))
        CosDue = cos(PiGre / sqrt(NodInt / ColRig))
        RhoJac2 = ((CosUno + RapDel ** 2 * CosDue) / (1 + RapDel ** 2)) ** 2
    else:
        RhoJac2 = 1 - (1 - 2 / FatSur) ** 2
    SurChe = 2
    for IndIte in range(1, MasIte+1):
        ErrMas = 0
        for IndSem in range(1, 3):
            SurChe = 1 / (1 - SurChe * RhoJac2 / 4)
            ColIni = IndSem
            for IndRig in range(1, NumRig+1):
                for IndCol in range(ColIni, NumCol+1, 2):
                    TopNod = MatTop[IndRig-1][IndCol-1]
                    if TopNod != 0 and TopNod != -1:
                        SomNum = 0
                        SomDen = 0
                        if IndRig > 1:
                            if MatTop[IndRig-2][IndCol-1] != 0:
                                CoeNor = ParCol[IndRig-2][IndCol-1] / RapDel
                                SomNum += CoeNor * FunInc[IndRig-2][IndCol-1]
                                SomDen += CoeNor
                        if IndCol > 1:
                            if MatTop[IndRig-1][IndCol-2] != 0:
                                CoeOve = ParRig[IndRig-1][IndCol-2] * RapDel
                                SomNum += CoeOve * FunInc[IndRig-1][IndCol-2]
                                SomDen += CoeOve
                        if IndCol < NumCol:
                            if MatTop[IndRig-1][IndCol] != 0:
                                CoeEst = ParRig[IndRig-1][IndCol-1] * RapDel
                                SomNum += CoeEst * FunInc[IndRig-1][IndCol]
                                SomDen += CoeEst
                        if IndRig < NumRig:
                            if MatTop[IndRig][IndCol-1] != 0:
                                CoeSud = ParCol[IndRig-1][IndCol-1] / RapDel
                                SomNum += CoeSud * FunInc[IndRig][IndCol-1]
                                SomDen += CoeSud
                        SomNum -= FluSor[IndRig-1][IndCol-1]
                        CorInc = FunInc[IndRig-1][IndCol-1] - SomNum / SomDen
                        if abs(CorInc) > ErrMas:
                            ErrMas = abs(CorInc)
                        FunInc[IndRig-1][IndCol-1] -= SurChe * CorInc
                ColIni = 3 - ColIni
        if ErrMas < TolCon:
            IndErr = 0
            return
    MesErr(IndErr, "INEell", 1)
    return FunInc, ParRig, ParCol


def MesErr(IndErr, NomSub, NumErr):
    if IndErr == 0:
        IndErr = NumErr
        return
    MsgBoxWrn = "Messaggio d'errore !"
    MsgBoxTxt = "INTERRUZIONE DELLE ELABORAZIONI\n" + " " * 22 + "errore n. " + \
        str(NumErr) + "\n" + " " * 13 + "nella routine: " + NomSub
    print(MsgBoxTxt + "\n\n" + MsgBoxWrn + "\n\nPremere INVIO per uscire...")
    exit()

def main():
    try:
        if sys.argv.__contains__("TipCar1") or sys.argv.__contains__("TipCar2"):
            Command_click1()
        else:
            print("Errore: nessun parametro passato")
    except Exception as e:
        print("Errore: " + str(e))

if __name__ == "__main__":
    main()