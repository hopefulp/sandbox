########################################################################
#
# $Header: /cx2/c2cap_c410/data/Cerius2-Resources/FORCE-FIELD/RCS/DREIDING2.21,v 1.2 1996/10/29 04:25:25 jan Exp $
#
# Author:S. Miller
#
# Purpose:
#
########################################################################
#
VERSION
 CERIUS2     1
END
#
HEADER

 ********************************************************************
                   DREIDING: A Generic Force Field

                               by

          S.L. Mayo, B.D. Olafson, W.A. Goddard III
           J. Phys. Chem. 1990, 94, 8897-8909

 ********************************************************************

END
#
PREFERENCES
 BONDS                                 T
 ANGLES                                T
 COULOMB                               T
 INVERSIONS                            F 
 TORSIONS                              F
 UREY_BRADLEY                          F
 HYDROGEN_BONDS                        F
 DIAGONAL_VAN_DER_WAALS                T
 OFF_DIAGONAL_VAN_DER_WAALS            T
 GENERATE_UNDEFINED_TERMS              F
 IGNORE_UNDEFINED_TERMS                F
 SHRINK_CH_BONDS                       F
 SINGLE_TORSION                        F
 SCALE_TORSIONS_ABOUT_COMMON_BOND      T
 EXOCYCLIC_TORSIONS_SCALE_FACTOR          1.00000
 SINGLE_INVERSION                      F
 H-BOND_METHOD                         SPLINE
 H-BOND_LIST                           T
 NON_BOND_LIST                         T
 DISTANCE_DEPENDENT_DIELETRIC_CONSTANT T
 COU_DIELETRIC_CONSTANT                   1.00000
 COU_EXCLUDE_1-2                       T
 COU_EXCLUDE_1-3                       T
 COU_EXCLUDE_1-4                       F
 COU_1-4_SCALE_FACTOR                     1.00000
 COU_METHOD                            EWALD
 EWALD_SUM_COU_ACCURACY                   0.00100
 EWALD_SUM_COU_ETA                        5.18900
 EWALD_SUM_COU_KCUT                       0.23862
 EWALD_SUM_COU_RCUT                      10.00000
 EWALD_SUM_COU_OPTIMIZE                NEVER
 VDW_RADII_COMBINATION_RULE            ARITHMETIC
 VDW_COMBINATION_RULE                  6TH-POWER
 VDW_INTER_CUT_OFF                       10.00000
 VDW_EXCLUDE_1-2                       T
 VDW_EXCLUDE_1-3                       T
 VDW_EXCLUDE_1-4                       F
 VDW_1-4_SCALE_FACTOR                     1.00000
 VDW_METHOD                            DIRECT    
 VDW_SPLINE_ON                            9.00000
 VDW_SPLINE_OFF                          14.00000
 EWALD_SUM_VDW_OPTIMIZE                NEVER
 EWALD_SUM_VDW_ACCURACY                   0.00100
 EWALD_SUM_VDW_ETA                        2.50000
 EWALD_SUM_VDW_KCUT                       0.50000
 EWALD_SUM_VDW_RCUT                       6.00000
 EWALD_SUM_VDW_REP_CUT                    6.00000
 FAST_EWALD_SUM_RATIO                    10.00000
 SLOW_EWALD_SUM_RATIO                     5.00000
 MINIMUM_IMAGE                         F
 ASSIGN_MASS                           T
 ASSIGN_CHARGE                         T
 ASSIGN_HYBRIDIZATION                  F
 ATOM_TYPE                             F 
 ATOM_TYPE_ALL                         F
 CALCULATE_BOND_ORDER                  F
END
#
ATOMTYPES
 H_         H       1.00800  0.0000   0   0   0
 H___A      H       1.00800  0.0000   0   0   0
 H___b      H       1.00800  0.0000   0   0   0
 B_3        B      10.81000  0.0000   3   0   0
 B_2        B      10.81000  0.0000   2   0   0
 C_34       C      16.04300  0.0000   3   4   0
 C_33       C      15.03500  0.0000   3   3   0
 C_32       C      14.02700  0.0000   3   2   0
 C_31       C      13.01900  0.0000   3   1   0
 C_3        C      12.01100  0.0000   3   0   0
 C_22       C      14.02700  0.0000   2   2   0
 C_21       C      13.01900  0.0000   2   1   0
 C_2        C      12.01100  0.0000   2   0   0
 C_R2       C      14.02700  0.0000   2   2   0
 C_R1       C      13.01900  0.0000   2   1   0
 C_R        C      12.01100  0.0000   2   0   0
 C_11       C      13.01900  0.0000   1   1   0
 C_1        C      12.01100  0.6512   3   0   0
 N_3        N      14.00670  0.0000   3   0   1
 N_2        N      14.00670  0.0000   2   0   1
 N_R        N      14.00670  0.0000   2   0   1
 N_1        N      14.00670  0.0000   1   0   1
 O_3        O      15.99940  0.0000   3   0   2
 O_2        O      15.99940 -0.3256   2   0   0
 O_R        O      15.99940  0.0000   2   0   1
 O_1        O      15.99940  0.0000   1   0   1
 F_         F      18.99840 -1.0000   0   0   3
 Al3        Al     26.98150  0.0000   3   0   0
 Si3        Si     28.08600  0.0000   3   0   0
 P_3        P      30.97380  0.0000   3   0   1
 S_3        S      32.06000  0.0000   3   0   2
 Cl         Cl     35.45300 -1.0000   0   0   3
 Ga3        Ga     69.72000  0.0000   3   0   0
 Ge3        Ge     72.59000  0.0000   3   0   0
 As3        As     74.92160  0.0000   3   0   1
 Se3        Se     78.96000  0.0000   3   0   2
 Br         Br     79.90400 -1.0000   0   0   3
 In3        In    114.82000  0.0000   3   0   0
 Sn3        Sn    118.69000  0.0000   3   0   0
 Sb3        Sb    121.75000  0.0000   3   0   1
 Te3        Te    127.60000  0.0000   3   0   2
 I_         I     126.90450 -1.0000   0   0   3
 Na         Na     22.99000  1.0000   0   0   0
 Ca         Ca     40.08000  2.0000   0   0   0
 Ti         Ti     47.90000  3.0000   0   0   0
 Fe         Fe     55.84700  3.0000   0   0   0
 Zn         Zn     65.37700  2.0000   0   0   0
 Tc         Tc     98.90620  2.0000   0   0   0
 Ru         Ru    101.07000  3.0000   0   0   0
END
#
DIAGONAL_VDW
 H_          LJ_6_12        3.1950  0.1520E-01
 H___A       LJ_6_12        3.1950  0.1000E-03
 H___b       LJ_6_12        3.1950  0.1520E-01
 B_3         LJ_6_12        4.0200  0.9500E-01
 B_2         LJ_6_12        4.0200  0.9500E-01
 C_34        LJ_6_12        4.2370  0.3016E+00
 C_33        LJ_6_12        4.1524  0.2500E+00
 C_32        LJ_6_12        4.0677  0.1984E+00
 C_31        LJ_6_12        3.9830  0.1467E+00
 C_3         LJ_6_12        3.8983  0.9510E-01
 C_22        LJ_6_12        4.0677  0.1984E+00
 C_21        LJ_6_12        3.9830  0.1467E+00
 C_2         LJ_6_12        3.8983  0.9510E-01
 C_R2        LJ_6_12        4.0677  0.1984E+00
 C_R1        LJ_6_12        4.2300  0.1356E+00
 C_R         LJ_6_12        3.8983  0.9510E-01
 C_11        LJ_6_12        3.9830  0.1467E+00
 C_1         LJ_6_12        3.0946  0.5589E-01
 N_3         LJ_6_12        3.6621  0.7740E-01
 N_2         LJ_6_12        3.6621  0.7740E-01
 N_R         LJ_6_12        3.6621  0.7740E-01
 N_1         LJ_6_12        3.6621  0.7740E-01
 O_3         LJ_6_12        3.4046  0.9570E-01
 O_2         LJ_6_12        3.4044  0.1600E-00
 O_R         LJ_6_12        3.4046  0.9570E-01
 O_1         LJ_6_12        3.4046  0.9570E-01
 F_          LJ_6_12        3.4720  0.7250E-01
 Al3         LJ_6_12        4.3900  0.3100E+00
 Si3         LJ_6_12        4.2700  0.3100E+00
 P_3         LJ_6_12        4.1500  0.3200E+00
 S_3         LJ_6_12        4.0300  0.3440E+00
 Cl          LJ_6_12        3.9503  0.2833E+00
 Ga3         LJ_6_12        4.3900  0.4000E+00
 Ge3         LJ_6_12        4.2700  0.4000E+00
 As3         LJ_6_12        4.1500  0.4100E+00
 Se3         LJ_6_12        4.0300  0.4300E+00
 Br          LJ_6_12        3.9500  0.3700E+00
 In3         LJ_6_12        4.5900  0.5500E+00
 Sn3         LJ_6_12        4.4700  0.5500E+00
 Sb3         LJ_6_12        4.3500  0.5500E+00
 Te3         LJ_6_12        4.2300  0.5700E+00
 I_          LJ_6_12        4.1500  0.5100E+00
 Na          LJ_6_12        3.1440  0.5000E+00
 Ca          LJ_6_12        3.4720  0.5000E-01
 Ti          LJ_6_12        4.5400  0.5500E-01
 Fe          LJ_6_12        4.5400  0.5500E-01
 Zn          LJ_6_12        4.5400  0.5500E-01
 Tc          LJ_6_12        4.5400  0.5500E-01
 Ru          LJ_6_12        4.5400  0.5500E-01
END
#
ATOM_TYPING_RULES
 H___A           H            0           0           0           1
 H___A           H            0           0           1           1
                 N            0           0           0           1
 H___A           H            0           0           1           1
                 O            0           0           0           1
 H___A           H            0           0           1           1
                 S            0           0           0           1
 H_              H            0           0           1           1
                 B            0           0           0           1
 H_              H            0           0           1           1
                 C            0           0           0           1
 H_              H            0           0           1           1
                 F            0           0           0           1
 H_              H            0           0           1           1
                 Al           0           0           0           1
 H_              H            0           0           1           1
                 Si           0           0           0           1
 H_              H            0           0           1           1
                 P            0           0           0           1
 H_              H            0           0           1           1
                 Cl           0           0           0           1
 H_              H            0           0           1           1
                 Ga           0           0           0           1
 H_              H            0           0           1           1
                 Ge           0           0           0           1
 H_              H            0           0           1           1
                 As           0           0           0           1
 H_              H            0           0           1           1
                 Se           0           0           0           1
 H_              H            0           0           1           1
                 Br           0           0           0           1
 H_              H            0           0           1           1
                 In           0           0           0           1
 H_              H            0           0           1           1
                 Sn           0           0           0           1
 H_              H            0           0           1           1
                 Sb           0           0           0           1
 H_              H            0           0           1           1
                 Te           0           0           0           1
 H_              H            0           0           1           1
                 I            0           0           0           1
 H___b           H            0           0           2           1
                 **           0           0           0           1
                 **           0           0           0           1
 B_3             B            0           0           4           1
                 **           0           0           0           1
                 **           0           0           0           1
                 **           0           0           0           1
                 **           0           0           0           1
 B_2             B            0           0           3           1
                 **           0           0           0           1
                 **           0           0           0           1
                 **           0           0           0           1
 B_2             B            0           0           1           1
                 N            0           0           0           1
 B_2             B            0           0           1           1
                 O            0           0           0           1
 C_3             C            3           0           0           1
 C_31            C            3           0           1           1
                 VH           0           0           0           1
 C_32            C            3           0           2           1
                 VH           0           0           0           1
                 VH           0           0           0           1
 C_33            C            3           0           3           1
                 VH           0           0           0           1
                 VH           0           0           0           1
                 VH           0           0           0           1
 C_34            C            3           0           4           1
                 VH           0           0           0           1
                 VH           0           0           0           1
                 VH           0           0           0           1
                 VH           0           0           0           1
 C_R2            C            2           0           3           1
                 **           2           0           2           1
                 VH           0           0           0           1
                 VH           0           0           0           1
                 **           2           0           0           1
                 **           2           0           0           1
 C_R2            C            2           0           3           1
                 **           2           0           2           1
                 VH           0           0           0           1
                 VH           0           0           0           1
                 **           2           0           0           1
                 O            3           0           0           1
 C_R1            C            2           0           2           1
                 **           2           0           2           1
                 VH           0           0           0           1
                 **           2           0           0           1
                 **           2           0           0           1
 C_R1            C            2           0           2           1
                 **           2           0           2           1
                 VH           0           0           0           1
                 **           2           0           0           1
                 O            3           0           0           1
 C_R1            C            2           0           3           1
                 **           2           0           0           1
                 **           2           0           0           1
                 VH           0           0           0           1
 C_R1            C            2           0           3           1
                 **           2           0           0           1
                 O            3           0           0           1
                 VH           0           0           0           1
 C_R             C            2           0           1           1
                 **           2           0           2           1
                 **           2           0           0           1
                 **           2           0           0           1
 C_R             C            2           0           1           1
                 **           2           0           2           1
                 **           2           0           0           1
                 O            3           0           0           1
 C_R             C            2           0           2           1
                 **           2           0           0           1
                 **           2           0           0           1
 C_R             C            2           0           2           1
                 **           2           0           0           1
                 O            3           0           0           1
 C_22            C            2           0           2           1
                 VH           0           0           0           1
                 VH           0           0           0           1
 C_21            C            2           0           1           1
                 VH           0           0           0           1
 C_2             C            2           0           0           1
 C_11            C            1           0           1           1
                 VH           0           0           0           1
 C_1             C            1           0           0           1
 N_3             N            3           0           0           1
 N_R             N            2           0           1           1
                 **           2           0           2           1
                 **           2           0           0           1
                 **           2           0           0           1
 N_R             N            2           0           1           1
                 **           2           0           2           1
                 **           2           0           0           1
                 O            3           0           0           1
 N_R             N            2           0           2           1
                 **           2           0           0           1
                 **           2           0           0           1
 N_R             N            2           0           2           1
                 **           2           0           0           1
                 O            3           0           0           1
 N_2             N            2           0           0           1
 N_2             N            0           0           1           1
                 B            0           0           0           1
 N_1             N            1           0           0           1
 O_3             O            3           0           0           1
 O_R             O            2           0           2           1
                 **           2           0           0           1
                 **           2           0           0           1
 O_R             O            3           0           2           1
                 **           2           0           0           1
                 **           0           0           0           1
 O_R             O            2           0           2           1
                 **           2           0           0           1
                 **           0           0           0           1
 O_2             O            3           0           2           1
                 VH           0           0           0           1
                 **           0           0           1           1
                 O            2           0           0           1
 O_2             O            2           0           0           1
 O_2             O            0           0           1           1
                 B            0           0           0           1
 O_1             O            1           0           0           1
 F_              F            0           0           0           1
 Al3             Al           3           0           0           1
 Si3             Si           3           0           0           1
 P_3             P            3           0           0           1
 S_3             S            3           0           0           1
 Cl              Cl           0           0           0           1
 Ga3             Ga           0           0           0           1
 Ge3             Ge           0           0           0           1
 As3             As           0           0           0           1
 Se3             Se           0           0           0           1
 Br              Br           0           0           0           1
 In3             In           0           0           0           1
 Sn3             Sn           0           0           0           1
 Sb3             Sb           0           0           0           1
 Te3             Te           0           0           0           1
 I_              I            0           0           0           1
 Na              Na           0           0           0           1
 Ca              Ca           0           0           0           1
 Ti              Ti           0           0           0           1
 Fe              Fe           0           0           0           1
 Zn              Zn           0           0           0           1
 Tc              Tc           0           0           0           1
 Ru              Ru           0           0           0           1
END
#
#
OFF_DIAGONAL_VDW
END
#
BOND_STRETCH
 H_       H_          HARMONIC     700.0000    0.6500
 H___A    H_          HARMONIC     700.0000    0.6500
 H___A    H___A       HARMONIC     700.0000    0.6500
 H___b    H_          HARMONIC     700.0000    0.8300
 H___b    H___A       HARMONIC     700.0000    0.8300
 H___b    H___b       HARMONIC     700.0000    1.0100
 B_3      H_          HARMONIC     700.0000    1.2000
 B_3      H___A       HARMONIC     700.0000    1.2000
 B_3      H___b       HARMONIC     700.0000    1.3800
 B_3      B_3         HARMONIC     700.0000    1.7500
 B_2      H_          HARMONIC     700.0000    1.1100
 B_2      H___A       HARMONIC     700.0000    1.1100
 B_2      H___b       HARMONIC     700.0000    1.2900
 B_2      B_3         HARMONIC     700.0000    1.6600
 B_2      B_2         HARMONIC    1400.0000    1.5700
 C_34     H_          HARMONIC     700.0000    1.0900
 C_34     H___A       HARMONIC     700.0000    1.0900
 C_34     H___b       HARMONIC     700.0000    1.2700
 C_34     B_3         HARMONIC     700.0000    1.6400
 C_34     B_2         HARMONIC     700.0000    1.5500
 C_34     C_34        HARMONIC     700.0000    1.5300
 C_33     H_          HARMONIC     700.0000    1.0900
 C_33     H___A       HARMONIC     700.0000    1.0900
 C_33     H___b       HARMONIC     700.0000    1.2700
 C_33     B_3         HARMONIC     700.0000    1.6400
 C_33     B_2         HARMONIC     700.0000    1.5500
 C_33     C_34        HARMONIC     700.0000    1.5300
 C_33     C_33        HARMONIC     700.0000    1.5300
 C_32     H_          HARMONIC     700.0000    1.0900
 C_32     H___A       HARMONIC     700.0000    1.0900
 C_32     H___b       HARMONIC     700.0000    1.2700
 C_32     B_3         HARMONIC     700.0000    1.6400
 C_32     B_2         HARMONIC     700.0000    1.5500
 C_32     C_34        HARMONIC     700.0000    1.5300
 C_32     C_33        HARMONIC     700.0000    1.5300
 C_32     C_32        HARMONIC     700.0000    1.5300
 C_31     H_          HARMONIC     700.0000    1.0900
 C_31     H___A       HARMONIC     700.0000    1.0900
 C_31     H___b       HARMONIC     700.0000    1.2700
 C_31     B_3         HARMONIC     700.0000    1.6400
 C_31     B_2         HARMONIC     700.0000    1.5500
 C_31     C_34        HARMONIC     700.0000    1.5300
 C_31     C_33        HARMONIC     700.0000    1.5300
 C_31     C_32        HARMONIC     700.0000    1.5300
 C_31     C_31        HARMONIC     700.0000    1.5300
 C_3      H_          HARMONIC     700.0000    1.0900
 C_3      H___A       HARMONIC     700.0000    1.0900
 C_3      H___b       HARMONIC     700.0000    1.2700
 C_3      B_3         HARMONIC     700.0000    1.6400
 C_3      B_2         HARMONIC     700.0000    1.5500
 C_3      C_34        HARMONIC     700.0000    1.5300
 C_3      C_33        HARMONIC     700.0000    1.5300
 C_3      C_32        HARMONIC     700.0000    1.5300
 C_3      C_31        HARMONIC     700.0000    1.5300
 C_3      C_3         HARMONIC     700.0000    1.5300
 C_22     H_          HARMONIC     700.0000    0.9900
 C_22     H___A       HARMONIC     700.0000    0.9900
 C_22     H___b       HARMONIC     700.0000    1.1700
 C_22     B_3         HARMONIC     700.0000    1.5400
 C_22     B_2         HARMONIC    1400.0000    1.4500
 C_22     C_34        HARMONIC     700.0000    1.4300
 C_22     C_33        HARMONIC     700.0000    1.4300
 C_22     C_32        HARMONIC     700.0000    1.4300
 C_22     C_31        HARMONIC     700.0000    1.4300
 C_22     C_3         HARMONIC     700.0000    1.4300
 C_22     C_22        HARMONIC    1400.0000    1.3300
 C_21     H_          HARMONIC     700.0000    0.9900
 C_21     H___A       HARMONIC     700.0000    0.9900
 C_21     H___b       HARMONIC     700.0000    1.1700
 C_21     B_3         HARMONIC     700.0000    1.5400
 C_21     B_2         HARMONIC    1400.0000    1.4500
 C_21     C_34        HARMONIC     700.0000    1.4300
 C_21     C_33        HARMONIC     700.0000    1.4300
 C_21     C_32        HARMONIC     700.0000    1.4300
 C_21     C_31        HARMONIC     700.0000    1.4300
 C_21     C_3         HARMONIC     700.0000    1.4300
 C_21     C_22        HARMONIC    1400.0000    1.3300
 C_21     C_21        HARMONIC    1400.0000    1.3300
 C_2      H_          HARMONIC     700.0000    0.9900
 C_2      H___A       HARMONIC     700.0000    0.9900
 C_2      H___b       HARMONIC     700.0000    1.1700
 C_2      B_3         HARMONIC     700.0000    1.5400
 C_2      B_2         HARMONIC    1400.0000    1.4500
 C_2      C_34        HARMONIC     700.0000    1.4300
 C_2      C_33        HARMONIC     700.0000    1.4300
 C_2      C_32        HARMONIC     700.0000    1.4300
 C_2      C_31        HARMONIC     700.0000    1.4300
 C_2      C_3         HARMONIC     700.0000    1.4300
 C_2      C_22        HARMONIC    1400.0000    1.3300
 C_2      C_21        HARMONIC    1400.0000    1.3300
 C_2      C_2         HARMONIC    1400.0000    1.3300
 C_R2     H_          HARMONIC     700.0000    1.0200
 C_R2     H___A       HARMONIC     700.0000    1.0200
 C_R2     H___b       HARMONIC     700.0000    1.2000
 C_R2     B_3         HARMONIC     700.0000    1.5700
 C_R2     B_2         HARMONIC    1400.0000    1.4800
 C_R2     C_34        HARMONIC     700.0000    1.4600
 C_R2     C_33        HARMONIC     700.0000    1.4600
 C_R2     C_32        HARMONIC     700.0000    1.4600
 C_R2     C_31        HARMONIC     700.0000    1.4600
 C_R2     C_3         HARMONIC     700.0000    1.4600
 C_R2     C_22        HARMONIC    1400.0000    1.3600
 C_R2     C_21        HARMONIC    1400.0000    1.3600
 C_R2     C_2         HARMONIC    1400.0000    1.3600
 C_R2     C_R2        HARMONIC    1050.0000    1.3900
 C_R1     H_          HARMONIC     700.0000    1.0200
 C_R1     H___A       HARMONIC     700.0000    1.0200
 C_R1     H___b       HARMONIC     700.0000    1.2000
 C_R1     B_3         HARMONIC     700.0000    1.5700
 C_R1     B_2         HARMONIC    1400.0000    1.4800
 C_R1     C_34        HARMONIC     700.0000    1.4600
 C_R1     C_33        HARMONIC     700.0000    1.4600
 C_R1     C_32        HARMONIC     700.0000    1.4600
 C_R1     C_31        HARMONIC     700.0000    1.4600
 C_R1     C_3         HARMONIC     700.0000    1.4600
 C_R1     C_22        HARMONIC    1400.0000    1.3600
 C_R1     C_21        HARMONIC    1400.0000    1.3600
 C_R1     C_2         HARMONIC    1400.0000    1.3600
 C_R1     C_R2        HARMONIC    1050.0000    1.3900
 C_R1     C_R1        HARMONIC    1050.0000    1.3900
 C_R      H_          HARMONIC     700.0000    1.0200
 C_R      H___A       HARMONIC     700.0000    1.0200
 C_R      H___b       HARMONIC     700.0000    1.2000
 C_R      B_3         HARMONIC     700.0000    1.5700
 C_R      B_2         HARMONIC    1400.0000    1.4800
 C_R      C_34        HARMONIC     700.0000    1.4600
 C_R      C_33        HARMONIC     700.0000    1.4600
 C_R      C_32        HARMONIC     700.0000    1.4600
 C_R      C_31        HARMONIC     700.0000    1.4600
 C_R      C_3         HARMONIC     700.0000    1.4600
 C_R      C_22        HARMONIC    1400.0000    1.3600
 C_R      C_21        HARMONIC    1400.0000    1.3600
 C_R      C_2         HARMONIC    1400.0000    1.3600
 C_R      C_R2        HARMONIC    1050.0000    1.3900
 C_R      C_R1        HARMONIC    1050.0000    1.3900
 C_R      C_R         HARMONIC    1050.0000    1.3900
 C_11     H_          HARMONIC     700.0000    0.9220
 C_11     H___A       HARMONIC     700.0000    0.9220
 C_11     H___b       HARMONIC     700.0000    1.1020
 C_11     B_3         HARMONIC     700.0000    1.4720
 C_11     B_2         HARMONIC     700.0000    1.3820
 C_11     C_34        HARMONIC     700.0000    1.3620
 C_11     C_33        HARMONIC     700.0000    1.3620
 C_11     C_32        HARMONIC     700.0000    1.3620
 C_11     C_31        HARMONIC     700.0000    1.3620
 C_11     C_3         HARMONIC     700.0000    1.3620
 C_11     C_22        HARMONIC     700.0000    1.2620
 C_11     C_21        HARMONIC     700.0000    1.2620
 C_11     C_2         HARMONIC     700.0000    1.2620
 C_11     C_R2        HARMONIC     700.0000    1.2920
 C_11     C_R1        HARMONIC     700.0000    1.2920
 C_11     C_R         HARMONIC     700.0000    1.2920
 C_11     C_11        HARMONIC    2100.0000    1.1940
 C_1      H_          HARMONIC     700.0000    0.9220
 C_1      H___A       HARMONIC     700.0000    0.9220
 C_1      H___b       HARMONIC     700.0000    1.1020
 C_1      B_3         HARMONIC     700.0000    1.4720
 C_1      B_2         HARMONIC     700.0000    1.3820
 C_1      C_34        HARMONIC     700.0000    1.3620
 C_1      C_33        HARMONIC     700.0000    1.3620
 C_1      C_32        HARMONIC     700.0000    1.3620
 C_1      C_31        HARMONIC     700.0000    1.3620
 C_1      C_3         HARMONIC     700.0000    1.3620
 C_1      C_22        HARMONIC     700.0000    1.2620
 C_1      C_21        HARMONIC     700.0000    1.2620
 C_1      C_2         HARMONIC     700.0000    1.2620
 C_1      C_R2        HARMONIC     700.0000    1.2920
 C_1      C_R1        HARMONIC     700.0000    1.2920
 C_1      C_R         HARMONIC     700.0000    1.2920
 C_1      C_11        HARMONIC    2100.0000    1.1940
 C_1      C_1         HARMONIC    2100.0000    1.1940
 N_3      H_          HARMONIC     700.0000    1.0220
 N_3      H___A       HARMONIC     700.0000    1.0220
 N_3      H___b       HARMONIC     700.0000    1.2020
 N_3      B_3         HARMONIC     700.0000    1.5720
 N_3      B_2         HARMONIC     700.0000    1.4820
 N_3      C_34        HARMONIC     700.0000    1.4620
 N_3      C_33        HARMONIC     700.0000    1.4620
 N_3      C_32        HARMONIC     700.0000    1.4620
 N_3      C_31        HARMONIC     700.0000    1.4620
 N_3      C_3         HARMONIC     700.0000    1.4620
 N_3      C_22        HARMONIC     700.0000    1.3620
 N_3      C_21        HARMONIC     700.0000    1.3620
 N_3      C_2         HARMONIC     700.0000    1.3620
 N_3      C_R2        HARMONIC     700.0000    1.3920
 N_3      C_R1        HARMONIC     700.0000    1.3920
 N_3      C_R         HARMONIC     700.0000    1.3920
 N_3      C_11        HARMONIC     700.0000    1.2940
 N_3      C_1         HARMONIC     700.0000    1.2940
 N_3      N_3         HARMONIC     700.0000    1.3940
 N_2      H_          HARMONIC     700.0000    0.9350
 N_2      H___A       HARMONIC     700.0000    0.9350
 N_2      H___b       HARMONIC     700.0000    1.1150
 N_2      B_3         HARMONIC     700.0000    1.4850
 N_2      B_2         HARMONIC    1400.0000    1.3950
 N_2      C_34        HARMONIC     700.0000    1.3750
 N_2      C_33        HARMONIC     700.0000    1.3750
 N_2      C_32        HARMONIC     700.0000    1.3750
 N_2      C_31        HARMONIC     700.0000    1.3750
 N_2      C_3         HARMONIC     700.0000    1.3750
 N_2      C_22        HARMONIC    1400.0000    1.2750
 N_2      C_21        HARMONIC    1400.0000    1.2750
 N_2      C_2         HARMONIC    1400.0000    1.2750
 N_2      C_R2        HARMONIC    1400.0000    1.3050
 N_2      C_R1        HARMONIC    1400.0000    1.3050
 N_2      C_R         HARMONIC    1400.0000    1.3050
 N_2      C_11        HARMONIC     700.0000    1.2070
 N_2      C_1         HARMONIC     700.0000    1.2070
 N_2      N_3         HARMONIC     700.0000    1.3070
 N_2      N_2         HARMONIC    1400.0000    1.2200
 N_R      H_          HARMONIC     700.0000    0.9700
 N_R      H___A       HARMONIC     700.0000    0.9700
 N_R      H___b       HARMONIC     700.0000    1.1500
 N_R      B_3         HARMONIC     700.0000    1.5200
 N_R      B_2         HARMONIC    1400.0000    1.4300
 N_R      C_34        HARMONIC     700.0000    1.4100
 N_R      C_33        HARMONIC     700.0000    1.4100
 N_R      C_32        HARMONIC     700.0000    1.4100
 N_R      C_31        HARMONIC     700.0000    1.4100
 N_R      C_3         HARMONIC     700.0000    1.4100
 N_R      C_22        HARMONIC    1400.0000    1.3100
 N_R      C_21        HARMONIC    1400.0000    1.3100
 N_R      C_2         HARMONIC    1400.0000    1.3100
 N_R      C_R2        HARMONIC    1050.0000    1.3400
 N_R      C_R1        HARMONIC    1050.0000    1.3400
 N_R      C_R         HARMONIC    1050.0000    1.3400
 N_R      C_11        HARMONIC     700.0000    1.2420
 N_R      C_1         HARMONIC     700.0000    1.2420
 N_R      N_3         HARMONIC     700.0000    1.3420
 N_R      N_2         HARMONIC    1400.0000    1.2550
 N_R      N_R         HARMONIC    1050.0000    1.2900
 N_1      H_          HARMONIC     700.0000    0.8760
 N_1      H___A       HARMONIC     700.0000    0.8760
 N_1      H___b       HARMONIC     700.0000    1.0560
 N_1      B_3         HARMONIC     700.0000    1.4260
 N_1      B_2         HARMONIC     700.0000    1.3360
 N_1      C_34        HARMONIC     700.0000    1.3160
 N_1      C_33        HARMONIC     700.0000    1.3160
 N_1      C_32        HARMONIC     700.0000    1.3160
 N_1      C_31        HARMONIC     700.0000    1.3160
 N_1      C_3         HARMONIC     700.0000    1.3160
 N_1      C_22        HARMONIC     700.0000    1.2160
 N_1      C_21        HARMONIC     700.0000    1.2160
 N_1      C_2         HARMONIC     700.0000    1.2160
 N_1      C_R2        HARMONIC     700.0000    1.2460
 N_1      C_R1        HARMONIC     700.0000    1.2460
 N_1      C_R         HARMONIC     700.0000    1.2460
 N_1      C_11        HARMONIC    2100.0000    1.1480
 N_1      C_1         HARMONIC    2100.0000    1.1480
 N_1      N_3         HARMONIC     700.0000    1.2480
 N_1      N_2         HARMONIC     700.0000    1.1610
 N_1      N_R         HARMONIC     700.0000    1.1960
 N_1      N_1         HARMONIC    2100.0000    1.1020
 O_3      H_          HARMONIC     700.0000    0.9800
 O_3      H___A       HARMONIC     700.0000    0.9800
 O_3      H___b       HARMONIC     700.0000    1.1600
 O_3      B_3         HARMONIC     700.0000    1.5300
 O_3      B_2         HARMONIC     700.0000    1.4400
 O_3      C_34        HARMONIC     700.0000    1.4200
 O_3      C_33        HARMONIC     700.0000    1.4200
 O_3      C_32        HARMONIC     700.0000    1.4200
 O_3      C_31        HARMONIC     700.0000    1.4200
 O_3      C_3         HARMONIC     700.0000    1.4200
 O_3      C_22        HARMONIC     700.0000    1.3200
 O_3      C_21        HARMONIC     700.0000    1.3200
 O_3      C_2         HARMONIC     700.0000    1.3200
 O_3      C_R2        HARMONIC     700.0000    1.3500
 O_3      C_R1        HARMONIC     700.0000    1.3500
 O_3      C_R         HARMONIC     700.0000    1.3500
 O_3      C_11        HARMONIC     700.0000    1.2520
 O_3      C_1         HARMONIC     700.0000    1.2520
 O_3      N_3         HARMONIC     700.0000    1.3520
 O_3      N_2         HARMONIC     700.0000    1.2650
 O_3      N_R         HARMONIC     700.0000    1.3000
 O_3      N_1         HARMONIC     700.0000    1.2060
 O_3      O_3         HARMONIC     700.0000    1.3100
 O_2      H_          HARMONIC     700.0000    0.8800
 O_2      H___A       HARMONIC     700.0000    0.8800
 O_2      H___b       HARMONIC     700.0000    1.0600
 O_2      B_3         HARMONIC     700.0000    1.4300
 O_2      B_2         HARMONIC    1400.0000    1.3400
 O_2      C_34        HARMONIC     700.0000    1.3200
 O_2      C_33        HARMONIC     700.0000    1.3200
 O_2      C_32        HARMONIC     700.0000    1.3200
 O_2      C_31        HARMONIC     700.0000    1.3200
 O_2      C_3         HARMONIC     700.0000    1.3200
 O_2      C_22        HARMONIC    1400.0000    1.2200
 O_2      C_21        HARMONIC    1400.0000    1.2200
 O_2      C_2         HARMONIC    1400.0000    1.2200
 O_2      C_R2        HARMONIC    1400.0000    1.2500
 O_2      C_R1        HARMONIC    1400.0000    1.2500
 O_2      C_R         HARMONIC    1400.0000    1.2500
 O_2      C_11        HARMONIC    1400.0000    1.1520
 O_2      C_1         HARMONIC    1283.38      1.1490
 O_2      N_3         HARMONIC     700.0000    1.2520
 O_2      N_2         HARMONIC    1400.0000    1.1650
 O_2      N_R         HARMONIC    1400.0000    1.2000
 O_2      N_1         HARMONIC     700.0000    1.1060
 O_2      O_3         HARMONIC     700.0000    1.2100
 O_2      O_2         HARMONIC    1400.0000    1.1100
 O_R      H_          HARMONIC     700.0000    0.9800
 O_R      H___A       HARMONIC     700.0000    0.9800
 O_R      H___b       HARMONIC     700.0000    1.1600
 O_R      B_3         HARMONIC     700.0000    1.5300
 O_R      B_2         HARMONIC    1400.0000    1.4400
 O_R      C_34        HARMONIC     700.0000    1.4200
 O_R      C_33        HARMONIC     700.0000    1.4200
 O_R      C_32        HARMONIC     700.0000    1.4200
 O_R      C_31        HARMONIC     700.0000    1.4200
 O_R      C_3         HARMONIC     700.0000    1.4200
 O_R      C_22        HARMONIC    1400.0000    1.3200
 O_R      C_21        HARMONIC    1400.0000    1.3200
 O_R      C_2         HARMONIC    1400.0000    1.3200
 O_R      C_R2        HARMONIC    1050.0000    1.3500
 O_R      C_R1        HARMONIC    1050.0000    1.3500
 O_R      C_R         HARMONIC    1050.0000    1.3500
 O_R      C_11        HARMONIC     700.0000    1.2520
 O_R      C_1         HARMONIC     700.0000    1.2520
 O_R      N_3         HARMONIC     700.0000    1.3520
 O_R      N_2         HARMONIC    1400.0000    1.2650
 O_R      N_R         HARMONIC    1050.0000    1.3000
 O_R      N_1         HARMONIC     700.0000    1.2060
 O_R      O_3         HARMONIC     700.0000    1.3100
 O_R      O_2         HARMONIC    1400.0000    1.2100
 O_R      O_R         HARMONIC    1050.0000    1.3100
 O_1      H_          HARMONIC     700.0000    0.8480
 O_1      H___A       HARMONIC     700.0000    0.8480
 O_1      H___b       HARMONIC     700.0000    1.0280
 O_1      B_3         HARMONIC     700.0000    1.3980
 O_1      B_2         HARMONIC     700.0000    1.3080
 O_1      C_34        HARMONIC     700.0000    1.2880
 O_1      C_33        HARMONIC     700.0000    1.2880
 O_1      C_32        HARMONIC     700.0000    1.2880
 O_1      C_31        HARMONIC     700.0000    1.2880
 O_1      C_3         HARMONIC     700.0000    1.2880
 O_1      C_22        HARMONIC     700.0000    1.1880
 O_1      C_21        HARMONIC     700.0000    1.1880
 O_1      C_2         HARMONIC     700.0000    1.1880
 O_1      C_R2        HARMONIC     700.0000    1.2180
 O_1      C_R1        HARMONIC     700.0000    1.2180
 O_1      C_R         HARMONIC     700.0000    1.2180
 O_1      C_11        HARMONIC    2100.0000    1.1200
 O_1      C_1         HARMONIC    2100.0000    1.1200
 O_1      N_3         HARMONIC     700.0000    1.2200
 O_1      N_2         HARMONIC     700.0000    1.1330
 O_1      N_R         HARMONIC     700.0000    1.1680
 O_1      N_1         HARMONIC    2100.0000    1.0740
 O_1      O_3         HARMONIC     700.0000    1.1780
 O_1      O_2         HARMONIC     700.0000    1.0780
 O_1      O_R         HARMONIC     700.0000    1.1780
 O_1      O_1         HARMONIC    2100.0000    1.0460
 F_       H_          HARMONIC     700.0000    0.9310
 F_       H___A       HARMONIC     700.0000    0.9310
 F_       H___b       HARMONIC     700.0000    1.1110
 F_       B_3         HARMONIC     700.0000    1.4810
 F_       B_2         HARMONIC     700.0000    1.3910
 F_       C_34        HARMONIC     700.0000    1.3710
 F_       C_33        HARMONIC     700.0000    1.3710
 F_       C_32        HARMONIC     700.0000    1.3710
 F_       C_31        HARMONIC     700.0000    1.3710
 F_       C_3         HARMONIC     700.0000    1.3710
 F_       C_22        HARMONIC     700.0000    1.2710
 F_       C_21        HARMONIC     700.0000    1.2710
 F_       C_2         HARMONIC     700.0000    1.2710
 F_       C_R2        HARMONIC     700.0000    1.3010
 F_       C_R1        HARMONIC     700.0000    1.3010
 F_       C_R         HARMONIC     700.0000    1.3010
 F_       C_11        HARMONIC     700.0000    1.2030
 F_       C_1         HARMONIC     700.0000    1.2030
 F_       N_3         HARMONIC     700.0000    1.3030
 F_       N_2         HARMONIC     700.0000    1.2160
 F_       N_R         HARMONIC     700.0000    1.2510
 F_       N_1         HARMONIC     700.0000    1.1570
 F_       O_3         HARMONIC     700.0000    1.2610
 F_       O_2         HARMONIC     700.0000    1.1610
 F_       O_R         HARMONIC     700.0000    1.2610
 F_       O_1         HARMONIC     700.0000    1.1290
 F_       F_          HARMONIC     700.0000    1.2120
 Al3      H_          HARMONIC     700.0000    1.3670
 Al3      H___A       HARMONIC     700.0000    1.3670
 Al3      H___b       HARMONIC     700.0000    1.5470
 Al3      B_3         HARMONIC     700.0000    1.9170
 Al3      B_2         HARMONIC     700.0000    1.8270
 Al3      C_34        HARMONIC     700.0000    1.8070
 Al3      C_33        HARMONIC     700.0000    1.8070
 Al3      C_32        HARMONIC     700.0000    1.8070
 Al3      C_31        HARMONIC     700.0000    1.8070
 Al3      C_3         HARMONIC     700.0000    1.8070
 Al3      C_22        HARMONIC     700.0000    1.7070
 Al3      C_21        HARMONIC     700.0000    1.7070
 Al3      C_2         HARMONIC     700.0000    1.7070
 Al3      C_R2        HARMONIC     700.0000    1.7370
 Al3      C_R1        HARMONIC     700.0000    1.7370
 Al3      C_R         HARMONIC     700.0000    1.7370
 Al3      C_11        HARMONIC     700.0000    1.6390
 Al3      C_1         HARMONIC     700.0000    1.6390
 Al3      N_3         HARMONIC     700.0000    1.7390
 Al3      N_2         HARMONIC     700.0000    1.6520
 Al3      N_R         HARMONIC     700.0000    1.6870
 Al3      N_1         HARMONIC     700.0000    1.5930
 Al3      O_3         HARMONIC     700.0000    1.6970
 Al3      O_2         HARMONIC     700.0000    1.5970
 Al3      O_R         HARMONIC     700.0000    1.6970
 Al3      O_1         HARMONIC     700.0000    1.5650
 Al3      F_          HARMONIC     700.0000    1.6480
 Al3      Al3         HARMONIC     700.0000    2.0840
 Si3      H_          HARMONIC     700.0000    1.2570
 Si3      H___A       HARMONIC     700.0000    1.2570
 Si3      H___b       HARMONIC     700.0000    1.4370
 Si3      B_3         HARMONIC     700.0000    1.8070
 Si3      B_2         HARMONIC     700.0000    1.7170
 Si3      C_34        HARMONIC     700.0000    1.6970
 Si3      C_33        HARMONIC     700.0000    1.6970
 Si3      C_32        HARMONIC     700.0000    1.6970
 Si3      C_31        HARMONIC     700.0000    1.6970
 Si3      C_3         HARMONIC     700.0000    1.6970
 Si3      C_22        HARMONIC     700.0000    1.5970
 Si3      C_21        HARMONIC     700.0000    1.5970
 Si3      C_2         HARMONIC     700.0000    1.5970
 Si3      C_R2        HARMONIC     700.0000    1.6270
 Si3      C_R1        HARMONIC     700.0000    1.6270
 Si3      C_R         HARMONIC     700.0000    1.6270
 Si3      C_11        HARMONIC     700.0000    1.5290
 Si3      C_1         HARMONIC     700.0000    1.5290
 Si3      N_3         HARMONIC     700.0000    1.6290
 Si3      N_2         HARMONIC     700.0000    1.5420
 Si3      N_R         HARMONIC     700.0000    1.5770
 Si3      N_1         HARMONIC     700.0000    1.4830
 Si3      O_3         HARMONIC     700.0000    1.5870
 Si3      O_2         HARMONIC     700.0000    1.4870
 Si3      O_R         HARMONIC     700.0000    1.5870
 Si3      O_1         HARMONIC     700.0000    1.4550
 Si3      F_          HARMONIC     700.0000    1.5380
 Si3      Al3         HARMONIC     700.0000    1.9740
 Si3      Si3         HARMONIC     700.0000    1.8640
 P_3      H_          HARMONIC     700.0000    1.2100
 P_3      H___A       HARMONIC     700.0000    1.2100
 P_3      H___b       HARMONIC     700.0000    1.3900
 P_3      B_3         HARMONIC     700.0000    1.7600
 P_3      B_2         HARMONIC     700.0000    1.6700
 P_3      C_34        HARMONIC     700.0000    1.6500
 P_3      C_33        HARMONIC     700.0000    1.6500
 P_3      C_32        HARMONIC     700.0000    1.6500
 P_3      C_31        HARMONIC     700.0000    1.6500
 P_3      C_3         HARMONIC     700.0000    1.6500
 P_3      C_22        HARMONIC     700.0000    1.5500
 P_3      C_21        HARMONIC     700.0000    1.5500
 P_3      C_2         HARMONIC     700.0000    1.5500
 P_3      C_R2        HARMONIC     700.0000    1.5800
 P_3      C_R1        HARMONIC     700.0000    1.5800
 P_3      C_R         HARMONIC     700.0000    1.5800
 P_3      C_11        HARMONIC     700.0000    1.4820
 P_3      C_1         HARMONIC     700.0000    1.4820
 P_3      N_3         HARMONIC     700.0000    1.5820
 P_3      N_2         HARMONIC     700.0000    1.4950
 P_3      N_R         HARMONIC     700.0000    1.5300
 P_3      N_1         HARMONIC     700.0000    1.4360
 P_3      O_3         HARMONIC     700.0000    1.5400
 P_3      O_2         HARMONIC     700.0000    1.4400
 P_3      O_R         HARMONIC     700.0000    1.5400
 P_3      O_1         HARMONIC     700.0000    1.4080
 P_3      F_          HARMONIC     700.0000    1.4910
 P_3      Al3         HARMONIC     700.0000    1.9270
 P_3      Si3         HARMONIC     700.0000    1.8170
 P_3      P_3         HARMONIC     700.0000    1.7700
 S_3      H_          HARMONIC     700.0000    1.3600
 S_3      H___A       HARMONIC     700.0000    1.3600
 S_3      H___b       HARMONIC     700.0000    1.5400
 S_3      B_3         HARMONIC     700.0000    1.9100
 S_3      B_2         HARMONIC     700.0000    1.8200
 S_3      C_34        HARMONIC     700.0000    1.8000
 S_3      C_33        HARMONIC     700.0000    1.8000
 S_3      C_32        HARMONIC     700.0000    1.8000
 S_3      C_31        HARMONIC     700.0000    1.8000
 S_3      C_3         HARMONIC     700.0000    1.8000
 S_3      C_22        HARMONIC     700.0000    1.7000
 S_3      C_21        HARMONIC     700.0000    1.7000
 S_3      C_2         HARMONIC     700.0000    1.7000
 S_3      C_R2        HARMONIC     700.0000    1.7300
 S_3      C_R1        HARMONIC     700.0000    1.7300
 S_3      C_R         HARMONIC     700.0000    1.7300
 S_3      C_11        HARMONIC     700.0000    1.6320
 S_3      C_1         HARMONIC     700.0000    1.6320
 S_3      N_3         HARMONIC     700.0000    1.7320
 S_3      N_2         HARMONIC     700.0000    1.6450
 S_3      N_R         HARMONIC     700.0000    1.6800
 S_3      N_1         HARMONIC     700.0000    1.5860
 S_3      O_3         HARMONIC     700.0000    1.6900
 S_3      O_2         HARMONIC     700.0000    1.5900
 S_3      O_R         HARMONIC     700.0000    1.6900
 S_3      O_1         HARMONIC     700.0000    1.5580
 S_3      F_          HARMONIC     700.0000    1.6410
 S_3      Al3         HARMONIC     700.0000    2.0770
 S_3      Si3         HARMONIC     700.0000    1.9670
 S_3      P_3         HARMONIC     700.0000    1.9200
 S_3      S_3         HARMONIC     700.0000    2.0700
 Cl       H_          HARMONIC     700.0000    1.3170
 Cl       H___A       HARMONIC     700.0000    1.3170
 Cl       H___b       HARMONIC     700.0000    1.4970
 Cl       B_3         HARMONIC     700.0000    1.8670
 Cl       B_2         HARMONIC     700.0000    1.7770
 Cl       C_34        HARMONIC     700.0000    1.7570
 Cl       C_33        HARMONIC     700.0000    1.7570
 Cl       C_32        HARMONIC     700.0000    1.7570
 Cl       C_31        HARMONIC     700.0000    1.7570
 Cl       C_3         HARMONIC     700.0000    1.7570
 Cl       C_22        HARMONIC     700.0000    1.6570
 Cl       C_21        HARMONIC     700.0000    1.6570
 Cl       C_2         HARMONIC     700.0000    1.6570
 Cl       C_R2        HARMONIC     700.0000    1.6870
 Cl       C_R1        HARMONIC     700.0000    1.6870
 Cl       C_R         HARMONIC     700.0000    1.6870
 Cl       C_11        HARMONIC     700.0000    1.5890
 Cl       C_1         HARMONIC     700.0000    1.5890
 Cl       N_3         HARMONIC     700.0000    1.6890
 Cl       N_2         HARMONIC     700.0000    1.6020
 Cl       N_R         HARMONIC     700.0000    1.6370
 Cl       N_1         HARMONIC     700.0000    1.5430
 Cl       O_3         HARMONIC     700.0000    1.6470
 Cl       O_2         HARMONIC     700.0000    1.5470
 Cl       O_R         HARMONIC     700.0000    1.6470
 Cl       O_1         HARMONIC     700.0000    1.5150
 Cl       F_          HARMONIC     700.0000    1.5980
 Cl       Al3         HARMONIC     700.0000    2.0340
 Cl       Si3         HARMONIC     700.0000    1.9240
 Cl       P_3         HARMONIC     700.0000    1.8770
 Cl       S_3         HARMONIC     700.0000    2.0270
 Cl       Cl          HARMONIC     700.0000    1.9840
 Ga3      H_          HARMONIC     700.0000    1.5300
 Ga3      H___A       HARMONIC     700.0000    1.5300
 Ga3      H___b       HARMONIC     700.0000    1.7100
 Ga3      B_3         HARMONIC     700.0000    2.0800
 Ga3      B_2         HARMONIC     700.0000    1.9900
 Ga3      C_34        HARMONIC     700.0000    1.9700
 Ga3      C_33        HARMONIC     700.0000    1.9700
 Ga3      C_32        HARMONIC     700.0000    1.9700
 Ga3      C_31        HARMONIC     700.0000    1.9700
 Ga3      C_3         HARMONIC     700.0000    1.9700
 Ga3      C_22        HARMONIC     700.0000    1.8700
 Ga3      C_21        HARMONIC     700.0000    1.8700
 Ga3      C_2         HARMONIC     700.0000    1.8700
 Ga3      C_R2        HARMONIC     700.0000    1.9000
 Ga3      C_R1        HARMONIC     700.0000    1.9000
 Ga3      C_R         HARMONIC     700.0000    1.9000
 Ga3      C_11        HARMONIC     700.0000    1.8020
 Ga3      C_1         HARMONIC     700.0000    1.8020
 Ga3      N_3         HARMONIC     700.0000    1.9020
 Ga3      N_2         HARMONIC     700.0000    1.8150
 Ga3      N_R         HARMONIC     700.0000    1.8500
 Ga3      N_1         HARMONIC     700.0000    1.7560
 Ga3      O_3         HARMONIC     700.0000    1.8600
 Ga3      O_2         HARMONIC     700.0000    1.7600
 Ga3      O_R         HARMONIC     700.0000    1.8600
 Ga3      O_1         HARMONIC     700.0000    1.7280
 Ga3      F_          HARMONIC     700.0000    1.8110
 Ga3      Al3         HARMONIC     700.0000    2.2470
 Ga3      Si3         HARMONIC     700.0000    2.1370
 Ga3      P_3         HARMONIC     700.0000    2.0900
 Ga3      S_3         HARMONIC     700.0000    2.2400
 Ga3      Cl          HARMONIC     700.0000    2.1970
 Ga3      Ga3         HARMONIC     700.0000    2.4100
 Ge3      H_          HARMONIC     700.0000    1.5300
 Ge3      H___A       HARMONIC     700.0000    1.5300
 Ge3      H___b       HARMONIC     700.0000    1.7100
 Ge3      B_3         HARMONIC     700.0000    2.0800
 Ge3      B_2         HARMONIC     700.0000    1.9900
 Ge3      C_34        HARMONIC     700.0000    1.9700
 Ge3      C_33        HARMONIC     700.0000    1.9700
 Ge3      C_32        HARMONIC     700.0000    1.9700
 Ge3      C_31        HARMONIC     700.0000    1.9700
 Ge3      C_3         HARMONIC     700.0000    1.9700
 Ge3      C_22        HARMONIC     700.0000    1.8700
 Ge3      C_21        HARMONIC     700.0000    1.8700
 Ge3      C_2         HARMONIC     700.0000    1.8700
 Ge3      C_R2        HARMONIC     700.0000    1.9000
 Ge3      C_R1        HARMONIC     700.0000    1.9000
 Ge3      C_R         HARMONIC     700.0000    1.9000
 Ge3      C_11        HARMONIC     700.0000    1.8020
 Ge3      C_1         HARMONIC     700.0000    1.8020
 Ge3      N_3         HARMONIC     700.0000    1.9020
 Ge3      N_2         HARMONIC     700.0000    1.8150
 Ge3      N_R         HARMONIC     700.0000    1.8500
 Ge3      N_1         HARMONIC     700.0000    1.7560
 Ge3      O_3         HARMONIC     700.0000    1.8600
 Ge3      O_2         HARMONIC     700.0000    1.7600
 Ge3      O_R         HARMONIC     700.0000    1.8600
 Ge3      O_1         HARMONIC     700.0000    1.7280
 Ge3      F_          HARMONIC     700.0000    1.8110
 Ge3      Al3         HARMONIC     700.0000    2.2470
 Ge3      Si3         HARMONIC     700.0000    2.1370
 Ge3      P_3         HARMONIC     700.0000    2.0900
 Ge3      S_3         HARMONIC     700.0000    2.2400
 Ge3      Cl          HARMONIC     700.0000    2.1970
 Ge3      Ga3         HARMONIC     700.0000    2.4100
 Ge3      Ge3         HARMONIC     700.0000    2.4100
 As3      H_          HARMONIC     700.0000    1.5300
 As3      H___A       HARMONIC     700.0000    1.5300
 As3      H___b       HARMONIC     700.0000    1.7100
 As3      B_3         HARMONIC     700.0000    2.0800
 As3      B_2         HARMONIC     700.0000    1.9900
 As3      C_34        HARMONIC     700.0000    1.9700
 As3      C_33        HARMONIC     700.0000    1.9700
 As3      C_32        HARMONIC     700.0000    1.9700
 As3      C_31        HARMONIC     700.0000    1.9700
 As3      C_3         HARMONIC     700.0000    1.9700
 As3      C_22        HARMONIC     700.0000    1.8700
 As3      C_21        HARMONIC     700.0000    1.8700
 As3      C_2         HARMONIC     700.0000    1.8700
 As3      C_R2        HARMONIC     700.0000    1.9000
 As3      C_R1        HARMONIC     700.0000    1.9000
 As3      C_R         HARMONIC     700.0000    1.9000
 As3      C_11        HARMONIC     700.0000    1.8020
 As3      C_1         HARMONIC     700.0000    1.8020
 As3      N_3         HARMONIC     700.0000    1.9020
 As3      N_2         HARMONIC     700.0000    1.8150
 As3      N_R         HARMONIC     700.0000    1.8500
 As3      N_1         HARMONIC     700.0000    1.7560
 As3      O_3         HARMONIC     700.0000    1.8600
 As3      O_2         HARMONIC     700.0000    1.7600
 As3      O_R         HARMONIC     700.0000    1.8600
 As3      O_1         HARMONIC     700.0000    1.7280
 As3      F_          HARMONIC     700.0000    1.8110
 As3      Al3         HARMONIC     700.0000    2.2470
 As3      Si3         HARMONIC     700.0000    2.1370
 As3      P_3         HARMONIC     700.0000    2.0900
 As3      S_3         HARMONIC     700.0000    2.2400
 As3      Cl          HARMONIC     700.0000    2.1970
 As3      Ga3         HARMONIC     700.0000    2.4100
 As3      Ge3         HARMONIC     700.0000    2.4100
 As3      As3         HARMONIC     700.0000    2.4100
 Se3      H_          HARMONIC     700.0000    1.5300
 Se3      H___A       HARMONIC     700.0000    1.5300
 Se3      H___b       HARMONIC     700.0000    1.7100
 Se3      B_3         HARMONIC     700.0000    2.0800
 Se3      B_2         HARMONIC     700.0000    1.9900
 Se3      C_34        HARMONIC     700.0000    1.9700
 Se3      C_33        HARMONIC     700.0000    1.9700
 Se3      C_32        HARMONIC     700.0000    1.9700
 Se3      C_31        HARMONIC     700.0000    1.9700
 Se3      C_3         HARMONIC     700.0000    1.9700
 Se3      C_22        HARMONIC     700.0000    1.8700
 Se3      C_21        HARMONIC     700.0000    1.8700
 Se3      C_2         HARMONIC     700.0000    1.8700
 Se3      C_R2        HARMONIC     700.0000    1.9000
 Se3      C_R1        HARMONIC     700.0000    1.9000
 Se3      C_R         HARMONIC     700.0000    1.9000
 Se3      C_11        HARMONIC     700.0000    1.8020
 Se3      C_1         HARMONIC     700.0000    1.8020
 Se3      N_3         HARMONIC     700.0000    1.9020
 Se3      N_2         HARMONIC     700.0000    1.8150
 Se3      N_R         HARMONIC     700.0000    1.8500
 Se3      N_1         HARMONIC     700.0000    1.7560
 Se3      O_3         HARMONIC     700.0000    1.8600
 Se3      O_2         HARMONIC     700.0000    1.7600
 Se3      O_R         HARMONIC     700.0000    1.8600
 Se3      O_1         HARMONIC     700.0000    1.7280
 Se3      F_          HARMONIC     700.0000    1.8110
 Se3      Al3         HARMONIC     700.0000    2.2470
 Se3      Si3         HARMONIC     700.0000    2.1370
 Se3      P_3         HARMONIC     700.0000    2.0900
 Se3      S_3         HARMONIC     700.0000    2.2400
 Se3      Cl          HARMONIC     700.0000    2.1970
 Se3      Ga3         HARMONIC     700.0000    2.4100
 Se3      Ge3         HARMONIC     700.0000    2.4100
 Se3      As3         HARMONIC     700.0000    2.4100
 Se3      Se3         HARMONIC     700.0000    2.4100
 Br       H_          HARMONIC     700.0000    1.4870
 Br       H___A       HARMONIC     700.0000    1.4870
 Br       H___b       HARMONIC     700.0000    1.6670
 Br       B_3         HARMONIC     700.0000    2.0370
 Br       B_2         HARMONIC     700.0000    1.9470
 Br       C_34        HARMONIC     700.0000    1.9270
 Br       C_33        HARMONIC     700.0000    1.9270
 Br       C_32        HARMONIC     700.0000    1.9270
 Br       C_31        HARMONIC     700.0000    1.9270
 Br       C_3         HARMONIC     700.0000    1.9270
 Br       C_22        HARMONIC     700.0000    1.8270
 Br       C_21        HARMONIC     700.0000    1.8270
 Br       C_2         HARMONIC     700.0000    1.8270
 Br       C_R2        HARMONIC     700.0000    1.8570
 Br       C_R1        HARMONIC     700.0000    1.8570
 Br       C_R         HARMONIC     700.0000    1.8570
 Br       C_11        HARMONIC     700.0000    1.7590
 Br       C_1         HARMONIC     700.0000    1.7590
 Br       N_3         HARMONIC     700.0000    1.8590
 Br       N_2         HARMONIC     700.0000    1.7720
 Br       N_R         HARMONIC     700.0000    1.8070
 Br       N_1         HARMONIC     700.0000    1.7130
 Br       O_3         HARMONIC     700.0000    1.8170
 Br       O_2         HARMONIC     700.0000    1.7170
 Br       O_R         HARMONIC     700.0000    1.8170
 Br       O_1         HARMONIC     700.0000    1.6850
 Br       F_          HARMONIC     700.0000    1.7680
 Br       Al3         HARMONIC     700.0000    2.2040
 Br       Si3         HARMONIC     700.0000    2.0940
 Br       P_3         HARMONIC     700.0000    2.0470
 Br       S_3         HARMONIC     700.0000    2.1970
 Br       Cl          HARMONIC     700.0000    2.1540
 Br       Ga3         HARMONIC     700.0000    2.3670
 Br       Ge3         HARMONIC     700.0000    2.3670
 Br       As3         HARMONIC     700.0000    2.3670
 Br       Se3         HARMONIC     700.0000    2.3670
 Br       Br          HARMONIC     700.0000    2.3240
 In3      H_          HARMONIC     700.0000    1.7100
 In3      H___A       HARMONIC     700.0000    1.7100
 In3      H___b       HARMONIC     700.0000    1.8900
 In3      B_3         HARMONIC     700.0000    2.2600
 In3      B_2         HARMONIC     700.0000    2.1700
 In3      C_34        HARMONIC     700.0000    2.1500
 In3      C_33        HARMONIC     700.0000    2.1500
 In3      C_32        HARMONIC     700.0000    2.1500
 In3      C_31        HARMONIC     700.0000    2.1500
 In3      C_3         HARMONIC     700.0000    2.1500
 In3      C_22        HARMONIC     700.0000    2.0500
 In3      C_21        HARMONIC     700.0000    2.0500
 In3      C_2         HARMONIC     700.0000    2.0500
 In3      C_R2        HARMONIC     700.0000    2.0800
 In3      C_R1        HARMONIC     700.0000    2.0800
 In3      C_R         HARMONIC     700.0000    2.0800
 In3      C_11        HARMONIC     700.0000    1.9820
 In3      C_1         HARMONIC     700.0000    1.9820
 In3      N_3         HARMONIC     700.0000    2.0820
 In3      N_2         HARMONIC     700.0000    1.9950
 In3      N_R         HARMONIC     700.0000    2.0300
 In3      N_1         HARMONIC     700.0000    1.9360
 In3      O_3         HARMONIC     700.0000    2.0400
 In3      O_2         HARMONIC     700.0000    1.9400
 In3      O_R         HARMONIC     700.0000    2.0400
 In3      O_1         HARMONIC     700.0000    1.9080
 In3      F_          HARMONIC     700.0000    1.9910
 In3      Al3         HARMONIC     700.0000    2.4270
 In3      Si3         HARMONIC     700.0000    2.3170
 In3      P_3         HARMONIC     700.0000    2.2700
 In3      S_3         HARMONIC     700.0000    2.4200
 In3      Cl          HARMONIC     700.0000    2.3770
 In3      Ga3         HARMONIC     700.0000    2.5900
 In3      Ge3         HARMONIC     700.0000    2.5900
 In3      As3         HARMONIC     700.0000    2.5900
 In3      Se3         HARMONIC     700.0000    2.5900
 In3      Br          HARMONIC     700.0000    2.5470
 In3      In3         HARMONIC     700.0000    2.7700
 Sn3      H_          HARMONIC     700.0000    1.6930
 Sn3      H___A       HARMONIC     700.0000    1.6930
 Sn3      H___b       HARMONIC     700.0000    1.8730
 Sn3      B_3         HARMONIC     700.0000    2.2430
 Sn3      B_2         HARMONIC     700.0000    2.1530
 Sn3      C_34        HARMONIC     700.0000    2.1330
 Sn3      C_33        HARMONIC     700.0000    2.1330
 Sn3      C_32        HARMONIC     700.0000    2.1330
 Sn3      C_31        HARMONIC     700.0000    2.1330
 Sn3      C_3         HARMONIC     700.0000    2.1330
 Sn3      C_22        HARMONIC     700.0000    2.0330
 Sn3      C_21        HARMONIC     700.0000    2.0330
 Sn3      C_2         HARMONIC     700.0000    2.0330
 Sn3      C_R2        HARMONIC     700.0000    2.0630
 Sn3      C_R1        HARMONIC     700.0000    2.0630
 Sn3      C_R         HARMONIC     700.0000    2.0630
 Sn3      C_11        HARMONIC     700.0000    1.9650
 Sn3      C_1         HARMONIC     700.0000    1.9650
 Sn3      N_3         HARMONIC     700.0000    2.0650
 Sn3      N_2         HARMONIC     700.0000    1.9780
 Sn3      N_R         HARMONIC     700.0000    2.0130
 Sn3      N_1         HARMONIC     700.0000    1.9190
 Sn3      O_3         HARMONIC     700.0000    2.0230
 Sn3      O_2         HARMONIC     700.0000    1.9230
 Sn3      O_R         HARMONIC     700.0000    2.0230
 Sn3      O_1         HARMONIC     700.0000    1.8910
 Sn3      F_          HARMONIC     700.0000    1.9740
 Sn3      Al3         HARMONIC     700.0000    2.4100
 Sn3      Si3         HARMONIC     700.0000    2.3000
 Sn3      P_3         HARMONIC     700.0000    2.2530
 Sn3      S_3         HARMONIC     700.0000    2.4030
 Sn3      Cl          HARMONIC     700.0000    2.3600
 Sn3      Ga3         HARMONIC     700.0000    2.5730
 Sn3      Ge3         HARMONIC     700.0000    2.5730
 Sn3      As3         HARMONIC     700.0000    2.5730
 Sn3      Se3         HARMONIC     700.0000    2.5730
 Sn3      Br          HARMONIC     700.0000    2.5300
 Sn3      In3         HARMONIC     700.0000    2.7530
 Sn3      Sn3         HARMONIC     700.0000    2.7360
 Sb3      H_          HARMONIC     700.0000    1.7520
 Sb3      H___A       HARMONIC     700.0000    1.7520
 Sb3      H___b       HARMONIC     700.0000    1.9320
 Sb3      B_3         HARMONIC     700.0000    2.3020
 Sb3      B_2         HARMONIC     700.0000    2.2120
 Sb3      C_34        HARMONIC     700.0000    2.1920
 Sb3      C_33        HARMONIC     700.0000    2.1920
 Sb3      C_32        HARMONIC     700.0000    2.1920
 Sb3      C_31        HARMONIC     700.0000    2.1920
 Sb3      C_3         HARMONIC     700.0000    2.1920
 Sb3      C_22        HARMONIC     700.0000    2.0920
 Sb3      C_21        HARMONIC     700.0000    2.0920
 Sb3      C_2         HARMONIC     700.0000    2.0920
 Sb3      C_R2        HARMONIC     700.0000    2.1220
 Sb3      C_R1        HARMONIC     700.0000    2.1220
 Sb3      C_R         HARMONIC     700.0000    2.1220
 Sb3      C_11        HARMONIC     700.0000    2.0240
 Sb3      C_1         HARMONIC     700.0000    2.0240
 Sb3      N_3         HARMONIC     700.0000    2.1240
 Sb3      N_2         HARMONIC     700.0000    2.0370
 Sb3      N_R         HARMONIC     700.0000    2.0720
 Sb3      N_1         HARMONIC     700.0000    1.9780
 Sb3      O_3         HARMONIC     700.0000    2.0820
 Sb3      O_2         HARMONIC     700.0000    1.9820
 Sb3      O_R         HARMONIC     700.0000    2.0820
 Sb3      O_1         HARMONIC     700.0000    1.9500
 Sb3      F_          HARMONIC     700.0000    2.0330
 Sb3      Al3         HARMONIC     700.0000    2.4690
 Sb3      Si3         HARMONIC     700.0000    2.3590
 Sb3      P_3         HARMONIC     700.0000    2.3120
 Sb3      S_3         HARMONIC     700.0000    2.4620
 Sb3      Cl          HARMONIC     700.0000    2.4190
 Sb3      Ga3         HARMONIC     700.0000    2.6320
 Sb3      Ge3         HARMONIC     700.0000    2.6320
 Sb3      As3         HARMONIC     700.0000    2.6320
 Sb3      Se3         HARMONIC     700.0000    2.6320
 Sb3      Br          HARMONIC     700.0000    2.5890
 Sb3      In3         HARMONIC     700.0000    2.8120
 Sb3      Sn3         HARMONIC     700.0000    2.7950
 Sb3      Sb3         HARMONIC     700.0000    2.8540
 Te3      H_          HARMONIC     700.0000    1.6000
 Te3      H___A       HARMONIC     700.0000    1.6000
 Te3      H___b       HARMONIC     700.0000    1.7800
 Te3      B_3         HARMONIC     700.0000    2.1500
 Te3      B_2         HARMONIC     700.0000    2.0600
 Te3      C_34        HARMONIC     700.0000    2.0400
 Te3      C_33        HARMONIC     700.0000    2.0400
 Te3      C_32        HARMONIC     700.0000    2.0400
 Te3      C_31        HARMONIC     700.0000    2.0400
 Te3      C_3         HARMONIC     700.0000    2.0400
 Te3      C_22        HARMONIC     700.0000    1.9400
 Te3      C_21        HARMONIC     700.0000    1.9400
 Te3      C_2         HARMONIC     700.0000    1.9400
 Te3      C_R2        HARMONIC     700.0000    1.9700
 Te3      C_R1        HARMONIC     700.0000    1.9700
 Te3      C_R         HARMONIC     700.0000    1.9700
 Te3      C_11        HARMONIC     700.0000    1.8720
 Te3      C_1         HARMONIC     700.0000    1.8720
 Te3      N_3         HARMONIC     700.0000    1.9720
 Te3      N_2         HARMONIC     700.0000    1.8850
 Te3      N_R         HARMONIC     700.0000    1.9200
 Te3      N_1         HARMONIC     700.0000    1.8260
 Te3      O_3         HARMONIC     700.0000    1.9300
 Te3      O_2         HARMONIC     700.0000    1.8300
 Te3      O_R         HARMONIC     700.0000    1.9300
 Te3      O_1         HARMONIC     700.0000    1.7980
 Te3      F_          HARMONIC     700.0000    1.8810
 Te3      Al3         HARMONIC     700.0000    2.3170
 Te3      Si3         HARMONIC     700.0000    2.2070
 Te3      P_3         HARMONIC     700.0000    2.1600
 Te3      S_3         HARMONIC     700.0000    2.3100
 Te3      Cl          HARMONIC     700.0000    2.2670
 Te3      Ga3         HARMONIC     700.0000    2.4800
 Te3      Ge3         HARMONIC     700.0000    2.4800
 Te3      As3         HARMONIC     700.0000    2.4800
 Te3      Se3         HARMONIC     700.0000    2.4800
 Te3      Br          HARMONIC     700.0000    2.4370
 Te3      In3         HARMONIC     700.0000    2.6600
 Te3      Sn3         HARMONIC     700.0000    2.6430
 Te3      Sb3         HARMONIC     700.0000    2.7020
 Te3      Te3         HARMONIC     700.0000    2.5500
 I_       H_          HARMONIC     700.0000    1.6800
 I_       H___A       HARMONIC     700.0000    1.6800
 I_       H___b       HARMONIC     700.0000    1.8600
 I_       B_3         HARMONIC     700.0000    2.2300
 I_       B_2         HARMONIC     700.0000    2.1400
 I_       C_34        HARMONIC     700.0000    2.1200
 I_       C_33        HARMONIC     700.0000    2.1200
 I_       C_32        HARMONIC     700.0000    2.1200
 I_       C_31        HARMONIC     700.0000    2.1200
 I_       C_3         HARMONIC     700.0000    2.1200
 I_       C_22        HARMONIC     700.0000    2.0200
 I_       C_21        HARMONIC     700.0000    2.0200
 I_       C_2         HARMONIC     700.0000    2.0200
 I_       C_R2        HARMONIC     700.0000    2.0500
 I_       C_R1        HARMONIC     700.0000    2.0500
 I_       C_R         HARMONIC     700.0000    2.0500
 I_       C_11        HARMONIC     700.0000    1.9520
 I_       C_1         HARMONIC     700.0000    1.9520
 I_       N_3         HARMONIC     700.0000    2.0520
 I_       N_2         HARMONIC     700.0000    1.9650
 I_       N_R         HARMONIC     700.0000    2.0000
 I_       N_1         HARMONIC     700.0000    1.9060
 I_       O_3         HARMONIC     700.0000    2.0100
 I_       O_2         HARMONIC     700.0000    1.9100
 I_       O_R         HARMONIC     700.0000    2.0100
 I_       O_1         HARMONIC     700.0000    1.8780
 I_       F_          HARMONIC     700.0000    1.9610
 I_       Al3         HARMONIC     700.0000    2.3970
 I_       Si3         HARMONIC     700.0000    2.2870
 I_       P_3         HARMONIC     700.0000    2.2400
 I_       S_3         HARMONIC     700.0000    2.3900
 I_       Cl          HARMONIC     700.0000    2.3470
 I_       Ga3         HARMONIC     700.0000    2.5600
 I_       Ge3         HARMONIC     700.0000    2.5600
 I_       As3         HARMONIC     700.0000    2.5600
 I_       Se3         HARMONIC     700.0000    2.5600
 I_       Br          HARMONIC     700.0000    2.5170
 I_       In3         HARMONIC     700.0000    2.7400
 I_       Sn3         HARMONIC     700.0000    2.7230
 I_       Sb3         HARMONIC     700.0000    2.7820
 I_       Te3         HARMONIC     700.0000    2.6300
 I_       I_          HARMONIC     700.0000    2.7100
 Na       H_          HARMONIC     700.0000    2.1800
 Na       H___A       HARMONIC     700.0000    2.1800
 Na       H___b       HARMONIC     700.0000    2.3600
 Na       B_3         HARMONIC     700.0000    2.7300
 Na       B_2         HARMONIC     700.0000    2.6400
 Na       C_34        HARMONIC     700.0000    2.6200
 Na       C_33        HARMONIC     700.0000    2.6200
 Na       C_32        HARMONIC     700.0000    2.6200
 Na       C_31        HARMONIC     700.0000    2.6200
 Na       C_3         HARMONIC     700.0000    2.6200
 Na       C_22        HARMONIC     700.0000    2.5200
 Na       C_21        HARMONIC     700.0000    2.5200
 Na       C_2         HARMONIC     700.0000    2.5200
 Na       C_R2        HARMONIC     700.0000    2.5500
 Na       C_R1        HARMONIC     700.0000    2.5500
 Na       C_R         HARMONIC     700.0000    2.5500
 Na       C_11        HARMONIC     700.0000    2.4520
 Na       C_1         HARMONIC     700.0000    2.4520
 Na       N_3         HARMONIC     700.0000    2.5520
 Na       N_2         HARMONIC     700.0000    2.4650
 Na       N_R         HARMONIC     700.0000    2.5000
 Na       N_1         HARMONIC     700.0000    2.4060
 Na       O_3         HARMONIC     700.0000    2.5100
 Na       O_2         HARMONIC     700.0000    2.4100
 Na       O_R         HARMONIC     700.0000    2.5100
 Na       O_1         HARMONIC     700.0000    2.3780
 Na       F_          HARMONIC     700.0000    2.4610
 Na       Al3         HARMONIC     700.0000    2.8970
 Na       Si3         HARMONIC     700.0000    2.7870
 Na       P_3         HARMONIC     700.0000    2.7400
 Na       S_3         HARMONIC     700.0000    2.8900
 Na       Cl          HARMONIC     700.0000    2.8470
 Na       Ga3         HARMONIC     700.0000    3.0600
 Na       Ge3         HARMONIC     700.0000    3.0600
 Na       As3         HARMONIC     700.0000    3.0600
 Na       Se3         HARMONIC     700.0000    3.0600
 Na       Br          HARMONIC     700.0000    3.0170
 Na       In3         HARMONIC     700.0000    3.2400
 Na       Sn3         HARMONIC     700.0000    3.2230
 Na       Sb3         HARMONIC     700.0000    3.2820
 Na       Te3         HARMONIC     700.0000    3.1300
 Na       I_          HARMONIC     700.0000    3.2100
 Na       Na          HARMONIC     700.0000    3.7100
 Ca       H_          HARMONIC     700.0000    2.2600
 Ca       H___A       HARMONIC     700.0000    2.2600
 Ca       H___b       HARMONIC     700.0000    2.4400
 Ca       B_3         HARMONIC     700.0000    2.8100
 Ca       B_2         HARMONIC     700.0000    2.7200
 Ca       C_34        HARMONIC     700.0000    2.7000
 Ca       C_33        HARMONIC     700.0000    2.7000
 Ca       C_32        HARMONIC     700.0000    2.7000
 Ca       C_31        HARMONIC     700.0000    2.7000
 Ca       C_3         HARMONIC     700.0000    2.7000
 Ca       C_22        HARMONIC     700.0000    2.6000
 Ca       C_21        HARMONIC     700.0000    2.6000
 Ca       C_2         HARMONIC     700.0000    2.6000
 Ca       C_R2        HARMONIC     700.0000    2.6300
 Ca       C_R1        HARMONIC     700.0000    2.6300
 Ca       C_R         HARMONIC     700.0000    2.6300
 Ca       C_11        HARMONIC     700.0000    2.5320
 Ca       C_1         HARMONIC     700.0000    2.5320
 Ca       N_3         HARMONIC     700.0000    2.6320
 Ca       N_2         HARMONIC     700.0000    2.5450
 Ca       N_R         HARMONIC     700.0000    2.5800
 Ca       N_1         HARMONIC     700.0000    2.4860
 Ca       O_3         HARMONIC     700.0000    2.5900
 Ca       O_2         HARMONIC     700.0000    2.4900
 Ca       O_R         HARMONIC     700.0000    2.5900
 Ca       O_1         HARMONIC     700.0000    2.4580
 Ca       F_          HARMONIC     700.0000    2.5410
 Ca       Al3         HARMONIC     700.0000    2.9770
 Ca       Si3         HARMONIC     700.0000    2.8670
 Ca       P_3         HARMONIC     700.0000    2.8200
 Ca       S_3         HARMONIC     700.0000    2.9700
 Ca       Cl          HARMONIC     700.0000    2.9270
 Ca       Ga3         HARMONIC     700.0000    3.1400
 Ca       Ge3         HARMONIC     700.0000    3.1400
 Ca       As3         HARMONIC     700.0000    3.1400
 Ca       Se3         HARMONIC     700.0000    3.1400
 Ca       Br          HARMONIC     700.0000    3.0970
 Ca       In3         HARMONIC     700.0000    3.3200
 Ca       Sn3         HARMONIC     700.0000    3.3030
 Ca       Sb3         HARMONIC     700.0000    3.3620
 Ca       Te3         HARMONIC     700.0000    3.2100
 Ca       I_          HARMONIC     700.0000    3.2900
 Ca       Na          HARMONIC     700.0000    3.7900
 Ca       Ca          HARMONIC     700.0000    3.8700
 Ti       H_          HARMONIC     700.0000    1.7500
 Ti       H___A       HARMONIC     700.0000    1.7500
 Ti       H___b       HARMONIC     700.0000    1.9300
 Ti       B_3         HARMONIC     700.0000    2.3000
 Ti       B_2         HARMONIC     700.0000    2.2100
 Ti       C_34        HARMONIC     700.0000    2.1900
 Ti       C_33        HARMONIC     700.0000    2.1900
 Ti       C_32        HARMONIC     700.0000    2.1900
 Ti       C_31        HARMONIC     700.0000    2.1900
 Ti       C_3         HARMONIC     700.0000    2.1900
 Ti       C_22        HARMONIC     700.0000    2.0900
 Ti       C_21        HARMONIC     700.0000    2.0900
 Ti       C_2         HARMONIC     700.0000    2.0900
 Ti       C_R2        HARMONIC     700.0000    2.1200
 Ti       C_R1        HARMONIC     700.0000    2.1200
 Ti       C_R         HARMONIC     700.0000    2.1200
 Ti       C_11        HARMONIC     700.0000    2.0220
 Ti       C_1         HARMONIC     700.0000    2.0220
 Ti       N_3         HARMONIC     700.0000    2.1220
 Ti       N_2         HARMONIC     700.0000    2.0350
 Ti       N_R         HARMONIC     700.0000    2.0700
 Ti       N_1         HARMONIC     700.0000    1.9760
 Ti       O_3         HARMONIC     700.0000    2.0800
 Ti       O_2         HARMONIC     700.0000    1.9800
 Ti       O_R         HARMONIC     700.0000    2.0800
 Ti       O_1         HARMONIC     700.0000    1.9480
 Ti       F_          HARMONIC     700.0000    2.0310
 Ti       Al3         HARMONIC     700.0000    2.4670
 Ti       Si3         HARMONIC     700.0000    2.3570
 Ti       P_3         HARMONIC     700.0000    2.3100
 Ti       S_3         HARMONIC     700.0000    2.4600
 Ti       Cl          HARMONIC     700.0000    2.4170
 Ti       Ga3         HARMONIC     700.0000    2.6300
 Ti       Ge3         HARMONIC     700.0000    2.6300
 Ti       As3         HARMONIC     700.0000    2.6300
 Ti       Se3         HARMONIC     700.0000    2.6300
 Ti       Br          HARMONIC     700.0000    2.5870
 Ti       In3         HARMONIC     700.0000    2.8100
 Ti       Sn3         HARMONIC     700.0000    2.7930
 Ti       Sb3         HARMONIC     700.0000    2.8520
 Ti       Te3         HARMONIC     700.0000    2.7000
 Ti       I_          HARMONIC     700.0000    2.7800
 Ti       Na          HARMONIC     700.0000    3.2800
 Ti       Ca          HARMONIC     700.0000    3.3600
 Ti       Ti          HARMONIC     700.0000    2.8500
 Fe       H_          HARMONIC     700.0000    1.6050
 Fe       H___A       HARMONIC     700.0000    1.6050
 Fe       H___b       HARMONIC     700.0000    1.7850
 Fe       B_3         HARMONIC     700.0000    2.1550
 Fe       B_2         HARMONIC     700.0000    2.0650
 Fe       C_34        HARMONIC     700.0000    2.0450
 Fe       C_33        HARMONIC     700.0000    2.0450
 Fe       C_32        HARMONIC     700.0000    2.0450
 Fe       C_31        HARMONIC     700.0000    2.0450
 Fe       C_3         HARMONIC     700.0000    2.0450
 Fe       C_22        HARMONIC     700.0000    1.9450
 Fe       C_21        HARMONIC     700.0000    1.9450
 Fe       C_2         HARMONIC     700.0000    1.9450
 Fe       C_R2        HARMONIC     700.0000    1.9750
 Fe       C_R1        HARMONIC     700.0000    1.9750
 Fe       C_R         HARMONIC     700.0000    1.9750
 Fe       C_11        HARMONIC     700.0000    1.8770
 Fe       C_1         HARMONIC     700.0000    1.8770
 Fe       N_3         HARMONIC     700.0000    1.9770
 Fe       N_2         HARMONIC     700.0000    1.8900
 Fe       N_R         HARMONIC     700.0000    1.9250
 Fe       N_1         HARMONIC     700.0000    1.8310
 Fe       O_3         HARMONIC     700.0000    1.9350
 Fe       O_2         HARMONIC     700.0000    1.8350
 Fe       O_R         HARMONIC     700.0000    1.9350
 Fe       O_1         HARMONIC     700.0000    1.8030
 Fe       F_          HARMONIC     700.0000    1.8860
 Fe       Al3         HARMONIC     700.0000    2.3220
 Fe       Si3         HARMONIC     700.0000    2.2120
 Fe       P_3         HARMONIC     700.0000    2.1650
 Fe       S_3         HARMONIC     700.0000    2.3150
 Fe       Cl          HARMONIC     700.0000    2.2720
 Fe       Ga3         HARMONIC     700.0000    2.4850
 Fe       Ge3         HARMONIC     700.0000    2.4850
 Fe       As3         HARMONIC     700.0000    2.4850
 Fe       Se3         HARMONIC     700.0000    2.4850
 Fe       Br          HARMONIC     700.0000    2.4420
 Fe       In3         HARMONIC     700.0000    2.6650
 Fe       Sn3         HARMONIC     700.0000    2.6480
 Fe       Sb3         HARMONIC     700.0000    2.7070
 Fe       Te3         HARMONIC     700.0000    2.5550
 Fe       I_          HARMONIC     700.0000    2.6350
 Fe       Na          HARMONIC     700.0000    3.1350
 Fe       Ca          HARMONIC     700.0000    3.2150
 Fe       Ti          HARMONIC     700.0000    2.7050
 Fe       Fe          HARMONIC     700.0000    2.5600
 Zn       H_          HARMONIC     700.0000    1.6500
 Zn       H___A       HARMONIC     700.0000    1.6500
 Zn       H___b       HARMONIC     700.0000    1.8300
 Zn       B_3         HARMONIC     700.0000    2.2000
 Zn       B_2         HARMONIC     700.0000    2.1100
 Zn       C_34        HARMONIC     700.0000    2.0900
 Zn       C_33        HARMONIC     700.0000    2.0900
 Zn       C_32        HARMONIC     700.0000    2.0900
 Zn       C_31        HARMONIC     700.0000    2.0900
 Zn       C_3         HARMONIC     700.0000    2.0900
 Zn       C_22        HARMONIC     700.0000    1.9900
 Zn       C_21        HARMONIC     700.0000    1.9900
 Zn       C_2         HARMONIC     700.0000    1.9900
 Zn       C_R2        HARMONIC     700.0000    2.0200
 Zn       C_R1        HARMONIC     700.0000    2.0200
 Zn       C_R         HARMONIC     700.0000    2.0200
 Zn       C_11        HARMONIC     700.0000    1.9220
 Zn       C_1         HARMONIC     700.0000    1.9220
 Zn       N_3         HARMONIC     700.0000    2.0220
 Zn       N_2         HARMONIC     700.0000    1.9350
 Zn       N_R         HARMONIC     700.0000    1.9700
 Zn       N_1         HARMONIC     700.0000    1.8760
 Zn       O_3         HARMONIC     700.0000    1.9800
 Zn       O_2         HARMONIC     700.0000    1.8800
 Zn       O_R         HARMONIC     700.0000    1.9800
 Zn       O_1         HARMONIC     700.0000    1.8480
 Zn       F_          HARMONIC     700.0000    1.9310
 Zn       Al3         HARMONIC     700.0000    2.3670
 Zn       Si3         HARMONIC     700.0000    2.2570
 Zn       P_3         HARMONIC     700.0000    2.2100
 Zn       S_3         HARMONIC     700.0000    2.3600
 Zn       Cl          HARMONIC     700.0000    2.3170
 Zn       Ga3         HARMONIC     700.0000    2.5300
 Zn       Ge3         HARMONIC     700.0000    2.5300
 Zn       As3         HARMONIC     700.0000    2.5300
 Zn       Se3         HARMONIC     700.0000    2.5300
 Zn       Br          HARMONIC     700.0000    2.4870
 Zn       In3         HARMONIC     700.0000    2.7100
 Zn       Sn3         HARMONIC     700.0000    2.6930
 Zn       Sb3         HARMONIC     700.0000    2.7520
 Zn       Te3         HARMONIC     700.0000    2.6000
 Zn       I_          HARMONIC     700.0000    2.6800
 Zn       Na          HARMONIC     700.0000    3.1800
 Zn       Ca          HARMONIC     700.0000    3.2600
 Zn       Ti          HARMONIC     700.0000    2.7500
 Zn       Fe          HARMONIC     700.0000    2.6050
 Zn       Zn          HARMONIC     700.0000    2.6500
 Tc       H_          HARMONIC     700.0000    1.6715
 Tc       H___A       HARMONIC     700.0000    1.6715
 Tc       H___b       HARMONIC     700.0000    1.8515
 Tc       B_3         HARMONIC     700.0000    2.2215
 Tc       B_2         HARMONIC     700.0000    2.1315
 Tc       C_34        HARMONIC     700.0000    2.1115
 Tc       C_33        HARMONIC     700.0000    2.1115
 Tc       C_32        HARMONIC     700.0000    2.1115
 Tc       C_31        HARMONIC     700.0000    2.1115
 Tc       C_3         HARMONIC     700.0000    2.1115
 Tc       C_22        HARMONIC     700.0000    2.0115
 Tc       C_21        HARMONIC     700.0000    2.0115
 Tc       C_2         HARMONIC     700.0000    2.0115
 Tc       C_R2        HARMONIC     700.0000    2.0415
 Tc       C_R1        HARMONIC     700.0000    2.0415
 Tc       C_R         HARMONIC     700.0000    2.0415
 Tc       C_11        HARMONIC     700.0000    1.9435
 Tc       C_1         HARMONIC     700.0000    1.9435
 Tc       N_3         HARMONIC     700.0000    2.0435
 Tc       N_2         HARMONIC     700.0000    1.9565
 Tc       N_R         HARMONIC     700.0000    1.9915
 Tc       N_1         HARMONIC     700.0000    1.8975
 Tc       O_3         HARMONIC     700.0000    2.0015
 Tc       O_2         HARMONIC     700.0000    1.9015
 Tc       O_R         HARMONIC     700.0000    2.0015
 Tc       O_1         HARMONIC     700.0000    1.8695
 Tc       F_          HARMONIC     700.0000    1.9525
 Tc       Al3         HARMONIC     700.0000    2.3885
 Tc       Si3         HARMONIC     700.0000    2.2785
 Tc       P_3         HARMONIC     700.0000    2.2315
 Tc       S_3         HARMONIC     700.0000    2.3815
 Tc       Cl          HARMONIC     700.0000    2.3385
 Tc       Ga3         HARMONIC     700.0000    2.5515
 Tc       Ge3         HARMONIC     700.0000    2.5515
 Tc       As3         HARMONIC     700.0000    2.5515
 Tc       Se3         HARMONIC     700.0000    2.5515
 Tc       Br          HARMONIC     700.0000    2.5085
 Tc       In3         HARMONIC     700.0000    2.7315
 Tc       Sn3         HARMONIC     700.0000    2.7145
 Tc       Sb3         HARMONIC     700.0000    2.7735
 Tc       Te3         HARMONIC     700.0000    2.6215
 Tc       I_          HARMONIC     700.0000    2.7015
 Tc       Na          HARMONIC     700.0000    3.2015
 Tc       Ca          HARMONIC     700.0000    3.2815
 Tc       Ti          HARMONIC     700.0000    2.7715
 Tc       Fe          HARMONIC     700.0000    2.6265
 Tc       Zn          HARMONIC     700.0000    2.6715
 Tc       Tc          HARMONIC     700.0000    2.6930
 Ru       H_          HARMONIC     700.0000    1.6500
 Ru       H___A       HARMONIC     700.0000    1.6500
 Ru       H___b       HARMONIC     700.0000    1.8300
 Ru       B_3         HARMONIC     700.0000    2.2000
 Ru       B_2         HARMONIC     700.0000    2.1100
 Ru       C_34        HARMONIC     700.0000    2.0900
 Ru       C_33        HARMONIC     700.0000    2.0900
 Ru       C_32        HARMONIC     700.0000    2.0900
 Ru       C_31        HARMONIC     700.0000    2.0900
 Ru       C_3         HARMONIC     700.0000    2.0900
 Ru       C_22        HARMONIC     700.0000    1.9900
 Ru       C_21        HARMONIC     700.0000    1.9900
 Ru       C_2         HARMONIC     700.0000    1.9900
 Ru       C_R2        HARMONIC     700.0000    2.0200
 Ru       C_R1        HARMONIC     700.0000    2.0200
 Ru       C_R         HARMONIC     700.0000    2.0200
 Ru       C_11        HARMONIC     700.0000    1.9220
 Ru       C_1         HARMONIC     700.0000    1.9220
 Ru       N_3         HARMONIC     700.0000    2.0220
 Ru       N_2         HARMONIC     700.0000    1.9350
 Ru       N_R         HARMONIC     700.0000    1.9700
 Ru       N_1         HARMONIC     700.0000    1.8760
 Ru       O_3         HARMONIC     700.0000    1.9800
 Ru       O_2         HARMONIC     700.0000    1.8800
 Ru       O_R         HARMONIC     700.0000    1.9800
 Ru       O_1         HARMONIC     700.0000    1.8480
 Ru       F_          HARMONIC     700.0000    1.9310
 Ru       Al3         HARMONIC     700.0000    2.3670
 Ru       Si3         HARMONIC     700.0000    2.2570
 Ru       P_3         HARMONIC     700.0000    2.2100
 Ru       S_3         HARMONIC     700.0000    2.3600
 Ru       Cl          HARMONIC     700.0000    2.3170
 Ru       Ga3         HARMONIC     700.0000    2.5300
 Ru       Ge3         HARMONIC     700.0000    2.5300
 Ru       As3         HARMONIC     700.0000    2.5300
 Ru       Se3         HARMONIC     700.0000    2.5300
 Ru       Br          HARMONIC     700.0000    2.4870
 Ru       In3         HARMONIC     700.0000    2.7100
 Ru       Sn3         HARMONIC     700.0000    2.6930
 Ru       Sb3         HARMONIC     700.0000    2.7520
 Ru       Te3         HARMONIC     700.0000    2.6000
 Ru       I_          HARMONIC     700.0000    2.6800
 Ru       Na          HARMONIC     700.0000    3.1800
 Ru       Ca          HARMONIC     700.0000    3.2600
 Ru       Ti          HARMONIC     700.0000    2.7500
 Ru       Fe          HARMONIC     700.0000    2.6050
 Ru       Zn          HARMONIC     700.0000    2.6500
 Ru       Tc          HARMONIC     700.0000    2.6715
 Ru       Ru          HARMONIC     700.0000    2.6500
END
#
ANGLE_BEND
 O_2      C_1      O_2         THETA_HARM    56.53    180.0000
 X        H___b    X           THETA_HARM   100.0000   90.0000
 X        B_3      X           THETA_HARM   100.0000  109.4710
 X        B_2      X           THETA_HARM   100.0000  120.0000
 X        C_34     X           THETA_HARM   100.0000  109.4710
 X        C_33     X           THETA_HARM   100.0000  109.4710
 X        C_32     X           THETA_HARM   100.0000  109.4710
 X        C_31     X           THETA_HARM   100.0000  109.4710
 X        C_3      X           THETA_HARM   100.0000  109.4710
 X        C_22     X           THETA_HARM   100.0000  120.0000
 X        C_21     X           THETA_HARM   100.0000  120.0000
 X        C_2      X           THETA_HARM   100.0000  120.0000
 X        C_R2     X           THETA_HARM   100.0000  120.0000
 X        C_R1     X           THETA_HARM   100.0000  120.0000
 X        C_R      X           THETA_HARM   100.0000  120.0000
 X        C_11     X           THETA_HARM   100.0000  180.0000
# X        C_1      X           THETA_HARM   100.0000  180.0000
 X        N_3      X           THETA_HARM   100.0000  106.7000
 X        N_2      X           THETA_HARM   100.0000  120.0000
 X        N_R      X           THETA_HARM   100.0000  120.0000
 X        N_1      X           THETA_HARM   100.0000  180.0000
 X        O_3      X           THETA_HARM   100.0000  104.5100
 X        O_2      X           THETA_HARM   100.0000  120.0000
 X        O_R      X           THETA_HARM   100.0000  120.0000
 X        Al3      X           THETA_HARM   100.0000  109.4710
 X        Si3      X           THETA_HARM   100.0000  109.4710
 X        P_3      X           THETA_HARM   100.0000   93.3000
 X        S_3      X           THETA_HARM   100.0000   92.1000
 X        Ga3      X           THETA_HARM   100.0000  109.4710
 X        Ge3      X           THETA_HARM   100.0000  109.4710
 X        As3      X           THETA_HARM   100.0000   92.1000
 X        Se3      X           THETA_HARM   100.0000   90.6000
 X        In3      X           THETA_HARM   100.0000  109.4710
 X        Sn3      X           THETA_HARM   100.0000  109.4710
 X        Sb3      X           THETA_HARM   100.0000   91.6000
 X        Te3      X           THETA_HARM   100.0000   90.3000
 X        Na       X           THETA_HARM   100.0000   90.0000
 X        Ca       X           THETA_HARM   100.0000   90.0000
 X        Ti       X           THETA_HARM   100.0000   90.0000
 X        Fe       X           THETA_HARM   100.0000   90.0000
 X        Zn       X           THETA_HARM   100.0000  109.4710
 X        Tc       X           THETA_HARM   100.0000   90.0000
 X        Ru       X           THETA_HARM   100.0000   90.0000
END
#
TORSIONS
 X        B_3      B_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        B_2      B_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 B_2      B_2      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_22     B_2      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_21     B_2      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_2      B_2      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R2     B_2      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R1     B_2      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R      B_2      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 N_2      B_2      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 N_R      B_2      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 O_2      B_2      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 O_R      B_2      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 X        B_2      B_2      X           DIHEDRAL      45.0000    2.0000    1.0000
 X        C_34     B_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        C_34     B_2      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        C_34     B_2      B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        C_34     B_2      C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        C_34     B_2      C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        C_34     B_2      C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        C_34     B_2      C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        C_34     B_2      C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        C_34     B_2      C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        C_34     B_2      N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        C_34     B_2      N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        C_34     B_2      O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        C_34     B_2      O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        C_34     C_34     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        C_33     B_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        C_33     B_2      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        C_33     B_2      B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        C_33     B_2      C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        C_33     B_2      C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        C_33     B_2      C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        C_33     B_2      C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        C_33     B_2      C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        C_33     B_2      C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        C_33     B_2      N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        C_33     B_2      N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        C_33     B_2      O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        C_33     B_2      O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        C_33     C_34     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        C_33     C_33     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        C_32     B_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        C_32     B_2      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        C_32     B_2      B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        C_32     B_2      C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        C_32     B_2      C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        C_32     B_2      C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        C_32     B_2      C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        C_32     B_2      C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        C_32     B_2      C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        C_32     B_2      N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        C_32     B_2      N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        C_32     B_2      O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        C_32     B_2      O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        C_32     C_34     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        C_32     C_33     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        C_32     C_32     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        C_31     B_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        C_31     B_2      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        C_31     B_2      B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        C_31     B_2      C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        C_31     B_2      C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        C_31     B_2      C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        C_31     B_2      C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        C_31     B_2      C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        C_31     B_2      C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        C_31     B_2      N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        C_31     B_2      N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        C_31     B_2      O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        C_31     B_2      O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        C_31     C_34     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        C_31     C_33     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        C_31     C_32     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        C_31     C_31     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        C_3      B_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        C_3      B_2      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        C_3      B_2      B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        C_3      B_2      C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        C_3      B_2      C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        C_3      B_2      C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        C_3      B_2      C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        C_3      B_2      C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        C_3      B_2      C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        C_3      B_2      N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        C_3      B_2      N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        C_3      B_2      O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        C_3      B_2      O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        C_3      C_34     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        C_3      C_33     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        C_3      C_32     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        C_3      C_31     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        C_3      C_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        C_22     B_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 B_2      C_22     B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_22     C_22     B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_21     C_22     B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_2      C_22     B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R2     C_22     B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R1     C_22     B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R      C_22     B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 N_2      C_22     B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 N_R      C_22     B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 O_2      C_22     B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 O_R      C_22     B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 X        C_22     B_2      X           DIHEDRAL      45.0000    2.0000    1.0000
 X        C_22     C_34     X           DIHEDRAL       2.0000    3.0000   -1.0000
 B_2      C_22     C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_22     C_22     C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_21     C_22     C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_2      C_22     C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R2     C_22     C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R1     C_22     C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R      C_22     C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_2      C_22     C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_R      C_22     C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_2      C_22     C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_R      C_22     C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 X        C_22     C_33     X           DIHEDRAL       2.0000    3.0000   -1.0000
 B_2      C_22     C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_22     C_22     C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_21     C_22     C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_2      C_22     C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R2     C_22     C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R1     C_22     C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R      C_22     C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_2      C_22     C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_R      C_22     C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_2      C_22     C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_R      C_22     C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 X        C_22     C_32     X           DIHEDRAL       2.0000    3.0000   -1.0000
 B_2      C_22     C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_22     C_22     C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_21     C_22     C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_2      C_22     C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R2     C_22     C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R1     C_22     C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R      C_22     C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_2      C_22     C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_R      C_22     C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_2      C_22     C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_R      C_22     C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 X        C_22     C_31     X           DIHEDRAL       2.0000    3.0000   -1.0000
 B_2      C_22     C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_22     C_22     C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_21     C_22     C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_2      C_22     C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R2     C_22     C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R1     C_22     C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R      C_22     C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_2      C_22     C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_R      C_22     C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_2      C_22     C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_R      C_22     C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 X        C_22     C_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 B_2      C_22     C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_22     C_22     C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_21     C_22     C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_2      C_22     C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R2     C_22     C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R1     C_22     C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R      C_22     C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 N_2      C_22     C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 N_R      C_22     C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 O_2      C_22     C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 O_R      C_22     C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 X        C_22     C_22     X           DIHEDRAL      45.0000    2.0000    1.0000
 X        C_21     B_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 B_2      C_21     B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_22     C_21     B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_21     C_21     B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_2      C_21     B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R2     C_21     B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R1     C_21     B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R      C_21     B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 N_2      C_21     B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 N_R      C_21     B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 O_2      C_21     B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 O_R      C_21     B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 X        C_21     B_2      X           DIHEDRAL      45.0000    2.0000    1.0000
 X        C_21     C_34     X           DIHEDRAL       2.0000    3.0000   -1.0000
 B_2      C_21     C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_22     C_21     C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_21     C_21     C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_2      C_21     C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R2     C_21     C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R1     C_21     C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R      C_21     C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_2      C_21     C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_R      C_21     C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_2      C_21     C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_R      C_21     C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 X        C_21     C_33     X           DIHEDRAL       2.0000    3.0000   -1.0000
 B_2      C_21     C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_22     C_21     C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_21     C_21     C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_2      C_21     C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R2     C_21     C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R1     C_21     C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R      C_21     C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_2      C_21     C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_R      C_21     C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_2      C_21     C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_R      C_21     C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 X        C_21     C_32     X           DIHEDRAL       2.0000    3.0000   -1.0000
 B_2      C_21     C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_22     C_21     C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_21     C_21     C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_2      C_21     C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R2     C_21     C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R1     C_21     C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R      C_21     C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_2      C_21     C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_R      C_21     C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_2      C_21     C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_R      C_21     C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 X        C_21     C_31     X           DIHEDRAL       2.0000    3.0000   -1.0000
 B_2      C_21     C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_22     C_21     C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_21     C_21     C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_2      C_21     C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R2     C_21     C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R1     C_21     C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R      C_21     C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_2      C_21     C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_R      C_21     C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_2      C_21     C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_R      C_21     C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 X        C_21     C_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 B_2      C_21     C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_22     C_21     C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_21     C_21     C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_2      C_21     C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R2     C_21     C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R1     C_21     C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R      C_21     C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 N_2      C_21     C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 N_R      C_21     C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 O_2      C_21     C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 O_R      C_21     C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 X        C_21     C_22     X           DIHEDRAL      45.0000    2.0000    1.0000
 X        C_21     C_21     X           DIHEDRAL      45.0000    2.0000    1.0000
 X        C_2      B_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 B_2      C_2      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_22     C_2      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_21     C_2      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_2      C_2      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R2     C_2      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R1     C_2      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R      C_2      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 N_2      C_2      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 N_R      C_2      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 O_2      C_2      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 O_R      C_2      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 X        C_2      B_2      X           DIHEDRAL      45.0000    2.0000    1.0000
 X        C_2      C_34     X           DIHEDRAL       2.0000    3.0000   -1.0000
 B_2      C_2      C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_22     C_2      C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_21     C_2      C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_2      C_2      C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R2     C_2      C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R1     C_2      C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R      C_2      C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_2      C_2      C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_R      C_2      C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_2      C_2      C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_R      C_2      C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 X        C_2      C_33     X           DIHEDRAL       2.0000    3.0000   -1.0000
 B_2      C_2      C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_22     C_2      C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_21     C_2      C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_2      C_2      C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R2     C_2      C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R1     C_2      C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R      C_2      C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_2      C_2      C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_R      C_2      C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_2      C_2      C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_R      C_2      C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 X        C_2      C_32     X           DIHEDRAL       2.0000    3.0000   -1.0000
 B_2      C_2      C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_22     C_2      C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_21     C_2      C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_2      C_2      C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R2     C_2      C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R1     C_2      C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R      C_2      C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_2      C_2      C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_R      C_2      C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_2      C_2      C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_R      C_2      C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 X        C_2      C_31     X           DIHEDRAL       2.0000    3.0000   -1.0000
 B_2      C_2      C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_22     C_2      C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_21     C_2      C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_2      C_2      C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R2     C_2      C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R1     C_2      C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R      C_2      C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_2      C_2      C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_R      C_2      C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_2      C_2      C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_R      C_2      C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 X        C_2      C_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 B_2      C_2      C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_22     C_2      C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_21     C_2      C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_2      C_2      C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R2     C_2      C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R1     C_2      C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R      C_2      C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 N_2      C_2      C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 N_R      C_2      C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 O_2      C_2      C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 O_R      C_2      C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 X        C_2      C_22     X           DIHEDRAL      45.0000    2.0000    1.0000
 X        C_2      C_21     X           DIHEDRAL      45.0000    2.0000    1.0000
 X        C_2      C_2      X           DIHEDRAL      45.0000    2.0000    1.0000
 X        C_R2     B_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 B_2      C_R2     B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_22     C_R2     B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_21     C_R2     B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_2      C_R2     B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R2     C_R2     B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R1     C_R2     B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R      C_R2     B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 N_2      C_R2     B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 N_R      C_R2     B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 O_2      C_R2     B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 O_R      C_R2     B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 X        C_R2     B_2      X           DIHEDRAL      25.0000    2.0000    1.0000
 X        C_R2     C_34     X           DIHEDRAL       2.0000    3.0000   -1.0000
 B_2      C_R2     C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_22     C_R2     C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_21     C_R2     C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_2      C_R2     C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R2     C_R2     C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R1     C_R2     C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R      C_R2     C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_2      C_R2     C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_R      C_R2     C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_2      C_R2     C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_R      C_R2     C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 X        C_R2     C_33     X           DIHEDRAL       2.0000    3.0000   -1.0000
 B_2      C_R2     C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_22     C_R2     C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_21     C_R2     C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_2      C_R2     C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R2     C_R2     C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R1     C_R2     C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R      C_R2     C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_2      C_R2     C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_R      C_R2     C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_2      C_R2     C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_R      C_R2     C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 X        C_R2     C_32     X           DIHEDRAL       2.0000    3.0000   -1.0000
 B_2      C_R2     C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_22     C_R2     C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_21     C_R2     C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_2      C_R2     C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R2     C_R2     C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R1     C_R2     C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R      C_R2     C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_2      C_R2     C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_R      C_R2     C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_2      C_R2     C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_R      C_R2     C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 X        C_R2     C_31     X           DIHEDRAL       2.0000    3.0000   -1.0000
 B_2      C_R2     C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_22     C_R2     C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_21     C_R2     C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_2      C_R2     C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R2     C_R2     C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R1     C_R2     C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R      C_R2     C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_2      C_R2     C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_R      C_R2     C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_2      C_R2     C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_R      C_R2     C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 X        C_R2     C_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 B_2      C_R2     C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_22     C_R2     C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_21     C_R2     C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_2      C_R2     C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R2     C_R2     C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R1     C_R2     C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R      C_R2     C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 N_2      C_R2     C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 N_R      C_R2     C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 O_2      C_R2     C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 O_R      C_R2     C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 X        C_R2     C_22     X           DIHEDRAL      25.0000    2.0000    1.0000
 X        C_R2     C_21     X           DIHEDRAL      25.0000    2.0000    1.0000
 X        C_R2     C_2      X           DIHEDRAL      25.0000    2.0000    1.0000
 X        C_R2     C_R2     X           DIHEDRAL      25.0000    2.0000    1.0000
 X        C_R1     B_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 B_2      C_R1     B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_22     C_R1     B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_21     C_R1     B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_2      C_R1     B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R2     C_R1     B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R1     C_R1     B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R      C_R1     B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 N_2      C_R1     B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 N_R      C_R1     B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 O_2      C_R1     B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 O_R      C_R1     B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 X        C_R1     B_2      X           DIHEDRAL      25.0000    2.0000    1.0000
 X        C_R1     C_34     X           DIHEDRAL       2.0000    3.0000   -1.0000
 B_2      C_R1     C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_22     C_R1     C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_21     C_R1     C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_2      C_R1     C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R2     C_R1     C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R1     C_R1     C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R      C_R1     C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_2      C_R1     C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_R      C_R1     C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_2      C_R1     C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_R      C_R1     C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 X        C_R1     C_33     X           DIHEDRAL       2.0000    3.0000   -1.0000
 B_2      C_R1     C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_22     C_R1     C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_21     C_R1     C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_2      C_R1     C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R2     C_R1     C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R1     C_R1     C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R      C_R1     C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_2      C_R1     C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_R      C_R1     C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_2      C_R1     C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_R      C_R1     C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 X        C_R1     C_32     X           DIHEDRAL       2.0000    3.0000   -1.0000
 B_2      C_R1     C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_22     C_R1     C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_21     C_R1     C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_2      C_R1     C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R2     C_R1     C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R1     C_R1     C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R      C_R1     C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_2      C_R1     C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_R      C_R1     C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_2      C_R1     C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_R      C_R1     C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 X        C_R1     C_31     X           DIHEDRAL       2.0000    3.0000   -1.0000
 B_2      C_R1     C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_22     C_R1     C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_21     C_R1     C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_2      C_R1     C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R2     C_R1     C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R1     C_R1     C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R      C_R1     C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_2      C_R1     C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_R      C_R1     C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_2      C_R1     C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_R      C_R1     C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 X        C_R1     C_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 B_2      C_R1     C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_22     C_R1     C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_21     C_R1     C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_2      C_R1     C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R2     C_R1     C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R1     C_R1     C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R      C_R1     C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 N_2      C_R1     C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 N_R      C_R1     C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 O_2      C_R1     C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 O_R      C_R1     C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 X        C_R1     C_22     X           DIHEDRAL      25.0000    2.0000    1.0000
 X        C_R1     C_21     X           DIHEDRAL      25.0000    2.0000    1.0000
 X        C_R1     C_2      X           DIHEDRAL      25.0000    2.0000    1.0000
 X        C_R1     C_R2     X           DIHEDRAL      25.0000    2.0000    1.0000
 X        C_R1     C_R1     X           DIHEDRAL      25.0000    2.0000    1.0000
 X        C_R      B_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 B_2      C_R      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_22     C_R      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_21     C_R      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_2      C_R      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R2     C_R      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R1     C_R      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R      C_R      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 N_2      C_R      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 N_R      C_R      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 O_2      C_R      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 O_R      C_R      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 X        C_R      B_2      X           DIHEDRAL      25.0000    2.0000    1.0000
 X        C_R      C_34     X           DIHEDRAL       2.0000    3.0000   -1.0000
 B_2      C_R      C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_22     C_R      C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_21     C_R      C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_2      C_R      C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R2     C_R      C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R1     C_R      C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R      C_R      C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_2      C_R      C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_R      C_R      C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_2      C_R      C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_R      C_R      C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 X        C_R      C_33     X           DIHEDRAL       2.0000    3.0000   -1.0000
 B_2      C_R      C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_22     C_R      C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_21     C_R      C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_2      C_R      C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R2     C_R      C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R1     C_R      C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R      C_R      C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_2      C_R      C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_R      C_R      C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_2      C_R      C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_R      C_R      C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 X        C_R      C_32     X           DIHEDRAL       2.0000    3.0000   -1.0000
 B_2      C_R      C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_22     C_R      C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_21     C_R      C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_2      C_R      C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R2     C_R      C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R1     C_R      C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R      C_R      C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_2      C_R      C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_R      C_R      C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_2      C_R      C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_R      C_R      C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 X        C_R      C_31     X           DIHEDRAL       2.0000    3.0000   -1.0000
 B_2      C_R      C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_22     C_R      C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_21     C_R      C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_2      C_R      C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R2     C_R      C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R1     C_R      C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R      C_R      C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_2      C_R      C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_R      C_R      C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_2      C_R      C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_R      C_R      C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 X        C_R      C_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 B_2      C_R      C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_22     C_R      C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_21     C_R      C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_2      C_R      C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R2     C_R      C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R1     C_R      C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R      C_R      C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 N_2      C_R      C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 N_R      C_R      C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 O_2      C_R      C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 O_R      C_R      C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 X        C_R      C_22     X           DIHEDRAL      25.0000    2.0000    1.0000
 X        C_R      C_21     X           DIHEDRAL      25.0000    2.0000    1.0000
 X        C_R      C_2      X           DIHEDRAL      25.0000    2.0000    1.0000
 X        C_R      C_R2     X           DIHEDRAL      25.0000    2.0000    1.0000
 X        C_R      C_R1     X           DIHEDRAL      25.0000    2.0000    1.0000
 X        C_R      C_R      X           DIHEDRAL      25.0000    2.0000    1.0000
 X        N_3      B_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        N_3      B_2      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        N_3      B_2      B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      B_2      C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      B_2      C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      B_2      C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      B_2      C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      B_2      C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      B_2      C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      B_2      N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      B_2      N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      B_2      O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      B_2      O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      C_34     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        N_3      C_33     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        N_3      C_32     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        N_3      C_31     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        N_3      C_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        N_3      C_22     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        N_3      C_22     B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      C_22     C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      C_22     C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      C_22     C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      C_22     C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      C_22     C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      C_22     C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      C_22     N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      C_22     N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      C_22     O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      C_22     O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      C_21     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        N_3      C_21     B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      C_21     C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      C_21     C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      C_21     C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      C_21     C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      C_21     C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      C_21     C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      C_21     N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      C_21     N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      C_21     O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      C_21     O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      C_2      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        N_3      C_2      B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      C_2      C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      C_2      C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      C_2      C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      C_2      C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      C_2      C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      C_2      C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      C_2      N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      C_2      N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      C_2      O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      C_2      O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      C_R2     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        N_3      C_R2     B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      C_R2     C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      C_R2     C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      C_R2     C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      C_R2     C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      C_R2     C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      C_R2     C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      C_R2     N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      C_R2     N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      C_R2     O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      C_R2     O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      C_R1     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        N_3      C_R1     B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      C_R1     C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      C_R1     C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      C_R1     C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      C_R1     C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      C_R1     C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      C_R1     C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      C_R1     N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      C_R1     N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      C_R1     O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      C_R1     O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      C_R      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        N_3      C_R      B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      C_R      C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      C_R      C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      C_R      C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      C_R      C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      C_R      C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      C_R      C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      C_R      N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      C_R      N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      C_R      O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      C_R      O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        N_3      N_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        N_2      B_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 B_2      N_2      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_22     N_2      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_21     N_2      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_2      N_2      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R2     N_2      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R1     N_2      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R      N_2      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 N_2      N_2      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 N_R      N_2      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 O_2      N_2      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 O_R      N_2      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 X        N_2      B_2      X           DIHEDRAL      45.0000    2.0000    1.0000
 X        N_2      C_34     X           DIHEDRAL       2.0000    3.0000   -1.0000
 B_2      N_2      C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_22     N_2      C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_21     N_2      C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_2      N_2      C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R2     N_2      C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R1     N_2      C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R      N_2      C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_2      N_2      C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_R      N_2      C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_2      N_2      C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_R      N_2      C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 X        N_2      C_33     X           DIHEDRAL       2.0000    3.0000   -1.0000
 B_2      N_2      C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_22     N_2      C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_21     N_2      C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_2      N_2      C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R2     N_2      C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R1     N_2      C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R      N_2      C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_2      N_2      C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_R      N_2      C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_2      N_2      C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_R      N_2      C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 X        N_2      C_32     X           DIHEDRAL       2.0000    3.0000   -1.0000
 B_2      N_2      C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_22     N_2      C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_21     N_2      C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_2      N_2      C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R2     N_2      C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R1     N_2      C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R      N_2      C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_2      N_2      C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_R      N_2      C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_2      N_2      C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_R      N_2      C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 X        N_2      C_31     X           DIHEDRAL       2.0000    3.0000   -1.0000
 B_2      N_2      C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_22     N_2      C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_21     N_2      C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_2      N_2      C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R2     N_2      C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R1     N_2      C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R      N_2      C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_2      N_2      C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_R      N_2      C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_2      N_2      C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_R      N_2      C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 X        N_2      C_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 B_2      N_2      C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_22     N_2      C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_21     N_2      C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_2      N_2      C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R2     N_2      C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R1     N_2      C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R      N_2      C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 N_2      N_2      C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 N_R      N_2      C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 O_2      N_2      C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 O_R      N_2      C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 X        N_2      C_22     X           DIHEDRAL      45.0000    2.0000    1.0000
 X        N_2      C_21     X           DIHEDRAL      45.0000    2.0000    1.0000
 X        N_2      C_2      X           DIHEDRAL      45.0000    2.0000    1.0000
 X        N_2      C_R2     X           DIHEDRAL      25.0000    2.0000    1.0000
 X        N_2      C_R1     X           DIHEDRAL      25.0000    2.0000    1.0000
 X        N_2      C_R      X           DIHEDRAL      25.0000    2.0000    1.0000
 X        N_2      N_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 B_2      N_2      N_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_22     N_2      N_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_21     N_2      N_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_2      N_2      N_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R2     N_2      N_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R1     N_2      N_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R      N_2      N_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 N_2      N_2      N_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 N_R      N_2      N_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 O_2      N_2      N_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 O_R      N_2      N_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 X        N_2      N_2      X           DIHEDRAL      45.0000    2.0000    1.0000
 X        N_R      B_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 B_2      N_R      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_22     N_R      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_21     N_R      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_2      N_R      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R2     N_R      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R1     N_R      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R      N_R      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 N_2      N_R      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 N_R      N_R      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 O_2      N_R      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 O_R      N_R      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 X        N_R      B_2      X           DIHEDRAL      25.0000    2.0000    1.0000
 X        N_R      C_34     X           DIHEDRAL       2.0000    3.0000   -1.0000
 B_2      N_R      C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_22     N_R      C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_21     N_R      C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_2      N_R      C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R2     N_R      C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R1     N_R      C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R      N_R      C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_2      N_R      C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_R      N_R      C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_2      N_R      C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_R      N_R      C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 X        N_R      C_33     X           DIHEDRAL       2.0000    3.0000   -1.0000
 B_2      N_R      C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_22     N_R      C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_21     N_R      C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_2      N_R      C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R2     N_R      C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R1     N_R      C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R      N_R      C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_2      N_R      C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_R      N_R      C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_2      N_R      C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_R      N_R      C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 X        N_R      C_32     X           DIHEDRAL       2.0000    3.0000   -1.0000
 B_2      N_R      C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_22     N_R      C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_21     N_R      C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_2      N_R      C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R2     N_R      C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R1     N_R      C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R      N_R      C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_2      N_R      C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_R      N_R      C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_2      N_R      C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_R      N_R      C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 X        N_R      C_31     X           DIHEDRAL       2.0000    3.0000   -1.0000
 B_2      N_R      C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_22     N_R      C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_21     N_R      C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_2      N_R      C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R2     N_R      C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R1     N_R      C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R      N_R      C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_2      N_R      C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_R      N_R      C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_2      N_R      C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_R      N_R      C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 X        N_R      C_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 B_2      N_R      C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_22     N_R      C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_21     N_R      C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_2      N_R      C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R2     N_R      C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R1     N_R      C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R      N_R      C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 N_2      N_R      C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 N_R      N_R      C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 O_2      N_R      C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 O_R      N_R      C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 X        N_R      C_22     X           DIHEDRAL      25.0000    2.0000    1.0000
 X        N_R      C_21     X           DIHEDRAL      25.0000    2.0000    1.0000
 X        N_R      C_2      X           DIHEDRAL      25.0000    2.0000    1.0000
 X        N_R      C_R2     X           DIHEDRAL      25.0000    2.0000    1.0000
 X        N_R      C_R1     X           DIHEDRAL      25.0000    2.0000    1.0000
 X        N_R      C_R      X           DIHEDRAL      25.0000    2.0000    1.0000
 X        N_R      N_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 B_2      N_R      N_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_22     N_R      N_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_21     N_R      N_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_2      N_R      N_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R2     N_R      N_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R1     N_R      N_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R      N_R      N_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 N_2      N_R      N_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 N_R      N_R      N_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 O_2      N_R      N_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 O_R      N_R      N_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 X        N_R      N_2      X           DIHEDRAL      25.0000    2.0000    1.0000
 X        N_R      N_R      X           DIHEDRAL      25.0000    2.0000    1.0000
 X        O_3      B_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        O_3      B_2      X           DIHEDRAL       2.0000    2.0000    1.0000
 X        O_3      C_34     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        O_3      C_33     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        O_3      C_32     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        O_3      C_31     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        O_3      C_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        O_3      C_22     X           DIHEDRAL       2.0000    2.0000    1.0000
 X        O_3      C_21     X           DIHEDRAL       2.0000    2.0000    1.0000
 X        O_3      C_2      X           DIHEDRAL       2.0000    2.0000    1.0000
 X        O_3      C_R2     X           DIHEDRAL       2.0000    2.0000    1.0000
 X        O_3      C_R1     X           DIHEDRAL       2.0000    2.0000    1.0000
 X        O_3      C_R      X           DIHEDRAL       2.0000    2.0000    1.0000
 X        O_3      N_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        O_3      N_2      X           DIHEDRAL       2.0000    2.0000    1.0000
 X        O_3      N_R      X           DIHEDRAL       2.0000    2.0000    1.0000
 X        O_3      O_3      X           DIHEDRAL       2.0000    2.0000   -1.0000
 X        O_2      B_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 B_2      O_2      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_22     O_2      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_21     O_2      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_2      O_2      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R2     O_2      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R1     O_2      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R      O_2      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 N_2      O_2      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 N_R      O_2      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 O_2      O_2      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 O_R      O_2      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 X        O_2      B_2      X           DIHEDRAL      45.0000    2.0000    1.0000
 X        O_2      C_34     X           DIHEDRAL       2.0000    3.0000   -1.0000
 B_2      O_2      C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_22     O_2      C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_21     O_2      C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_2      O_2      C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R2     O_2      C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R1     O_2      C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R      O_2      C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_2      O_2      C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_R      O_2      C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_2      O_2      C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_R      O_2      C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 X        O_2      C_33     X           DIHEDRAL       2.0000    3.0000   -1.0000
 B_2      O_2      C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_22     O_2      C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_21     O_2      C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_2      O_2      C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R2     O_2      C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R1     O_2      C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R      O_2      C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_2      O_2      C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_R      O_2      C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_2      O_2      C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_R      O_2      C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 X        O_2      C_32     X           DIHEDRAL       2.0000    3.0000   -1.0000
 B_2      O_2      C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_22     O_2      C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_21     O_2      C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_2      O_2      C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R2     O_2      C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R1     O_2      C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R      O_2      C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_2      O_2      C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_R      O_2      C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_2      O_2      C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_R      O_2      C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 X        O_2      C_31     X           DIHEDRAL       2.0000    3.0000   -1.0000
 B_2      O_2      C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_22     O_2      C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_21     O_2      C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_2      O_2      C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R2     O_2      C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R1     O_2      C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R      O_2      C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_2      O_2      C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_R      O_2      C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_2      O_2      C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_R      O_2      C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 X        O_2      C_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 B_2      O_2      C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_22     O_2      C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_21     O_2      C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_2      O_2      C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R2     O_2      C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R1     O_2      C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R      O_2      C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 N_2      O_2      C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 N_R      O_2      C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 O_2      O_2      C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 O_R      O_2      C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 X        O_2      C_22     X           DIHEDRAL      45.0000    2.0000    1.0000
 X        O_2      C_21     X           DIHEDRAL      45.0000    2.0000    1.0000
 X        O_2      C_2      X           DIHEDRAL      45.0000    2.0000    1.0000
 X        O_2      C_R2     X           DIHEDRAL      25.0000    2.0000    1.0000
 X        O_2      C_R1     X           DIHEDRAL      25.0000    2.0000    1.0000
 X        O_2      C_R      X           DIHEDRAL      25.0000    2.0000    1.0000
 X        O_2      N_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 B_2      O_2      N_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_22     O_2      N_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_21     O_2      N_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_2      O_2      N_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R2     O_2      N_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R1     O_2      N_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R      O_2      N_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 N_2      O_2      N_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 N_R      O_2      N_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 O_2      O_2      N_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 O_R      O_2      N_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 X        O_2      N_2      X           DIHEDRAL      45.0000    2.0000    1.0000
 X        O_2      N_R      X           DIHEDRAL      25.0000    2.0000    1.0000
 X        O_2      O_3      X           DIHEDRAL       2.0000    2.0000    1.0000
 X        O_2      O_2      X           DIHEDRAL      45.0000    2.0000    1.0000
 X        O_R      B_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 B_2      O_R      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_22     O_R      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_21     O_R      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_2      O_R      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R2     O_R      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R1     O_R      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R      O_R      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 N_2      O_R      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 N_R      O_R      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 O_2      O_R      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 O_R      O_R      B_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 X        O_R      B_2      X           DIHEDRAL      25.0000    2.0000    1.0000
 X        O_R      C_34     X           DIHEDRAL       2.0000    3.0000   -1.0000
 B_2      O_R      C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_22     O_R      C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_21     O_R      C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_2      O_R      C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R2     O_R      C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R1     O_R      C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R      O_R      C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_2      O_R      C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_R      O_R      C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_2      O_R      C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_R      O_R      C_34     X           DIHEDRAL       1.0000    6.0000    1.0000
 X        O_R      C_33     X           DIHEDRAL       2.0000    3.0000   -1.0000
 B_2      O_R      C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_22     O_R      C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_21     O_R      C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_2      O_R      C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R2     O_R      C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R1     O_R      C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R      O_R      C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_2      O_R      C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_R      O_R      C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_2      O_R      C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_R      O_R      C_33     X           DIHEDRAL       1.0000    6.0000    1.0000
 X        O_R      C_32     X           DIHEDRAL       2.0000    3.0000   -1.0000
 B_2      O_R      C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_22     O_R      C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_21     O_R      C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_2      O_R      C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R2     O_R      C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R1     O_R      C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R      O_R      C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_2      O_R      C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_R      O_R      C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_2      O_R      C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_R      O_R      C_32     X           DIHEDRAL       1.0000    6.0000    1.0000
 X        O_R      C_31     X           DIHEDRAL       2.0000    3.0000   -1.0000
 B_2      O_R      C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_22     O_R      C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_21     O_R      C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_2      O_R      C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R2     O_R      C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R1     O_R      C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R      O_R      C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_2      O_R      C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 N_R      O_R      C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_2      O_R      C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 O_R      O_R      C_31     X           DIHEDRAL       1.0000    6.0000    1.0000
 X        O_R      C_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 B_2      O_R      C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_22     O_R      C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_21     O_R      C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_2      O_R      C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R2     O_R      C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R1     O_R      C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R      O_R      C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 N_2      O_R      C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 N_R      O_R      C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 O_2      O_R      C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 O_R      O_R      C_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 X        O_R      C_22     X           DIHEDRAL      25.0000    2.0000    1.0000
 X        O_R      C_21     X           DIHEDRAL      25.0000    2.0000    1.0000
 X        O_R      C_2      X           DIHEDRAL      25.0000    2.0000    1.0000
 X        O_R      C_R2     X           DIHEDRAL      25.0000    2.0000    1.0000
 X        O_R      C_R1     X           DIHEDRAL      25.0000    2.0000    1.0000
 X        O_R      C_R      X           DIHEDRAL      25.0000    2.0000    1.0000
 X        O_R      N_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 B_2      O_R      N_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_22     O_R      N_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_21     O_R      N_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_2      O_R      N_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R2     O_R      N_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R1     O_R      N_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 C_R      O_R      N_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 N_2      O_R      N_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 N_R      O_R      N_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 O_2      O_R      N_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 O_R      O_R      N_3      X           DIHEDRAL       1.0000    6.0000    1.0000
 X        O_R      N_2      X           DIHEDRAL      25.0000    2.0000    1.0000
 X        O_R      N_R      X           DIHEDRAL      25.0000    2.0000    1.0000
 X        O_R      O_3      X           DIHEDRAL       2.0000    2.0000    1.0000
 X        O_R      O_2      X           DIHEDRAL      25.0000    2.0000    1.0000
 X        O_R      O_R      X           DIHEDRAL      25.0000    2.0000    1.0000
 X        Al3      B_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Al3      B_2      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Al3      B_2      B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      B_2      C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      B_2      C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      B_2      C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      B_2      C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      B_2      C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      B_2      C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      B_2      N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      B_2      N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      B_2      O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      B_2      O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      C_34     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Al3      C_33     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Al3      C_32     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Al3      C_31     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Al3      C_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Al3      C_22     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Al3      C_22     B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      C_22     C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      C_22     C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      C_22     C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      C_22     C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      C_22     C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      C_22     C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      C_22     N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      C_22     N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      C_22     O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      C_22     O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      C_21     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Al3      C_21     B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      C_21     C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      C_21     C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      C_21     C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      C_21     C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      C_21     C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      C_21     C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      C_21     N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      C_21     N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      C_21     O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      C_21     O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      C_2      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Al3      C_2      B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      C_2      C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      C_2      C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      C_2      C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      C_2      C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      C_2      C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      C_2      C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      C_2      N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      C_2      N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      C_2      O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      C_2      O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      C_R2     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Al3      C_R2     B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      C_R2     C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      C_R2     C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      C_R2     C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      C_R2     C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      C_R2     C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      C_R2     C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      C_R2     N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      C_R2     N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      C_R2     O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      C_R2     O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      C_R1     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Al3      C_R1     B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      C_R1     C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      C_R1     C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      C_R1     C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      C_R1     C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      C_R1     C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      C_R1     C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      C_R1     N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      C_R1     N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      C_R1     O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      C_R1     O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      C_R      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Al3      C_R      B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      C_R      C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      C_R      C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      C_R      C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      C_R      C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      C_R      C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      C_R      C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      C_R      N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      C_R      N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      C_R      O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      C_R      O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      N_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Al3      N_2      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Al3      N_2      B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      N_2      C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      N_2      C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      N_2      C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      N_2      C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      N_2      C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      N_2      C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      N_2      N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      N_2      N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      N_2      O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      N_2      O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      N_R      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Al3      N_R      B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      N_R      C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      N_R      C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      N_R      C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      N_R      C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      N_R      C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      N_R      C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      N_R      N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      N_R      N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      N_R      O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      N_R      O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      O_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Al3      O_2      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Al3      O_2      B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      O_2      C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      O_2      C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      O_2      C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      O_2      C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      O_2      C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      O_2      C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      O_2      N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      O_2      N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      O_2      O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      O_2      O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      O_R      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Al3      O_R      B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      O_R      C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      O_R      C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      O_R      C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      O_R      C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      O_R      C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      O_R      C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      O_R      N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      O_R      N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      O_R      O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      O_R      O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Al3      Al3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Si3      B_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Si3      B_2      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Si3      B_2      B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      B_2      C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      B_2      C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      B_2      C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      B_2      C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      B_2      C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      B_2      C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      B_2      N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      B_2      N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      B_2      O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      B_2      O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      C_34     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Si3      C_33     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Si3      C_32     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Si3      C_31     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Si3      C_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Si3      C_22     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Si3      C_22     B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      C_22     C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      C_22     C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      C_22     C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      C_22     C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      C_22     C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      C_22     C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      C_22     N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      C_22     N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      C_22     O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      C_22     O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      C_21     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Si3      C_21     B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      C_21     C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      C_21     C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      C_21     C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      C_21     C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      C_21     C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      C_21     C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      C_21     N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      C_21     N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      C_21     O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      C_21     O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      C_2      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Si3      C_2      B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      C_2      C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      C_2      C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      C_2      C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      C_2      C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      C_2      C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      C_2      C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      C_2      N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      C_2      N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      C_2      O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      C_2      O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      C_R2     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Si3      C_R2     B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      C_R2     C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      C_R2     C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      C_R2     C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      C_R2     C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      C_R2     C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      C_R2     C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      C_R2     N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      C_R2     N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      C_R2     O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      C_R2     O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      C_R1     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Si3      C_R1     B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      C_R1     C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      C_R1     C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      C_R1     C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      C_R1     C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      C_R1     C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      C_R1     C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      C_R1     N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      C_R1     N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      C_R1     O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      C_R1     O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      C_R      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Si3      C_R      B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      C_R      C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      C_R      C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      C_R      C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      C_R      C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      C_R      C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      C_R      C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      C_R      N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      C_R      N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      C_R      O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      C_R      O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      N_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Si3      N_2      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Si3      N_2      B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      N_2      C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      N_2      C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      N_2      C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      N_2      C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      N_2      C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      N_2      C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      N_2      N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      N_2      N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      N_2      O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      N_2      O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      N_R      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Si3      N_R      B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      N_R      C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      N_R      C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      N_R      C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      N_R      C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      N_R      C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      N_R      C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      N_R      N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      N_R      N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      N_R      O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      N_R      O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      O_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Si3      O_2      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Si3      O_2      B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      O_2      C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      O_2      C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      O_2      C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      O_2      C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      O_2      C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      O_2      C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      O_2      N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      O_2      N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      O_2      O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      O_2      O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      O_R      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Si3      O_R      B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      O_R      C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      O_R      C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      O_R      C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      O_R      C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      O_R      C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      O_R      C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      O_R      N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      O_R      N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      O_R      O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      O_R      O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Si3      Al3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Si3      Si3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        P_3      B_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        P_3      B_2      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        P_3      B_2      B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      B_2      C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      B_2      C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      B_2      C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      B_2      C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      B_2      C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      B_2      C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      B_2      N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      B_2      N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      B_2      O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      B_2      O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      C_34     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        P_3      C_33     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        P_3      C_32     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        P_3      C_31     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        P_3      C_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        P_3      C_22     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        P_3      C_22     B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      C_22     C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      C_22     C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      C_22     C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      C_22     C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      C_22     C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      C_22     C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      C_22     N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      C_22     N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      C_22     O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      C_22     O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      C_21     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        P_3      C_21     B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      C_21     C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      C_21     C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      C_21     C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      C_21     C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      C_21     C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      C_21     C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      C_21     N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      C_21     N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      C_21     O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      C_21     O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      C_2      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        P_3      C_2      B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      C_2      C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      C_2      C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      C_2      C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      C_2      C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      C_2      C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      C_2      C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      C_2      N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      C_2      N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      C_2      O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      C_2      O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      C_R2     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        P_3      C_R2     B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      C_R2     C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      C_R2     C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      C_R2     C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      C_R2     C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      C_R2     C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      C_R2     C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      C_R2     N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      C_R2     N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      C_R2     O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      C_R2     O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      C_R1     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        P_3      C_R1     B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      C_R1     C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      C_R1     C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      C_R1     C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      C_R1     C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      C_R1     C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      C_R1     C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      C_R1     N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      C_R1     N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      C_R1     O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      C_R1     O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      C_R      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        P_3      C_R      B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      C_R      C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      C_R      C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      C_R      C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      C_R      C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      C_R      C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      C_R      C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      C_R      N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      C_R      N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      C_R      O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      C_R      O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      N_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        P_3      N_2      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        P_3      N_2      B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      N_2      C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      N_2      C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      N_2      C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      N_2      C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      N_2      C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      N_2      C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      N_2      N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      N_2      N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      N_2      O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      N_2      O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      N_R      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        P_3      N_R      B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      N_R      C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      N_R      C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      N_R      C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      N_R      C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      N_R      C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      N_R      C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      N_R      N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      N_R      N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      N_R      O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      N_R      O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      O_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        P_3      O_2      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        P_3      O_2      B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      O_2      C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      O_2      C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      O_2      C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      O_2      C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      O_2      C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      O_2      C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      O_2      N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      O_2      N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      O_2      O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      O_2      O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      O_R      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        P_3      O_R      B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      O_R      C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      O_R      C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      O_R      C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      O_R      C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      O_R      C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      O_R      C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      O_R      N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      O_R      N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      O_R      O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      O_R      O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        P_3      Al3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        P_3      Si3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        P_3      P_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        S_3      B_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        S_3      B_2      X           DIHEDRAL       2.0000    2.0000    1.0000
 X        S_3      C_34     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        S_3      C_33     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        S_3      C_32     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        S_3      C_31     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        S_3      C_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        S_3      C_22     X           DIHEDRAL       2.0000    2.0000    1.0000
 X        S_3      C_21     X           DIHEDRAL       2.0000    2.0000    1.0000
 X        S_3      C_2      X           DIHEDRAL       2.0000    2.0000    1.0000
 X        S_3      C_R2     X           DIHEDRAL       2.0000    2.0000    1.0000
 X        S_3      C_R1     X           DIHEDRAL       2.0000    2.0000    1.0000
 X        S_3      C_R      X           DIHEDRAL       2.0000    2.0000    1.0000
 X        S_3      N_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        S_3      N_2      X           DIHEDRAL       2.0000    2.0000    1.0000
 X        S_3      N_R      X           DIHEDRAL       2.0000    2.0000    1.0000
 X        S_3      O_3      X           DIHEDRAL       2.0000    2.0000   -1.0000
 X        S_3      O_2      X           DIHEDRAL       2.0000    2.0000    1.0000
 X        S_3      O_R      X           DIHEDRAL       2.0000    2.0000    1.0000
 X        S_3      Al3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        S_3      Si3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        S_3      P_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        S_3      S_3      X           DIHEDRAL       2.0000    2.0000   -1.0000
 X        Ga3      B_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Ga3      B_2      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Ga3      B_2      B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      B_2      C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      B_2      C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      B_2      C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      B_2      C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      B_2      C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      B_2      C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      B_2      N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      B_2      N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      B_2      O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      B_2      O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      C_34     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Ga3      C_33     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Ga3      C_32     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Ga3      C_31     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Ga3      C_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Ga3      C_22     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Ga3      C_22     B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      C_22     C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      C_22     C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      C_22     C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      C_22     C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      C_22     C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      C_22     C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      C_22     N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      C_22     N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      C_22     O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      C_22     O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      C_21     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Ga3      C_21     B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      C_21     C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      C_21     C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      C_21     C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      C_21     C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      C_21     C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      C_21     C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      C_21     N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      C_21     N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      C_21     O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      C_21     O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      C_2      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Ga3      C_2      B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      C_2      C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      C_2      C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      C_2      C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      C_2      C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      C_2      C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      C_2      C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      C_2      N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      C_2      N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      C_2      O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      C_2      O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      C_R2     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Ga3      C_R2     B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      C_R2     C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      C_R2     C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      C_R2     C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      C_R2     C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      C_R2     C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      C_R2     C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      C_R2     N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      C_R2     N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      C_R2     O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      C_R2     O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      C_R1     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Ga3      C_R1     B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      C_R1     C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      C_R1     C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      C_R1     C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      C_R1     C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      C_R1     C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      C_R1     C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      C_R1     N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      C_R1     N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      C_R1     O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      C_R1     O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      C_R      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Ga3      C_R      B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      C_R      C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      C_R      C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      C_R      C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      C_R      C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      C_R      C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      C_R      C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      C_R      N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      C_R      N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      C_R      O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      C_R      O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      N_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Ga3      N_2      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Ga3      N_2      B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      N_2      C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      N_2      C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      N_2      C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      N_2      C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      N_2      C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      N_2      C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      N_2      N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      N_2      N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      N_2      O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      N_2      O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      N_R      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Ga3      N_R      B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      N_R      C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      N_R      C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      N_R      C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      N_R      C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      N_R      C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      N_R      C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      N_R      N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      N_R      N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      N_R      O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      N_R      O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      O_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Ga3      O_2      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Ga3      O_2      B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      O_2      C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      O_2      C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      O_2      C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      O_2      C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      O_2      C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      O_2      C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      O_2      N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      O_2      N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      O_2      O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      O_2      O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      O_R      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Ga3      O_R      B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      O_R      C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      O_R      C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      O_R      C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      O_R      C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      O_R      C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      O_R      C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      O_R      N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      O_R      N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      O_R      O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      O_R      O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ga3      Al3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Ga3      Si3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Ga3      P_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Ga3      S_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Ga3      Ga3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Ge3      B_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Ge3      B_2      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Ge3      B_2      B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      B_2      C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      B_2      C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      B_2      C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      B_2      C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      B_2      C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      B_2      C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      B_2      N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      B_2      N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      B_2      O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      B_2      O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      C_34     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Ge3      C_33     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Ge3      C_32     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Ge3      C_31     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Ge3      C_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Ge3      C_22     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Ge3      C_22     B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      C_22     C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      C_22     C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      C_22     C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      C_22     C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      C_22     C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      C_22     C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      C_22     N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      C_22     N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      C_22     O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      C_22     O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      C_21     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Ge3      C_21     B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      C_21     C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      C_21     C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      C_21     C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      C_21     C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      C_21     C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      C_21     C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      C_21     N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      C_21     N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      C_21     O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      C_21     O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      C_2      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Ge3      C_2      B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      C_2      C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      C_2      C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      C_2      C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      C_2      C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      C_2      C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      C_2      C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      C_2      N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      C_2      N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      C_2      O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      C_2      O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      C_R2     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Ge3      C_R2     B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      C_R2     C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      C_R2     C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      C_R2     C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      C_R2     C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      C_R2     C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      C_R2     C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      C_R2     N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      C_R2     N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      C_R2     O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      C_R2     O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      C_R1     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Ge3      C_R1     B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      C_R1     C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      C_R1     C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      C_R1     C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      C_R1     C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      C_R1     C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      C_R1     C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      C_R1     N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      C_R1     N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      C_R1     O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      C_R1     O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      C_R      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Ge3      C_R      B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      C_R      C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      C_R      C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      C_R      C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      C_R      C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      C_R      C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      C_R      C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      C_R      N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      C_R      N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      C_R      O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      C_R      O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      N_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Ge3      N_2      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Ge3      N_2      B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      N_2      C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      N_2      C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      N_2      C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      N_2      C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      N_2      C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      N_2      C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      N_2      N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      N_2      N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      N_2      O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      N_2      O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      N_R      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Ge3      N_R      B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      N_R      C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      N_R      C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      N_R      C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      N_R      C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      N_R      C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      N_R      C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      N_R      N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      N_R      N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      N_R      O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      N_R      O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      O_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Ge3      O_2      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Ge3      O_2      B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      O_2      C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      O_2      C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      O_2      C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      O_2      C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      O_2      C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      O_2      C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      O_2      N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      O_2      N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      O_2      O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      O_2      O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      O_R      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Ge3      O_R      B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      O_R      C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      O_R      C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      O_R      C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      O_R      C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      O_R      C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      O_R      C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      O_R      N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      O_R      N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      O_R      O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      O_R      O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Ge3      Al3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Ge3      Si3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Ge3      P_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Ge3      S_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Ge3      Ga3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Ge3      Ge3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        As3      B_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        As3      B_2      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        As3      B_2      B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      B_2      C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      B_2      C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      B_2      C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      B_2      C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      B_2      C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      B_2      C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      B_2      N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      B_2      N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      B_2      O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      B_2      O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      C_34     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        As3      C_33     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        As3      C_32     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        As3      C_31     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        As3      C_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        As3      C_22     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        As3      C_22     B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      C_22     C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      C_22     C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      C_22     C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      C_22     C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      C_22     C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      C_22     C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      C_22     N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      C_22     N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      C_22     O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      C_22     O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      C_21     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        As3      C_21     B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      C_21     C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      C_21     C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      C_21     C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      C_21     C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      C_21     C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      C_21     C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      C_21     N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      C_21     N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      C_21     O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      C_21     O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      C_2      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        As3      C_2      B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      C_2      C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      C_2      C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      C_2      C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      C_2      C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      C_2      C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      C_2      C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      C_2      N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      C_2      N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      C_2      O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      C_2      O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      C_R2     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        As3      C_R2     B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      C_R2     C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      C_R2     C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      C_R2     C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      C_R2     C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      C_R2     C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      C_R2     C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      C_R2     N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      C_R2     N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      C_R2     O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      C_R2     O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      C_R1     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        As3      C_R1     B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      C_R1     C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      C_R1     C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      C_R1     C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      C_R1     C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      C_R1     C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      C_R1     C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      C_R1     N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      C_R1     N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      C_R1     O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      C_R1     O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      C_R      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        As3      C_R      B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      C_R      C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      C_R      C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      C_R      C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      C_R      C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      C_R      C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      C_R      C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      C_R      N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      C_R      N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      C_R      O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      C_R      O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      N_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        As3      N_2      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        As3      N_2      B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      N_2      C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      N_2      C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      N_2      C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      N_2      C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      N_2      C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      N_2      C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      N_2      N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      N_2      N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      N_2      O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      N_2      O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      N_R      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        As3      N_R      B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      N_R      C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      N_R      C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      N_R      C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      N_R      C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      N_R      C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      N_R      C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      N_R      N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      N_R      N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      N_R      O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      N_R      O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      O_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        As3      O_2      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        As3      O_2      B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      O_2      C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      O_2      C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      O_2      C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      O_2      C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      O_2      C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      O_2      C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      O_2      N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      O_2      N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      O_2      O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      O_2      O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      O_R      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        As3      O_R      B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      O_R      C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      O_R      C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      O_R      C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      O_R      C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      O_R      C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      O_R      C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      O_R      N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      O_R      N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      O_R      O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      O_R      O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        As3      Al3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        As3      Si3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        As3      P_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        As3      S_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        As3      Ga3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        As3      Ge3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        As3      As3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Se3      B_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Se3      B_2      X           DIHEDRAL       2.0000    2.0000    1.0000
 X        Se3      C_34     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Se3      C_33     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Se3      C_32     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Se3      C_31     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Se3      C_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Se3      C_22     X           DIHEDRAL       2.0000    2.0000    1.0000
 X        Se3      C_21     X           DIHEDRAL       2.0000    2.0000    1.0000
 X        Se3      C_2      X           DIHEDRAL       2.0000    2.0000    1.0000
 X        Se3      C_R2     X           DIHEDRAL       2.0000    2.0000    1.0000
 X        Se3      C_R1     X           DIHEDRAL       2.0000    2.0000    1.0000
 X        Se3      C_R      X           DIHEDRAL       2.0000    2.0000    1.0000
 X        Se3      N_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Se3      N_2      X           DIHEDRAL       2.0000    2.0000    1.0000
 X        Se3      N_R      X           DIHEDRAL       2.0000    2.0000    1.0000
 X        Se3      O_3      X           DIHEDRAL       2.0000    2.0000   -1.0000
 X        Se3      O_2      X           DIHEDRAL       2.0000    2.0000    1.0000
 X        Se3      O_R      X           DIHEDRAL       2.0000    2.0000    1.0000
 X        Se3      Al3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Se3      Si3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Se3      P_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Se3      S_3      X           DIHEDRAL       2.0000    2.0000   -1.0000
 X        Se3      Ga3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Se3      Ge3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Se3      As3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Se3      Se3      X           DIHEDRAL       2.0000    2.0000   -1.0000
 X        In3      B_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        In3      B_2      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        In3      B_2      B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      B_2      C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      B_2      C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      B_2      C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      B_2      C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      B_2      C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      B_2      C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      B_2      N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      B_2      N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      B_2      O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      B_2      O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      C_34     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        In3      C_33     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        In3      C_32     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        In3      C_31     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        In3      C_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        In3      C_22     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        In3      C_22     B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      C_22     C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      C_22     C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      C_22     C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      C_22     C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      C_22     C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      C_22     C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      C_22     N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      C_22     N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      C_22     O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      C_22     O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      C_21     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        In3      C_21     B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      C_21     C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      C_21     C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      C_21     C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      C_21     C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      C_21     C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      C_21     C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      C_21     N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      C_21     N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      C_21     O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      C_21     O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      C_2      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        In3      C_2      B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      C_2      C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      C_2      C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      C_2      C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      C_2      C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      C_2      C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      C_2      C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      C_2      N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      C_2      N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      C_2      O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      C_2      O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      C_R2     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        In3      C_R2     B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      C_R2     C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      C_R2     C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      C_R2     C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      C_R2     C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      C_R2     C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      C_R2     C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      C_R2     N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      C_R2     N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      C_R2     O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      C_R2     O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      C_R1     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        In3      C_R1     B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      C_R1     C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      C_R1     C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      C_R1     C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      C_R1     C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      C_R1     C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      C_R1     C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      C_R1     N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      C_R1     N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      C_R1     O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      C_R1     O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      C_R      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        In3      C_R      B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      C_R      C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      C_R      C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      C_R      C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      C_R      C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      C_R      C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      C_R      C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      C_R      N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      C_R      N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      C_R      O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      C_R      O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      N_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        In3      N_2      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        In3      N_2      B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      N_2      C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      N_2      C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      N_2      C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      N_2      C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      N_2      C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      N_2      C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      N_2      N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      N_2      N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      N_2      O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      N_2      O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      N_R      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        In3      N_R      B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      N_R      C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      N_R      C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      N_R      C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      N_R      C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      N_R      C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      N_R      C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      N_R      N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      N_R      N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      N_R      O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      N_R      O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      O_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        In3      O_2      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        In3      O_2      B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      O_2      C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      O_2      C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      O_2      C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      O_2      C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      O_2      C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      O_2      C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      O_2      N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      O_2      N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      O_2      O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      O_2      O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      O_R      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        In3      O_R      B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      O_R      C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      O_R      C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      O_R      C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      O_R      C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      O_R      C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      O_R      C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      O_R      N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      O_R      N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      O_R      O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      O_R      O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        In3      Al3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        In3      Si3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        In3      P_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        In3      S_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        In3      Ga3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        In3      Ge3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        In3      As3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        In3      Se3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        In3      In3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Sn3      B_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Sn3      B_2      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Sn3      B_2      B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      B_2      C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      B_2      C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      B_2      C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      B_2      C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      B_2      C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      B_2      C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      B_2      N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      B_2      N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      B_2      O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      B_2      O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      C_34     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Sn3      C_33     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Sn3      C_32     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Sn3      C_31     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Sn3      C_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Sn3      C_22     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Sn3      C_22     B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      C_22     C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      C_22     C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      C_22     C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      C_22     C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      C_22     C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      C_22     C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      C_22     N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      C_22     N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      C_22     O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      C_22     O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      C_21     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Sn3      C_21     B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      C_21     C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      C_21     C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      C_21     C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      C_21     C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      C_21     C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      C_21     C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      C_21     N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      C_21     N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      C_21     O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      C_21     O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      C_2      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Sn3      C_2      B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      C_2      C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      C_2      C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      C_2      C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      C_2      C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      C_2      C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      C_2      C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      C_2      N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      C_2      N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      C_2      O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      C_2      O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      C_R2     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Sn3      C_R2     B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      C_R2     C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      C_R2     C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      C_R2     C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      C_R2     C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      C_R2     C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      C_R2     C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      C_R2     N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      C_R2     N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      C_R2     O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      C_R2     O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      C_R1     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Sn3      C_R1     B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      C_R1     C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      C_R1     C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      C_R1     C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      C_R1     C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      C_R1     C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      C_R1     C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      C_R1     N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      C_R1     N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      C_R1     O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      C_R1     O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      C_R      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Sn3      C_R      B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      C_R      C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      C_R      C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      C_R      C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      C_R      C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      C_R      C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      C_R      C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      C_R      N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      C_R      N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      C_R      O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      C_R      O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      N_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Sn3      N_2      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Sn3      N_2      B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      N_2      C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      N_2      C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      N_2      C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      N_2      C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      N_2      C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      N_2      C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      N_2      N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      N_2      N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      N_2      O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      N_2      O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      N_R      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Sn3      N_R      B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      N_R      C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      N_R      C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      N_R      C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      N_R      C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      N_R      C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      N_R      C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      N_R      N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      N_R      N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      N_R      O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      N_R      O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      O_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Sn3      O_2      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Sn3      O_2      B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      O_2      C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      O_2      C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      O_2      C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      O_2      C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      O_2      C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      O_2      C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      O_2      N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      O_2      N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      O_2      O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      O_2      O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      O_R      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Sn3      O_R      B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      O_R      C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      O_R      C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      O_R      C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      O_R      C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      O_R      C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      O_R      C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      O_R      N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      O_R      N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      O_R      O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      O_R      O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sn3      Al3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Sn3      Si3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Sn3      P_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Sn3      S_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Sn3      Ga3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Sn3      Ge3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Sn3      As3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Sn3      Se3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Sn3      In3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Sn3      Sn3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Sb3      B_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Sb3      B_2      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Sb3      B_2      B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      B_2      C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      B_2      C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      B_2      C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      B_2      C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      B_2      C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      B_2      C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      B_2      N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      B_2      N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      B_2      O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      B_2      O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      C_34     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Sb3      C_33     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Sb3      C_32     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Sb3      C_31     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Sb3      C_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Sb3      C_22     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Sb3      C_22     B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      C_22     C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      C_22     C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      C_22     C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      C_22     C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      C_22     C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      C_22     C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      C_22     N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      C_22     N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      C_22     O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      C_22     O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      C_21     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Sb3      C_21     B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      C_21     C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      C_21     C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      C_21     C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      C_21     C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      C_21     C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      C_21     C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      C_21     N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      C_21     N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      C_21     O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      C_21     O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      C_2      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Sb3      C_2      B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      C_2      C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      C_2      C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      C_2      C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      C_2      C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      C_2      C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      C_2      C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      C_2      N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      C_2      N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      C_2      O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      C_2      O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      C_R2     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Sb3      C_R2     B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      C_R2     C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      C_R2     C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      C_R2     C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      C_R2     C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      C_R2     C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      C_R2     C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      C_R2     N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      C_R2     N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      C_R2     O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      C_R2     O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      C_R1     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Sb3      C_R1     B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      C_R1     C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      C_R1     C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      C_R1     C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      C_R1     C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      C_R1     C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      C_R1     C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      C_R1     N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      C_R1     N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      C_R1     O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      C_R1     O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      C_R      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Sb3      C_R      B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      C_R      C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      C_R      C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      C_R      C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      C_R      C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      C_R      C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      C_R      C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      C_R      N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      C_R      N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      C_R      O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      C_R      O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      N_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Sb3      N_2      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Sb3      N_2      B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      N_2      C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      N_2      C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      N_2      C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      N_2      C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      N_2      C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      N_2      C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      N_2      N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      N_2      N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      N_2      O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      N_2      O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      N_R      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Sb3      N_R      B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      N_R      C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      N_R      C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      N_R      C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      N_R      C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      N_R      C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      N_R      C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      N_R      N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      N_R      N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      N_R      O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      N_R      O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      O_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Sb3      O_2      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Sb3      O_2      B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      O_2      C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      O_2      C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      O_2      C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      O_2      C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      O_2      C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      O_2      C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      O_2      N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      O_2      N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      O_2      O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      O_2      O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      O_R      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Sb3      O_R      B_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      O_R      C_22        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      O_R      C_21        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      O_R      C_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      O_R      C_R2        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      O_R      C_R1        DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      O_R      C_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      O_R      N_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      O_R      N_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      O_R      O_2         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      O_R      O_R         DIHEDRAL       1.0000    6.0000    1.0000
 X        Sb3      Al3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Sb3      Si3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Sb3      P_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Sb3      S_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Sb3      Ga3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Sb3      Ge3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Sb3      As3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Sb3      Se3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Sb3      In3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Sb3      Sn3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Sb3      Sb3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Te3      B_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Te3      B_2      X           DIHEDRAL       2.0000    2.0000    1.0000
 X        Te3      C_34     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Te3      C_33     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Te3      C_32     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Te3      C_31     X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Te3      C_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Te3      C_22     X           DIHEDRAL       2.0000    2.0000    1.0000
 X        Te3      C_21     X           DIHEDRAL       2.0000    2.0000    1.0000
 X        Te3      C_2      X           DIHEDRAL       2.0000    2.0000    1.0000
 X        Te3      C_R2     X           DIHEDRAL       2.0000    2.0000    1.0000
 X        Te3      C_R1     X           DIHEDRAL       2.0000    2.0000    1.0000
 X        Te3      C_R      X           DIHEDRAL       2.0000    2.0000    1.0000
 X        Te3      N_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Te3      N_2      X           DIHEDRAL       2.0000    2.0000    1.0000
 X        Te3      N_R      X           DIHEDRAL       2.0000    2.0000    1.0000
 X        Te3      O_3      X           DIHEDRAL       2.0000    2.0000   -1.0000
 X        Te3      O_2      X           DIHEDRAL       2.0000    2.0000    1.0000
 X        Te3      O_R      X           DIHEDRAL       2.0000    2.0000    1.0000
 X        Te3      Al3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Te3      Si3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Te3      P_3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Te3      S_3      X           DIHEDRAL       2.0000    2.0000   -1.0000
 X        Te3      Ga3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Te3      Ge3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Te3      As3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Te3      Se3      X           DIHEDRAL       2.0000    2.0000   -1.0000
 X        Te3      In3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Te3      Sn3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Te3      Sb3      X           DIHEDRAL       2.0000    3.0000   -1.0000
 X        Te3      Te3      X           DIHEDRAL       2.0000    2.0000   -1.0000
END
#
INVERSIONS
#
# Added explicit IGNORE for N_3, P_3, S_3, As3 and Sb3 [Jan, 29Oct96]
#
 B_2      X        X        X           UMBRELLA      40.0000    0.0000
 C_31     X        X        X           UMBRELLA      40.0000   54.7360
 C_22     X        X        X           UMBRELLA      40.0000    0.0000
 C_21     X        X        X           UMBRELLA      40.0000    0.0000
 C_2      X        X        X           UMBRELLA      40.0000    0.0000
 C_R2     X        X        X           UMBRELLA      40.0000    0.0000
 C_R1     X        X        X           UMBRELLA      40.0000    0.0000
 C_R      X        X        X           UMBRELLA      40.0000    0.0000
 N_3      X        X        X           IGNORE
 N_2      X        X        X           UMBRELLA      40.0000    0.0000
 N_R      X        X        X           UMBRELLA      40.0000    0.0000
 O_2      X        X        X           UMBRELLA      40.0000    0.0000
 O_R      X        X        X           UMBRELLA      40.0000    0.0000
 P_3      X        X        X           IGNORE
 S_3      X        X        X           IGNORE
 As3      X        X        X           IGNORE
 Sb3      X        X        X           IGNORE
END
#
UREY_BRADLEY
END
#
HYDROGEN_BONDS
 X        X           LJ_12_10     0.4000E+01    2.7500
END
#
GENERATOR
END