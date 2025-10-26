'''
fdf     siesta input files
    STRUCT.fdf
    BASIS.fdf
    KPT.fdf
    RUN.fdf

'''
ordered_fdf_keys=['SystemName','SystemLabel','NumberOfAtoms', 'NumberOfSpecies', 'ChemicalSpeciesLabel',\
                    'LatticeConstant', 'LatticeVectors', 'AtomicCoordinatesFormat', 'AtomicCoordinatesAndAtomicSpecies'\
                    'kgrid_Monkhorst_Pack', 'PAO.BasisType', 'PAO.BasisSize', 'PAO.EnergyShift', 'PAO.SplitNorm'\
                    'XC.functional', 'XC.authors', 'MeshCutoff', 'MaxSCFIterations', 'SCF.DM.Converge', 'SCF.DM.Tolerance',\
                    'SCF.Mixer.Method', 'SCF.Mixer.Weight', 'SCF.Mixer.History', 'DM.UseSaveDM', 'NeglNonOverlapInt', 'SolutionMethod', 'ElectronicTemperature',\
                    'MD.TypeOfRun', ''  ]