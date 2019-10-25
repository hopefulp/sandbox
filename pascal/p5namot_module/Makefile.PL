use ExtUtils::MakeMaker;
     WriteMakefile(
             'NAME'    => 'p5namot',                 # Name of package
             'LIBS'    => ['-L/usr/X11R6/lib -L/exec/python/lib -L/opt/intel/intel_fc_80/lib -L/ul/tpascal/programs/lib -lreadline -lncurses -lutil -ldl -lpthread -lpng -lm -lXm -lXt -lXext -lX11 -lSM -lICE -lgsl -lgslcblas'],                    # Name of custom libraries
             'OBJECT'  => 'p5namot_wrap.o '  # Object files
     );
