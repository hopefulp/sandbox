class template_fit:
    
    def __init__(self):
        return

    def initialize(self, pd, ref):
        self.pd = pd
        self.ref = ref
        self.generate_reference()
        self.get_weights()
        self.cycle = 0
        return

    def __call__(self):
        self.msd, self.a_msd = self.calc_msd()
        self.cycle += 1
        return self.msd, [self.a_msd]
