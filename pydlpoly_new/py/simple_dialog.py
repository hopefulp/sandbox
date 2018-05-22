import pexpect

class simple_dialog:
    
    def __init__(self, protocoll, log):
        self.log = log
        self.logfile = open(log, "w")
        self.protocoll = protocoll
        self.questions = []
        self.answers   = []
        self.verbose = False
        self.timeout = 30
        # setup
        self.parse_protocoll()
        
    def parse_protocoll(self):
        if self.verbose: print ("Parsing protocoll file %s" % self.protocoll)
        f = open(self.protocoll, "r")
        stop = False
        while not stop:
            ql = f.readline()
            if len(ql) > 0:
                qa = f.readline()
                if ql[0:1] != ">":
                    raise SyntaxError, ql
                ql = ql[2:-1]
                if qa[0:1] != "<":
                    raise SyntaxError, qa
                qa = qa[2:-1]
                self.questions.append(ql)
                self.answers.append(qa)
            else:
                stop = True
        f.close()
        return
        
    def talk_to(self, program):
        self.p = pexpect.spawn(program, timeout=self.timeout)
        if self.verbose: print "STARTING TO TALK"
        self.p.logfile = self.logfile
        for i,q in enumerate(self.questions):
            a = self.answers[i]
            if self.verbose: print "question asked : %s" % q
            ret = self.p.expect_exact(q)
            if ret != 0:
                raise IOError, "Unexpected Question!"
            if self.verbose: print "saying         : %s" % a
            self.p.sendline(a)
            self.logfile.flush()
        return
        
    def __del__(self):
        self.logfile.close()
        return