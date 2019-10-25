import commands, os, re, string, sys


class PBS(dict):
    def __init__(self):
        dict.__init__(self);

        status, output = self.getstatusoutput('qstat -B')
        if status == 0:
            output = output.split('\n')
            fields = output[2].split()
            self.server = fields[0]
        else:
            print output
            raise Exception, "qstat -B did not respond"

    def __str__(self):
        s = ''
        for job in self:
            s += str(job) + " " + str(self[job]) + '\n'
        return s

    def getstatusoutput(self, cmd):
        status, output = commands.getstatusoutput(cmd)

        return status, output

    '''
    def findnodes(self):
        status, output = self.getstatusoutput('pbsnodes -x')
        if status is 0:
            import xmltodict
            return xmltodict.parse(output)['Data']['Node']
    '''

    def parse_nodes(self):
        nodelist = {}

        status, output = self.getstatusoutput('pbsnodes -a')
        output = output.split('\n')
        nodename = ''

        hostre = re.compile('^\S')
        for line in output:
            if hostre.match(line):
                nodename = line
                node = {}
                continue

            if line == '':
                nodelist[nodename] = node

            line.strip()
            if line is not '':
                tag, value = line.split(' = ')
                node[tag.strip()] = value.strip()

        return nodelist

    def parse_qstat(self):
        status, output = self.getstatusoutput('qstat -f')
        output = output.split('\n')

        if not output:
            return

        linere = re.compile('^    [a-zA-Z]')
        linecont = re.compile('^\t')
        jobidre = re.compile('^Job Id:')
        jobid = 0

        for line in output:
            if jobidre.search(line):
                job = {}
                tag, value = line.split(':')
                job[tag.strip()] = value.strip()
            elif linere.search(line):
                tag, value = line.split(' = ')
                tag = tag.strip()
                value = value.strip()
                job[tag] = value
            elif linecont.search(line):
                job[tag] += line.strip()
            elif not line:
                self[job['Job Id']] = job
            else:
                print('no match', line)

    def show_queue(self):
        pass

    def show_node(self, nodename):
        pass
