### Element-wise math

def _sum(list, value):
    sum=[]
    for x in list:
        sum.append(x+value)

    return sum        

def l_smalleq(a, b):
    if not len(a) == len(b):
        print "Error: different length a, b lists"
        exit(1)
    for x, y in zip(a, b):
        if x <= y:
            pass
        else:
            return False
    return True            
