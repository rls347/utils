
def pert75():
    cases = ['aug11','aug17','feb23']
    modeldirs = []
    for case in cases:
        modeldirs.append(case+'-control')
        for i in range(1,25):
            modeldirs.append(case+'-pert'+str(i))
    return modeldirs

def case25():
    cases = ['aug11','aug17','feb23']
    output = [[],[],[]]
    for i in range(3):
        output[i].append(cases[i]+'-control')
        for f in range(1,25):
            output[i].append(cases[i]+'-pert'+str(f))
    return output

