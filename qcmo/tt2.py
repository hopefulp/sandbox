def find_level_5base(qc_list, ihomo,imo_homo, tag_hl, matom):
    level_base=[]
    if tag_hl == 'homo':
        locc = qc_list[:ihomo+1]
        locc.reverse()
    elif tag_hl == 'lumo':
        locc = qc_list[ihomo:]

    for qc in locc:     # homo homo-1 homo-2 ...
        buffer = heapq.nlargest(1, qc.bas_dic.items(), lambda i: i[1])
        for mbasis in buffer:           # 1d
            if re.search(matom, mbasis[0]):
                ind = qc.imo - imo_homo
                level_base.append([ind,mbasis[0]])
                print "%s-%2d level [1st base:coeff] = %-7s %5.2f in %20s" % (tag_hl,ind, mbasis[0], mbasis[1], whereami())
                break
    return level_base

