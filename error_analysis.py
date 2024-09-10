e4_no_char = {}

def to_list(lst_str):
    lst_str = lst_str.strip()[1:-1].split(',')
    return [s.strip()[1:-1] for s in lst_str if len(s) > 3]



with open('e4.txt', 'r') as file:
    for line in file:
        idx = line.index(' ')
        fam = line[:idx]
        lst = to_list(line[idx + 1:])
        e4_no_char[fam] = lst
    file.close()

with open('error.txt', 'w') as err_file:
    with open('e3.txt', 'r') as file:
        err_file.write("'e3 miss' contains the set of systems that did not have a characteristic domain only with error 1e-3.\n")
        err_file.write("'e4 miss' contains the set of systems that did not have a characteristic domain only with error 1e-4.\n")
        err_file.write("'both miss' contains the set of systems that did not have a characteristic domain in either query.\n\n")
        for line in file:
            idx = line.index(' ')
            fam = line[:idx]
            e3_no_char = to_list(line[idx + 1:])
            if fam not in e4_no_char:
                err_file.write(">" + fam + "\n no hit in e4.\n\n")
            elif set(e4_no_char[fam]) != set(e3_no_char):
                res = " "
                e3_miss = set(e3_no_char) - set(e4_no_char[fam])
                e4_miss = set(e4_no_char[fam]) - set(e3_no_char) 
                both_miss = set(e4_no_char[fam]).intersection(set(e3_no_char))
                if len(e3_miss) > 0:
                    res += '\ne3 miss ' + str(e3_miss) + ' '
                if len(e4_miss) > 0:
                    res += '\ne4 miss ' + str(e4_miss) + ' '
                if len(both_miss) > 0:
                    res += '\nboth miss ' + str(both_miss)
                err_file.write(">" + fam + res + '\n\n')
        file.close()

    err_file.close()