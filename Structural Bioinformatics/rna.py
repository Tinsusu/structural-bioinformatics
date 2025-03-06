import numpy as np
import math
import matplotlib.pyplot as plt


def init_matrix(rna, l):
    ''' initialize a new matrix to the length of seq '''
    m = [[0 for i in range(l)] for j in range(l)]
    return m

def fill_matrix(m,h,l,scores,x,seq):   # max of all
    ''' find max(down, left, diag, bifurcation term) then fill in the matrix'''
    for j0 in range(h + 1, l):  # the diagonal of the matrix to loop over (h = min_loop)
        for i in range(0, l - j0):  # the entry on the diagonal to fill
            j = i + j0
            # rule 1) i,j paired // diagnal
            if seq[i] + seq[j] in scores and (x[i] + x[j]) / 2 < 0.5:
                m[i][j] = m[i + 1][j - 1] + scores[seq[i] + seq[j]]
            # rule 2) i unpaired // down
            if m[i + 1][j] > m[i][j]:
                m[i][j] = m[i + 1][j]
            # rule 3) j unpaired // left side
            if m[i][j - 1] > m[i][j]:
                m[i][j] = m[i][j - 1]
            # rule 4) bifurcation k
            for k in range(i + 1 + h, j - 1 - h):
                if m[i][k] + m[k + 1][j] > m[i][j]:
                    m[i][j] = m[i][k] + m[k + 1][j]
    return m


def backtracking(m,l,x,seq,h):
    ''' Tracking back to find the path and construct dot bracket notation'''
    structure = ['.' for i in range(l)]
    stack = []
    stack.append((0, l - 1))
    while len(stack) > 0:
        top = stack.pop(),
        i = top[0][0]
        j = top[0][1]
        if i >= j:
            continue
        elif m[i + 1][j] == m[i][j]:
            stack.append((i + 1, j))
        elif m[i][j - 1] == m[i][j]:
            stack.append((i, j - 1))
        elif seq[i] + seq[j] in scores  and  m[i + 1][j - 1] + scores[seq[i] + seq[j]] == m[i][j]:
            # record basepair i,j
            structure[i] = "("
            structure[j] = ")"
            stack.append((i + 1, j - 1))
        else:
            for k in range(i + 1 + h, j - 1 - h):
                if m[i][k] + m[k + 1][j] == m[i][j]:
                    stack.append((k + 1, j))
                    stack.append((i, k))
                    break

    # print out the output
    print(seq, l, "\n", ''.join(structure), "\n", "Score:", m[0][l - 1])


###############################################

def identifyBP(seq):
    "Identify sequence"
    stack = []
    pairs= []
    for index,x in enumerate(seq):
        if x == "(":
            stack.append(index)
        if x == ')':
            pairs.append((stack.pop(),index))
    return pairs
def count_BP(list_1, list_2):
    ''' Count base pair distance'''
    my_list = [(a,b) for (a,b) in list_2 for (c,d) in list_1  if ((a==c) and (b==d))]
    match = []
    mismatch= []
    for i in list_1:
        if i in my_list:
            match.append(i)
        else:
            mismatch.append(i)
    return len(mismatch)

def hamming(list1,list2):
    ''' calculate hamming distance'''
    return len(set(list_1)^ set(list_2))

if __name__ == "__main__":
    scores = {'AU': 1, 'UA': 1, 'GU': 1, 'UG': 1, 'GC': 1, 'CG': 1}
    min_loop = 3
    seq = 'AUCUAUAUAGUAUAAAAGUAUAUUUGACUUCCAAUCAUAAGGUCUAUUAAUUAAUAGUAUAGAUA'
    x = [0.364523623524748, 0.153212539431237, 0.138223862760304, 0.0376615921718191, 0.47365746932932,
         0.208168769364413, 0.0438007380100043, 0.584052272120477, 0.523366820420774, 0.244384944817625,
         0.337766067947902, 0.229236246431118, 0.0543621087131161, 0.573269366827951, 0.626209832342058,
         0.885664923189271, 0.874654380735221, 0.178273558325159, 0.110588265614956, 0.328783306912104,
         0.302229944068083, 0.983005549967608, 0.18902671734743, 0.161354366564478, 0.0117015350480063,
         0.416886995342026, 0.296000971040533, 0.87426871101091, 0.779523667221019, 0.505478192280712,
         0.997012621986956, 0.937335885951516, 0.72353937499946, 0.77760672608629, 0.0212545744246935,
         0.415371477241511, 0.00704423417969657, 0.009896413159872, 0.449106772508896, 0.714199932809056,
         0.808192409750502, 0.807601730834051, 0.683493953670032, 0.359688962617657, 0.409624211083889,
         0.0571881337868447, 0.127831234783626, 0.160595408745504, 0.97077686365972, 0.526461146450785,
         0.790212430455052, 0.9431686060156, 0.285808346263147, 0.201925877544305, 0.135556224168921, 0.270655971845496,
         0.146443052144535, 0.136682348881294, 0.191701828815317, 0.291398671601957, 0.48478737563047,
         0.196880703148715, 0.0062362593518781, 0.359099085471509, 0.653317753803888, ]
    l = len(seq)
    m = init_matrix(seq,l)
    fill_matrix(m,min_loop,l,scores,x,seq)
    print(fill_matrix(m,min_loop,l,scores,x,seq))
    backtracking(m,l,x,seq,min_loop)
    print("\n")

    # CASE1
    NC = "(((((((((((((.(...)))))((((.((..(((.((...))..))))))))))).)))))))."
    WC = "((.(.((((((((.(...)))(((((...)..))..)))...)(((((....))))))))))))."
    list_1 = identifyBP(NC)
    list_2 = identifyBP(WC)
    print("Basepare distance: ", count_BP(list_1, list_2))
    print("Hamming distance:", hamming(list_1,list_2))

