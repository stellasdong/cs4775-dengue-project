# dict = {1: {1: 0.0, 2: 7.0, 3: 8.0, 4: 11.0, 5: 13.0, 6: 16.0, 7: 13.0, 8: 17.0}, 
#         2: {1: 7.0, 2: 0.0, 3: 5.0, 4: 8.0, 5: 10.0, 6: 13.0, 7: 10.0, 8: 14.0}, 
#         3: {1: 8.0, 2: 5.0, 3: 0.0, 4: 7.0, 5: 10.0, 6: 10.0, 7: 7.0, 8: 11.0},
#         4: {1: 11.0, 2: 8.0, 3: 5.0, 4: 0.0, 5: 8.0, 6: 11.0, 7: 8.0, 8: 12.0},
#         5: {1: 13.0, 2: 10.0, 3: 7.0, 4: 8.0, 5: 0.0, 6: 5.0, 7: 6.0, 8: 10.0},
#         6: {1: 16.0, 2: 13.0, 3: 10.0, 4: 11.0, 5: 5.0, 6: 0.0, 7: 9.0, 8: 13.0},
#         7: {1: 13.0, 2: 10.0, 3: 7.0, 4: 8.0, 5: 6.0, 6: 9.0, 7: 0.0, 8: 8.0},
#         8: {1: 17.0, 2: 14.0, 3: 11.0, 4: 12.0, 5: 10.0, 6: 13.0, 7: 8.0, 8: 0.0}}

# d_sum = {}
# for row in dict:
#     for col in dict[row]:
#         if col not in d_sum:
#             d_sum[col] = 0
#         d_sum[col] += dict[row][col]
# # print(d_sum)

# def find_Q(D, d_sum):
#     Q = {}
#     m = None
#     m_index = (0, 0)
#     for row in D:
#         dict = {}
#         for col in D[row]:
#             if row > col:
#                 dict[col] = (len(D) - 2)*D[row][col] - d_sum[row] - d_sum[col]
#                 if m is None or dict[col] < m:
#                     m = dict[col]
#                     m_index = (col, row)
#         if row != 1:
#             Q[row] = dict
#     # print(Q)

#     return m, m_index

# # print(find_Q(dict, d_sum))

# x = 1
# y = 2
# new_row = {}
# new_row[x] = (0.5*dict[x][y])+(d_sum[x]-d_sum[y])/(2*(len(dict)-2))
# new_row[y] = (0.5*dict[x][y])+(d_sum[y]-d_sum[x])/(2*(len(dict)-2))
# for k in range(1, len(dict)+1):
#     if k != x and k != y:
#         new_row[k] = 0.5*(dict[x][k]+dict[y][k]-dict[x][y])

# # print(new_row)
# def add_row(uD, new_row, z):
#     uD[z] = new_row
#     for row in uD:
#         if row != z:
#             uD[row][z] = new_row[row]
#         if row == z:
#             uD[row][z] = 0

# add_row(dict, new_row, 9)

# def remove_rows(uD, x, y):
#     temp_D = uD
#     del temp_D[x]
#     del temp_D[y]
#     print(temp_D)
#     for row in temp_D:
#         del temp_D[row][x]
#         del temp_D[row][y]
#     return temp_D

# temp_D = remove_rows(dict, 1, 2)
# print(dict)
from copy import deepcopy

dict = {1: {1: 0.0, 2: 7.0, 3: 8.0, 4: 11.0, 5: 13.0, 6: 16.0, 7: 13.0, 8: 17.0}, 
        2: {1: 7.0, 2: 0.0, 3: 5.0, 4: 8.0, 5: 10.0, 6: 13.0, 7: 10.0, 8: 14.0}, 
        3: {1: 8.0, 2: 5.0, 3: 0.0, 4: 5.0, 5: 7.0, 6: 10.0, 7: 7.0, 8: 11.0},
        4: {1: 11.0, 2: 8.0, 3: 5.0, 4: 0.0, 5: 8.0, 6: 11.0, 7: 8.0, 8: 12.0},
        5: {1: 13.0, 2: 10.0, 3: 7.0, 4: 8.0, 5: 0.0, 6: 5.0, 7: 6.0, 8: 10.0},
        6: {1: 16.0, 2: 13.0, 3: 10.0, 4: 11.0, 5: 5.0, 6: 0.0, 7: 9.0, 8: 13.0},
        7: {1: 13.0, 2: 10.0, 3: 7.0, 4: 8.0, 5: 6.0, 6: 9.0, 7: 0.0, 8: 8.0},
        8: {1: 17.0, 2: 14.0, 3: 11.0, 4: 12.0, 5: 10.0, 6: 13.0, 7: 8.0, 8: 0.0}}

def neighbor_join(D, og):
    ''' Complete this function. '''
    uD = deepcopy(D)
    temp_D = deepcopy(D)
    E = []
    removed = []

    for z in range(len(D)+1, len(D)+(len(D)-1)):
        n = len(temp_D)
        d_sum = column_sums(temp_D)
        if n > 3:
            x, y = find_Q(temp_D, d_sum)
            E += [(z, x), (z, y)]
            removed += [x, y]
            new_row = {}
            new_row[x] = (0.5*uD[x][y])+(d_sum[x]-d_sum[y])/(2*(n-2))
            new_row[y] = (0.5*uD[x][y])+(d_sum[y]-d_sum[x])/(2*(n-2))
            # print(removed)
            for k in range(1, z):
                if k not in removed:
                    # print(x, y, k)
                    # print(uD)
                    new_row[k] = 0.5*(uD[x][k]+uD[y][k]-uD[x][y])
            
            add_row(uD, new_row, z)
            add_row(temp_D, new_row, z)
            # print(x, y)
            # print(uD)
            temp_D = remove_rows(temp_D, x, y)
        if n == 3:
            x = list(temp_D.keys())[0]
            y = list(temp_D.keys())[1]
            k = list(temp_D.keys())[2]
            E += [(x, z), (y, z), (k, z)]
            new_row = {}
            new_row[x] = (0.5*uD[x][y])+(d_sum[x]-d_sum[y])/(2*(n-2))
            new_row[y] = uD[x][y]-new_row[x]
            new_row[k] = 0.5*(uD[x][k]+uD[y][k]-uD[x][y])
            add_row(uD, new_row, z)
            add_row(temp_D, new_row, z)
    fake_root = len(uD)+1
    node = 0
    dist = 0
    # print(uD)
    for x, y in E:
        if y == og:
            node = x
            dist = uD[x][y]
            E.remove((x, y))

    E += [(fake_root, og), (node, fake_root)]
    uD[fake_root] = {og: dist/2.0, node: dist/2.0}
        
    return E, uD, fake_root
#wrong at: 12, 11, 14

def remove_rows(uD, x, y):
    temp_D = deepcopy(uD)
    del temp_D[x]
    del temp_D[y]
    for row in temp_D:
        if x in temp_D[row]:
            del temp_D[row][x]
        if y in temp_D[row]:
            del temp_D[row][y]
    return temp_D

def add_row(uD, new_row, z):
    uD[z] = new_row
    for row in new_row:
        if row != z:
            uD[row][z] = new_row[row]
    uD[z][z] = 0
    # print(uD)

def column_sums(D):
    d_sum = {}
    for row in D:
        for col in D[row]:
            if col not in d_sum:
                d_sum[col] = 0
            d_sum[col] += D[row][col]
    return d_sum

def find_Q(D, d_sum):
    Q = {}
    m = None
    m_index = (0, 0)
    for row in D:
        dict = {}
        for col in D[row]:
            if row > col:
                dict[col] = (len(D) - 2)*D[row][col] - d_sum[row] - d_sum[col]
                if m is None or dict[col] < m:
                    m = dict[col]
                    m_index = (col, row)
        if row != 1:
            Q[row] = dict
        if len(D) == 3:
            m_index = (list(D.keys())[0], list(D.keys())[2])
    print(Q)
    return m_index[0], m_index[1]

print(neighbor_join(dict, 8))