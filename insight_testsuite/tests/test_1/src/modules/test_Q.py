
import query
# sys.path.append('../PyC_API_exercises')
import sort_map
import sys
sys.path.append('../') 
import events


N1=10
K = 5
S1 = list()
for idx in range(K):
	S1 = S1 + list(range(N1))

# sort_map.randomize(S1)

# SS1=list(S1)
# L1=sort_map.sort(SS1)
# BS1 = sort_map.build_bin(SS1)

# idx=-1

# for i in SS1:
#     idx+=1
#     if(idx >0):
#         B = (int(SS1[idx]) - int(SS1[idx-1]))
#         if (B<0):
#         	print("ERROR")


def string_number_compound_sort(string_list, number_list):
    """Return a list of sort order keys of type list[int].prescriber_records
    
    When given a unique list of strings, and numeric list, function returns
    the sort order  with primary order dictated by number_list and secondary
    sort order dictated by string_list
    """
   

    print(len(string_list), len(number_list))
    print("number_list = \n{0}".format(number_list))
    print("string_list = \n{0}".format(string_list))
    # A is unique, L is Not. Sort according to L then A
    # Map A elements to ordinal numbers
    string_orderedMap = sorted(enumerate(string_list), key= lambda x: x[1],reverse=False)
    print("string_orderedMap = {}".format(string_orderedMap))
    string_ordinal = [0]*len(string_list)
    for i in range(len(string_orderedMap)):
        string_ordinal[string_orderedMap[i][0]] = i
    print("string_ordinal = {}".format(string_ordinal))
    string_map_dict = dict()
    print("\nstring_map_dict = {}".format(string_map_dict))
    for k in range(len(string_list)):
        string_map_dict[str(string_ordinal[k])] = k
    print("\nstring_map_dict = {}".format(string_map_dict))

    str_o = string_ordinal

    ss_num = list(number_list)
    sm_num = events.SORT(ss_num)
    sm_num.sortbymedian()

    ss_str_o = list(str_o)
    sm_str_o = events.SORT(ss_str_o)
    sm_str_o.sortbymedian()
    sort_str_num = events.SORT.compound_sort(sm_num, sm_str_o)
    strings = []
    for i in range(len(string_list)):
        strings.append(string_list[string_map_dict[str(sm_str_o.original[sm_num.lookup[i]])]])
        # print(A_letter[i])
    for i in range(len(string_list)):
        print("number_list[{0:2}] {1:8}= ,string_ordinal[{0:2}] = {2:8}".format(i, sm_num.series[i], strings[i]))
    sort_order = []
    for i in range(len(string_list)):
        sort_order.append(string_map_dict[str(sm_str_o.original[sm_num.lookup[i]])])
    return sort_order

num_list = [6, 3,2, 5,6, 5, 5, 3, 3 ,3, 9, 4]
string_list = ['z', 'r', 'c', 'e' , 'b','t','f', 'p', 'q', 'x','y', 'm']
sort_order = string_number_compound_sort(string_list, num_list)


# letter = A[int(A_map_dict[str(smAo.original[smL.lookup[i]])])]
# print(A[A_map_dict[str(smAo.original[smL.lookup[i]])]])
# for i in range(len(Ao)):
#     a_str = A[A_map_dict[str(smAo.original[smL.lookup[i]])]]
#     print("Li[{0:2}] {2:8}= ,Ao[{0:2}] = {2:8}".\
#     	format(i, smL.series[i], a_str))

# sort_map.randomize(S1)

# SS1 = list(S1)
# SM1 = events.SORT(SS1)
# SM1.sortbymedian()

# N2 = 50
# K = 1
# S2 = list()
# for idx in range(K):
#     S2 = S2 + list(range(N2))

# sort_map.randomize(S2)
# SS2 = list(S2)
# SM2 = events.SORT(SS2)
# SM2.sortbymedian()

# S12 = events.SORT.compound_sort(SM1, SM2)
# V1 = [8, 1, 4, 5, 7]
# V2 = [2.1, 3.1]

# for i in range(len(S1)):
#     # print("SM1.lookup[{0}] = {1}".format(i, SM1.lookup[i]))
#     print("SS1[{0:2}] {1:8}= ,SS2[{0:2}] = {2:8}".format(i, SM1.series[i], SM2.original[SM1.lookup[i]]))
# relational = 'E'


# R = query.find(V1,S1,SS1,L1,BS1,V2,S2,SS2,L2,BS2)
# R1 = query.find(relational, V1,S1,SS1,L1,BS1)
# R2 = query.find(relational,V2,S2,SS2,L2,BS2, R1)