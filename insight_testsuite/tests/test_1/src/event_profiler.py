import events



    
def string_number_compound_sort(string_list, number_list):
    """Return a list of sort order keys of type list[int].prescriber_records
    
    When given a unique list of strings, and numeric list, function returns
    the sort order  with primary order dictated by number_list and secondary
    sort order dictated by string_list
    
    Usage:
        S = ['g', 'cv', 's', 'f', 'u', 'd', 'wed','er']
        N = [1,    4,   5,   5,   5,   4,    3,    1]
        SO = string_number_compound_sort(S, N)

        print("sort_order = {0}".format(SO))
        S = [S[i] for i in SO]
        N = [N[i] for i in SO]
        print("S = {0}".format(S))
        print("N = {0}".format(N))
    """
    
    # Create a string-number map
    string_orderedMap = \
        sorted(enumerate(string_list), key= lambda x: x[1],reverse=True)
    string_ordinal = [0]*len(string_list)
    for i in range(len(string_orderedMap)):
        string_ordinal[string_orderedMap[i][0]] = i
    string_map_dict = dict()
    for k in range(len(string_list)):
        string_map_dict[str(string_ordinal[k])] = k

    str_o = string_ordinal

    # Set up SORT objects based uopn string lis and number list inputs
    ss_num = list(number_list)
    sm_num = events.SORT(ss_num)
    sm_num.sortbymedian()

    ss_str_o = list(str_o)
    sm_str_o = events.SORT(ss_str_o)
    sm_str_o.sortbymedian()
    events.SORT.compound_sort(sm_num, sm_str_o)
    strings = []
    
    # Use the string_orderedMap to help keep track of input and \
    # new sort order
    for i in range(len(string_list)):
        strings.append(string_list[string_map_dict[str(sm_str_o.original[sm_num.lookup[i]])]])

    for i in range(len(string_list)):
        print("number_list[{0:2}] {1:8}= ,string_ordinal[{0:2}] = {2:8}".format(i, sm_num.series[i], strings[i]))
    
    sort_order = []
    for i in range(len(string_list)):
        sort_order.append(string_map_dict[str(sm_str_o.original[sm_num.lookup[i]])])
    return sort_order

def analyze_block(a_stream):
    a_stream.make_db()
    H1B = a_stream['H1B']
                                # Get SORT objects
    CS = H1B['CASE_STATUS']     # query result for CS.records is initially None
    SN = H1B['SOC_NAME']
    ES = H1B['WORKSITE_STATE']
                                # query result for CERTIFIED records is now CS.records
    CS == [hash('CERTIFIED')]
    cert_records = CS.records   
                                # now slice the SORT objects to contain only cert_records
    SN.slice(cert_records)
    ES.slice(cert_records)
    
    occupations = events.DB(H1B.hash_table)
    occupations.build({'SOC_NAME': SN.original, 'WORKSITE_STATE': ES.original})
    occs = occupations.contract('SOC_NAME')
    states = occupations.contract('WORKSITE_STATE')
 
    return occs, states

def compile(totals, totals_name):
    """takes a contracted DB object and a column name. Computes 
    the number tally indexed by totals_name (i.e., 'WORK_STATE').original
    Also keeps track of those applications that are certifified. Returns
    convenient totals lists and sort_order (compound alpha-numeric sort)
    for use in ranked data presentation. total_ordered_list and 
    totals_number are lists that must be indexed by each element
    of sort_order for correct ranking presentation.
    """
    totals_hash = totals.hash_table[totals_name]
    totals_lookup = totals['UNIQUE_COUNT'].lookup
    totals_ordered = [totals[totals_name].original[i] for i in totals_lookup]
    
    # Convert totals from hash values, pretty up before multi-sort
    length = len(totals_ordered)
    totals_ordered_list = [totals_hash[totals_ordered[i]] for i in range(length)]
    totals_ordered_list = [elem.strip('\"') for elem in totals_ordered_list]
    totals_number = [totals['UNIQUE_COUNT'].original[i] for i in totals_lookup]
    
    # Compute the total of certified applications
    total_certified = sum(totals_number)
    print("\ntotal_certified = {0}".format(total_certified))
    
    sort_order = string_number_compound_sort(totals_ordered_list, totals_number)
    
    return totals_ordered_list, total_certified, totals_number, sort_order

def main():
    
    # SCALABILITY STEP 1: 
    # Partition or deploy data files to slave nodes and execute
    # analyze_block at each slave node
    a_stream = events.STREAM(2016)
    occs, states = analyze_block(a_stream)
    
    # SCALABILITY STEP 2: 
    # DEPLOY
    # Collect all analyzed data in the form of occs, and states 
    # which are DB objects

    occs_ordered_list, total_certified, \
            occs_number, occs_sort_order = compile(occs, 'SOC_NAME')

    states_ordered_list, total_certified, \
            states_number, states_sort_order = compile(states, 'WORKSITE_STATE')    
   
    # SCALABILITY STEP 3:
    # COLLECT
    # Now build a new master DB and iteratively concatenate each
    # DB formed from slave node output (two compile statements above)
    # Use DB.concat()

    # SCALABILITY STEP 4:
    # CONCATENATE
    # Similar to procedures in analyze_block() contract the fully
    # concatenated DB formed of data from all slave nodes. This data will
    # be of N*m length where N is the number of nodes and m is the number
    # of unique values within each H1B category (i.e, 50 states).
    # The concatenated DB will now have only unique elements and the
    # number of occurences of each

    # Write TOP OCCUPATIONS file
    new_fo = open("./output/top_10_occupations.txt",'w')
    new_fo.write("TOP_OCCUPATIONS;NUMBER_CERTIFIED_APPLICATIONS;PERCENTAGE\n")
    for i in reversed(occs_sort_order):
        file_line = occs_ordered_list[i] + ';'\
         + str(occs_number[i]) + \
            ';' + str('{0:1.1f}%'.format(100*occs_number[i]/total_certified)) + '\n'
        new_fo.write(file_line)

    
    # Write TOP STATES file
    new_fo = open("./output/top_10_states.txt",'w')
    new_fo.write("TOP_STATES;NUMBER_CERTIFIED_APPLICATIONS;PERCENTAGE\n")
    for i in reversed(states_sort_order):
        file_line = states_ordered_list[i] + ';'\
         + str(states_number[i]) + \
            ';' + str('{0:1.1f}%'.format(100*states_number[i]/total_certified)) + '\n'
        new_fo.write(file_line)

   
    new_fo.close()
    a_stream.close()


    # UNIT_TESTS ... 
    def UT_string_number_compound_sort():
        S = ['g', 'cv', 's', 'f', 'u', 'd', 'wed','er']
        N = [1,    4,   5,   5,   5,   4,    3,    1]
        SO = string_number_compound_sort(S, N)

        print("sort_order = {0}".format(SO))
        S = [S[i] for i in SO]
        N = [N[i] for i in SO]
        print("S = {0}".format(S))
        print("N = {0}".format(N))
    
    # Execute Unit Tests
    UT_string_number_compound_sort()

if __name__ == "__main__":
    main()
print('\n...goodbye...')
