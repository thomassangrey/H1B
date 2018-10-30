Input_Files:
(in /input)
        h1b_input.csv
Python Files:
(in /src)
        event_profiler.py
        events.py
C files:
(in src/modules/query, src/modules/sort, src/modules/uitilities)

Output files:
(in /output)
        top_10_occupations.txt
        top_10_states.txt

Instructions:
i) type ./run.sh from the root directory of the installation folder
2) or type ./run_test.sh from the insight_testsuite root folder

run.sh script should build the C-API extensions developed by this author (T.Sangrey)

__________

Results:

i)  A correct PASS was obtained fromt the tests suite
ii) This project utilized a custom database project in the works since mid-summer. I've found it useful for sessionization type problems.

ii) events.py uses built-ins developed in C using Py-C-API extension tools. Its executable is built in place when you run ./run.sh

iv) The relational database is able to do O[n*log(n)] sorting, multi-sorting, single queries, batch-queries of many types and across files. Queries are of O[log(n)]. Queries are of the set relational type (range comparisons, unions, intersections, complements, etc.) or of point relational type (single value comparisons)

    A summary of relational database queries is given below:

    # set comparators
    v = events.SORT.v        # set intersection
    u = events.SORT.u        # set union
    c = events.SORT.c        # intersect complement
    lte = events.SORT.lte    # less than or equal to
    gte = events.SORT.gte    # greater than or equal to
    lt = events.SORT.lt      # less than
    gt = events.SORT.gt      # greater than 
    
    # point comparators
    n = events.SORT.n        # point comparison: nearest to
    # <, >, <=, >=, ==       # for point comparisons

vi) Scalability: To scale the H1B project, split master files if necessary into smaller files. Most likely separate files will already exist. This implementation as written in event_profiler.py starts out by analyzing data from one year (one file). It can deploy each file-related analysis to each node and then merge the returned databases called OCCS and STATES. Tools for merging and then contracting databases are given in events.py under the DB class. After deployment of each task, an important question is whether the returned analysis called STATES and OCCS, in this case, are in theselves smaller or more efficiently handled at the master node. The answer is yes, the size of each file could be hundreds of thousands of applications, however since slave nodes contract each sub database (STATES, OCCS), Their size is capped around 50 rows for STATES and a few hundred OCCS according to the number of unique SOC_CODES.


