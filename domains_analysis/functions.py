# Jaccard Index
def get_borders(windows, tads):
    borders = []
    for i in range(0, len(windows)):
        if windows[1].values[i] in tads[1].values:
            borders.append(1)
        elif windows[2].values[i] in tads[2].values:
            borders.append(1)
        else:
            borders.append(0)
    return(borders)
    
def jaccard_index(vector1, vector2):
    vec1_true = [v1 for (v1, v2) in zip(vector1, vector2) if v1!=0 or v2!=0]
    vec2_true = [v2 for (v1, v2) in zip(vector1, vector2) if v1!=0 or v2!=0]

    if len(vec1_true)==0 and len(vec2_true)==0:
        return(0)
    else:
        intersection = 0
        for (v1, v2) in zip(vec1_true, vec2_true):
            if v1 == 1 and v2 == 1:
                intersection += 1
        ji = intersection / (len(vec1_true))
        return(ji)

# MoC
def overlap(a, b):
    return(max(0, min(a[1], b[1]) - max(a[0], b[0])))
    
def calculate_moc(tads1, tads2):
    score = 0
    for i in range(0, len(tads1)):
        for j in range(0, len(tads2)):
            F = overlap(tads1[i], tads2[j])
            P = tads1[i][1]-tads1[i][0]
            Q = tads2[j][1]-tads2[j][0]
            score += (F/P * F/Q)

    moc = 1/(np.sqrt(len(tads1)*len(tads2))-1) * (score-1)
    return(moc)

# Genome coverage with domains
def get_coverage(tool_name):
    covs = []
    
    for f in os.listdir('../final_cids/' + tool_name + '/'):

        cids = pd.read_csv('../final_cids/' + tool_name + '/' + f, sep = '\s+', header = None)
        cids = cids.drop([0], axis = 1).as_matrix()
        temp_tuple = [list(x) for x in cids]
        
        temp_tuple.sort(key=lambda interval: interval[0])
        merged = [temp_tuple[0]]
        for current in temp_tuple:
            previous = merged[-1]
            if current[0] <= previous[1]:
                previous[1] = max(previous[1], current[1])
            else:
                merged.append(current)
        
        cids_sizes = [x[1]-x[0] for x in merged]

        if 'Bacillus' in f:
            genome_size = 4033459
        elif 'Caulobacter' in f:
            genome_size = 4042929
        elif 'Escherichia' in f:
            genome_size = 4641652
        elif 'Mycoplasma' in f:
            genome_size = 816394
        elif 'Sulfolobus' in f:
            genome_size = 2225959
            
        coverage = np.sum(cids_sizes) / genome_size
        if coverage > 1:
            coverage = 1
        covs.append(coverage) 
        
    return(covs)

# Shared domain boundaries
def get_borders(windows, tads):
    borders = []
    for i in range(0, len(windows)):
        if windows[1].values[i] in tads[1].values:
            borders.append(1)
        elif windows[2].values[i] in tads[2].values:
            borders.append(1)
        else:
            borders.append(0)
    return(borders)
    
def count_borders(algo_main):
    
    all_borders = []
    for (rep, win) in zip(reps, wins):
        
        windows = pd.read_csv(win + '_windows.bed', sep = '\t', header = None)

        if not os.stat('../final_cids/' + algo_main + '/' + rep).st_size == 0:
            tads_main = pd.read_csv('../final_cids/' + algo_main + '/' + rep, sep = '\s+', header = None)
            borders1 = get_borders(windows, tads_main)
            borders1 = [x[0] for x in enumerate(borders1) if x[1]!=0]
        else:
            borders1 = []
        borders_count = [0]*len(borders1)

        for algo in tools:
            if algo != algo_main and not os.stat('../final_cids/' + algo + '/' + rep).st_size == 0:
                tads = pd.read_csv('../final_cids/' + algo + '/' + rep, sep = '\s+', header = None)
                borders2 = get_borders(windows, tads)
                borders2 = [x[0] for x in enumerate(borders2) if x[1]!=0]
                
                for i in range(0, len(borders1)):
                    if borders1[i] in borders2:
                        borders_count[i] += 1
        
        all_borders.append(borders_count)

    all_borders = [item for sublist in all_borders for item in sublist]
    return(all_borders)