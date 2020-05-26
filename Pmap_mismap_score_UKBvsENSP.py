# df has ID col and protein sequences
# pos df has ID col and dictionary of positions
# merge col is shared col between two files


from ast import literal_eval # for mismap_score func
def mismap_score(df, posdf, mergecol):
    score_file = []
    # merge two df
    mer = pd.merge(df, posdf, how='inner', on=[mergecol])
    # QC merge
    print('shape ORIGINAL  df: ', df.shape)
    print('shape merge: ',mer.shape)
    print('shape ORIGINAL KEY ref: ', posdf.shape)
    print()
    for index, row in mer.iterrows():
        newrow = []
        ukb_id = row['ID']
        #ensp = row['ENSP']
        #ensp_len = row['Length']
        #ukb_len = row['UniProt_length']
        #C_abun = row['C_abundance']
        #K_abun = row['K_abundance']
        identical = row['identical_seq']
        len18 = row['Length']
        tot_tar = row['labeled_pos_count']
        pep = row['proSequence']
        posCK = row['pos_dict']
        # col order
        # newrow.append(ensp)
        newrow.append(ukb_id)
        newrow.append(identical)
        newrow.append(len18)
        #newrow.append(ukb_len)
        #newrow.append(ensp_len)
        #newrow.append(C_abun)
        #newrow.append(K_abun)
        newrow.append(tot_tar)

        # positions found counter
        count = 0
        foundC = 0
        foundK = 0
        missedC = 0
        missedK = 0

        # evaluate as dictionary
        #python_dict = literal_eval(posCK) 
        python_dict = posCK # already was dict object
        resultfound = []
        resultnotfound = []
        
        # loop thru position dict
        for key in python_dict:
            k = int(key)
            # convert pos # to index 0-based
            i = k - 1
            # if index within ensembl/or ukb protein string
            if i < int(len18):
                # get AA from peptide string
                AA = pep[i]
                # get AA from pos dic
                checker = python_dict[key]
                posid = checker + str(k)
                # match, increase count
                if AA == checker:
                    count += 1
                    resultfound.append(posid)
                    if AA == 'C':
                        foundC += 1
                    if AA == 'K':
                        foundK += 1
                if AA != checker:
                    resultnotfound.append(posid)
                    if AA == 'C':
                        missedC += 1
                    if AA == 'K':
                        missedK += 1

        # fraction correct pos
        cor_per = count/tot_tar
        miss = tot_tar - count
        # fraction missed
        miss_per = miss/tot_tar

        # add columns
        newrow.append(count) # found count
        newrow.append(miss) # missed count
        newrow.append(round(cor_per, 2))
        newrow.append(round(miss_per, 2))
        resfound = ', '.join(resultfound)
        resnotfound = ', '.join(resultnotfound)
        newrow.append(resfound)
        newrow.append(resnotfound)
        newrow.append(posCK)
        newrow.append(foundC)
        newrow.append(missedC)
        newrow.append(foundK)
        newrow.append(missedK)
        score_file.append(newrow)

    # turn list of list into dataframe
    df = pd.DataFrame(score_file)
    df.columns = ['UniProtID', 'identical_seq','len_uniprot18','total_targets','found_count','missed_count','correct_frac','missed_frac', 'positions_found', 'positions_notfound','pos_dict','found_C_count','missed_C_count','found_K_count','missed_K_count']
    #df.columns = ['UniProtID', 'identical_seq','len_uniprot18','len_ensp', 'C_abun_ukb','K_abun_ukb','total_targets','found_count','missed_count','correct_frac','missed_frac']
    return df